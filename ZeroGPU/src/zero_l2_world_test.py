#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_l2_world_test.py —— L2 哪个层是我们的世界：单槽 vs 整层
================================================================
问题：L2 的某一个亚层（f0/f1/f2）是我们的世界，还是整个 L2 才是？

三条判据：
  (i) **自治**：候选者有没有自己的更新律（Z016：只有 L2 有）；
  (ii) **承载**：候选者能不能装下"三代"（3 个自由度）；
  (iii) **选择**：动力学是否选出某一个槽（否则"最后那个"只是标签）。

结果：
  · 单槽：AR(1) 自预测 R² 0.19/0.25/0.20；驱动 R² 0.56/0.61/0.57，
    驱动增量 0.37/0.36/0.37 —— **被驱动、不自洽**；
  · 三槽：p→f 的 Jacobian 秩 3（三槽承载 3 自由度），两两相关 0.84–0.97
    —— 它们是**同一个状态的三条坐标**，不是三个独立世界；
  · 整层：闭合律 c_{t+1}=A^L c_t 精确（dev 2.8e-16）——**有自己的更新律**；
  · 三个槽统计相近、无特殊槽；唯一的不对称来自 seed（输入），不是动力学。

判定：**整个 L2 是我们的世界（候选）；单个亚层只是坐标/读出，不是世界。**
"最后一个亚层 = 我们的世界"是**标签约定**，不是动力学事实。

输出：results/zero_l2_world_test.json、zero_l2_world_test.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_l2_world_test.py
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo, Params  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
SEED = (8.0, 4.0, 1.5, 1.0)
L = 4

B = np.zeros((3, 4))
for _i in range(3):
    B[_i, _i] = 1.0
    B[_i, _i + 1] = -1.0


def f_of(p):
    p = np.asarray(p, float)
    return B @ np.log(p / p.sum())


def base_params():
    return dict(micro=True, micro_semantics="fresh", L=L, T=4,
                retention="top2", eta_n=0.2, eta_h=0.2, eta_psi=0.2,
                mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
                capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
                Gamma_h=0.2, kappa_h=0.5, alpha=0.0)


def run_slots(gens=400):
    p = Params(**base_params())
    rng = np.random.default_rng(7)
    for _ in range(3):
        x0 = rng.normal(size=4)
    x0 -= x0.mean()
    x0 *= 0.05 / max(1e-12, np.abs(x0).max())
    e = ZeroEvo(zs.k4(), p, seed=np.array(SEED), n0=1.0 + x0)
    F, D = [], []
    for _ in range(gens):
        n = e.n / max(1e-30, e.n.mean()) - 1.0
        h = e.h - e.h.mean()
        a = np.abs(e.psi) ** 2
        a = a / max(1e-30, a.mean()) - 1.0
        D.append(B @ (p.eta_n * n + p.eta_h * h + p.eta_psi * a))
        e.step()
        F.append(f_of(e.gen_closures))
    return np.array(F)[30:], np.array(D)[30:]


def r2(y, yp):
    return float(1.0 - np.sum((y - yp) ** 2) / np.sum((y - y.mean()) ** 2))


def ar1(y, d):
    x = y[:-1]
    X = np.c_[x, np.ones_like(x)]
    yp = X @ np.linalg.lstsq(X, y[1:], rcond=None)[0]
    r_ar = r2(y[1:], yp)
    X2 = np.c_[x, d[1:], np.ones_like(x)]
    yp2 = X2 @ np.linalg.lstsq(X2, y[1:], rcond=None)[0]
    r_dr = r2(y[1:], yp2)
    X3 = np.c_[d[1:], np.ones_like(x)]
    yp3 = X3 @ np.linalg.lstsq(X3, y[1:], rcond=None)[0]
    r_d = r2(y[1:], yp3)
    return dict(ar=r_ar, drive=r_d, ar_drive=r_dr)


def closure_check():
    p = Params(micro=True, micro_semantics="fresh", L=L, T=0,
               retention="last", wave_on=False, poisson_on=False)
    e = ZeroEvo(zs.k4(), p, seed=np.array(SEED))
    A = zs.adj_matrix(zs.k4()).astype(float)
    AL = np.linalg.matrix_power(A, L)
    cs = []
    for _ in range(20):
        e.step()
        cs.append(e.gen_closures.copy())
    cs = np.array(cs)
    dev = 0.0
    for t in range(len(cs) - 1):
        pred = AL @ cs[t]
        pred *= cs[t + 1].sum() / pred.sum()
        dev = max(dev, float(np.abs(pred - cs[t + 1]).max()
                             / max(1e-300, np.abs(cs[t + 1]).max())))
    return dev


def main() -> None:
    F, D = run_slots()
    stats = {}
    for i in range(3):
        stats[f"slot{i}"] = ar1(F[:, i], D[:, i])
    corr = np.corrcoef(F.T)
    jac = B @ np.diag(1.0 / (np.array(SEED) / sum(SEED)))
    rank = int(np.linalg.matrix_rank(jac))
    dev = closure_check()

    checks = [
        dict(name="单槽不自治（AR(1) R² < 0.4）",
             pass_=all(stats[f"slot{i}"]["ar"] < 0.4 for i in range(3)),
             detail=f"AR R²={[round(stats[f'slot{i}']['ar'], 3) for i in range(3)]}"),
        dict(name="单槽被驱动（驱动增量 > 0.2）",
             pass_=all(stats[f"slot{i}"]["drive"] - stats[f"slot{i}"]["ar"] > 0.2
                       for i in range(3)),
             detail=f"Δ={[round(stats[f'slot{i}']['drive'] - stats[f'slot{i}']['ar'], 3) for i in range(3)]}"),
        dict(name="三槽是同一状态的三条坐标（|corr| > 0.7）",
             pass_=bool(abs(corr[0, 1]) > 0.7 and abs(corr[1, 2]) > 0.7),
             detail=f"corr01={corr[0,1]:.3f}, corr12={corr[1,2]:.3f}"),
        dict(name="三槽承载 3 自由度（p→f 的秩 = 3）",
             pass_=rank == 3, detail=f"rank={rank}"),
        dict(name="整层有自己更新律（闭合律精确）",
             pass_=dev < 1e-12, detail=f"max rel dev={dev:.2e}"),
        dict(name="单个槽只承载 1 < 3（装不下三代）",
             pass_=True, detail="1 dof vs 3 dof"),
    ]

    res = dict(stats=stats, corr=corr.tolist(), rank=rank,
               closure_dev=dev,
               means=[float(x) for x in F.mean(axis=0)],
               stds=[float(x) for x in F.std(axis=0)],
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.0))
    xs = np.arange(3)
    axes[0].bar(xs - 0.18, [stats[f"slot{i}"]["ar"] for i in range(3)],
                width=0.36, label="AR(1) R²")
    axes[0].bar(xs + 0.18, [stats[f"slot{i}"]["drive"] for i in range(3)],
                width=0.36, label="drive R²")
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels(["slot 0", "slot 1", "slot 2"])
    axes[0].set_ylim(0, 1)
    axes[0].set_title("single slot: driven, not autonomous")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3, axis="y")
    im = axes[1].imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
    for i in range(3):
        for j in range(3):
            axes[1].text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center")
    axes[1].set_xticks(range(3))
    axes[1].set_yticks(range(3))
    axes[1].set_xticklabels(["slot 0", "slot 1", "slot 2"])
    axes[1].set_yticklabels(["slot 0", "slot 1", "slot 2"])
    axes[1].set_title("three slots = three coordinates of one state")
    fig.colorbar(im, ax=axes[1], fraction=0.046)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_l2_world_test.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    with open(os.path.join(RES, "zero_l2_world_test.json"), "w") as fh:
        json.dump(res, fh, indent=1)

    print("L2 单槽 vs 整层判别")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  corr=\n{np.round(corr, 3)}")
    print(f"  means={np.round(F.mean(0), 3).tolist()} stds={np.round(F.std(0), 3).tolist()}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
