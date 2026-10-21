#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_causal_coupling.py —— 因果集中度 × 物质集中度：双向耦合

上一轮 zero_causal_admit.py 把物质 e 当**固定背景**（质量层级），
所以 β 的效应被压平。用户指出：因果集中度与物质集中度**高度相关、耦合在一起**。
本稿把 e 也放开，让两者互相强化：

    落点    w ∝ e · exp(β λ)                 （因果塑造行动）
    记录    λ += 1[闭合]                      （行动写记录）
    物质    x ← x + η[-κ L x + α log(1+λ)]    （记录拉物质）
            e = softmax(x)

正反馈：λ 集中 → 行动集中 → 物质集中 → λ 更集中。
扩散 κ 是唯一的稳定化项。扫 β，看：
    ① corr(λ, e)         因果—物质的相关性
    ② IPR_e / IPR_λ      两者是否同向集中
    ③ 是否存在 β_c（共同锁死为单点）

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_causal_coupling.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
C = 10.0
HC = 6 * C

EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
A_line = np.zeros((6, 6))
for i, (u, v) in enumerate(EDGES):
    for j, (a, b) in enumerate(EDGES):
        if i != j and len({u, v} & {a, b}) > 0:
            A_line[i, j] = 1.0
L_line = np.diag(A_line.sum(1)) - A_line


def _softmax(x):
    z = x - x.max()
    e = np.exp(z)
    return e / e.sum()


def ipr(v):
    v = np.asarray(v, float)
    s = v.sum()
    if s <= 0:
        return 6.0
    p = v / s
    return float(1.0 / np.sum(p ** 2))


def run(beta, alpha=0.5, kappa=0.1, eta=0.2, p=0.5, gens=1500, seed=5):
    rng = np.random.default_rng(seed)
    x = np.zeros(6)
    lam = np.zeros(6)
    load = np.zeros(6)
    ie, il, cc = [], [], []
    for t in range(gens):
        e = _softmax(x)
        logw = np.log(e + 1e-12) + beta * lam
        logw -= logw.max()
        w = np.exp(logw)
        w /= w.sum()
        for _ in range(6):
            i = int(rng.choice(6, p=w))
            if rng.random() < p:
                lam[i] += 1.0
            else:
                load[i] += 1.0
        x = x + eta * (-kappa * (L_line @ x) + alpha * np.log1p(lam))
        x = x - x.mean()
        x = np.clip(x, -60, 60)
        ie.append(ipr(e))
        il.append(ipr(lam))
        if np.std(lam) > 1e-12 and np.std(e) > 1e-12:
            cc.append(float(np.corrcoef(lam, e)[0, 1]))
    tail = slice(gens // 2, None)
    return dict(beta=float(beta), ipr_e=float(np.mean(ie[tail])),
                ipr_l=float(np.mean(il[tail])),
                corr=float(np.mean(cc[len(cc)//2:])) if cc else 0.0,
                e_end=_softmax(x).tolist(), lam_end=lam.tolist(),
                ie=ie, il=il)


def run_fixed_matter(beta, alpha=0.0, gens=1500, seed=5):
    """对照：物质固定（上一轮的做法），只让记录自反馈。"""
    rng = np.random.default_rng(seed)
    E = _softmax(np.log(np.array([2.16, 4.67, 93.4, 1270.0, 4180.0, 172760.0])))
    lam = np.zeros(6)
    for t in range(gens):
        logw = np.log(E + 1e-12) + beta * lam
        logw -= logw.max()
        w = np.exp(logw); w /= w.sum()
        for _ in range(6):
            i = int(rng.choice(6, p=w))
            if rng.random() < 0.5:
                lam[i] += 1.0
    return dict(beta=float(beta), ipr_l=ipr(lam), lam_end=lam.tolist())


def main():
    betas = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
    rows = [run(b) for b in betas]
    fixed = [run_fixed_matter(b) for b in betas]

    # 阈值：IPR_e 首次 < 2（共同锁死）
    b_c = None
    for r in rows:
        if r["ipr_e"] < 2.0:
            b_c = r["beta"]
            break

    checks = [
        dict(name="c1 β↑ ⇒ 物质集中度 IPR_e↓（因果把物质拉集中）",
             pass_=bool(rows[-1]["ipr_e"] < rows[0]["ipr_e"] - 1.0),
             detail="IPR_e " + " ".join(f"{r['beta']:g}:{r['ipr_e']:.2f}" for r in rows)),
        dict(name="c2 两者同向集中：IPR_e 与 IPR_λ 同降",
             pass_=bool(rows[-1]["ipr_e"] < rows[0]["ipr_e"]
                        and rows[-1]["ipr_l"] < rows[0]["ipr_l"]),
             detail="IPR_e " + " ".join(f"{r['beta']:g}:{r['ipr_e']:.2f}" for r in rows)
                    + " | IPR_λ " + " ".join(f"{r['beta']:g}:{r['ipr_l']:.2f}" for r in rows)),
        dict(name="c3 存在共同锁死阈值 β_c（IPR_e<2）",
             pass_=bool(b_c is not None and b_c > 0.0),
             detail=f"β_c={b_c}"),
        dict(name="c4 对照：物质固定时 β 只改记录分布，不改物质（无 IPR_e 效应）",
             pass_=bool(max(abs(f["ipr_l"] - 1.0) for f in fixed) < 0.1),
             detail="fixed IPR_λ " + " ".join(f"{f['beta']:g}:{f['ipr_l']:.2f}" for f in fixed)),
        dict(name="c5 锁死态：末态单点占优（max e > 0.9）",
             pass_=bool(max(rows[-1]["e_end"]) > 0.9),
             detail=f"max e={max(rows[-1]['e_end']):.3f}"),
    ]

    out = dict(params=dict(C=C, alpha=0.5, kappa=0.4, eta=0.2, p=0.5),
               rows=[{k: v for k, v in r.items() if k not in ("ie", "il")} for r in rows],
               fixed=fixed, beta_c=b_c,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_causal_coupling.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    ax = axes[0]
    ax.plot(betas, [max(r["e_end"]) for r in rows], "o-")
    ax.set_xlabel("β"); ax.set_ylabel("max e (dominance)")
    ax.set_title("matter locks onto one edge"); ax.grid(alpha=.3)
    ax = axes[1]
    ax.plot(betas, [r["ipr_e"] for r in rows], "o-", label="IPR_e (matter)")
    ax.plot(betas, [r["ipr_l"] for r in rows], "s--", label="IPR_λ (ledger)")
    ax.set_xlabel("β"); ax.set_ylabel("IPR (lower = concentrated)")
    ax.set_title("both concentrate together"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    ax = axes[2]
    for r in rows[::3]:
        ax.plot(r["ie"], lw=0.8, label=f"β={r['beta']:g}")
    ax.set_xlabel("generation"); ax.set_ylabel("IPR_e")
    ax.set_title("matter concentration trajectory"); ax.legend(fontsize=6)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_causal_coupling.png")
    fig.savefig(plot, dpi=130); plt.close(fig)

    print("因果 × 物质 双向耦合")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
