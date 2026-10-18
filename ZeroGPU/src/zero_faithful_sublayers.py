#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""忠实引擎里的 16/8/4 分级毁灭（第一版：内容 n 的亚亚层簿记）。

结构：4 个 K4 顶点，每个顶点 7 个亚亚层通道 =
  底 4（全纤维）+ 中 2（一个二元因子）+ 顶 1（两个因子都固定）。
  → 底层 16、中层 8、顶层 4。

毁灭（每 T=4，在 e.step() 之前）：
  · 底 4 个通道全灭；
  · 中 2 个通道只保留“最高占据顶点”的那一对，其余 3 个顶点清零；
  · 顶 1 个通道不动。
  → 顶点总量 = 该顶点 7 个通道之和；把 e.n 设为这个缩减后的总量；
     e.step() 之后，通道按固定 profile 重设为 e.n/7（脉冲式毁灭）。

对照：同一引擎的全局毁灭（T=4 清 W）作为 baseline。
输出：results/zero_faithful_sublayers.json、results/zero_faithful_sublayers.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo  # noqa: E402
from zero_evo_phase import base_params, init_random  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
T = 4
CK = []
# 每个顶点的通道权重 profile：底4、中2、顶1；按观测粗分组标定
# 目标 shares 底/中/顶 = 0.1724/0.2759/0.5517
PROFILE = np.array([0.1724 / 4] * 4 + [0.2759 / 2] * 2 + [0.5517 / 1])


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def make_full(seed=1, T_destroy=0):
    adj = zs.cycle(16)
    p = base_params(micro=True, wave_on=True, poisson_on=True, c_h=1.7,
                    g_src=1.0, kappa_h=0.5, Gamma_h=0.2, kappa_p=0.1,
                    gamma=0.2, mu=0.5, K=2.0)
    p.T = T_destroy
    e = ZeroEvo(adj, p, n0=1.0 + init_random(len(adj), seed))
    return e


def init_comp(e):
    # 每个顶点 7 个通道：0..3 底、4..5 中、6 顶；按 PROFILE 分配
    return e.n[:, None] * PROFILE[None, :]


def shares(comp):
    p = comp.sum(axis=0)
    tot = p.sum()
    return (np.array([p[0:4].sum(), p[4:6].sum(), p[6:7].sum()]) / tot
            if tot > 0 else np.zeros(3))


def main():
    steps = 400
    # 分级：T=0（关全局），手动作
    e = make_full(seed=1, T_destroy=0)
    comp = init_comp(e)
    seed_frac = 0.05
    rows = []
    for gen in range(1, steps + 1):
        destroyed = (gen % T == 0)
        pre_top = float(comp[:, 6].sum())
        if destroyed:
            top_v = int(np.argmax(e.n))
            comp[:, 0:4] = 0.0                 # 底全灭
            mask = np.ones(e.V, dtype=bool)
            mask[top_v] = False
            comp[mask, 4:6] = 0.0              # 中只留最高顶点的一对
            # 顶 comp[:,6] 不动
            n_red = comp.sum(axis=1)
            e.n = np.maximum(n_red, 0.0)
            post_d = shares(comp)
            post_d_top_abs = float(comp[:, 6].sum())
        else:
            post_d = None
            post_d_top_abs = pre_top
        e.step()
        # 保留毁灭记忆 + L1 seed 重新填充：按 comp+seed 重分配到新 n
        base = comp + seed_frac * PROFILE[None, :]
        base = np.maximum(base, 0.0)
        rs = base.sum(axis=1, keepdims=True)
        base = np.divide(base, rs, out=np.zeros_like(base), where=rs > 0)
        comp = base * e.n[:, None]
        post_f = shares(comp)
        s = e.snapshot()
        rows.append({"gen": gen, "destroyed": destroyed,
                     "post_destroy": None if post_d is None else post_d.tolist(),
                     "post_destroy_top_abs": post_d_top_abs,
                     "pre_top_abs": pre_top,
                     "post_flow": post_f.tolist(),
                     "n_std": s["n_std"], "q": s["q"],
                     "h_energy": s["h_energy"]})
        if gen <= 8 or destroyed:
            print(f"gen={gen:>3} D={int(destroyed)} "
                  f"post_d={None if post_d is None else np.round(post_d,4)} "
                  f"post_f={np.round(post_f,4)} n_std={s['n_std']:.6g}")

    # baseline：全局毁灭
    e2 = make_full(seed=1, T_destroy=T)
    base = []
    for gen in range(1, steps + 1):
        e2.step()
        s = e2.snapshot()
        base.append({"gen": gen, "n_std": s["n_std"], "q": s["q"],
                     "h_energy": s["h_energy"]})

    print("=" * 72)
    print("判据")
    print("=" * 72)
    dd = [r for r in rows if r["destroyed"]]
    ck("毁灭事件每 T=4", len(dd) == steps // T, f"n={len(dd)}")
    ck("毁灭后底份额 = 0", all(abs(r["post_destroy"][0]) < 1e-15 for r in dd), "")
    ck("毁灭后中份额 > 0（只留最高顶点的一对）",
       all(r["post_destroy"][1] > 0 for r in dd), f"min={min(r['post_destroy'][1] for r in dd):.4f}")
    ck("毁灭后顶绝对量不变（4 个都不动）",
       all(abs(r["post_destroy_top_abs"] - r["pre_top_abs"]) < 1e-12 for r in dd), "")
    ck("无 NaN / 负 n", all(np.all(np.isfinite([r["n_std"], r["q"], r["h_energy"]])) for r in rows)
       and all(r["n_std"] >= 0 for r in rows), "")
    # 分级后的 n_std 是否有周期 4（自相关/FFT 粗判）
    ns = np.array([r["n_std"] for r in rows])
    fft = np.abs(np.fft.rfft(ns - ns.mean()))
    k4 = fft[100] if len(fft) > 100 else 0.0  # 400 点，周期 4 -> 第 100 个 bin
    ck("n_std 有周期 4 分量（FFT bin 100 > 中位数）", k4 > np.median(fft[1:]), f"bin100={k4:.3g}")
    tail = rows[-20:]
    mf = np.mean([r["post_flow"] for r in tail], axis=0)
    ck("尾段层级 底<中<顶", mf[0] < mf[1] < mf[2], f"{np.round(mf,4)}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    g = np.array([r["gen"] for r in rows])
    pd = np.array([r["post_destroy"] if r["post_destroy"] is not None else [np.nan]*3 for r in rows])
    ax[0].plot(g, pd[:, 0], label="bottom (post-destroy)", lw=1)
    ax[0].plot(g, pd[:, 1], label="middle (post-destroy)", lw=1)
    ax[0].plot(g, pd[:, 2], label="top (post-destroy)", lw=1)
    ax[0].set_title("(a) graded destruction"); ax[0].legend(fontsize=8)
    ax[1].plot(g, [r["n_std"] for r in rows], label="graded", lw=1)
    ax[1].plot(g, [r["n_std"] for r in base], label="global", lw=1)
    ax[1].set_title("(b) n_std"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout()
    fp = os.path.join(RES, "zero_faithful_sublayers.png")
    fig.savefig(fp, dpi=130)

    out = {"steps": steps, "T": T, "rows": rows, "baseline": base,
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "16/8/4 integer fiber; no new dimensionful anchor"}
    with open(os.path.join(RES, "zero_faithful_sublayers.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_faithful_sublayers.json"), "和", fp)


if __name__ == "__main__":
    main()
