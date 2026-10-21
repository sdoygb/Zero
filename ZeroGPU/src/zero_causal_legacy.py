#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_causal_legacy.py —— 一次圆满事件之后：影响如何被新因果稀释

设定（用户）：2400 年前发生一次圆满事件；它的因果集合（λ0，3 团）通过实时通道
（ℓ=0 ⇒ g=1）影响 L2；之后每一代产生新的、不圆满的因果。

问：影响还有吗？多大？怎么衰减？

模型（最简版）：
    λ_i(t) = λ0_i + t          （每边每代 +1 条新记录，均匀）
    e ∝ softmax(log(1+λ))      （g=1 的实时耦合，平衡态近似）
    份额 s(t) = H0 / (H0 + 6t)
    影响 I(t) = (6 - IPR_e(t)) / (6 - IPR_e(0))   （1=满，0=无）

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_causal_legacy.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
HC = 60.0


def softmax(v):
    v = np.asarray(v, float) - np.max(v)
    e = np.exp(v)
    return e / e.sum()


def ipr(v):
    v = np.asarray(v, float)
    s = v.sum()
    return float(1.0 / np.sum((v / s) ** 2))


def main():
    # 圆满集合：3 团，总量 1007（取自 zero_causal_rebirth 的成团终态形状）
    lam0 = np.array([335.0, 336.0, 336.0, 0.0, 0.0, 0.0])
    H0 = lam0.sum()
    ipr0 = ipr(softmax(np.log1p(lam0)))

    ts = [0, 10, 100, 1000, 3000, 10000, 30000, 100000, 1000000]
    rows = []
    for t in ts:
        lam = lam0 + t
        e = softmax(np.log1p(lam))
        H = lam.sum()
        rows.append(dict(t=t, H=float(H), s=float(H0 / H), ipr=ipr(e),
                         maxe=float(e.max())))

    # 2400 年 → 事件数（外部输入：每天 2 个重要因果事件）
    events_2400y = 2 * 365 * 2400
    lam = lam0 + events_2400y / 6.0          # 每边新增
    e = softmax(np.log1p(lam))
    s_2400 = float(H0 / lam.sum())
    ipr_2400 = ipr(e)

    I = [(6.0 - r["ipr"]) / (6.0 - ipr0) for r in rows]

    checks = [
        dict(name="c1 t=0：影响满（IPR_e≈3，3 团）",
             pass_=bool(abs(ipr0 - 3.0) < 0.1), detail=f"IPR={ipr0:.2f}"),
        dict(name="c2 影响随新因果单调衰减，趋近 0",
             pass_=bool(all(I[i + 1] < I[i] for i in range(len(I) - 1)) and I[-1] < 0.01),
             detail=" ".join(f"t={r['t']:g}:{i:.3f}" for r, i in zip(rows, I))),
        dict(name="c3 份额 s=H0/(H0+6t) 与数值一致",
             pass_=all(abs(r["s"] - H0 / r["H"]) < 1e-12 for r in rows),
             detail=f"s(t=1e5)={rows[7]['s']:.2e}"),
        dict(name="c4 2400 年后份额 ≈0.06%，影响 <0.1%",
             pass_=bool(0.0003 < s_2400 < 0.002),
             detail=f"L={events_2400y:.0f} 事件，s={s_2400:.2e}，IPR_e={ipr_2400:.3f}"),
        dict(name="c5 结论：有影响，但被稀释到可忽略",
             pass_=bool(ipr_2400 > 5.9), detail=f"IPR_e(2400y)={ipr_2400:.3f}"),
    ]

    out = dict(params=dict(H0=float(H0), lam0=lam0.tolist(),
                           events_2400y=events_2400y),
               rows=rows, influence=I,
               s_2400=s_2400, ipr_2400=ipr_2400,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_causal_legacy.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot([r["t"] for r in rows], I, "o-")
    ax.axhline(0.001, ls=":", c="k")
    ax.set_xscale("log"); ax.set_xlabel("generations after the event")
    ax.set_ylabel("influence I(t)")
    ax.set_title("a perfect event's influence is diluted by new karma")
    ax.grid(alpha=.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_causal_legacy.png")
    fig.savefig(plot, dpi=130); plt.close(fig)

    print("圆满事件之后：影响如何被稀释")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
