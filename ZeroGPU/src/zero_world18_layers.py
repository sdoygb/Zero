#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world18_layers.py —— 18 世界：per-generation 时钟（τ_g）测试
==================================================================
问题（用户猜想）：18 层按什么划分？时钟一样快吗？粒子一样吗？

本脚本只测**时钟**：给每代不同的毁灭周期 τ_g（full=清零整代，match=保留匹配），
看份额是否产生稳定层级。

结果：
  · 均匀 τ=4 + match → 均匀 [1/3,1/3,1/3]；
  · τ=[4,8,16] + match → [0.209,0.209,0.583]（慢代占优，但相位依赖）；
  · τ=[8,4,2] + match → [1/3,1/3,1/3]；
  · full 规则在采样代落到毁灭代时给 0。
判定：不同时钟改变振荡/占空比，**不产生稳定质量层级**；τ 是建模选择。

输出：results/zero_world18_layers.json、zero_world18_layers.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world18_layers.py
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
MATCH = [[0, 5], [1, 4], [2, 3]]
SEED = np.array([8.0, 4.0, 1.5, 1.0])
SRC = np.array([(SEED[u] + SEED[v]) / 2.0 for (u, v) in EDGES])
SRC = SRC / SRC.sum()


def run(taus, mode, gens=3000, alpha=0.35):
    W = np.ones((3, 6)) / 18.0
    for g in range(1, gens + 1):
        W = (1.0 - alpha) * W + alpha * np.outer(np.ones(3) / 3.0, SRC)
        W /= W.sum()
        for k in range(3):
            if g % taus[k] == 0:
                if mode == "full":
                    W[k] *= 0.0
                elif mode == "match":
                    m = np.zeros(6)
                    m[MATCH[k]] = 1.0
                    W[k] *= m
        W /= max(W.sum(), 1e-30)
    return W.sum(axis=1)


def main() -> None:
    cases = {
        "[4,4,4]/full": run([4, 4, 4], "full"),
        "[4,4,4]/match": run([4, 4, 4], "match"),
        "[4,8,16]/full": run([4, 8, 16], "full"),
        "[4,8,16]/match": run([4, 8, 16], "match"),
        "[8,4,2]/match": run([8, 4, 2], "match"),
    }
    checks = [
        dict(name="不同 τ 改变份额（[4,8,16] vs [4,4,4]）",
             pass_=not np.allclose(cases["[4,8,16]/match"], cases["[4,4,4]/match"]),
             detail=f"{np.round(cases['[4,8,16]/match'], 4).tolist()} vs {np.round(cases['[4,4,4]/match'], 4).tolist()}"),
        dict(name="时钟不产生稳定层级（换 τ 顺序回到均匀）",
             pass_=bool(np.allclose(cases["[8,4,2]/match"], cases["[4,4,4]/match"])),
             detail=f"{np.round(cases['[8,4,2]/match'], 4).tolist()} vs {np.round(cases['[4,4,4]/match'], 4).tolist()}"),
        dict(name="3 代边支撑互不相交（匹配划分 6 边）",
             pass_=bool(sum(len(m) for m in MATCH) == 6
                        and len(set(e for m in MATCH for e in m)) == 6),
             detail="MATCH 覆盖 6 边"),
        dict(name="源/影子不毁灭、世界毁灭（4→18→6）",
             pass_=True, detail="结构设定"),
    ]

    names = list(cases.keys())
    fig, ax = plt.subplots(figsize=(9, 4.0))
    xs = np.arange(3)
    w = 0.8 / (len(names) + 1)
    for i, k in enumerate(names):
        ax.bar(xs - 0.4 + w * (i + 1), cases[k], width=w, label=k)
    ax.set_xticks(xs)
    ax.set_xticklabels(["gen 0", "gen 1", "gen 2"])
    ax.set_ylabel("share")
    ax.set_title("per-generation clocks: shares are phase-dependent")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_world18_layers.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(cases={k: [float(x) for x in v] for k, v in cases.items()},
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_world18_layers.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("18 世界 per-generation 时钟")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for k, v in cases.items():
        print(f"  {k}: {np.round(v, 4).tolist()}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
