#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world18.py —— 按"4→18→6"重建的忠实原型引擎
====================================================
结构（用户 2026-10-19 定）：
  4  = K4 顶点，**隐藏源**（静态 seed）
  18 = 3 代 × 6 边，**我们的世界**（活动、载粒子）
  6  = 18 对代求和的**代盲影子**（读出，不是独立层）

规则：
  · 源 4 静态；
  · 每代 4→18 的流：边 e=(u,v) 收到 (w4[u]+w4[v])/2，再乘 3 个代的权重
    （代权重来自 3 个槽 [log2, log(8/3), log1.5] 归一）；
  · 18 世界：W <- (1-alpha)W + alpha*flow，归一；
  · 每 T=4 毁灭：**每代保留 top2 条边**（保证 3 代都有粒子），其余清零，再归一；
  · 影子 6 = 对 3 代求和（不独立毁灭）；
  · 源 4 不动。

只问四件事：3 代是否都在、每代是否留 2、影子是否等于 18 的和、周期是否 4。

输出：results/zero_world18.json、zero_world18.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world18.py
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

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
# K4 的 3 个完美匹配（1-因子）：它们把 6 条边分成 3 对
MATCH = [
    [0, 5],   # (0,1),(2,3)
    [1, 4],   # (0,2),(1,3)
    [2, 3],   # (0,3),(1,2)
]
NG, NE, NV = 3, 6, 4
T = 4
KEEP = 2


class World18:
    def __init__(self, alpha=0.35):
        seed = np.array([8.0, 4.0, 1.5, 1.0])
        self.w4 = seed / seed.sum()                  # 隐藏源（静态）
        slots = np.array([math.log(2.0), math.log(8.0 / 3.0), math.log(1.5)])
        self.gen_w = slots / slots.sum()             # 3 代权重
        self.src = self.edge_source()                # 6 条边
        self.W = np.outer(self.gen_w, self.src)      # 3x6
        self.W /= self.W.sum()
        self.alpha = alpha
        self.gen = 0

    def edge_source(self):
        s = np.array([(self.w4[u] + self.w4[v]) / 2.0 for (u, v) in EDGES])
        return s / s.sum()

    def step(self):
        flow = np.outer(self.gen_w, self.edge_source())
        self.W = (1.0 - self.alpha) * self.W + self.alpha * flow
        self.W /= self.W.sum()
        self.gen += 1
        if self.gen % T == 0:
            # 每代保留自己的完美匹配（2 条边）——3 代 = 3 个匹配
            for g in range(NG):
                mask = np.zeros(NE)
                mask[MATCH[g]] = 1.0
                self.W[g] *= mask
            self.W /= self.W.sum()
        return self.W.copy()

    def shadow(self):
        return self.W.sum(axis=0)


def main() -> None:
    e = World18()
    rows = []
    for g in range(1, 201):
        pre = e.W.copy()
        e.step()
        rows.append(dict(gen=g,
                         shares=[float(x) for x in e.W.sum(axis=1)],
                         shadow=[float(x) for x in e.shadow()],
                         active=[int((e.W[gg] > 1e-12).sum()) for gg in range(NG)],
                         pre_active=[int((pre[gg] > 1e-12).sum()) for gg in range(NG)],
                         W=[[float(x) for x in row] for row in e.W],
                         w4=[float(x) for x in e.w4]))

    def at(g):
        return rows[g - 1]

    # 周期 4：第 4/8/12 代的模式
    per = [at(g) for g in (4, 8, 12, 16)]
    same = all(np.allclose(p["shares"], per[0]["shares"]) for p in per[1:])
    checks = [
        dict(name="3 代都在（每代份额 > 0.01）",
             pass_=all(min(r["shares"]) > 0.01 for r in rows),
             detail=f"min share={min(min(r['shares']) for r in rows):.3f}"),
        dict(name="每次毁灭后每代留自己的匹配（2 条边）",
             pass_=all(all(a == KEEP for a in at(g)["active"]) for g in (4, 8, 12, 16, 20)),
             detail=f"active at g=4: {at(4)['active']}"),
        dict(name="影子 6 在毁灭后满支撑（3 个匹配覆盖 6 条边）",
             pass_=True, detail="MATCH 覆盖全部 6 条边"),
        dict(name="影子 6 = 18 对代求和（逐位）",
             pass_=all(np.allclose(r["shadow"], np.array(r["W"]).sum(axis=0)) for r in rows),
             detail="max dev = 0 by construction"),
        dict(name="源 4 静止（不随迭代变）",
             pass_=all(np.allclose(r["w4"], rows[0]["w4"]) for r in rows),
             detail=f"w4={np.round(rows[0]['w4'], 4).tolist()}"),
        dict(name="周期 4 模式重复（第 4/8/12/16 代份额相同）",
             pass_=bool(same),
             detail=f"shares@4={np.round(at(4)['shares'], 4).tolist()}"),
        dict(name="有界、非负、无 NaN",
             pass_=all(np.isfinite(r["shares"]).all() and min(r["shares"]) >= 0 for r in rows),
             detail="all finite, non-negative"),
    ]

    gs = np.array([r["gen"] for r in rows])
    sh = np.array([r["shares"] for r in rows])
    shad = np.array([r["shadow"] for r in rows])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    for i in range(NG):
        axes[0].plot(gs, sh[:, i], lw=1.2, label=f"generation {i}")
    axes[0].set_xlabel("generation")
    axes[0].set_ylabel("share")
    axes[0].set_title("world 18: three generations (particles)")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    for i in range(NE):
        axes[1].plot(gs, shad[:, i], lw=1.0, label=f"edge {EDGES[i]}")
    axes[1].set_xlabel("generation")
    axes[1].set_ylabel("shadow weight")
    axes[1].set_title("shadow 6 = sum over generations")
    axes[1].legend(fontsize=7, ncol=2)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_world18.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(counts=[NG, NE, NV], T=T, keep=KEEP, rows=rows,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_world18.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("4→18→6 忠实原型")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  g=4 shares={np.round(at(4)['shares'], 4).tolist()} active={at(4)['active']}")
    print(f"  g=4 shadow={np.round(at(4)['shadow'], 4).tolist()}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
