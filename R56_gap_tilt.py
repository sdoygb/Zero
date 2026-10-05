#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R56 · 量子线：Z* 下投影 + 层高线性权重（能隙型等比衰减）

构造：
  (1) Z* 成员 = 长度 L 的闭合 ±1 路径（Σw=0）
  (2) 下投影  = 每条闭合路径的全部前缀 σ（演化层状态）
  (3) 原生权重 = 挂在 σ 上的闭合路径数 W(σ)（补全数，纯计数）
  (4) 能隙项   = 按层高的**线性**代价：ω(σ) ∝ W(σ) · r^{h(σ)}
                 r = 每爬一层的衰减因子，r = e^{-γ}
                 依据：层高每加一层付固定代价 ⇒ 等比衰减 ⇒ 唯一能过门的形状

分划：按层高 h 分层（h = 前缀的高度层级，规范不变的层指标）

判据：S_max(ω) > 2
"""
from __future__ import annotations

import json
import math
import os
from collections import Counter
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
MU = (math.sqrt(5), (5 - math.sqrt(5)) / 2, (5 - math.sqrt(5)) / 2)


def smax(ws):
    v = sorted(ws, reverse=True)[:3]
    s = sum(v)
    if s <= 0:
        return 0.0
    v = [x / s for x in v]
    return sum(a * b for a, b in zip(v, MU))


def prefix_classes(L):
    """返回 {(k,s): W}，W = 补全数。"""
    out = {}
    for k in range(L + 1):
        for s in range(-k, k + 1, 2):
            rem = L - k
            if (rem - s) % 2:
                continue
            a = (rem - s) // 2
            if a < 0 or a > rem:
                continue
            out[(k, s)] = comb(rem, a)
    return out


def height_of_prefix(k, s):
    """前缀的高度层级：把前缀当作一条从 0 出发的路径，其高度 = 当前高度的绝对值
    加上内部振幅。这里取规范不变的层指标 = 该前缀可达到的振幅下界估计：
    用 |s| 与"剩余可爬高度"的组合。取最简：h = |s|。"""
    return abs(s)


def height_of_prefix_full(L, k, s):
    """更准确：前缀 σ 的层指标 = 它作为"半条路径"的高度份额。
    用 |s| 作为层高（离开零的距离），这是规范不变且 O(1) 的量。"""
    return abs(s)


def analyse(L, r):
    pc = prefix_classes(L)
    blk = Counter()
    for (k, s), W in pc.items():
        h = height_of_prefix_full(L, k, s)
        # W 是整数计数；r^h 是等比衰减
        val = W * (r ** h)
        blk[h] += val
    ws = list(blk.values())
    tot = sum(ws)
    if tot <= 0:
        return None
    ws = [x / tot for x in ws]
    return {"L": L, "r": r, "classes": len(ws),
            "max_share": max(ws), "S_max": smax(ws),
            "q": sum(x * x for x in ws)}


if __name__ == "__main__":
    OUT = {}
    print("等比衰减扫描：ω ∝ W(σ)·r^{|s|}")
    print()
    for L in [8, 12, 16, 20]:
        print("=" * 68)
        print("L = %d" % L)
        print("=" * 68)
        print("   r      最大类占比   S_max     量子?")
        rows = []
        for i in range(0, 41):
            r = 1.0 - i * 0.025
            if r <= 0:
                continue
            d = analyse(L, r)
            if d is None:
                continue
            rows.append(d)
            mark = ""
            if d["S_max"] > 2:
                mark = "  <== 量子"
            print("  %.3f   %.5f     %.6f %s"
                  % (r, d["max_share"], d["S_max"], mark))
        OUT["L%d" % L] = rows
        ok = [d for d in rows if d["S_max"] > 2]
        if ok:
            first = min(ok, key=lambda d: -d["r"])
            best = max(rows, key=lambda d: d["S_max"])
            print("  -> 首次越门 r=%.3f（最大类 %.4f）" % (first["r"], first["max_share"]))
            print("  -> 最强 S_max=%.6f 在 r=%.3f" % (best["S_max"], best["r"]))
            print("  -> 对应 gamma = -ln r = %.4f" % -math.log(first["r"]))
            OUT.setdefault("threshold", {})["L%d" % L] = {
                "first_r": first["r"], "gamma": -math.log(first["r"]),
                "best": best}
        else:
            print("  -> 未越门（最大 S_max=%.6f）" % max(d["S_max"] for d in rows))
        print()
    with open(os.path.join(HERE, "R56_gap_tilt_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("-> R56_gap_tilt_results.json")
