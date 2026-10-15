#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z047 两个猜想的数值核验（动态规划，L=4..32）—— 按"猜想先算、再形式化"的纪律。

待核项（覆盖账把 Z047 的 `a_h = 1/b_h` 标为 A?：中间一步未展开）：
  (i)   b_h 恒偶   （每层词数）
  (ii)  b_H = 2    （顶层恰一个 ± 轨道）
  (iii) a_h = 1/b_h 的候选读法：两半等权 ⟹ Bloch 半径 = 1/b_h

层高定义（与 Z001/Z079 一致）：h(w) = max(M⁺, M⁻)，M⁺/M⁻ = 前缀和的最大/最小（取负）。
枚举 2^L 不可行（L=32 是 4×10^9），故用 DP：状态 (当前和, 最大前缀, 最小前缀) 计数。

用法：/usr/bin/python3 z047_layer_parity.py
输出：results/z047_layer_parity.json
"""
from __future__ import annotations

import json
import math
import pathlib
from collections import Counter

HERE = pathlib.Path(__file__).resolve()
OUT = HERE.parents[1] / "results" / "z047_layer_parity.json"


def layer_distribution(L: int) -> Counter:
    """平衡词（长度 L、和为 0）的层高分布；DP 状态 = (和, 最大前缀, 最小前缀)。"""
    dp = {(0, 0, 0): 1}
    for _ in range(L):
        ndp: dict = {}
        for (s, hi, lo), cnt in dp.items():
            for step in (1, -1):
                s2 = s + step
                key = (s2, max(hi, s2), min(lo, s2))
                ndp[key] = ndp.get(key, 0) + cnt
        dp = ndp
    c: Counter = Counter()
    for (s, hi, lo), cnt in dp.items():
        if s == 0:
            c[max(hi, -lo)] += cnt
    return c


def run(lmax: int = 32):
    rows = []
    for L in range(4, lmax + 1, 2):
        c = layer_distribution(L)
        total = sum(c.values())
        top = max(c)
        rows.append({
            "L": L,
            "words": total,
            "choose_L_half": math.comb(L, L // 2),
            "layers": {str(h): c[h] for h in sorted(c)},
            "all_even": all(v % 2 == 0 for v in c.values()),
            "top_shell_b_H": c[top],
            "a_h_candidate": {str(h): 1.0 / c[h] for h in sorted(c)},
        })
        print(f"L={L:2d} |W|={total:9d} C(L,L/2)={math.comb(L, L//2):9d} "
              f"b_h={dict(sorted(c.items()))} 全偶={rows[-1]['all_even']} b_H={c[top]}")
    return rows


def main():
    rows = run()
    summary = {
        "note": "Z047 数值核验：b_h 恒偶、b_H=2、a_h=1/b_h 候选值（DP 精确计数）",
        "L_range": [rows[0]["L"], rows[-1]["L"]],
        "count_matches_choose": all(r["words"] == r["choose_L_half"] for r in rows),
        "all_layers_even": all(r["all_even"] for r in rows),
        "all_top_shell_two": all(r["top_shell_b_H"] == 2 for r in rows),
        "rows": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print()
    print(f"|W_L| = C(L,L/2)  全部相符 : {summary['count_matches_choose']}")
    print(f"b_h 恒偶（L=4..{rows[-1]['L']}）    : {summary['all_layers_even']}")
    print(f"b_H = 2 全部成立        : {summary['all_top_shell_two']}")
    print(f"→ results/{OUT.name}")


if __name__ == "__main__":
    main()
