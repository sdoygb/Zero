#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z3_check.py —— 【I2a 维数漂移判定】的核验
========================================
独立实断言：
  F1  q_eff = ln|V|/ln diam 单调上升，且与探针数据逐项一致
  F2  节点数 = 项链计数 K(T)；直径 ≈ T²/16（拟合检查）
  F3  q_eff 的解析发散：K ~ 2^T/T 指数 vs diam ~ T² 二次 ⇒ q_eff → ∞
  F4  闭合图**非正则**（故非顶点传递）——BFT 前提不满足
  F5  平均度无界（doubling 失效的第二面证据）
  F6  文档的判定与更正在位（固有性三面证据；Z2 §5 的错误归属已更正；I2a 仍开放）
"""
import io
import itertools
import math
import os
import sys
from collections import Counter

import json

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = ("v" if ok else "x") if (level == "ind" or LEDGER_MODE == "A" or not ok) else "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


DOC = io.open(os.path.join(HERE, "Z3_i2a_dimension_drift_verdict.md"), encoding="utf-8").read()
PROBE = json.load(io.open(os.path.join(HERE, "simulations",
                                       "zero_sum_geometry_probe_results.json"), encoding="utf-8"))
FP = PROBE["summary"]["fixed_period_graphs"]


def canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))


def modes(T):
    return sorted({canon(w) for w in itertools.product((1, -1), repeat=T) if sum(w) == 0})


def neigh(w):
    T = len(w)
    out = set()
    for i in range(T):
        j = (i + 1) % T
        if w[i] != w[j]:
            v = list(w)
            v[i], v[j] = v[j], v[i]
            c = canon(tuple(v))
            if c != w:
                out.add(c)
    return out


# ---------------------------------------------------------------- F1
head("F1  q_eff = ln|V| / ln diam（探针数据）")
REF = {8: (10, 4), 10: (26, 6), 12: (80, 9), 14: (246, 12), 16: (810, 16), 18: (2704, 20)}
REF_Q = {8: 1.661, 10: 1.818, 12: 1.994, 14: 2.216, 16: 2.415, 18: 2.638}
qs = []
for T, (n, d) in sorted(REF.items()):
    pn, pd = int(FP[str(T)]["nodes"]), int(FP[str(T)]["diameter"])
    check("T=%-3d 探针给出 节点=%d 直径=%d" % (T, pn, pd), (pn, pd) == (n, d))
    q = math.log(n) / math.log(d)
    qs.append(q)
    check("T=%-3d q_eff=%.3f（表列 %.3f）" % (T, q, REF_Q[T]), abs(q - REF_Q[T]) < 5e-4)
check("q_eff 单调上升（无平台）",
      all(qs[i] < qs[i + 1] for i in range(len(qs) - 1)),
      "%s" % [round(x, 3) for x in qs])

# ---------------------------------------------------------------- F2
head("F2  节点数 = K(T)；直径 ≈ T²/16")
for T in (8, 10, 12, 14):
    check("T=%-3d 节点数 = |modes|" % T, len(modes(T)) == REF[T][0],
          "实算 %d" % len(modes(T)))
dev = [abs(REF[T][1] - T * T / 16) for T in REF]
check("直径与 T²/16 的偏差 ≤ 0.3（六点）", max(dev) <= 0.3,
      "最大偏差 %.2f" % max(dev))

# ---------------------------------------------------------------- F3
head("F3  q_eff 的解析发散（K 指数 vs diam 二次）")
def q_analytic(T):
    return math.log(2 ** T / T) / math.log(T * T / 16)


seq = [q_analytic(T) for T in (20, 30, 50, 100)]
check("解析 q_eff 随 T 递增", all(seq[i] < seq[i + 1] for i in range(len(seq) - 1)),
      "%s" % [round(x, 2) for x in seq])
check("T=30 处 q_eff > 4", q_analytic(30) > 4, "实算 %.2f" % q_analytic(30))
check("T=100 处 q_eff > 9", q_analytic(100) > 9, "实算 %.2f" % q_analytic(100))
check("文档写明 q_eff → ∞ 与两条趋势",
      "q_{\\rm eff}" in DOC and "\\infty" in DOC and "2^{T}" in DOC)

# ---------------------------------------------------------------- F4
head("F4  闭合图非正则（故非顶点传递）")
for T in (6, 8, 10):
    M = modes(T)
    dist = dict(sorted(Counter(len(neigh(w)) for w in M).items()))
    check("T=%-3d 度数分布非单点（非正则）" % T, len(dist) > 1, "%s" % dist)
check("文档写明非正则 ⇒ BFT 前提不满足",
      "非正则" in DOC and "顶点传递" in DOC)

# ---------------------------------------------------------------- F5
head("F5  平均度无界（doubling 失效）")
degs = {}
for T in (8, 10, 12, 14):
    M = modes(T)
    E = sum(len(neigh(w)) for w in M) // 2
    degs[T] = 2 * E / len(M)
check("平均度单调上升", all(degs[a] < degs[b] for a, b in zip((8, 10, 12), (10, 12, 14))),
      "%s" % {k: round(v, 2) for k, v in degs.items()})
check("T: 8→14 平均度增长 > 2 倍", degs[14] / degs[8] > 2.0,
      "%.2f 倍" % (degs[14] / degs[8]))
check("文档把'度无界'列为第二面证据", "doubling" in DOC and "度无界" in DOC)

# ---------------------------------------------------------------- F6
head("F6  判定与自我更正在位")
check("写明对 G_T 是固有的（三面证据）", "固有" in DOC and "三面" in DOC)
check("写明'不是截断假象'", "不是截断假象" in DOC)
check("写明 Z2 §5 的归属错误已更正",
      "错误归属" in DOC and "Z2_zero_to_gr_direct_route.md" in DOC)
check("写明 I2a 仍然开放", "I2a" in DOC and "仍然开放" in DOC)
check("区分 G_T（位形空间）与 Γ（物理空间）",
      "位形" in DOC and "物理空间" in DOC)
check("引用两篇文献（CDT 谱维数 / BFT 标度极限）",
      "1411.7712" in DOC and "1203.5624" in DOC)
check("旧理论三条查证在位（D122 / K30-FAZ / cc 的维数定义）",
      "D122" in DOC and "1312.7856" in DOC and "谱维数" in DOC)
check("给出下一步命题 P1/P2/P3", "P1" in DOC and "P2" in DOC and "P3" in DOC)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
