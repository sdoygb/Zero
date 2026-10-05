#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sum_closure_graph_theorems_check.py —— Zero 文章③【闭合图定理】的核验
=========================================================================
独立实断言（**从零重算**闭合图，不 import 任何 Zero 程序）：
  F1  邻接规则 = 交换一对相邻异号步（即 Z1 定理 1 的补偿移动）
  F2  节点数 = 项链计数 K(T)，与文章表逐项一致
  F3  平均度与文章表逐项一致，且随 T 单调增长（无界趋势）
  F4  直径单调增长
  F5  文章写明"无截断 ⇒ 平均度无界 ⇒ 非流形 ⇒ 无稳定谱维数"
  F6  文章写明它对 Z4 终端款与 D=4 no-go 的两处后果
"""
import io
import itertools
import os
import sys
from collections import deque
from math import comb, gcd

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


DOC = io.open(os.path.join(HERE, "zero_sum_closure_graph_theorems.md"),
              encoding="utf-8").read()


def canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def K_closed(T):
    g = gcd(T, T // 2)
    return sum(phi(d) * comb(T // d, T // (2 * d)) for d in range(1, g + 1) if g % d == 0) // T


def modes(T):
    return sorted({canon(w) for w in itertools.product((1, -1), repeat=T) if sum(w) == 0})


def neigh(w):
    T = len(w)
    out = set()
    for i in range(T):
        j = (i + 1) % T
        if w[i] != w[j]:                              # 相邻异号
            v = list(w)
            v[i], v[j] = v[j], v[i]                   # 交换 ⇒ 补偿移动
            c = canon(tuple(v))
            if c != w:
                out.add(c)
    return out


# ---------------------------------------------------------------- F1
head("F1  邻接规则 = 交换相邻异号步（补偿移动）")
w0 = canon((1, -1, 1, -1))
nb = neigh(w0)
check("相邻异号交换给出非空邻域", len(nb) > 0, "%s 的邻居 %d 个" % (list(w0), len(nb)))
check("每个邻居的零和与周期都不变（补偿移动保荷）",
      all(sum(x) == 0 and len(x) == len(w0) for x in nb))
check("邻接对称（无向图）", all(w0 in neigh(x) for x in nb))

# ---------------------------------------------------------------- F2
head("F2  节点数 = K(T)（文章 §2 表）")
REF_NODES = {8: 10, 10: 26, 12: 80, 14: 246}
for T, ref in sorted(REF_NODES.items()):
    M = modes(T)
    check("T=%-3d 节点数 %-4d = K(闭式) %-4d = 文章 %-4d" % (T, len(M), K_closed(T), ref),
          len(M) == K_closed(T) == ref)

# ---------------------------------------------------------------- F3
head("F3  平均度（文章 §3 表）")
REF_DEG = {8: 3.60, 10: 4.92, 12: 6.15, 14: 7.39}
degs = []
for T, ref in sorted(REF_DEG.items()):
    M = modes(T)
    E = sum(len(neigh(w)) for w in M) // 2
    d = 2 * E / len(M)
    degs.append(d)
    check("T=%-3d 平均度 %.2f = 文章 %.2f" % (T, d, ref), abs(d - ref) < 0.01,
          "边数 %d" % E)
check("平均度随 T 单调增长（无界趋势）",
      all(degs[i] < degs[i + 1] for i in range(len(degs) - 1)), "序列 %s" % [round(x, 2) for x in degs])
check("T: 8→14 平均度增长 > 2 倍", degs[-1] / degs[0] > 2.0, "%.2f 倍" % (degs[-1] / degs[0]))

# ---------------------------------------------------------------- F4
head("F4  直径单调增长")
REF_DIAM = {8: 4, 10: 6, 12: 9, 14: 12}


def diameter(T):
    M = modes(T)
    idx = {w: i for i, w in enumerate(M)}
    adj = {w: neigh(w) for w in M}
    far = 0
    for s in M:
        dist = {s: 0}
        q = deque([s])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if y not in dist:
                    dist[y] = dist[x] + 1
                    q.append(y)
        far = max(far, max(dist.values()))
    return far


for T, ref in sorted(REF_DIAM.items()):
    d = diameter(T)
    check("T=%-3d 直径 %-3d = 文章 %-3d" % (T, d, ref), d == ref)

# ---------------------------------------------------------------- F5
head("F5  文章的定理链在位")
check("写明『无截断 ⇒ 平均度无界 ⇒ 非流形 ⇒ 无稳定谱维数』",
      "平均度无界" in DOC and "不是流形的离散化" in DOC and "无稳定谱维数" in DOC)
check("写明邻接即 Z1 定理 1 的补偿移动", "Z1 定理 1 的补偿移动" in DOC)
check("引用程序自述（no stable four-dimensional plateau）",
      "no stable four-dimensional plateau" in DOC)

# ---------------------------------------------------------------- F6
head("F6  两处后果在位")
check("推论 1：Z4 终端款（条件：闭环图须给出稳定谱维数）",
      "Z4 的终端款" in DOC and "稳定谱维数" in DOC)
check("推论 2：D=4 no-go 的可执行确认", "G89" in DOC or "no-go" in DOC)
check("诚实边界写明『无界』含解析论证而非无穷极限的严格证明",
      "解析论证" in DOC)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
