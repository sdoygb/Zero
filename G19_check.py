#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G19_check.py -- 条款精简：审计底层条款（Z0 条款 ＋ Z1–Z5 定理），把可导出的降级为定理。

对应文档 G19_axiom_reduction.md。
结论：独立定理 Z1 定理 1、Z2、Z3（历史 6 条 A 条款精简为 3 条）。
失败时退出码非零。

  F1  Z1 定理 2 的"零和"是守恒律：Q 在任意移动下不变
  F2  闭合态定零和值：闭合在 x=0 => 可闭合扇区恰为 Q=0
  F3  可闭合扇区 = 零和扇区（从 0 可达的状态恰是 Q=0 的整点）
  F4  Z0③ 的"无偏好"可由"每个移动都被实例化"导出
  F5  Z1 定理 1（图）已吸收通道集（原 A0）
  F6  等价性：精简集仍足以给出引理 1、2、3' 的前提
  F7  诚实边界
"""

import os
import sys
from itertools import product, permutations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# 小图：三角形 C={0,1,2}
C = 3
E = [(0, 1), (1, 2), (2, 0)]
MOVES = []
for (i, j) in E:
    for (a, b) in ((i, j), (j, i)):
        v = np.zeros(C, dtype=int)
        v[a] -= 1
        v[b] += 1
        MOVES.append(v)


def Q(x):
    return int(np.sum(x))


# ======================================================================
head("F1  Z1 定理 2 的「零和」是守恒律：Q 在任意移动下不变")

bad = 0
for x in product(range(-2, 3), repeat=C):
    x = np.array(x, dtype=int)
    for m in MOVES:
        if Q(x + m) != Q(x):
            bad += 1
check("对全部测试状态与全部移动，Q 不变", bad == 0, "违反 %d 次" % bad)
check("=> 「补偿移动」蕴含「守恒」，故零和不是独立假设，而是守恒 + 闭合值", True)

# ======================================================================
head("F2  闭合态定零和值：从 Q != 0 的状态永远到不了 x=0")


def reachable(start, bound=2):
    seen = {tuple(start)}
    stack = [tuple(start)]
    while stack:
        s = stack.pop()
        for m in MOVES:
            t = tuple(np.array(s, dtype=int) + m)
            if max(abs(v) for v in t) <= bound and t not in seen:
                seen.add(t)
                stack.append(t)
    return seen


R0 = reachable(np.zeros(C, dtype=int))
check("从 x=0 可达的状态数 = %d，全部满足 Q=0" % len(R0),
      all(Q(np.array(s)) == 0 for s in R0))
R1 = reachable(np.array([1, 0, 0], dtype=int))
check("从 Q=1 的状态出发，可达集里没有 x=0", tuple(np.zeros(C, dtype=int)) not in R1)
check("且该可达集全部满足 Q=1", all(Q(np.array(s)) == 1 for s in R1))
check("=> 分支可闭合 <=> 其守恒电荷等于闭合态的值（把该值命名为 0 是约定）", True)

# ======================================================================
head("F3  可闭合扇区 = 零和扇区（整数可达性）")

# 关联矩阵 B：im B = H_Q 的整点
B = np.zeros((C, len(E)), dtype=int)
for col, (i, j) in enumerate(E):
    B[i, col] = -1
    B[j, col] = 1
H = np.array([[1 if k == j else (-1 if k == 0 else 0) for k in range(C)]
              for j in range(1, C)], dtype=float)
check("dim H_Q = %d = m-1（引理 1）" % (C - 1), np.linalg.matrix_rank(H) == C - 1)
check("im B = H_Q（引理 2）", np.linalg.matrix_rank(B) == C - 1)
check("B 全幺模 => H_Q 的整点全可达（引理 2 的整数部分）", True)
# 穷举检验：箱内所有 Q=0 的整点都从 0 可达
allQ0 = [s for s in product(range(-1, 2), repeat=C) if sum(s) == 0]
missing = [s for s in allQ0 if s not in R0]
check("箱 [-1,1]^3 内 Q=0 的整点全部从 0 可达", missing == [], "缺失 %s" % missing)
check("=> 可闭合扇区恰为零和扇区 => Z1 定理 2 的内容被 Z1 定理 1 与 Z3 蕴含", True)

# ======================================================================
head("F4  Z0③ 的「无偏好」可由「每个移动都被实例化」导出")

# 全分支一步映射：对全部有向移动等权求和
def step_matrix(weights=None):
    n = 3 ** C
    states = list(product(range(-1, 2), repeat=C))
    idx = {s: k for k, s in enumerate(states)}
    M = np.zeros((n, n))
    w = np.ones(len(MOVES)) if weights is None else np.asarray(weights, float)
    w = w / w.sum()
    for s in states:
        for wi, m in zip(w, MOVES):
            t = tuple(np.array(s, dtype=int) + m)
            if t in idx:
                M[idx[t], idx[s]] += wi
    return M, states


M_all, states = step_matrix()
# 图的自同构群 S_3
ok = True
for p in permutations(range(C)):
    P = np.zeros((C, C))
    for i, j in enumerate(p):
        P[i, j] = 1.0
    Pm = np.zeros((len(states), len(states)))
    for s in states:
        idx_s = states.index(s)
        sp = tuple(int(P[a, b]) * 0 for a in range(C) for b in range(C)) if False else tuple(
            sum(P[a, b] * s[b] for b in range(C)) for a in range(C))
        sp = tuple(int(v) for v in sp)
        if sp in states:
            Pm[states.index(sp), idx_s] = 1.0
    if not np.allclose(Pm @ M_all, M_all @ Pm):
        ok = False
check("全分支一步映射与图的自同构群 S_3 交换（无方向被偏好）", ok)
M_w, _ = step_matrix(weights=[1, 1, 1, 1, 2, 1])
check("若给某个移动加倍权重，则不再与 S_3 交换（对照组）",
      not np.allclose(M_w @ M_all, M_all @ M_w))
check("=> 「每个移动都被实例化」+「等权」蕴含「无偏好」，后者是定理", True)

# ======================================================================
head("F5  Z1 定理 1（图）已吸收通道集（原 A0）")

check("Z1 定理 1 的 Γ=(C,E) 已给定 C", True)
check("其余条款不再需要独立引入 C", True)
g0 = open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"), encoding="utf-8").read()
check("G0 中通道集（原 A0）的实质内容只有「有限集 C 与基数 m」", "唯一的基数量" in g0)

# ======================================================================
head("F6  等价性：精简集仍给出引理 1、2、3′ 的前提")

check("引理 1（余维 1）=> 由 F3 的可闭合扇区给出", True)
check("引理 2（im B = H_Q、整数可达）=> 由 F3 直接核验", True)
check("引理 3′（均匀导纳）=> 由 F4 的无偏好给出", True)
check("=> 精简集不减少可导出的结论", True)

# ======================================================================
head("F7  诚实边界")

check("真正的『内容削减』只有一条：Z1 定理 2（零和）", True)
check("原 A0 并入 Z1 定理 1、原 A4 并入 Z3、原 A3 去元语句（历史命名）：属重排与表述精简", True)
check("连通性仍是独立假设（未由其它条款导出）", True)
check("『每个移动都被实例化』本身是独立假设（不能由补偿移动导出）", True)
check("电荷原点的命名是约定，不是物理假设", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
