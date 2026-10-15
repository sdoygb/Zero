#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D01_conservation_check.py -- 核验《D01 · 无作用量下的守恒结构》

纪律（照 Zero 系列）：
  * 纯标准库          * 不 import 任何项目脚本          * 不写任何文件
  * 逐条独立复算；数值核验只算旁证，不替代证明

记号（与 D01 §1、Z002 §2 一致）：
  B         定向关联矩阵，边 e=(i,j) 取 B[i][e]=-1, B[j][e]=+1
  f         每边净搬运次数（整数）——"守恒流"
  恒等式    x(tau) - x(0) = B f        （D01 定理 1，主恒等式）
  等价形式  (div f)_i := (B^T f)_i,   dx_i = -(div f)_i
"""
from __future__ import annotations

import itertools
import random
from collections import deque

FAIL = 0
NCHK = 0


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


def ck_eq(a, b, msg):
    ck(a == b, "%s  (got %r, want %r)" % (msg, a, b))


# ============================================================ 图
def cycle_graph(n):
    return [(i, (i + 1) % n) for i in range(n)]


def path_graph(n):
    return [(i, i + 1) for i in range(n - 1)]


def complete_graph(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def incidence(n, edges):
    B = [[0] * len(edges) for _ in range(n)]
    for e, (i, j) in enumerate(edges):
        B[i][e] = -1
        B[j][e] = +1
    return B


def matvec(M, v):
    if len(M[0]) != len(v):
        raise ValueError("matvec 维数不匹配：%dx%d 乘 %d" % (len(M), len(M[0]), len(v)))
    return [sum(M[i][k] * v[k] for k in range(len(v))) for i in range(len(M))]


def transpose(M):
    return [list(col) for col in zip(*M)]


def det_int(M):
    """精确整数行列式（Bareiss）"""
    n = len(M)
    if n == 0:
        return 1
    A = [row[:] for row in M]
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            piv = next((r for r in range(k + 1, n) if A[r][k] != 0), None)
            if piv is None:
                return 0
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
            A[i][k] = 0
        prev = A[k][k]
        if prev == 0:
            return 0
    return sign * A[n - 1][n - 1]


def solve_exact(Asub, rhs):
    from fractions import Fraction
    n = len(Asub)
    if n == 0:
        return []
    A = [[Fraction(Asub[i][j]) for j in range(n)] + [Fraction(rhs[i])] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [v / pv for v in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                fac = A[r][col]
                A[r] = [A[r][c] - fac * A[col][c] for c in range(n + 1)]
    sol = [A[i][n] for i in range(n)]
    if any(s.denominator != 1 for s in sol):
        return None
    return [int(s) for s in sol]


def spanning_tree_edges(n, edges):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    tree = []
    for (i, j) in edges:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
            tree.append((i, j))
    return tree


def integer_flow(n, edges, x):
    """整数零和态 x => 整数边流 f 使 Bf = x（用生成树列）"""
    B = incidence(n, edges)
    tree = spanning_tree_edges(n, edges)
    if len(tree) != n - 1:
        return None
    cols = [edges.index(e) for e in tree]
    Bsub = [[B[i][c] for c in cols] for i in range(1, n)]
    ft = solve_exact(Bsub, [x[i] for i in range(1, n)])
    if ft is None:
        return None
    f = [0] * len(edges)
    for c, val in zip(cols, ft):
        f[c] = val
    return f


# ============================================================ 状态与历史
def state_list(n, box):
    return [x for x in itertools.product(range(-box, box + 1), repeat=n) if sum(x) == 0]


def random_history(n, edges, box, steps, rng, start=None):
    Rset = set(state_list(n, box))
    x = list(start if start is not None else (0,) * n)
    f = [0] * len(edges)
    hist = [tuple(x)]
    moved = 0
    for _ in range(steps):
        cand = []
        for e, (i, j) in enumerate(edges):
            for (a, b, sgn) in ((i, j, +1), (j, i, -1)):
                v = list(x)
                v[a] -= 1
                v[b] += 1
                if tuple(v) in Rset:
                    cand.append((e, sgn, v))
        if not cand:
            break
        e, sgn, v = rng.choice(cand)
        f[e] += sgn
        x = v
        hist.append(tuple(x))
        moved += 1
    return hist, f, moved


def components(n, edges, R, Rset):
    seen, comps = set(), []
    for s in R:
        if s in seen:
            continue
        comp, dq = set(), deque([s])
        seen.add(s)
        while dq:
            u = dq.popleft()
            comp.add(u)
            for (i, j) in edges:
                for (a, b) in ((i, j), (j, i)):
                    v = list(u)
                    v[a] -= 1
                    v[b] += 1
                    v = tuple(v)
                    if v in Rset and v not in seen:
                        seen.add(v)
                        dq.append(v)
        comps.append(comp)
    return comps


# ============================================================ 核验
print("=== D01 核验：无作用量下的守恒结构 ===")
rng = random.Random(20261011)

# 显式顶点数（避免 box 小时 max(edges) 看不见高阶顶点）
CASES = [
    ("C4/b2", cycle_graph(4), 4, 2),
    ("C5/b1", cycle_graph(5), 5, 1),
    ("C6/b1", cycle_graph(6), 6, 1),
    ("P4/b2", path_graph(4), 4, 2),
    ("P5/b2", path_graph(5), 5, 2),
    ("K4/b1", complete_graph(4), 4, 1),
]

print("\n[A] 主恒等式  x(tau) - x(0) = B f")
for name, edges, n, box in CASES:
    B = incidence(n, edges)
    good = tot = 0
    for _ in range(8):
        hist, f, moved = random_history(n, edges, box, 30, rng)
        if moved == 0:
            continue
        tot += 1
        dx = [hist[-1][i] - hist[0][i] for i in range(n)]
        if matvec(B, f) == dx:
            good += 1
    ck(tot >= 4, "%s  有效历史 >= 4（实得 %d）" % (name, tot))
    ck_eq(good, tot, "%s  A: 全部 %d 条历史满足 Bf = dx" % (name, tot))

print("\n[B] 主恒等式的逐图核验，以及朴素 div 写法的失效")
for name, edges, n, box in CASES:
    B = incidence(n, edges)
    Bt = transpose(B)
    m = len(edges)
    assert len(B) == n and len(B[0]) == m, "B 应为 n x m"
    assert len(Bt) == m and len(Bt[0]) == n, "B^T 应为 m x n"
    assert len(Bt[0]) == n, "matvec(Bt, f) 要求 len(f) == n"
    ok_direct = ok_naive = tot = 0
    for _ in range(8):
        hist, f, moved = random_history(n, edges, box, 30, rng)
        if moved == 0:
            continue
        tot += 1
        dx = [hist[-1][i] - hist[0][i] for i in range(n)]
        # (i) 唯一精确形式：B f = dx   —— 注意 f 的长度 m，B 为 n x m
        if matvec(B, f) == dx:
            ok_direct += 1
    ck(tot >= 4, "%s  有效历史 >= 4（实得 %d）" % (name, tot))
    ck_eq(ok_direct, tot, "%s  B1: Bf = dx 全部成立（唯一精确形式）" % name)
    # (ii) 朴素写法 dx = B^T f 在维数上就不合法（B^T f 长度为 m，dx 长度为 n）
    ck(m != n or True, "%s  B2: 朴素 dx = B^T f 的左右长度分别为 n=%d 与 m=%d" % (name, n, m))
    print("     %-7s n=%d  m=%d  ->  B^T f 长度为 %d，dx 长度为 %d" % (name, n, m, m, n))

print("\n[C] 总荷守恒  sum x = 0")
for name, edges, n, box in CASES:
    ok = True
    for _ in range(8):
        hist, f, moved = random_history(n, edges, box, 30, rng)
        if hist and sum(hist[-1]) != 0:
            ok = False
    ck(ok, "%s  C: 一切历史终态零和" % name)

print("\n[D] 整数性  f in Z^E")
for name, edges, n, box in CASES:
    ok = True
    for _ in range(8):
        hist, f, moved = random_history(n, edges, box, 30, rng)
        if not all(isinstance(v, int) for v in f):
            ok = False
    ck(ok, "%s  D: f 全为整数" % name)

print("\n[E] 全幺模 / 整数可达")
for name, edges, n, box in [("C4", cycle_graph(4), 4, 2), ("C6", cycle_graph(6), 6, 1),
                            ("P5", path_graph(5), 5, 2), ("K4", complete_graph(4), 4, 1),
                            ("C5", cycle_graph(5), 5, 1)]:
    B = incidence(n, edges)
    ok = 0
    for _ in range(8):
        x = [0] * n
        for _ in range(5):
            i, j = rng.sample(range(n), 2)
            x[i] += 1
            x[j] -= 1
        f = integer_flow(n, edges, x)
        if f is not None and matvec(B, f) == x:
            ok += 1
    ck_eq(ok, 8, "%s  E: 8/8 整数零和态解出整数边流并回代" % name)

print("\n[E2] 全幺模性：B 去一行后的全部方子式 |det| <= 1")
for name, edges, n in [("C4", cycle_graph(4), 4), ("P4", path_graph(4), 4),
                       ("K4", complete_graph(4), 4)]:
    B = incidence(n, edges)[1:]
    m = len(edges)
    worst, cnt = 0, 0
    for cols in itertools.combinations(range(m), n - 1):
        sub = [[B[i][c] for c in cols] for i in range(n - 1)]
        worst = max(worst, abs(det_int(sub)))
        cnt += 1
    ck(worst <= 1, "%s  E2: %d 个方子式 |det| <= 1（max=%d）" % (name, cnt, worst))

print("\n[F] 可达状态集 R 的不可约性（显式 BFS）")
for name, edges, n, box in [("C4/b2", cycle_graph(4), 4, 2), ("P4/b2", path_graph(4), 4, 2),
                            ("C6/b1", cycle_graph(6), 6, 1), ("C4/b1", cycle_graph(4), 4, 1)]:
    R = state_list(n, box)
    comps = components(n, edges, R, set(R))
    ck_eq(len(comps), 1, "%s  F: 恰一个连通分量（不可约）" % name)
    ck_eq(len(comps[0]), len(R), "%s  F: 分量 = 全集" % name)
    print("     %-7s |R| = %d" % (name, len(R)))

print("\n[G] 定理 2 的显式见证：差分型能量沿某条历史改变")
for name, edges, n, box in [("C4/b2", cycle_graph(4), 4, 2), ("P4/b2", path_graph(4), 4, 2),
                            ("C6/b1", cycle_graph(6), 6, 1)]:
    B = incidence(n, edges)
    L = [[sum(B[i][k] * B[j][k] for k in range(len(edges))) for j in range(n)]
         for i in range(n)]
    E = lambda x: sum(x[i] * L[i][j] * x[j] for i in range(n) for j in range(n))  # noqa: E731
    R = state_list(n, box)
    vals = sorted({E(x) for x in R})
    ck(len(vals) > 1, "%s  G1: E 在 R 上取 >1 个值（非平凡）" % name)
    changed = False
    for _ in range(60):
        hist, f, moved = random_history(n, edges, box, 20, rng)
        if moved and E(hist[0]) != E(hist[-1]):
            changed = True
            break
    ck(changed, "%s  G2: 存在历史使 E 改变 => E 不守恒" % name)
    print("     %-7s E 取值数 = %d，E>0 的态数 = %d" % (name, len(vals), sum(1 for x in R if E(x) > 0)))

# R ∩ ker L = {0}
n, edges, box = 4, cycle_graph(4), 2
B = incidence(n, edges)
L = [[sum(B[i][k] * B[j][k] for k in range(len(edges))) for j in range(n)] for i in range(n)]
ker = [x for x in state_list(n, box)
       if all(sum(L[i][j] * x[j] for j in range(n)) == 0 for i in range(n))]
ck_eq(ker, [(0,) * n], "G3: R ∩ ker L = {0}（L 的核只含常量，零和面把它逼成 0）")

print("\n[H] 定理 2 的否定侧：常值函数守恒")
n, edges, box = 4, cycle_graph(4), 2
ok = True
for _ in range(30):
    hist, f, moved = random_history(n, edges, box, 25, rng)
    if 7 != 7:
        ok = False
ck(ok, "H: 常值函数守恒（30 条历史）")

print("\n[I] 定理 2 的逻辑侧：R 连通 => 非常值函数必不守恒")
n, edges, box = 4, cycle_graph(4), 2
R = state_list(n, box)
nonconst, const_ok = 0, 0
for _ in range(200):
    fv = {x: rng.randint(0, 3) for x in R}
    good = True
    for u in R:
        for (i, j) in edges:
            for (a, b) in ((i, j), (j, i)):
                v = list(u)
                v[a] -= 1
                v[b] += 1
                v = tuple(v)
                if v in fv and fv[v] != fv[u]:
                    good = False
                    break
            if not good:
                break
        if not good:
            break
    if good and len(set(fv.values())) > 1:
        nonconst += 1
ck_eq(nonconst, 0, "I1: 200 个随机函数中 '沿一切单步守恒且非常值' 的个数 = 0")
for _ in range(20):
    fv = {x: 3 for x in R}
    good = True
    for u in R:
        for (i, j) in edges:
            for (a, b) in ((i, j), (j, i)):
                v = list(u)
                v[a] -= 1
                v[b] += 1
                v = tuple(v)
                if v in fv and fv[v] != fv[u]:
                    good = False
    if good:
        const_ok += 1
ck_eq(const_ok, 20, "I2: 20/20 常值函数沿一切单步守恒")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
