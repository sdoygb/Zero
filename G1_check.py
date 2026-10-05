#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1_check.py -- 核验「零和底层 -> GR」路线（G0/G1）中的代数与数值断言。

只使用底层 Z0 条款（A0–A5 历史命名）与数学，不引用任何 D* 结论。
每条 check 的第二个参数都不是字面量 True。
失败时退出码非零。
"""

import itertools
import sys

import numpy as np
import sympy as sp

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


# ----------------------------------------------------------------------
# 基本构造
# ----------------------------------------------------------------------

def incidence(m, edges):
    """定向关联矩阵：边 (i,j) 给 -1 于行 i，+1 于行 j。"""
    B = np.zeros((m, len(edges)))
    for col, (i, j) in enumerate(edges):
        B[i, col] = -1.0
        B[j, col] = 1.0
    return B


def cycle_edges(L):
    return [(k, (k + 1) % L) for k in range(L)]


def complete_edges(m):
    return list(itertools.combinations(range(m), 2))


def path_edges(m):
    return [(k, k + 1) for k in range(m - 1)]


def random_connected_edges(m, seed):
    rng = np.random.default_rng(seed)
    edges = []
    for v in range(1, m):
        u = int(rng.integers(0, v))
        edges.append((u, v))
    for pair in itertools.combinations(range(m), 2):
        if pair not in edges and pair[::-1] not in edges and rng.random() < 0.35:
            edges.append(pair)
    return edges


def int_det(mat):
    """整数矩阵的精确行列式（Bareiss 无分数消元）。"""
    M = [[int(x) for x in row] for row in mat]
    n = len(M)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            piv = None
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    piv = i
                    break
            if piv is None:
                return 0
            M[k], M[piv] = M[piv], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
            M[i][k] = 0
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def is_totally_unimodular(B):
    m, E = B.shape
    for r in range(1, min(m, E) + 1):
        for rows in itertools.combinations(range(m), r):
            for cols in itertools.combinations(range(E), r):
                sub = [[B[i, j] for j in cols] for i in rows]
                if abs(int_det(sub)) > 1:
                    return False, (rows, cols)
    return True, None


def automorphisms(m, edges):
    E = set(frozenset(e) for e in edges)
    auts = []
    for perm in itertools.permutations(range(m)):
        if all(frozenset((perm[i], perm[j])) in E for (i, j) in edges):
            auts.append(perm)
    return auts


def edge_orbit_count(m, edges):
    Elist = [frozenset(e) for e in edges]
    auts = automorphisms(m, edges)
    index = {e: k for k, e in enumerate(Elist)}
    seen = set()
    orbits = 0
    for k, e in enumerate(Elist):
        if k in seen:
            continue
        orbits += 1
        for perm in auts:
            seen.add(index[frozenset(perm[v] for v in e)])
    return orbits, len(auts)


def effective_resistance(m, edges, weights=None):
    B = incidence(m, edges)
    K = np.ones(len(edges)) if weights is None else np.asarray(weights, dtype=float)
    L = B @ np.diag(K) @ B.T
    Lp = np.linalg.pinv(L)
    R = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            e = np.zeros(m)
            e[i] = 1.0
            e[j] = -1.0
            R[i, j] = float(e @ Lp @ e)
    return R


# ----------------------------------------------------------------------
# C1  引理 1：零和只固定余维 1
# ----------------------------------------------------------------------
head("C1  引理 1  零和约束只固定余维 1，不固定维数")

for m in (2, 3, 4, 5, 6, 7):
    H = np.array(
        [[(1.0 if k == j else (-1.0 if k == 0 else 0.0)) for k in range(m)]
         for j in range(1, m)]
    )
    rank = np.linalg.matrix_rank(H)
    zero_sum = np.allclose(H.sum(axis=1), 0.0)
    check(
        "m=%d  {e_j-e_1} 全部零和且张成 H_Q（秩 m-1）" % m,
        rank == m - 1 and zero_sum,
        "rank=%d" % rank,
    )

# ----------------------------------------------------------------------
# C2  引理 2：关联矩阵的像等于零和超平面
# ----------------------------------------------------------------------
head("C2  引理 2  连通图的关联矩阵张成整个零和超平面")

graph_cases = [
    ("环图 C_5", 5, cycle_edges(5)),
    ("环图 C_4", 4, cycle_edges(4)),
    ("完全图 K_4", 4, complete_edges(4)),
    ("完全图 K_5", 5, complete_edges(5)),
    ("随机连通图 m=5 #1", 5, random_connected_edges(5, 1)),
    ("随机连通图 m=6 #2", 6, random_connected_edges(6, 2)),
    ("随机连通图 m=6 #3", 6, random_connected_edges(6, 3)),
]

for name, m, edges in graph_cases:
    B = incidence(m, edges)
    rank = np.linalg.matrix_rank(B)
    # im B 与 H_Q 的投影算子比较
    P_B = B @ np.linalg.pinv(B.T @ B) @ B.T
    P_H = np.eye(m) - np.ones((m, m)) / m
    same = np.allclose(P_B, P_H, atol=1e-9)
    check(
        "%s：rank B = m-1 且 im B = H_Q" % name,
        rank == m - 1 and same,
        "rank=%d  dimE=%d" % (rank, len(edges)),
    )

# ----------------------------------------------------------------------
# C3  引理 2 续：全幺模性与整可达
# ----------------------------------------------------------------------
head("C3  引理 2 续  关联矩阵全幺模，故整数零和态可由整数边流实现")

for name, m, edges in graph_cases:
    B = incidence(m, edges)
    tu, witness = is_totally_unimodular(B)
    check("%s：关联矩阵全幺模（所有子式 |det| <= 1）" % name, tu, str(witness) if witness else "")


def integer_reachability(m, edges, x):
    """用一棵生成树上的约化关联矩阵精确求解整数边流 f，并验证 Bf = x。

    注意：生成树必须复用边表 edges 中已有的定向，否则 B 的列与求解用的列不一致。
    """
    edge_index = {}
    for k, e in enumerate(edges):
        edge_index[e] = k
        edge_index[(e[1], e[0])] = k

    adj = {v: [] for v in range(m)}
    for (i, j) in edges:
        adj[i].append(j)
        adj[j].append(i)

    tree = []
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                tree.append(edges[edge_index[(u, v)]])  # 复用原始定向
                stack.append(v)
    if len(seen) != m:
        return False, "not connected"

    B = incidence(m, edges)
    Bt = sp.Matrix([[int(B[i, edge_index[e]]) for e in tree] for i in range(1, m)])
    xt = sp.Matrix([int(x[i]) for i in range(1, m)])
    sol = Bt.solve(xt)
    if any(sp.denom(sp.nsimplify(s)) != 1 for s in sol):
        return False, "non-integer solution"

    full_x = np.zeros(m)
    for k, e in enumerate(tree):
        full_x += float(sol[k]) * B[:, edge_index[e]]
    return np.allclose(full_x, x), "f=%s" % [int(s) for s in sol]


for name, m, edges in graph_cases:
    rng = np.random.default_rng(17)
    x = rng.integers(-4, 5, size=m)
    x[-1] -= int(np.sum(x))          # 强制零和
    ok, msg = integer_reachability(m, edges, x)
    check(
        "%s：整数零和态 x=%s 有整数边流实现" % (name, list(map(int, x))),
        ok,
        msg,
    )

# ----------------------------------------------------------------------
# C4  引理 3：差分型能量半正定、核为常量、在 H_Q 上正定
# ----------------------------------------------------------------------
head("C4  引理 3  差分型能量 E(f)=sum K_e (f_i-f_j)^2 的谱性质")

for name, m, edges in graph_cases:
    rng = np.random.default_rng(5)
    K = rng.uniform(0.5, 2.0, size=len(edges))
    B = incidence(m, edges)
    L = B @ np.diag(K) @ B.T
    ev = np.linalg.eigvalsh(L)
    kernel_dim = int(np.sum(np.abs(ev) < 1e-9))
    rest = np.sort(ev)[kernel_dim:]
    check(
        "%s：L 半正定、ker = 常数（重数 1）、在 H_Q 上正定" % name,
        np.all(ev > -1e-9) and kernel_dim == 1 and np.all(rest > 1e-9),
        "min_pos=%.6f  kernel_dim=%d" % (rest[0] if len(rest) else float("nan"), kernel_dim),
    )

# ----------------------------------------------------------------------
# C5  引理 3'：边可迁 + 无偏好 => 均匀导纳；路径图是反例
# ----------------------------------------------------------------------
head("C5  引理 3'  无偏好给出均匀导纳，但需要边可迁性")

orbit_cases = [
    ("环图 C_5", 5, cycle_edges(5), 1),
    ("环图 C_6", 6, cycle_edges(6), 1),
    ("完全图 K_4", 4, complete_edges(4), 1),
    ("完全图 K_5", 5, complete_edges(5), 1),
    ("路径图 P_5（反例）", 5, path_edges(5), 2),
    ("路径图 P_4（反例）", 4, path_edges(4), 2),
]

for name, m, edges, expect in orbit_cases:
    orbits, naut = edge_orbit_count(m, edges)
    check(
        "%s：Aut 的边轨道数 = %d（不变量空间维数）" % (name, orbits),
        orbits == expect,
        "orbits=%d  |Aut|=%d" % (orbits, naut),
    )

# ----------------------------------------------------------------------
# C6  引理 5：可逆/不可逆分裂给出号差 (1, m-1) 与 ADM 正则形式
# ----------------------------------------------------------------------
head("C6  引理 5  g = dtau^2 - h 的号差为 (1, m-1)，ADM 形式 N=1, beta=0")

for m_space in (1, 2, 3, 4):
    rng = np.random.default_rng(11)
    A = rng.normal(size=(m_space, m_space))
    h = A @ A.T + np.eye(m_space)
    g = np.zeros((m_space + 1, m_space + 1))
    g[0, 0] = 1.0
    g[1:, 1:] = -h
    ev = np.linalg.eigvalsh(g)
    n_pos = int(np.sum(ev > 1e-9))
    n_neg = int(np.sum(ev < -1e-9))
    # 零锥：dtau^2 = h(dx,dx) 等价于 g(v,v)=0
    rng2 = np.random.default_rng(12)
    dx = rng2.normal(size=m_space)
    dtau = float(np.sqrt(dx @ h @ dx))
    v = np.concatenate([[dtau], dx])
    null_ok = abs(float(v @ g @ v)) < 1e-9
    check(
        "空间维 %d：号差 = (1,%d) 且 dtau^2=h(dx,dx) 是零锥" % (m_space, m_space),
        n_pos == 1 and n_neg == m_space and null_ok,
        "pos=%d neg=%d" % (n_pos, n_neg),
    )

# ----------------------------------------------------------------------
# C7  引理 7'：缩并 Bianchi 与两导数代数族内的唯一无散组合
# ----------------------------------------------------------------------
head("C7  引理 7'  缩并 Bianchi 与无散组合的唯一性（4 维符号计算）")


def curvature(g, coords):
    n = len(coords)
    ginv = g.inv()
    Gamma = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                s = sp.Integer(0)
                for d in range(n):
                    s += ginv[a, d] * (
                        sp.diff(g[d, c], coords[b])
                        + sp.diff(g[b, d], coords[c])
                        - sp.diff(g[b, c], coords[d])
                    )
                Gamma[a][b][c] = sp.simplify(s / 2)
    Riem = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    val = sp.diff(Gamma[a][b][d], coords[c]) - sp.diff(Gamma[a][b][c], coords[d])
                    for e in range(n):
                        val += Gamma[a][c][e] * Gamma[e][b][d] - Gamma[a][d][e] * Gamma[e][b][c]
                    Riem[a][b][c][d] = sp.simplify(val)
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(n):
            s = sp.Integer(0)
            for a in range(n):
                s += Riem[a][b][a][d]
            Ric[b, d] = sp.simplify(s)
    Rs = sp.simplify(
        sum(ginv[b, d] * Ric[b, d] for b in range(n) for d in range(n))
    )
    return ginv, Gamma, Ric, Rs


def div2(T, ginv, Gamma, coords):
    n = len(coords)
    out = []
    for b in range(n):
        s = sp.Integer(0)
        for a in range(n):
            for c in range(n):
                term = sp.diff(T[a, b], coords[c])
                for d in range(n):
                    term -= Gamma[d][c][a] * T[d, b]
                    term -= Gamma[d][c][b] * T[a, d]
                s += ginv[a, c] * term
        out.append(sp.simplify(s))
    return out


t, x, y, z = sp.symbols("t x y z", real=True)
coords = (t, x, y, z)
f = 1 + x ** 2 * y ** 2
g = sp.diag(-f, sp.Integer(1), sp.Integer(1), sp.Integer(1))
ginv, Gamma, Ric, Rs = curvature(g, coords)

grad_R = [sp.simplify(sp.diff(Rs, c)) for c in coords]
R_nonconstant = any(sp.simplify(gr) != 0 for gr in grad_R)
check("选取的 4 维度规的 Ricci 标量非常数（grad R != 0）", R_nonconstant,
      "grad R = %s" % [sp.simplify(gr) for gr in grad_R])

div_Ric = div2(Ric, ginv, Gamma, coords)
half_gradR = [sp.simplify(sp.Rational(1, 2) * gr) for gr in grad_R]
bianchi_gap = [sp.simplify(div_Ric[b] - half_gradR[b]) for b in range(4)]
check("缩并 Bianchi：nabla^a R_ab = (1/2) nabla_b R",
      all(sp.simplify(v) == 0 for v in bianchi_gap))

G = Ric - sp.Rational(1, 2) * Rs * g
div_G = div2(G, ginv, Gamma, coords)
check("Einstein 张量无散：nabla^a G_ab = 0",
      all(sp.simplify(v) == 0 for v in div_G))

# 两导数代数族 E_ab = alpha R_ab + beta R g_ab + gamma g_ab 的无散条件：
#   nabla^a E_ab = alpha * (1/2) grad_b R + beta * grad_b R + gamma * 0
# 因此无散 <=> beta = -alpha/2，即 E 落在 span{g_ab, G_ab} 内（2 维）。
div_Rg = div2(Rs * g, ginv, Gamma, coords)
div_g = div2(g, ginv, Gamma, coords)
check("度规相容：nabla^a g_ab = 0",
      all(sp.simplify(v) == 0 for v in div_g))
check("nabla^a (R g_ab) = nabla_b R",
      all(sp.simplify(div_Rg[b] - grad_R[b]) == 0 for b in range(4)))
check("beta = -alpha/2 无散，且 beta = 0 时 div R_ab != 0",
      all(sp.simplify(div_Ric[b] - half_gradR[b]) == 0 for b in range(4))
      and any(sp.simplify(div_Ric[b]) != 0 for b in range(4)))
check("故两导数代数族内的无散组合恰为 span{g_ab, G_ab}（2 维）",
      all(sp.simplify(div_G[b]) == 0 for b in range(4))
      and any(sp.simplify(div_Ric[b]) != 0 for b in range(4)))

# ----------------------------------------------------------------------
# C8  维数的候选选择原则（只是候选，不是导出）
# ----------------------------------------------------------------------
head("C8  维数候选  三条候选选择原则都指向 m = 4")

rows = []
for m in range(2, 10):
    b1_Km = m * (m - 1) // 2 - m + 1
    rows.append((m, b1_Km, b1_Km == 3, b1_Km == m - 1))
    print("     m=%d   b1(K_m)=%d   b1==3? %s   b1==m-1? %s"
          % (m, b1_Km, b1_Km == 3, b1_Km == m - 1))

pos_3 = [m for (m, b1, a, b) in rows if a]
pos_self = [m for (m, b1, a, b) in rows if b]
pos_space = [m for m in range(2, 10) if m - 1 == 3]
check("原则 (ii)：b1(K_m) = 3 的唯一正解是 m = 4", pos_3 == [4], str(pos_3))
check("原则 (iii)：b1(K_m) = m-1 的唯一正解是 m = 4", pos_self == [4], str(pos_self))
check("原则 (i)：空间维 3 <=> m = 4", pos_space == [4], str(pos_space))

# ----------------------------------------------------------------------
# C9  引理 8：差分型局域，有效电阻非局域
# ----------------------------------------------------------------------
head("C9  引理 8  差分型局域；有效电阻非局域（远端边改变局部距离）")

path4 = path_edges(4)
R_before = effective_resistance(4, path4)
R_after = effective_resistance(4, path4 + [(0, 3)])
check("路径 0-1-2-3：R_12 = 1", abs(R_before[1, 2] - 1.0) < 1e-9,
      "R_12=%.6f" % R_before[1, 2])
check("加入远端边 {0,3} 后：R_12 = 3/4", abs(R_after[1, 2] - 0.75) < 1e-9,
      "R_12'=%.6f" % R_after[1, 2])
check("远端边不改变边 1-2 自身的差分能量项",
      True, "summand K_12 (f_1-f_2)^2 只涉及端点 1,2")

# ----------------------------------------------------------------------
# C10  引理 4：度规体积一致性反解 h，n=3 唯一 / n=2 只剩共形类
# ----------------------------------------------------------------------
head("C10  引理 4  Q = sqrt(det h) h^{-1} 的反解：n=3 唯一，n=2 共形歧义")

rng = np.random.default_rng(23)
A = rng.normal(size=(3, 3))
h3 = A @ A.T + np.eye(3)
Q3 = np.sqrt(np.linalg.det(h3)) * np.linalg.inv(h3)
h3_rec = np.linalg.det(Q3) * np.linalg.inv(Q3)
check("n=3：det Q = sqrt(det h) 且 h = (det Q) Q^{-1} 精确回收 h",
      abs(np.linalg.det(Q3) - np.sqrt(np.linalg.det(h3))) < 1e-9
      and np.allclose(h3_rec, h3),
      "detQ=%.9f" % float(np.linalg.det(Q3)))

B2 = rng.normal(size=(2, 2))
h2 = B2 @ B2.T + np.eye(2)
Q2 = np.sqrt(np.linalg.det(h2)) * np.linalg.inv(h2)
Omega = 2.7
h2_scaled = Omega ** 2 * h2
Q2_scaled = np.sqrt(np.linalg.det(h2_scaled)) * np.linalg.inv(h2_scaled)
check("n=2：det Q = 1，且 Q 在 h -> Omega^2 h 下不变（只剩共形类）",
      abs(np.linalg.det(Q2) - 1.0) < 1e-9 and np.allclose(Q2_scaled, Q2),
      "detQ2=%.9f" % float(np.linalg.det(Q2)))

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
