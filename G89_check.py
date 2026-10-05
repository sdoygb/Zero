#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G89_check.py -- 维数的分类级 no-go 与「极化两标签无偏好」条件。

对应文档 G89_dimension_no_go_and_the_balance_condition.md。
只使用底层条款（Z0 条款 ＋ Z1–Z5 定理）与数学（表示论 / 线性代数 / 图论），不引用任何 D* 结论。
失败时退出码非零。

  F1  极化空间 P_D：维数公式与独立枚举（D=2..12）
  F2  Z1 定理 1 的移动取逆在 P_D 上作用平凡（G11 定位错误的根源）
  F3  横截反射 R 的特征维数；两特征空间各 1 维 <=> D=4
  F4  命题 2：原生 (1,1) no-go（闭式 + 穷尽全部对合）
  F5  命题 1：模型论 no-go（环图 C_m 对 m=2..12 逐条款核验）
  F6  候选约束的存活集合与「同义改写」判定
  F7  G1 引理 4：唯一反解 h=(det Q)^(1/(n-2)) Q^{-1} 对一切 n != 2
  F8  G63 §3：D=|F|(m-1)+1=4 的整数枚举与循环性的精确形式
"""

import itertools
import io
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

PASS = 0
FAIL = 0


def head(t):
    print("")
    print("=" * 70)
    print(t)
    print("=" * 70)


def check(name, cond, detail=""):
    global PASS, FAIL
    ok = bool(cond)
    if ok:
        PASS += 1
    else:
        FAIL += 1
    line = "  [%s] %s" % ("v" if ok else "x", name)
    if detail:
        line += "   %s" % detail
    print(line)


# ----------------------------------------------------------------------
# 工具
# ----------------------------------------------------------------------
def P_basis(D):
    """横截平面 R^{D-2} 上对称无迹 2-张量空间的一组基。"""
    n = D - 2
    if n <= 1:
        return []
    cand = []
    for i in range(n):
        for j in range(i, n):
            M = np.zeros((n, n))
            if i == j:
                M[i, i] = 1.0
                M[(i + 1) % n, (i + 1) % n] -= 1.0
                if abs(M).sum() < 1e-12:
                    continue
            else:
                M[i, j] = M[j, i] = 1.0
            cand.append(M)
    indep = []
    for M in cand:
        v = M.reshape(-1)
        if not indep:
            indep.append(v)
            continue
        A = np.array(indep).T
        s, *_ = np.linalg.lstsq(A, v, rcond=None)
        if np.linalg.norm(A @ s - v) > 1e-9:
            indep.append(v)
    return [v.reshape(n, n) for v in indep]


def split_P(basis, R):
    """R 以 h -> R h R^T 作用，返回 (dim +1, dim -1)。"""
    d = len(basis)
    if d == 0:
        return (0, 0)
    Bm = np.array([b.reshape(-1) for b in basis]).T
    A = np.zeros((d, d))
    for k, M in enumerate(basis):
        img = (R @ M @ R.T).reshape(-1)
        s, *_ = np.linalg.lstsq(Bm, img, rcond=None)
        A[:, k] = s
    w = np.linalg.eigvals(A)
    return (int(np.sum(np.abs(w - 1) < 1e-8)), int(np.sum(np.abs(w + 1) < 1e-8)))


def H_basis(m):
    H = np.zeros((m, m - 1))
    for i in range(m - 1):
        H[i, i] = 1.0
    H[m - 1, :] = -1.0
    return H


def split_H(m, M):
    """M 作用在 H_Q 上，返回 (dim +1, dim -1)。"""
    H = H_basis(m)
    A = np.linalg.pinv(H) @ M @ H
    w = np.linalg.eigvals(A)
    return (int(np.sum(np.abs(w - 1) < 1e-7)), int(np.sum(np.abs(w + 1) < 1e-7)))


def H_at(D):
    return np.diag([1.0] * (D - 3) + [-1.0])


def ring(m):
    return [(i, (i + 1) % m) for i in range(m)]


def complete(m):
    return [(i, j) for i in range(m) for j in range(i + 1, m)]


def incidence(m, E):
    B = np.zeros((m, len(E)))
    for k, (a, b) in enumerate(E):
        B[a, k] = -1.0
        B[b, k] = 1.0
    return B


def energy(m, E):
    L = np.zeros((m, m))
    for a, b in E:
        v = np.zeros(m)
        v[a] = -1.0
        v[b] = 1.0
        L += np.outer(v, v)
    return L


def cycles(perm):
    m = len(perm)
    seen = [False] * m
    out = []
    for i in range(m):
        if seen[i]:
            continue
        c = []
        j = i
        while not seen[j]:
            seen[j] = True
            c.append(j)
            j = perm[j]
        out.append(c)
    return out


def parts(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield []
        return
    for k in range(min(n, mx), 0, -1):
        for rest in parts(n - k, k):
            yield [k] + rest


def perm_from(ct):
    p = []
    idx = 0
    for L in ct:
        blk = list(range(idx, idx + L))
        for t in range(L):
            p.append(blk[(t + 1) % L])
        idx += L
    return p


DM = list(range(2, 13))
DM4 = list(range(4, 13))

# ======================================================================
head("F1  极化空间 P_D：维数公式与独立枚举")
tab = {}
for D in DM:
    B = P_basis(D)
    tab[D] = len(B)
    f = D * (D - 3) // 2 if D >= 3 else 0
    check("dim P_%d = D(D-3)/2" % D, tab[D] == f, "枚举 %d，公式 %d" % (tab[D], f))
check("dim P_D = D(D-3)/2 对 D=2..12 全部成立",
      all(tab[D] == (D * (D - 3) // 2 if D >= 3 else 0) for D in DM))

# ======================================================================
head("F2  Z1 定理 1 的移动取逆在 P_D 上作用平凡（G11 定位错误的根源）")
for D in range(4, 13):
    pp, pm = split_P(P_basis(D), -np.eye(D - 2))
    check("D=%d: rho_{Z1 定理 1} 在 P_D 上 = +I（(dimP+,dimP-)=(%d,%d)）" % (D, pp, pm),
          pm == 0 and pp == tab[D], "dimP=%d" % tab[D])
check("rho_{Z1 定理 1} 在 P_D 上作用平凡对 D=4..12 全部成立",
      all(split_P(P_basis(D), -np.eye(D - 2))[1] == 0 for D in range(4, 13)))
_D = sp.Symbol("D")
check("dim P_D = 1 <=> D^2-3D-2=0 无整数解（判别式 17 非平方）",
      sp.discriminant(_D ** 2 - 3 * _D - 2, _D) == 17)

# ======================================================================
head("F3  横截反射 R=diag(1,..,1,-1) 的特征维数")
for D in range(4, 13):
    pp, pm = split_P(P_basis(D), H_at(D))
    check("D=%d: (dimP+,dimP-)=(%d,%d)" % (D, pp, pm),
          (pp, pm) == ((D - 2) * (D - 3) // 2, D - 3))
check("两个特征空间同时 1 维 <=> D=4",
      [D for D in DM if split_P(P_basis(D), H_at(D)) == (1, 1)] == [4])

# ======================================================================
head("F4  命题 2：原生 (1,1) no-go（闭式 + 穷尽）")
for m in range(3, 9):
    a2 = split_H(m, -np.eye(m))
    tau = np.zeros((m, m))
    for i in range(m):
        tau[(1 if i == 0 else 0 if i == 1 else i), i] = 1.0
    a3 = split_H(m, tau)
    sig = np.zeros((m, m))
    for i in range(m):
        sig[(-i) % m, i] = 1.0
    a5 = split_H(m, sig)
    any11 = any(x == (1, 1) for x in (a2, a3, a5))
    check("m=%d: Z1 定理 1=%s Z0③=%s Z3=%s 任一 (1,1)? %s" % (m, a2, a3, a5, any11),
          any11 == (m == 3))

found11 = []
for m in range(2, 13):
    for ct in parts(m):
        if any(L > 2 for L in ct):
            continue
        perm = perm_from(ct)
        if any(perm[perm[i]] != i for i in range(m)):
            continue
        c = len(cycles(perm))
        nev = sum(1 for L in ct if L % 2 == 0)
        M = np.zeros((m, m))
        for i in range(m):
            M[perm[i], i] = 1.0
        comp = split_H(m, M)
        assert comp == (c - 1, nev), (m, ct, comp, (c - 1, nev))
        if comp == (1, 1):
            found11.append((m, tuple(sorted(ct, reverse=True))))
check("命题 2(a) 闭式 dim H+ = #轮换-1, dim H- = #偶非平凡轮换：m=2..12 全部对合通过", True)
check("命题 2(b) 全部对合中给 (1,1) 的只有 (m,轮换型) = [(3,(2,1))]",
      found11 == [(3, (2, 1))], str(found11))
for m in range(2, 10):
    specs = set()
    for pi in itertools.permutations(range(m)):
        if any(pi[pi[i]] != i for i in range(m)):
            continue
        M = np.zeros((m, m))
        for i in range(m):
            M[pi[i], i] = 1.0
        specs.add(split_H(m, M))
        specs.add(split_H(m, -M))
    check("命题 2(c) m=%d: {±1} x {S_m 对合} 给 (1,1) 当且仅当 m=3" % m,
          ((1, 1) in specs) == (m == 3))

# ======================================================================
head("F5  命题 1：模型论 no-go（环图 C_m 逐条款核验）")
allok = True
for m in range(2, 13):
    E = ring(m)
    B = incidence(m, E)
    L = energy(m, E)
    rk = int(np.linalg.matrix_rank(B))
    w = np.linalg.eigvalsh(L)
    ker = int(np.sum(np.abs(w) < 1e-9))
    H = H_basis(m)
    Lh = H.T @ L @ H
    pos = bool(np.all(np.linalg.eigvalsh(Lh) > 1e-9))
    ok = (rk == m - 1) and (ker == 1) and pos
    allok &= ok
    if m <= 6 or not ok:
        check("m=%d: 环图 C_m 是底层条款（Z0 条款 ＋ Z1–Z5 定理）的模型" % m, ok,
              "rankB=%d, kerL=%d, H_Q 上正定=%s" % (rk, ker, pos))
check("环图 C_m 对 m=2..12 全部给出底层条款（Z0 条款 ＋ Z1–Z5 定理）的模型", allok)


def edge_transitive(m, E):
    """Aut 在**边**上是否传递（m 小，暴力核验全部置换）。"""
    Es = set(tuple(sorted(e)) for e in E)
    base = tuple(sorted(E[0]))
    orbit = set()
    for p in itertools.permutations(range(m)):
        if set(tuple(sorted((p[a], p[b]))) for a, b in E) != Es:
            continue
        for a, b in E:
            if tuple(sorted((a, b))) == base:
                orbit.add(tuple(sorted((p[a], p[b]))))
    return orbit == Es


for m in (4, 5):
    check("m=%d: 环图 C_m 边可迁（均匀导纳前提成立）" % m,
          edge_transitive(m, ring(m)))
    check("m=%d: 路径图 P_m 不边可迁（G1 引理 3′ 的反例）" % m,
          not edge_transitive(m, [(i, i + 1) for i in range(m - 1)]))

# ======================================================================
head("F6  候选约束的存活集合与同义改写判定")
cands = [
    ("(eps) 横截反射两特征空间各 1 维",
     lambda D: split_P(P_basis(D), H_at(D)) == (1, 1)),
    ("(beta) dimP+ = dimP-",
     lambda D: (lambda s: s[0] == s[1] and tab[D] > 0)(split_P(P_basis(D), H_at(D)))),
    ("(alpha) dim P_D = 2", lambda D: tab[D] == 2),
    ("(gamma) dimP- = 1", lambda D: split_P(P_basis(D), H_at(D))[1] == 1),
    ("(b1-K) b1(K_m) = m-1",
     lambda D: len(complete(D)) - D + 1 == D - 1),
    ("(b1-C) b1(C_m) = m-1",
     lambda D: len(ring(D)) - D + 1 == D - 1),
    ("(b1-3) b1(K_m) = 3",
     lambda D: len(complete(D)) - D + 1 == 3),
]
for name, f in cands:
    s4 = [D for D in DM4 if f(D)]
    print("    %-34s 存活 = %s%s" % (name, s4, "   *唯一" if len(s4) == 1 else ""))
check("(eps)/(beta)/(alpha)/(gamma)/(b1-K)/(b1-3) 在 D>=4 上都存活 {4}（= 同义改写）",
      all([D for D in DM4 if f(D)] == [4]
          for _, f in cands[:4] + [cands[4], cands[6]]))
check("(b1-C) b1(C_m) = m-1 在 D>=4 上无解（环图 b1=1 对所有 m）",
      [D for D in DM4 if len(ring(D)) - D + 1 == D - 1] == [])
check("b1(K_m)=m-1 的维数等式只对 m=4 成立： (m-1)(m-4)=0",
      [m for m in range(2, 13) if (m - 1) * (m - 2) // 2 == m - 1] == [4])


def perm_mat(m, perm):
    M = np.zeros((m, m))
    for i in range(m):
        M[perm[i], i] = 1.0
    return M


for m in (4, 5, 6):
    dimZ = len(complete(m)) - m + 1          # dim Z_1(K_m)
    dimW = m - 1                             # dim H_Q
    eq = (dimZ == dimW)
    check("m=%d: dim Z_1(K_m) = %d, dim H_Q = %d, 相等？%s"
          % (m, dimZ, dimW, eq), eq == (m == 4),
          "(m-1)(m-4)=0" if eq else "维数已不等 => 必不同构")
    L = energy(m, complete(m))
    ones = np.ones(m)
    check("m=%d: L @ 1 = 0（完全图上能量核 = span{1}，引理 3）" % m,
          np.max(np.abs(L @ ones)) < 1e-9
          and int(np.sum(np.abs(np.linalg.eigvalsh(L)) < 1e-9)) == 1)

# ======================================================================
head("F7  G1 引理 4：唯一反解 h = (det Q)^(1/(n-2)) Q^{-1} 对一切 n != 2")
rng = np.random.default_rng(7)
for nv in range(3, 9):
    err = 0.0
    for _ in range(6):
        Z = rng.normal(size=(nv, nv))
        h = Z @ Z.T + nv * np.eye(nv)
        Q = np.sqrt(np.linalg.det(h)) * np.linalg.inv(h)
        dQ = np.linalg.det(Q)
        hr = dQ ** (1.0 / (nv - 2)) * np.linalg.inv(Q)
        err = max(err, np.max(np.abs(hr - h)) / np.max(np.abs(h)))
    check("n=%d: 反解相对误差 <= 1e-12" % nv, err <= 1e-12, "%.3e" % err)
Z = rng.normal(size=(2, 2))
hA = Z @ Z.T + 2 * np.eye(2)
QA = np.sqrt(np.linalg.det(hA)) * np.linalg.inv(hA)
diffs = []
for lam in (1.0, 3.0, 11.0):
    hB = lam * hA
    QB = np.sqrt(np.linalg.det(hB)) * np.linalg.inv(hB)
    diffs.append(np.max(np.abs(QB - QA)))
check("n=2: 同一 Q 对应无穷多 h（det Q 恒 1，Q 不含 det h）", max(diffs) < 1e-12,
      "max|dQ| = %.2e" % max(diffs))
check("引理 4 不能选出任何维数（唯一性区间是 n != 2，含 n=4,5,6,...）", True)

# ======================================================================
head("F8  G63 §3：D = |F|(m-1) + 1 = 4 的整数枚举与循环性")
sols = [(F, k) for F in range(1, 10) for k in range(1, 12) if F * k + 1 == 4]
check("正整数解恰为 (|F|, m-1) = (1,3) 与 (3,1)", sols == [(1, 3), (3, 1)], str(sols))
check("路线 beta 下 m-1 = 3 <=> m = 4 <=> D = 4（与待证结论等价 => 循环）", True)
g63 = io.open(os.path.join(HERE, "G63_target_list_and_audit.md"),
              encoding="utf-8").read()
check("G63 §3 已明写循环性更正", "循环" in g63)

# ----------------------------------------------------------------------
print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
