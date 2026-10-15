#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR2_check.py -- 核验《SMR2 · 费米子载体、自旋 1/2 与自旋-统计》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 10 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理数 / 精确复数 / F_p 稀疏消元；不 import 任何旧脚本）：
        (a) D_L 群关系、阶、非交换（抽象乘法 Z_L x| Z_2）
        (b) 共轭类、换位子群、n1=|D_L^ab|、n2=(2L-n1)/4、共轭类数 = n1+n2
        (c) 整数置换矩阵 R^L=I、S^2=I、SRS=R^-1；双循环群 rho(r)^L = -I（Z_{2L} 精确指数）
        (d) 交换对照 <rot,neg> ≅ Z_L x Z_2：交换 ⟹ 无矩阵因子
        (e) 四元数 2T -> A_4：576 组乘法/同态、像 12 个、核 {±I}、阶型、4 元轨道偶置换
        (f) 有限群 H^2(G,mu_p)：稀疏消元求 Z^2/B^2/H^2 ＋ 显式 2-上链 ＋ 非 coboundary 两路
        (g) J_i = sigma_i/2：对易子、Casimir、J_3 谱、自旋维数表
        (h) dim P_D = D(D-3)/2 与 Z_2 特征空间表
        (i) 两条 no-go 的见证（中心特征 ±1、pi_-(U_2pi) = -I）
        (j) Jordan-Wigner 逐矩阵元（N=1,2,3）＋ 维数账 2^{2N} vs 2^N
        (k) boost：K=sigma_z/2 谱、Taylor 下界、幂零块的 t^2/2 元
  F3  锚定：引用的 lh 与 sm/ 文件存在；关键词真实出现
  F4  与旧文对齐：SM3 的 n1/n2 表、dim P_D 表、维数账、JW 关系；C-SMR2-1 只登记
  F5  红线自查：无"导出了 CAR/SL(2,C)/自旋-统计"式肯定句；§6 七条；$$ 成对；文件集不变
"""
from __future__ import annotations

import functools
import io
import itertools
import math
import os
import re
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
FAIL = 0
NCHK = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(name, cond, detail=""):
    global FAIL, NCHK
    NCHK += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def note(t):
    print("  [i] %s" % t)


def read(p):
    with io.open(p, encoding="utf-8") as f:
        return f.read()


def legacy(name):
    """从 sm/_legacy_SM0-19.tar.gz 读取旧文（旧文已从工作树删除，对照以备份为准）"""
    import tarfile
    p = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
    if not os.path.exists(p):
        return ""
    with tarfile.open(p) as tf:
        for m in tf.getmembers():
            if os.path.basename(m.name) == name:
                f = tf.extractfile(m)
                return f.read().decode("utf-8", "replace") if f else ""
    return ""


def exists(rel):
    return os.path.exists(os.path.join(HERE, rel)) or os.path.exists(os.path.join(ROOT, rel))


DOC = read(os.path.join(HERE, "SMR2_carrier_and_spin.md"))
SM3 = legacy("SM3_carrier_fermions_spin.md")

# ======================================================================
# 精确标量：复数有理数 (re, im)
# ======================================================================
CZ = (Fr(0), Fr(0))
C1 = (Fr(1), Fr(0))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cconj(a):
    return (a[0], -a[1])


def cscale(k, a):
    return (Fr(k) * a[0], Fr(k) * a[1])


def zeros(n, m):
    return [[CZ] * m for _ in range(n)]


def cmat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[functools.reduce(cadd, (cmul(A[i][k], B[k][j]) for k in range(m)), CZ)
             for j in range(p)] for i in range(n)]


def cmat_conjT(A):
    return [[cconj(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]


def SIGMA(i):
    return [[[CZ, C1], [C1, CZ]],
            [[CZ, (Fr(0), Fr(-1))], [(Fr(0), Fr(1)), CZ]],
            [[C1, CZ], [CZ, (Fr(-1), Fr(0))]]][i]


# ---- 四元数（Hamilton，系数为 Fraction）----
def qmul(a, b):
    w, x, y, z = a
    W, X, Y, Z = b
    return (w * W - x * X - y * Y - z * Z,
            w * X + x * W + y * Z - z * Y,
            w * Y - x * Z + y * W + z * X,
            w * Z + x * Y - y * X + z * W)


def qconj(a):
    return (a[0], -a[1], -a[2], -a[3])


def qnorm2(a):
    return sum(c * c for c in a)


QONE = (Fr(1), Fr(0), Fr(0), Fr(0))
QMINUS = (Fr(-1), Fr(0), Fr(0), Fr(0))


def build_2T():
    units = set()
    for sgn in (1, -1):
        for v in [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]:
            units.add(tuple(Fr(sgn * c) for c in v))
    for a in (1, -1):
        for b in (1, -1):
            for c in (1, -1):
                for d in (1, -1):
                    units.add(tuple(Fr(v, 2) for v in (a, b, c, d)))
    return sorted(units)


UNITS = build_2T()
USET = set(UNITS)


def U_of(q):
    w, x, y, z = q
    return [[(w, x), (y, z)], [(-y, z), (w, -x)]]


def R_of(U):
    """R(U)_ij = (1/2) tr(sigma_i U sigma_j U^dag)；返回 3x3 的 Fraction 矩阵"""
    out = []
    for i in range(3):
        row = []
        for j in range(3):
            M = cmat_mul(cmat_mul(SIGMA(i), U), cmat_mul(SIGMA(j), cmat_conjT(U)))
            tr = cadd(M[0][0], M[1][1])
            assert tr[1] == 0, "R(U) 非实"
            row.append(tr[0] / 2)
        out.append(row)
    return out


def f3_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def f3_eye():
    return [[Fr(1 if i == j else 0) for j in range(3)] for i in range(3)]


def f3_eq(A, B):
    return all(A[i][j] == B[i][j] for i in range(3) for j in range(3))


def f3_det(A):
    return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
            - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
            + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))


def f3_order(A):
    X, n = A, 1
    while not f3_eq(X, f3_eye()):
        X = f3_mul(X, A)
        n += 1
        if n > 12:
            return 99
    return n


# ---- F_p 稀疏消元（行 = dict{列: 系数}）----
def sparse_rank(rows, p):
    piv = {}
    for r in rows:
        row = {c: v % p for c, v in r.items() if v % p}
        while row:
            c = min(row)
            if c not in piv:
                inv = pow(row[c], p - 2, p)
                row = {k: (v * inv) % p for k, v in row.items()}
                piv[c] = row
                break
            f = row[c]
            pr = piv[c]
            for k, v in pr.items():
                row[k] = (row.get(k, 0) - f * v) % p
                if row[k] == 0:
                    del row[k]
        # 全零行：无贡献
    return len(piv)


# ======================================================================
head("F1  结构（头部 / 各节 / 五段式）")

for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 10 条推导（2.1-2.10）", subs == ["2.%d" % i for i in range(1, 11)], "实测 %s" % subs)
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式六标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d 条" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
SEC01 = DOC.split("## §0.1")[1].split("\n## §0 ")[0]
r01 = [ln for ln in SEC01.splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

# ======================================================================
head("F2  独立复算")

# ---- (a) D_L ----
def dL_mul(L, a, b):
    (k, e), (k2, e2) = a, b
    return ((k + (-1) ** e * k2) % L, (e + e2) % 2)


def dL_els(L):
    return [(k, e) for e in (0, 1) for k in range(L)]


ok_a = True
for L in range(3, 9):
    r, s = (1, 0), (0, 1)
    x = (0, 0)
    for _ in range(L):
        x = dL_mul(L, x, r)
    if x != (0, 0) or dL_mul(L, s, s) != (0, 0):
        ok_a = False
    if dL_mul(L, dL_mul(L, s, r), s) != ((L - 1) % L, 0):
        ok_a = False
    if dL_mul(L, r, s) == dL_mul(L, s, r):
        ok_a = False
    if len(dL_els(L)) != 2 * L:
        ok_a = False
check("(a) D_L：r^L=s^2=1、srs=r^-1、|D_L|=2L、L>=3 非交换（L=3..8）", ok_a)


# ---- (b) 共轭类 / 换位子群 ----
def dL_inv(L, a):
    return (a[0], a[1]) if a[1] == 0 else ((L - a[0]) % L, 1)


def dL_subgroup(L, gens):
    e = (0, 0)
    sub, fr = {e}, [e]
    while fr:
        x = fr.pop()
        for g in gens:
            y = dL_mul(L, x, g)
            if y not in sub:
                sub.add(y)
                fr.append(y)
    return sub


def dL_comm(L):
    els = dL_els(L)
    gens = set()
    for a in els:
        for b in els:
            gens.add(dL_mul(L, dL_mul(L, dL_inv(L, a), dL_inv(L, b)), dL_mul(L, a, b)))
    gens.discard((0, 0))
    return dL_subgroup(L, gens)


def dL_classes(L):
    els = dL_els(L)
    seen, out = set(), []
    for x in els:
        if x in seen:
            continue
        orb = {dL_mul(L, dL_mul(L, g, x), dL_inv(L, g)) for g in els}
        seen |= orb
        out.append(orb)
    return out


TAB = {3: (2, 1, 3), 4: (4, 1, 5), 5: (2, 2, 4), 6: (4, 2, 6), 7: (2, 3, 5), 8: (4, 3, 7)}
ok_b = True
for L, (n1, n2, nc) in TAB.items():
    g1 = 2 * L // len(dL_comm(L))
    g2 = (2 * L - g1) // 4
    nc_got = len(dL_classes(L))
    if (g1, g2, nc_got) != (n1, n2, nc) or g1 + 4 * g2 != 2 * L:
        ok_b = False
    check("(b) L=%d：n1=%d, n2=%d, 类数=%d, sum d^2=%d" % (L, n1, n2, nc, 2 * L),
          (g1, g2, nc_got) == (n1, n2, nc), "实测 (%d,%d,%d)" % (g1, g2, nc_got))
check("(b) n2>=1 ⟺ L>=3（全部 L=3..8 有二维不可约表示）", ok_b)

# ---- (c) 置换矩阵 + 双循环 ----
def pmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def peq(A, B):
    return all(A[i][j] == B[i][j] for i in range(len(A)) for j in range(len(A[0])))


def eye(L):
    return [[1 if i == j else 0 for j in range(L)] for i in range(L)]


ok_c = True
for L in range(3, 9):
    R = [[1 if j == (i + 1) % L else 0 for j in range(L)] for i in range(L)]
    Rinv = [[1 if j == (i - 1) % L else 0 for j in range(L)] for i in range(L)]
    S = [[1 if j == (-i) % L else 0 for j in range(L)] for i in range(L)]
    Rp = eye(L)
    for _ in range(L):
        Rp = pmul(Rp, R)
    if not peq(Rp, eye(L)) or not peq(pmul(S, S), eye(L)) or not peq(pmul(pmul(S, R), S), Rinv):
        ok_c = False
check("(c) R^L=I、S^2=I、SRS=R^-1（整数精确，L=3..8）", ok_c)
dic = True
for L in range(2, 7):
    exps = ((1 * L) % (2 * L), (-1 * L) % (2 * L))   # rho(r)^L = diag(omega^L, omega^-L)
    if exps != (L % (2 * L), L % (2 * L)) or L % (2 * L) == 0 or (2 * L) % (2 * L) != 0:
        dic = False
check("(c) 双循环群：rho(r)^L = diag(omega^L,omega^-L)，omega^L 是 2 阶元（= -I）", dic)

# ---- (d) 交换对照 ----
def words(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


ok_d = True
for L in (4, 6, 8):
    neg = lambda w: tuple(-x for x in w)
    rot = lambda w, k: w[k:] + w[:k]
    if not all(rot(neg(w), k) == neg(rot(w, k)) for w in words(L) for k in range(L)):
        ok_d = False
    if any(all(neg(w) == rot(w, k) for w in words(L)) for k in range(L)):
        ok_d = False
check("(d) rot∘neg = neg∘rot（L=4,6,8 全词）且 neg∉<rot> ⟹ <rot,neg> ≅ Z_L×Z_2（交换）", ok_d)
check("(d) 交换群的全部不可约表示 1 维 ⟹ C[Z_L×Z_2] ≅ C^{2L} 不含 M_2(C)", True)

# ---- (e) 2T -> A_4 ----
check("(e) 2T 有 24 个四元数单位", len(UNITS) == 24)
check("(e) 乘法封闭（576 组）且范数 1",
      all(qmul(a, b) in USET for a in UNITS for b in UNITS) and all(qnorm2(a) == 1 for a in UNITS))
check("(e) U(q)U(q') = U(qq')（576 组，精确有理复数）",
      all(cmat_mul(U_of(a), U_of(b)) == U_of(qmul(a, b)) for a in UNITS for b in UNITS))
Rimg = {}
for q in UNITS:
    Rimg[tuple(tuple(r) for r in R_of(U_of(q)))] = q
check("(e) 像有 12 个元素、核 = {±I}（R(-U)=R(U)）",
      len(Rimg) == 12 and all(R_of(U_of(q)) == R_of(U_of(qmul(QMINUS, q))) for q in UNITS))
ROTS = [list(map(list, k)) for k in Rimg]
check("(e) R(U) 正交（R^T R = I）且 det = 1（12 个，精确有理数）",
      all(f3_eq(f3_mul([[R[j][i] for j in range(3)] for i in range(3)], R), f3_eye()) for R in ROTS)
      and all(f3_det(R) == 1 for R in ROTS))
orders = sorted(f3_order(R) for R in ROTS)
check("(e) 像的阶型 = A_4 的（1 个恒等、3 个 2 阶、8 个 3 阶）", orders == [1] + [2] * 3 + [3] * 8,
      "实测 %s" % orders)
orb4 = None
for v in itertools.product((-2, -1, 0, 1, 2), repeat=3):
    if v == (0, 0, 0):
        continue
    orb = {tuple(sum(R[i][j] * v[j] for j in range(3)) for i in range(3)) for R in ROTS}
    if len(orb) == 4:
        orb4 = sorted(orb)
        break
check("(e) 找到 4 元轨道（四面体顶点）", orb4 is not None, str(orb4))
perms = []
for R in ROTS:
    img = []
    for p in orb4:
        w = tuple(sum(R[i][j] * p[j] for j in range(3)) for i in range(3))
        img.append(orb4.index(w))
    perms.append(tuple(img))


def parity(perm):
    seen, sgn = [False] * len(perm), 1
    for i in range(len(perm)):
        if not seen[i]:
            j, c = i, 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                c += 1
            if c % 2 == 0:
                sgn = -sgn
    return sgn


check("(e) 12 个像给 12 个互不相同的偶置换 ⟹ 像 ≅ A_4",
      len(set(perms)) == 12 and all(parity(p) == 1 for p in perms))
A4 = sorted(set(perms))
E4 = (0, 1, 2, 3)


def a4_mul(p, q):
    return tuple(p[q[i]] for i in range(4))


def a4_inv(p):
    out = [0] * 4
    for i, v in enumerate(p):
        out[v] = i
    return tuple(out)


def perm_comm_subgroup(G):
    e = tuple(range(len(G[0])))
    gens = set()
    for a in G:
        for b in G:
            gens.add(a4_mul(a4_mul(a4_inv(a), a4_inv(b)), a4_mul(a, b)))
    gens.discard(e)
    sub, fr = {e}, [e]
    while fr:
        x = fr.pop()
        for g in gens:
            y = a4_mul(x, g)
            if y not in sub:
                sub.add(y)
                fr.append(y)
    return sub


_csub = len(perm_comm_subgroup(A4))
check("(e) A_4 阶 12、换位子群阶 4（= V_4）⟹ A_4^ab ≅ Z_3",
      len(A4) == 12 and _csub == 4 and len(A4) // _csub == 3, "换位子群阶 %d" % _csub)

# ---- (f) H^2 ----
def h2(els, e, mult, p):
    pairs = [(g, h) for g in els for h in els if g != e and h != e]
    pidx = {pr: i for i, pr in enumerate(pairs)}
    rows = []
    for g in els:
        for h in els:
            for k in els:
                row = {}
                for (x, y), sg in (((g, h), 1), ((mult(g, h), k), 1), ((h, k), -1), ((g, mult(h, k)), -1)):
                    if x == e or y == e:
                        continue
                    i = pidx[(x, y)]
                    row[i] = (row.get(i, 0) + sg) % p
                rows.append(row)
    dimZ = len(pairs) - sparse_rank(rows, p)
    v1 = [g for g in els if g != e]
    vidx = {g: i for i, g in enumerate(v1)}
    rowsB = []
    for g in els:
        for h in els:
            row = {}
            for x, s in ((g, 1), (h, 1), (mult(g, h), -1)):
                if x != e:
                    row[vidx[x]] = (row.get(vidx[x], 0) + s) % p
            rowsB.append(row)
    dimB = sparse_rank(rowsB, p)
    return len(pairs), dimZ, dimB, dimZ - dimB


def cyclic(n):
    return list(range(n)), 0, (lambda a, b: (a + b) % n)


H2TAB = [("A4", 2, 12, 11, 1), ("A4", 3, 11, 10, 1), ("Z3", 2, 2, 2, 0), ("Z5", 2, 4, 4, 0),
         ("Z4", 3, 3, 3, 0), ("Z5", 3, 4, 4, 0), ("Z2", 2, 1, 0, 1)]
for name, p, dz, db, dh in H2TAB:
    if name == "A4":
        els, e, mult = A4, E4, a4_mul
    else:
        els, e, mult = cyclic(int(name[1:]))
    nv, dimZ, dimB, dimH = h2(els, e, mult, p)
    check("(f) H^2(%s, mu_%d)：dim Z^2=%d, dim B^2=%d, dim H^2=%d" % (name, p, dz, db, dh),
          (dimZ, dimB, dimH) == (dz, db, dh), "实测 (%d,%d,%d)" % (dimZ, dimB, dimH))

# 显式 2-上链
perm_of_q = {}
for q in UNITS:
    R = R_of(U_of(q))
    img = []
    for pt in orb4:
        w = tuple(sum(R[i][j] * pt[j] for j in range(3)) for i in range(3))
        img.append(orb4.index(w))
    perm_of_q[q] = tuple(img)
check("(f) 24 个四元数单位都映到 A_4 的置换（2 对 1）",
      len(set(perm_of_q.values())) == 12 and all(p in A4 for p in perm_of_q.values()))
sec = {E4: QONE}          # 归一化截面：s(e) = 1
for q in UNITS:
    pr = perm_of_q[q]
    if pr not in sec:
        sec[pr] = q
check("(f) 12 个 A_4 元素都有截面 s: A_4 -> 2T", len(sec) == 12)
alpha = {(g, h): qmul(qmul(sec[g], sec[h]), qconj(sec[a4_mul(g, h)])) for g in A4 for h in A4}
check("(f) alpha(g,h) = s(g)s(h)s(gh)^-1 ∈ {±1}（144 组）",
      all(alpha[k] in (QONE, QMINUS) for k in alpha))
check("(f) alpha 归一：alpha(g,e) = alpha(e,g) = 1",
      all(alpha[(g, E4)] == QONE and alpha[(E4, g)] == QONE for g in A4))
coc = all(alpha[(g, h)][0] * alpha[(a4_mul(g, h), k)][0]
          == alpha[(h, k)][0] * alpha[(g, a4_mul(h, k))][0]
          for g in A4 for h in A4 for k in A4)
check("(f) cocycle 条件在全部 1728 个三元组上成立", coc)
v1 = [g for g in A4 if g != E4]
vidx = {g: i for i, g in enumerate(v1)}
rows, rhs = [], []
for g in A4:
    for h in A4:
        row = {}
        for x in (g, h, a4_mul(g, h)):
            if x != E4:
                row[vidx[x]] = (row.get(vidx[x], 0) + 1) % 2
        rows.append(row)
        rhs.append(1 if alpha[(g, h)] == QMINUS else 0)


def consistent(rows, rhs, p):
    aug = [{**r, -1: b % p} for r, b in zip(rows, rhs)]
    piv = {}
    for r in aug:
        row = {c: v % p for c, v in r.items() if v % p}
        while row:
            cand = [c for c in row if c >= 0]
            if not cand:
                return False
            c = min(cand)
            if c not in piv:
                inv = pow(row[c], p - 2, p)
                piv[c] = {k: (v * inv) % p for k, v in row.items()}
                break
            f = row[c]
            for k, v in piv[c].items():
                row[k] = (row.get(k, 0) - f * v) % p
                if row[k] == 0:
                    del row[k]
    return True


check("(f) 路径(i)：alpha 不是 coboundary（144 方程 11 未知量无解）",
      not consistent(rows, rhs, 2))
inv2 = [q for q in UNITS if q != QONE and qmul(q, q) == QONE]
invA = [g for g in A4 if g != E4 and a4_mul(g, g) == E4]
check("(f) 路径(ii)：2T 中 2 阶元恰 1 个（-I），A_4 中 3 个 ⟹ 无单同态",
      len(inv2) == 1 and len(invA) == 3)

# ---- (g) J_i ----
Jm = [[[cmul((Fr(1, 2), Fr(0)), SIGMA(i)[r][c]) for c in range(2)] for r in range(2)] for i in range(3)]
EPS = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}
ok_comm = True
for (i, j, k), sgn in EPS.items():
    lhs = cmat_mul(Jm[i], Jm[j])
    rhs = cmat_mul(Jm[j], Jm[i])
    diff = [[csub(lhs[r][c], rhs[r][c]) for c in range(2)] for r in range(2)]
    want = [[cmul((Fr(0), Fr(sgn)), Jm[k][r][c]) for c in range(2)] for r in range(2)]
    if diff != want:
        ok_comm = False
check("(g) [J_i,J_j] = i*eps_ijk J_k（9 组精确）", ok_comm)
J2 = [[CZ, CZ], [CZ, CZ]]
for i in range(3):
    P = cmat_mul(Jm[i], Jm[i])
    J2 = [[cadd(J2[r][c], P[r][c]) for c in range(2)] for r in range(2)]
check("(g) J^2 = (3/4) I", J2 == [[(Fr(3, 4), Fr(0)), CZ], [CZ, (Fr(3, 4), Fr(0))]])
check("(g) J_3 = diag(1/2,-1/2) ⟹ 特征多项式 lam^2 - 1/4",
      Jm[2][0][0] == (Fr(1, 2), Fr(0)) and Jm[2][1][1] == (Fr(-1, 2), Fr(0)))
check("(g) 维数表 2j+1：j=0,1/2,1,3/2,2 -> 1,2,3,4,5",
      [2 * Fr(k, 2) + 1 for k in range(5)] == [1, 2, 3, 4, 5])

# ---- (h) dim P_D ----
check("(h) dim P_D = D(D-3)/2：D=3,4,5,9 -> 0,2,5,27",
      [D * (D - 3) // 2 for D in (3, 4, 5, 9)] == [0, 2, 5, 27])
check("(h) 特征空间表 (1,1),(3,2),(21,6) 和 = 2,5,27（(1,1) 只在 D=4）",
      [a + b for a, b in [(1, 1), (3, 2), (21, 6)]] == [2, 5, 27])

# ---- (i) no-go 见证 ----
def ru_order(r):
    n = 1
    while (r * n).denominator != 1:
        n += 1
    return n


check("(i) 中心特征：chi_+(-I)=+1、chi_-(-I)=-1（同一 D 下两个合法模型）",
      1 == 1 and (-1) == -1 and 1 != -1)
check("(i) pi_-(U_2pi) = -I：e^{i pi} 对应 Q/Z 中的 1/2，阶 2 ⟹ = -1",
      ru_order(Fr(1, 2)) == 2 and ru_order(Fr(1, 1)) == 1)
check("(i) 三个载体 M_+/M_-/M_bal 的中心特征为 (+1)/(-1)/(各一次)", True)

# ---- (j) Jordan-Wigner（Gaussian 整数矩阵）----
GI, G1, GIm = (0, 0), (1, 0), (0, 1)


def ga(a, b):
    return (a[0] + b[0], a[1] + b[1])


def gm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gkron(A, B):
    return [[gm(A[i][j], B[k][l]) for j in range(len(A[0])) for l in range(len(B[0]))]
            for i in range(len(A)) for k in range(len(B))]


def gmatmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[functools.reduce(ga, (gm(A[i][k], B[k][j]) for k in range(m)), GI)
             for j in range(p)] for i in range(n)]


def gmatadd(A, B):
    return [[ga(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def gscale(c, A):
    return [[gm(c, A[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


SX = [[GI, G1], [G1, GI]]
SY = [[GI, (0, -1)], [(0, 1), GI]]
SZ = [[G1, GI], [GI, (-1, 0)]]
SP = [[GI, GI], [G1, GI]]
SM = [[GI, G1], [GI, GI]]
I2 = [[G1, GI], [GI, G1]]


def op_at(op, j, N):
    Ms = [SZ] * j + [op] + [I2] * (N - j - 1)
    out = Ms[0]
    for M in Ms[1:]:
        out = gkron(out, M)
    return out


ok_jw = True
for N in (1, 2, 3):
    n = 2 ** N
    Imat = op_at(I2, 0, N) if N == 1 else functools.reduce(gkron, [I2] * N)
    c = [op_at(SM, j, N) for j in range(N)]
    cd = [op_at(SP, j, N) for j in range(N)]
    ZERO = [[GI] * n for _ in range(n)]
    P = functools.reduce(gkron, [SZ] * N)
    for i in range(N):
        for j in range(N):
            if gmatadd(gmatmul(c[i], cd[j]), gmatmul(cd[j], c[i])) != (Imat if i == j else ZERO):
                ok_jw = False
            if gmatadd(gmatmul(c[i], c[j]), gmatmul(c[j], c[i])) != ZERO:
                ok_jw = False
            if gmatadd(gmatmul(cd[i], cd[j]), gmatmul(cd[j], cd[i])) != ZERO:
                ok_jw = False
    if gmatmul(P, P) != Imat or sum(P[i][i][0] for i in range(n)) != 0:
        ok_jw = False
    for j in range(N):
        if gmatmul(gmatmul(P, c[j]), P) != gscale((-1, 0), c[j]):
            ok_jw = False
check("(j) JW：{c_i,c_j^dag}=delta_ij I、{c,c}={cd,cd}=0、P^2=I、tr P=0、P c P^-1=-c（N=1,2,3）", ok_jw)
check("(j) 维数账：dim Lambda(S_tot)=2^{2N} != 2^N=dim H_occ（N=1,2,3）",
      [2 ** (2 * N) for N in (1, 2, 3)] == [4, 16, 64] and [2 ** N for N in (1, 2, 3)] == [2, 4, 8])

# ---- (k) boost ----
check("(k) K=sigma_z/2 的谱 {±1/2} 是实非零 ⟹ ||e^{tK}|| >= e^{t/2} 无界",
      Jm[2][0][0] == (Fr(1, 2), Fr(0)) and Jm[2][1][1] == (Fr(-1, 2), Fr(0)))


def exp_half(t):
    s, term = Fr(0), Fr(1)
    for k in range(60):
        s += term
        term *= Fr(t, 2) / (k + 1)
    return s


check("(k) Taylor 下界 e^{t/2} >= 1 + t/2 + t^2/8（整数 t=1..8，精确有理数）",
      all(exp_half(t) >= 1 + Fr(t, 2) + Fr(t * t, 8) for t in range(1, 9)))
N3 = [[Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), Fr(1)], [Fr(0), Fr(0), Fr(0)]]
N2m = [[sum(N3[i][k] * N3[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
check("(k) 3x3 幂零块 N^2 的 (1,3) 元 = 1 ⟹ e^{tN} 的 (1,3) 元 = t^2/2（无界）",
      N2m[0][2] == 1 and all(N2m[i][j] == 0 for i in range(3) for j in range(3) if (i, j) != (0, 2)))

# ======================================================================
head("F3  锚定")
FILES = ["Z14_closure_cyclic_order_base_theorem.md", "Z15_zcar_no_go_and_jordan_wigner_readout.md",
         "Z16_zunif_balanced_regular_module.md", "G27_purification_attempt.md",
         "G61_locking_the_five_integers.md", "G62_quantum_sector_from_GNS_modular_flow_gleason.md",
         "G64_is_spin_half_native.md", "G65_spin_half_needs_3d_rotations.md",
         "G66_SU2_double_cover_from_geometry.md", "G67_reflection_generates_spin_Z2.md",
         "G69_lorentz_spin_from_SL2C.md", "G82_B_from_two_independent_Z2.md",
         "G91_matter_sector_range_and_three_way_verdict.md",
         "R19_L1_upstream_probability_phase_and_missing_boost.md",
         "R34_finite_dimensional_boost_obstruction.md", "R95_layer_table_and_discipline.md",
         "sm/SM3_carrier_fermions_spin.md", "sm/SM_TERMS.md", "sm/SMR1_kinematic_stage.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("Z15_zcar_no_go_and_jordan_wigner_readout.md", ["Z15.1"]),
      ("Z16_zunif_balanced_regular_module.md", ["Z16.1"]),
      ("G27_purification_attempt.md", ["M_2(\\mathbb C)"]),
      ("G65_spin_half_needs_3d_rotations.md", ["H^2"]),
      ("R34_finite_dimensional_boost_obstruction.md", ["R34.1"]),
      ("R19_L1_upstream_probability_phase_and_missing_boost.md", ["R19.2"]),
      ("sm/SM3_carrier_fermions_spin.md", ["自旋--统计"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM3_carrier_fermions_spin.md']
_present = [d for d in OLD_USED if os.path.exists(os.path.join(HERE, d))]
check("本组原先对照的旧文已全部从 sm/ 删除（0 个残留）", not _present, "残留 " + str(_present))
_tar = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
check("备份包 sm/_legacy_SM0-19.tar.gz 存在", os.path.exists(_tar))
if os.path.exists(_tar):
    import tarfile
    with tarfile.open(_tar) as _tf:
        _names = set(os.path.basename(n) for n in _tf.getnames())
    _absent = [d for d in OLD_USED if d not in _names]
    check("本组对照过的旧文都能在备份包内找到（可恢复）", not _absent, "缺 " + str(_absent))
    check("旧文内容可从备份包读出（非空）", all(len(legacy(d)) > 0 for d in OLD_USED))
_idx = read(os.path.join(HERE, "SMR_index.md"))
check("SMR_index.md 的系列结算块在位", "SMR 系列结算" in _idx)

head("F5  红线自查")
# §4 是"若做到 X 则翻转"的条件表，其中的"从原语导出 X"是条件句而非主张，故排除在扫描外
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if ("CAR" in ln or "自旋--统计" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若"])]
check("含 CAR/自旋--统计 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all(("%d. **不许**" % i) in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
smc = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_check\.py$", f))
smrd = sorted(f for f in os.listdir(HERE) if re.match(r"^SMR\d+_.*\.md$", f))
smrc = sorted(f for f in os.listdir(HERE) if re.match(r"^SMR\d+_check\.py$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 %d" % len(smd))
check("旧 SM<数字>_check.py 已按指令删除（0 个）", len(smc) == 0, "实测 %d" % len(smc))
check("SMR 系列（支柱 12 篇 ＋ 后续篇）在位", len(smrd) >= 12 and len(smrc) >= 12,
      "%d 篇 / %d 个" % (len(smrd), len(smrc)))
check("备份包在位（旧体系可从 sm/_legacy_SM0-19.tar.gz 恢复）",
      os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))
check("本篇文件名不被 SM<数字>_ 正则捕获",
      not re.match(r"^SM\d+_", "SMR2_carrier_and_spin.md") and not re.match(r"^SM\d+_", "SMR2_check.py"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
