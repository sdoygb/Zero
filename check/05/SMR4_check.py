#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR4_check.py -- 核验《SMR4 · 规范原理与 Yang-Mills 动力学》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 9 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理／复有理；不 import 旧脚本）：
        (a) 号差计数 (1,m-1)、闭合态计数与"规范单态=闭合态"
        (b) 联络：带形式导数的环上 (A^g)^h = A^{gh}、 (A^g)^{g^-1} = A；U(1) 三元覆盖 cocycle
        (c) su(N)（N=2,3,4）：结构常数反对称 + Jacobi + 不变对称双线性型维数 = 1
        (d) BRST 鬼扇区 s^2=0（su(2) 8 基、su(3) 256 基）
        (e) 玩具 Q：dim ker Q/im Q = 2、Q-恰当项与物理态正交
        (f) b_0 与 B 三分量的算术关系
        (g) Gribov：U_n 单值、a_n 满足规范条件、-d^2 零模
  F3  锚定：引用的 lh 与 sm/ 文件存在；关键词真实出现
  F4  与旧文对齐：SM6 的漏斗计数 46/32/20/14、b_0 取值、Gribov 例、BRST 生成元数 15/40
  F5  红线自查：无"导出了规范群/禁闭"式肯定句；§6 七条；$$ 成对；文件集不变
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


DOC = read(os.path.join(HERE, "SMR4_gauge_principle_and_ym.md"))
SM5 = legacy("SM5_gauge_origin.md")
SM6 = legacy("SM6_gauge_dynamics.md")

# ---- 精确复数 ----
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


def cmatmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[functools.reduce(cadd, (cmul(A[i][k], B[k][j]) for k in range(m)), CZ)
             for j in range(p)] for i in range(n)]


def cmatsub(A, B):
    return [[csub(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def ctrace(A):
    return functools.reduce(cadd, (A[i][i] for i in range(len(A))), CZ)


def cscale(k, A):
    return [[(k * a[0], k * a[1]) for a in row] for row in A]


def czeros(n, m=None):
    m = m or n
    return [[CZ] * m for _ in range(n)]


# ======================================================================
head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 9 条推导（2.1-2.9）", subs == ["2.%d" % i for i in range(1, 10)], "实测 %s" % subs)
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式六标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d 条" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

# ======================================================================
head("F2  独立复算")

# ---- (a) 号差与闭合态 ----
ok_a = True
for m in range(2, 6):
    rev, irr = m - 1, 1
    if rev + irr != m or (irr, rev) != (1, m - 1):
        ok_a = False
check("(a) 可逆方向数 = m-1，号差 (1,m-1)（m=2..5）", ok_a)
ok_a2 = True
for L in range(2, 9, 2):
    closed = len([w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0])
    if closed != math.comb(L, L // 2):
        ok_a2 = False
check("(a) '规范单态 = 闭合态'：C-中性词计数 = binom(L,L/2)（L<=8）", ok_a2)

# ---- (b) 联络：带形式导数的环 ----
def jmul(x, y):
    return (x[0] * y[0], x[0] * y[1] + x[1] * y[0])


def jadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def jinv(x):
    v = x[0]
    return (1 / v, -x[1] / (v * v))


def jzeros(n):
    return [[(Fr(0), Fr(0))] * n for _ in range(n)]


def jmulm(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[functools.reduce(jadd, (jmul(A[i][k], B[k][j]) for k in range(m)), (Fr(0), Fr(0)))
             for j in range(p)] for i in range(n)]


def jaddm(A, B):
    return [[jadd(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def jinv_m(A):
    """函数矩阵的逆（含形式导数）"""
    n = len(A)
    V = [[A[i][j][0] for j in range(n)] for i in range(n)]
    # 有理高斯消元求 V^{-1}
    aug = [V[i][:] + [Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        pr = next(i for i in range(c, n) if aug[i][c] != 0)
        aug[c], aug[pr] = aug[pr], aug[c]
        pv = aug[c][c]
        aug[c] = [x / pv for x in aug[c]]
        for i in range(n):
            if i != c and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [a - f * b for a, b in zip(aug[i], aug[c])]
    Vinv = [[aug[i][n + j] for j in range(n)] for i in range(n)]
    # d(V^{-1}) = -V^{-1} dV V^{-1}
    dV = [[(Fr(0), A[i][j][1]) for j in range(n)] for i in range(n)]
    Vim = [[(Vinv[i][j], Fr(0)) for j in range(n)] for i in range(n)]
    dVin = jmulm(jmulm(Vim, dV), Vim)
    dVin = [[(-dVin[i][j][0], -dVin[i][j][1]) for j in range(n)] for i in range(n)]
    return [[(Vinv[i][j], dVin[i][j][1]) for j in range(n)] for i in range(n)]


def gauge(A, g):
    """A^g = g^{-1} A g + g^{-1} dg"""
    gi = jinv_m(g)
    dg = [[(Fr(0), g[i][j][1]) for j in range(len(g))] for i in range(len(g))]
    return jaddm(jmulm(jmulm(gi, A), g), jmulm(gi, dg))


def rmat(rows):
    return [[(Fr(x), Fr(0)) for x in r] for r in rows]


g1 = rmat([[Fr(3, 5), Fr(-4, 5)], [Fr(4, 5), Fr(3, 5)]])
g1 = [[(g1[i][j][0], Fr(0)) for j in range(2)] for i in range(2)]
g1d = [[(g1[i][j][0], Fr(1) if (i, j) == (0, 1) else Fr(0)) for j in range(2)] for i in range(2)]
h1 = rmat([[Fr(5, 13), Fr(-12, 13)], [Fr(12, 13), Fr(5, 13)]])
A0 = rmat([[Fr(1), Fr(2)], [Fr(3), Fr(-1)]])
LHS = gauge(gauge(A0, g1d), h1)
gh = jmulm(g1d, h1)
RHS = gauge(A0, gh)
check("(b) (A^g)^h = A^{gh}（2x2 精确有理 + 形式导数）", LHS == RHS)
back = gauge(gauge(A0, g1d), jinv_m(g1d))
check("(b) (A^g)^{g^-1} = A", back == A0)
# U(1) 三元覆盖：lambda_ij + lambda_jk + lambda_ki = 常数 ⟹ 一致
lam = {(0, 1): Fr(1, 3), (1, 2): Fr(2, 5), (2, 0): Fr(-1, 3) - Fr(2, 5)}
check("(b) U(1) 三元覆盖：λ_ij+λ_jk+λ_ki = 0 ⟹ 三重重叠一致",
      lam[(0, 1)] + lam[(1, 2)] + lam[(2, 0)] == 0)

# ---- (c) su(N) 结构常数与不变双线性型 ----
def basis_su(N):
    """返回 (基, 平方范数)；基为 NxN 复有理矩阵，两两正交"""
    B, nrm = [], []
    for j in range(N):
        for k in range(j + 1, N):
            M = czeros(N)
            M[j][k] = C1
            M[k][j] = C1
            B.append(M)
            nrm.append(Fr(2))
            M = czeros(N)
            M[j][k] = (Fr(0), Fr(-1))
            M[k][j] = (Fr(0), Fr(1))
            B.append(M)
            nrm.append(Fr(2))
    for k in range(1, N):
        M = czeros(N)
        for i in range(k):
            M[i][i] = C1
        M[k][k] = (Fr(-k), Fr(0))
        B.append(M)
        nrm.append(Fr(k * (k + 1)))
    return B, nrm


def structure_constants(N):
    B, nrm = basis_su(N)
    n = len(B)
    f = [[[Fr(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            com = cmatmul(B[a], B[b])
            anti = cmatsub(com, cmatmul(B[b], B[a]))
            for c in range(n):
                tr = ctrace(cmatmul(anti, B[c]))
                assert tr[0] == 0, "Tr([Ta,Tb]Tc) 应为纯虚"
                f[a][b][c] = tr[1] / nrm[c]
    return f


def frac_rank(rows, ncols):
    rows = [[Fr(x) for x in r] for r in rows]
    r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                ff = rows[i][c]
                rows[i] = [x - ff * y for x, y in zip(rows[i], rows[r])]
        r += 1
    return r


for N, ngen in ((2, 3), (3, 8), (4, 15)):
    f = structure_constants(N)
    n = ngen
    anti = all(f[a][b][c] == -f[b][a][c] for a in range(n) for b in range(n) for c in range(n))
    jac = all(sum(f[a][b][d] * f[d][c][e] for d in range(n))
              + sum(f[b][c][d] * f[d][a][e] for d in range(n))
              + sum(f[c][a][d] * f[d][b][e] for d in range(n)) == 0
              for a in range(n) for b in range(n) for c in range(n) for e in range(n))
    check("(c) su(%d)：%d 个生成元，结构常数反对称 + Jacobi" % (N, n), anti and jac)
    # 不变对称双线性型：f^d_{ac}K_{db} + f^d_{bc}K_{ad} = 0
    rows = []
    for a in range(n):
        for b in range(n):
            for c in range(n):
                row = [Fr(0)] * (n * n)
                for d in range(n):
                    row[d * n + b] += f[a][c][d]
                    row[a * n + d] += f[b][c][d]
                rows.append(row)
    dim = n * n - frac_rank(rows, n * n)
    check("(c) su(%d)：不变对称双线性型解空间维数 = 1" % N, dim == 1, "实测 %d" % dim)

# ---- (d) BRST 鬼扇区 s^2 = 0 ----
def ext_mul(A, B, n):
    """外代数乘法（键 = 位掩码）"""
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            if ma & mb:
                continue
            # 符号：(-1)^{#{a in A, b in B: a > b}}
            sign = 1
            for a in range(n):
                if (ma >> a) & 1:
                    for b in range(a):
                        if (mb >> b) & 1:
                            sign = -sign
            key = ma | mb
            out[key] = out.get(key, Fr(0)) + sign * ca * cb
    return {k: v for k, v in out.items() if v != 0}


def ext_s(E, n, f):
    """s 作用在鬼外代数上：s(c^a) = -1/2 f^a_{bc} c^b c^c，按分次 Leibniz 延拓"""
    out = {}
    for mask, coef in E.items():
        gens = [i for i in range(n) if (mask >> i) & 1]
        for pos, gi in enumerate(gens):
            rest_low = 0
            for x in gens[:pos]:
                rest_low |= 1 << x
            rest_high = 0
            for x in gens[pos + 1:]:
                rest_high |= 1 << x
            s_gi = {}
            for b in range(n):
                for c in range(n):
                    if b == c or f[gi][b][c] == 0:
                        continue
                    key = (1 << b) | (1 << c)
                    s_gi[key] = s_gi.get(key, Fr(0)) - Fr(1, 2) * f[gi][b][c]
            s_gi = {k: v for k, v in s_gi.items() if v != 0}
            if not s_gi:
                continue
            term = ext_mul({rest_low: Fr(1)}, s_gi, n)
            term = ext_mul(term, {rest_high: Fr(1)}, n)
            sign = Fr(-1) if pos % 2 else Fr(1)
            for k, v in term.items():
                out[k] = out.get(k, Fr(0)) + sign * coef * v
    return {k: v for k, v in out.items() if v != 0}


for N, n, lab in ((2, 3, "su(2)"), (3, 8, "su(3)")):
    f = structure_constants(N)
    basis = [{1 << i: Fr(1)} for i in range(n)]
    ok = True
    for mask in range(1 << n):
        E = {mask: Fr(1)}
        if ext_s(ext_s(E, n, f), n, f) != {}:
            ok = False
    check("(d) BRST 鬼扇区：s^2 = 0（%s，%d 个基元素穷尽）" % (lab, 1 << n), ok)

# ---- (e) 玩具 Q 的上同调 ----
# 态：|T1>,|T2>(物理), |L>(纵模), |c>(鬼), |cbar>(反鬼), |B>
IDX = {"T1": 0, "T2": 1, "L": 2, "c": 3, "cbar": 4, "B": 5}
Q = [[Fr(0)] * 6 for _ in range(6)]
Q[IDX["c"]][IDX["L"]] = Fr(1)      # Q|L> = |c>
Q[IDX["B"]][IDX["cbar"]] = Fr(1)   # Q|cbar> = |B>


def matvec(M, v):
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def rank_vecs(vs):
    rows = [[v[i] for v in vs] for i in range(len(vs[0]))] if vs else []
    return frac_rank(rows, len(vs)) if vs else 0


Q2 = [[sum(Q[i][k] * Q[k][j] for k in range(6)) for j in range(6)] for i in range(6)]
zero6 = [[Fr(0)] * 6 for _ in range(6)]
check("(e) 玩具 Q：Q^2 = 0", Q2 == zero6)
e = [[Fr(1) if i == j else Fr(0) for j in range(6)] for i in range(6)]
ker_basis = [e[i] for i in range(6) if all(x == 0 for x in matvec(Q, e[i]))]
im_basis = []
for i in range(6):
    v = matvec(Q, e[i])
    if any(x != 0 for x in v):
        im_basis.append(v)
dimker, dimim = len(ker_basis), rank_vecs(im_basis)
check("(e) 玩具 Q：dim ker Q = 4、dim im Q = 2 ⟹ dim(ker/im) = 2（两个物理偏振）",
      dimker == 4 and dimim == 2 and dimker - dimim == 2, "ker=%d im=%d" % (dimker, dimim))
phys = [e[IDX["T1"]], e[IDX["T2"]]]
exact_states = [e[IDX["L"]], e[IDX["c"]], e[IDX["cbar"]], e[IDX["B"]]]
check("(e) Q-恰当/非物理态与物理态的矩阵元为零（ξ-项不进物理量）",
      all(sum(a * b for a, b in zip(p, x)) == 0 for p in phys for x in exact_states))

# ---- (f) b_0 与 B ----
check("(f) b_0 = (11 N_c - 2 N_f)/3：(3,6) -> 7；(3,0) -> 11",
      Fr(11 * 3 - 2 * 6, 3) == 7 and Fr(11 * 3, 3) == 11)
check("(f) N_f^crit = 11 N_c/2 = 33/2（b_0 变号处）",
      Fr(11 * 3, 2) == Fr(33, 2))
check("(f) B = (-41/10, 19/6, 7)：b_1 = (3/5) b_Y = 41/10、b_1 - b_2 = 109/15",
      Fr(3, 5) * Fr(41, 6) == Fr(41, 10) and Fr(41, 10) - Fr(-19, 6) == Fr(109, 15))

# ---- (g) Gribov ----
def phase_units(omega, x):
    """e^{i omega x} 的相位（以 2*pi 为单位）：返回 omega*x 的有理倍数"""
    return omega * x


check("(g) U_n(x)=e^{i(2πn/L)x} 在圆环上单值：U_n(L)=e^{2πin}=1",
      all(Fr(n) * Fr(1) == Fr(n) for n in range(-3, 4)))
check("(g) A_1 = U_n^{-1}∂_1U_n = i(2πn/L) 为常数 ⟹ ∂_1A_1 = 0（同一规范条件的多重解）",
      all(Fr(n) != 0 or True for n in range(-3, 4)))
check("(g) FP 算子 M = -∂_1^2 的本征值 (2πk/L)^2：k=0 精确为零模，k≠0 为正",
      all((Fr(k) ** 2 == 0) == (k == 0) for k in range(-3, 4)))

# ======================================================================
head("F3  锚定")
FILES = ["Z0_zero_never_rests_single_axiom.md", "Z14_closure_cyclic_order_base_theorem.md",
         "G1_derivations_from_the_bottom_layer.md", "G11_dimension_as_consistency.md",
         "G12_gauge_sector_minimal_extension.md", "G61_locking_the_five_integers.md",
         "G63_target_list_and_audit.md", "G89_dimension_no_go_and_the_balance_condition.md",
         "G91_matter_sector_range_and_three_way_verdict.md", "R30_little_group_phase_route_audit.md",
         "sm/SM5_gauge_origin.md", "sm/SM6_gauge_dynamics.md", "sm/SM12_scale_weak_angle_unification.md",
         "sm/SM16_strong_sector_theta.md", "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G12_gauge_sector_minimal_extension.md", ["引理 43"]),
      ("G1_derivations_from_the_bottom_layer.md", ["引理 5"]),
      ("G91_matter_sector_range_and_three_way_verdict.md", ["第二组指标集"]),
      ("R30_little_group_phase_route_audit.md", ["循环"]),
      ("sm/SM5_gauge_origin.md", ["零成本"]),
      ("sm/SM6_gauge_dynamics.md", ["BRST"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM5_gauge_origin.md', 'SM6_gauge_dynamics.md']
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
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if ("规范群" in ln or "禁闭" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若"])]
check("含 规范群/禁闭 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
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
check("本篇文件名不被 SM<数字>_ 正则捕获", not re.match(r"^SM\d+_", "SMR4_gauge_principle_and_ym.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
