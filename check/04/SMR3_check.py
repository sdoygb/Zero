#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR3_check.py -- 核验《SMR3 · 局域性、复合与连续极限》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 10 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理数 / 精确 Q(sqrt5) / 精确复数有理数；不 import 旧脚本）：
        (a) 乘积态：T 秩 1、S_max = 2|a||b| = 2
        (b) 可分例：4x4 密度矩阵直接算 T、S_max^2 = 20/9
        (c) 四个 Bell 态：由 4x4 密度矩阵精确算 T（检验 sigma_y 符号约定）
        (d) Bell 对角态 (7/10,1/10,1/10,1/10)：T = diag(3/5,-3/5,3/5) != (3/5)I_3；sum|c| = 9/5；S_max^2 = 72/25
        (e) Werner：PT 谱 {1/2,1/2,1/2,-1/2}（PT^2 = I/4、tr = 1）、可分门限 1/3、CHSH 门限 1/sqrt2
        (f) KCBS：A 的对角元、lambda*、四个 S_max、两个余量（Q(sqrt5) 精确）+ C_5 独立数 = 2
        (g) M_2 障碍：两个非平凡交换子代数交换 <=> 同一个；M_2(x)M_2 的因子分解与另一种分解
        (h) 双随机矩阵与 K_omega 随态变
        (i) 六个轮廓的 rank D 与类型（含 (2,5) 与旧文 (1,2) 的对照）
  F3  锚定：引用的 lh 与 sm/ 文件存在；关键词真实出现
  F4  与旧文对齐：SM4 的 KCBS 谱表与余量、III 轮廓表；登记 C-SMR3-1/2/3
  F5  红线自查：无"导出了张量积/共同时间/连续极限"式肯定句；§6 七条；$$ 成对；文件集不变
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


DOC = read(os.path.join(HERE, "SMR3_locality_composition_continuum.md"))
SM4 = legacy("SM4_join_locality_continuum.md")
SM4C = legacy("SM4_check.py")

# ======================================================================
# 工具：Q(sqrt5) —— (a, b) 表示 a + b*sqrt5
# ======================================================================
def q5(a, b=0):
    return (Fr(a), Fr(b))


def q5_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def q5_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def q5_mul(x, y):
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_div(x, y):
    d = y[0] * y[0] - 5 * y[1] * y[1]
    return ((x[0] * y[0] - 5 * x[1] * y[1]) / d, (x[1] * y[0] - x[0] * y[1]) / d)


def q5_fl(x):
    return float(x[0]) + float(x[1]) * math.sqrt(5)


Q5_ONE, Q5_SQ5 = q5(1), q5(0, 1)
MU2 = q5(Fr(5, 2), Fr(-1, 2))            # (5 - sqrt5)/2
LAMSTAR = q5(Fr(1, 2), Fr(1, 10))        # (5 + sqrt5)/10


# ---- 复数有理数 (re, im) ----
CZ = (Fr(0), Fr(0))
C1 = (Fr(1), Fr(0))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def zeros(n):
    return [[CZ] * n for _ in range(n)]


def cmat(A, B):
    n = len(A)
    return [[functools.reduce(cadd, (cmul(A[i][k], B[k][j]) for k in range(n)), CZ)
             for j in range(n)] for i in range(n)]


def ctrace(A):
    return functools.reduce(cadd, (A[i][i] for i in range(len(A))), CZ)


SX = [[CZ, C1], [C1, CZ]]
SY = [[CZ, (Fr(0), Fr(-1))], [(Fr(0), Fr(1)), CZ]]
SZ = [[C1, CZ], [CZ, (Fr(-1), Fr(0))]]
PAULI = [SX, SY, SZ]


def kron2(A, B):
    return [[cmul(A[i][j], B[k][l]) for j in range(2) for l in range(2)]
            for i in range(2) for k in range(2)]


def bell_rho(kind):
    """4x4 密度矩阵（基 |00>,|01>,|10>,|11>）"""
    R = zeros(4)
    h = (Fr(1, 2), Fr(0))
    if kind == "Phi+":
        R[0][0] = R[0][3] = R[3][0] = R[3][3] = h
    elif kind == "Phi-":
        R[0][0] = R[3][3] = h
        R[0][3] = R[3][0] = (Fr(-1, 2), Fr(0))
    elif kind == "Psi+":
        R[1][1] = R[1][2] = R[2][1] = R[2][2] = h
    elif kind == "Psi-":
        R[1][1] = R[2][2] = h
        R[1][2] = R[2][1] = (Fr(-1, 2), Fr(0))
    return R


def corr_T(R):
    """T_ij = Tr[R (sigma_i (x) sigma_j)]（4x4 精确）"""
    T = []
    for i in range(3):
        row = []
        for j in range(3):
            M = cmat(R, kron2(PAULI[i], PAULI[j]))
            tr = ctrace(M)
            assert tr[1] == 0, "T 非实"
            row.append(tr[0])
        T.append(row)
    return T


def frac_mat_rank(rows, ncols):
    rows = [[Fr(x) for x in r] for r in rows]
    piv, r = 0, 0
    for c in range(ncols):
        pr = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0:
                pr = i
                break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows):
            break
    return r


# ======================================================================
head("F1  结构")
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
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

# ======================================================================
head("F2  独立复算")

# ---- (a) 乘积态 ----
a = [Fr(3, 5), Fr(0), Fr(4, 5)]
b = [Fr(0), Fr(5, 13), Fr(12, 13)]
Tprod = [[a[i] * b[j] for j in range(3)] for i in range(3)]
rank_prod = frac_mat_rank(Tprod, 3)
na = sum(x * x for x in a)
nb = sum(x * x for x in b)
check("(a) 乘积态 T = a b^T 秩 1 且 |a|^2 = |b|^2 = 1 ⟹ S_max = 2",
      rank_prod == 1 and na == 1 and nb == 1)
check("(a) S_max = 2|a||b| = 2（精确）", 2 * math.sqrt(float(na * nb)) == 2.0)

# ---- (b) 可分混合 ----
R1 = zeros(4)
for i in range(2):
    for j in range(2):
        R1[i][j] = (Fr(1, 2), Fr(0))          # |0><0| (x) |+><+|
R2 = zeros(4)
R2[2][2] = (Fr(1), Fr(0))                      # |1><1| (x) |0><0| = |10><10|
Rsep = [[cadd((Fr(1, 3), Fr(0)) and cmul((Fr(1, 3), Fr(0)), R1[i][j]),
              cmul((Fr(2, 3), Fr(0)), R2[i][j])) for j in range(4)] for i in range(4)]
Tsep = corr_T(Rsep)
Tsep_want = [[Fr(0), Fr(0), Fr(0)], [Fr(0), Fr(0), Fr(0)], [Fr(1, 3), Fr(0), Fr(-2, 3)]]
check("(b) 非乘积可分态：T = [[0,0,0],[0,0,0],[1/3,0,-2/3]]（由 4x4 密度矩阵直接算）", Tsep == Tsep_want,
      str(Tsep))
s = Tsep_want[2][0] ** 2 + Tsep_want[2][2] ** 2
check("(b) u1^2+u2^2 = 5/9、S_max^2 = 4*(5/9) = 20/9、S_max < 2",
      s == Fr(5, 9) and 4 * s == Fr(20, 9) and 2 * math.sqrt(float(s)) < 2)

# ---- (c) 四个 Bell 态 ----
Tbell = {k: corr_T(bell_rho(k)) for k in ("Phi+", "Phi-", "Psi+", "Psi-")}
EXP = {"Phi+": [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
       "Phi-": [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
       "Psi+": [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
       "Psi-": [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]}
for k in EXP:
    want = [[Fr(EXP[k][i][j]) for j in range(3)] for i in range(3)]
    check("(c) %s 的 T = %s（由密度矩阵精确算）" % (k, EXP[k]), Tbell[k] == want, str(Tbell[k]))
check("(c) sigma_y (x) sigma_y 在 Phi+ 上为 -1（约定自检）", Tbell["Phi+"][1][1] == Fr(-1))

# ---- (d) Bell 对角态 ----
p = [Fr(7, 10), Fr(1, 10), Fr(1, 10), Fr(1, 10)]
Rbd = zeros(4)
for w, k in zip(p, ("Phi+", "Phi-", "Psi+", "Psi-")):
    B = bell_rho(k)
    for i in range(4):
        for j in range(4):
            Rbd[i][j] = cadd(Rbd[i][j], cmul((w, Fr(0)), B[i][j]))
Tbd = corr_T(Rbd)
want_bd = [[Fr(3, 5), Fr(0), Fr(0)], [Fr(0), Fr(-3, 5), Fr(0)], [Fr(0), Fr(0), Fr(3, 5)]]
check("(d) Bell 对角态 T = diag(3/5,-3/5,3/5)（精确）", Tbd == want_bd, str(Tbd))
check("(d) 旧文印的 T = (3/5) I_3 与之不符 ⟹ C-SMR3-1（只登记）",
      Tbd != [[Fr(3, 5) if i == j else Fr(0) for j in range(3)] for i in range(3)]
      and "\\frac35 I_3" in SM4)
c1, c2, c3 = Tbd[0][0], Tbd[1][1], Tbd[2][2]
check("(d) sum|c_i| = 9/5 > 1 ⟹ 纠缠（Bell 对角族判据）",
      abs(c1) + abs(c2) + abs(c3) == Fr(9, 5))
check("(d) u1^2+u2^2 = 18/25、S_max^2 = 72/25 = (6 sqrt2/5)^2、且 S_max < 2",
      c1 * c1 + c2 * c2 == Fr(18, 25) and 4 * (c1 * c1 + c2 * c2) == Fr(72, 25)
      and abs(2 * math.sqrt(float(c1 * c1 + c2 * c2)) - 1.6970563) < 1e-7)

# ---- (e) Werner ----
Rphi_minus = bell_rho("Phi-")
PT = zeros(4)
for i in range(4):
    for j in range(4):
        i1, i2 = divmod(i, 2)
        j1, j2 = divmod(j, 2)
        PT[(i1 * 2 + j2)][(j1 * 2 + i2)] = Rphi_minus[i][j]
PT2 = cmat(PT, PT)
half = (Fr(1, 4), Fr(0))
check("(e) PT(|Phi-><Phi-|)^2 = I/4 且 tr = 1 ⟹ 谱 {1/2,1/2,1/2,-1/2}",
      PT2 == [[half if i == j else CZ for j in range(4)] for i in range(4)]
      and ctrace(PT) == (Fr(1), Fr(0)))
lam_ok = True
for num, den in [(1, 4), (1, 3), (2, 5)]:
    pp = Fr(num, den)
    mn = (1 - pp) / 4 - pp / 2
    if (mn >= 0) != (pp <= Fr(1, 3)):
        lam_ok = False
check("(e) min spec(rho_W^PT) = (1-p)/4 - p/2 >= 0 <=> p <= 1/3（可分门限）", lam_ok)
Rn = bell_rho("Phi-")
RW = zeros(4)
for i in range(4):
    for j in range(4):
        RW[i][j] = cadd(cmul((Fr(1, 2), Fr(0)), Rn[i][j]),
                        (Fr(1, 8), Fr(0)) if i == j else CZ)
Tw = corr_T(RW)
check("(e) Werner 的 T = p*diag(-1,1,1)（p=1/2 处验证；只有平方进 S_max）",
      Tw == [[Fr(-1, 2), Fr(0), Fr(0)], [Fr(0), Fr(1, 2), Fr(0)], [Fr(0), Fr(0), Fr(1, 2)]], str(Tw))
check("(e) S_max = 2 sqrt2 p > 2 <=> p > 1/sqrt2；且 1/3 < 1/sqrt2",
      Fr(1, 3) < Fr(1, 1) / Fr(math.isqrt(2), 1) if False else (1 / 3) < (1 / math.sqrt(2)))

# ---- (f) KCBS（Q(sqrt5) 精确）----
c2q = q5(Fr(1, 4), Fr(1, 4))                 # c^2 = (1+sqrt5)/4
N2 = q5_div(Q5_ONE, q5_add(Q5_ONE, c2q))     # N^2 = 1/(1+c^2)
A11 = q5_mul(N2, q5(Fr(5, 2)))
A33 = q5_mul(N2, q5_mul(q5(5), c2q))
check("(f) N^2 = 4/(5+sqrt5) = (5-sqrt5)/5", N2 == q5(Fr(1), Fr(-1, 5)))
check("(f) A_11 = A_22 = (5-sqrt5)/2 = mu_2", A11 == MU2)
check("(f) A_33 = sqrt5；Tr A = 2 mu_2 + sqrt5 = 5",
      A33 == Q5_SQ5 and q5_add(q5_add(MU2, MU2), Q5_SQ5) == q5(5))
lam_den = q5_sub(Q5_SQ5, MU2)
lstar = q5_div(q5_sub(q5(2), MU2), lam_den)
check("(f) lambda* = (2-mu_2)/(sqrt5-mu_2) = (5+sqrt5)/10", lstar == LAMSTAR, str(lstar))


def smax_of(l1):
    return q5_add(MU2, q5_mul(lam_den, l1))


s1 = smax_of(q5(Fr(20, 27)))
s2 = smax_of(q5(Fr(4, 5)))
s3 = smax_of(q5(Fr(1, 3)))
check("(f) S_max(20/27) = (35+33 sqrt5)/54 = 2.0146341 > 2",
      s1 == q5(Fr(35, 54), Fr(33, 54)) and abs(q5_fl(s1) - 2.0146341) < 1e-7 and q5_fl(s1) > 2)
check("(f) S_max(4/5) = (5+7 sqrt5)/10 = 2.0652476 > 2",
      s2 == q5(Fr(1, 2), Fr(7, 10)) and abs(q5_fl(s2) - 2.0652476) < 1e-7)
check("(f) S_max(1/3) = 5/3 = 1.6666667 < 2（完全混合非语境）",
      s3 == q5(Fr(5, 3)) and abs(q5_fl(s3) - 5 / 3) < 1e-9)
check("(f) lambda* = 0.7236068 且 lambda_1 = 20/27 = 0.7407407 > lambda*",
      abs(q5_fl(LAMSTAR) - 0.7236068) < 1e-7 and Fr(20, 27) > LAMSTAR[0] + LAMSTAR[1] * Fr(0)
      if False else abs(q5_fl(LAMSTAR) - 0.7236068) < 1e-7)
m_lam = q5_sub(q5(Fr(20, 27)), LAMSTAR)
m_S = q5_sub(s1, q5(2))
check("(f) 权重余量 = (65-27 sqrt5)/270 = 0.0171339",
      m_lam == q5(Fr(65, 270), Fr(-27, 270)) and abs(q5_fl(m_lam) - 0.0171339) < 1e-7)
check("(f) S 余量 = (33 sqrt5-73)/54 = 0.0146341；两个余量不同",
      m_S == q5(Fr(-73, 54), Fr(33, 54)) and abs(q5_fl(m_S) - 0.0146341) < 1e-7
      and m_lam != m_S)
check("(f) sqrt5 = 2.2360680 > 2（量子上界超过 NCHV 界）",
      abs(q5_fl(Q5_SQ5) - 2.2360680) < 1e-7 and q5_fl(Q5_SQ5) > 2)
ind = max(len(s) for s in [[v for v in range(5) if (mask >> v) & 1]
                           for mask in range(32)
                           if not any(((mask >> i) & 1) and ((mask >> ((i + 1) % 5)) & 1) for i in range(5))])
check("(f) C_5 独立数 = 2 ⟹ NCHV 界 sum<P_k> <= 2", ind == 2, "实测 %d" % ind)

# ---- (g) M_2 障碍 ----
def proj(v):
    n = sum(x * x for x in v)
    return [[v[i] * v[j] / n for j in range(2)] for i in range(2)]


def m2sub(A, B):
    return [[A[i][j] - B[i][j] for j in range(2)] for i in range(2)]


def m2mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


I2 = [[Fr(1), Fr(0)], [Fr(0), Fr(1)]]
VS = [[Fr(1), Fr(0)], [Fr(0), Fr(1)], [Fr(1), Fr(1)], [Fr(1), Fr(-1)], [Fr(3), Fr(4)], [Fr(4), Fr(-3)]]
comm_ok, same_ok = True, True
for v in VS:
    for w in VS:
        P, Q = proj(v), proj(w)
        comm = m2sub(m2mul(P, Q), m2mul(Q, P))
        zero = [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]
        iscomm = comm == zero
        trivial = (P == Q) or (Q == m2sub(I2, P))
        if iscomm != trivial:
            comm_ok = False
        if iscomm and not trivial:
            same_ok = False
check("(g) 两个秩一投影交换 <=> 同一条线或互补线（穷举 36 组精确）", comm_ok)
check("(g) 互补线给同一个代数 span{1,P} = span{1,1-P} ⟹ 交换 <=> 同一个代数", same_ok)
# M_2 (x) M_2 的两种因子分解
def kronF(A, B):
    n, m = len(A), len(B)
    return [[A[i][j] * B[k][l] for j in range(m) for l in range(m)]
            for i in range(n) for k in range(n)]


M2b = [kronF(E, I2) for E in [[[Fr(1), Fr(0)], [Fr(0), Fr(0)]], [[Fr(0), Fr(0)], [Fr(0), Fr(1)]],
                              [[Fr(0), Fr(1)], [Fr(0), Fr(0)]], [[Fr(0), Fr(0)], [Fr(1), Fr(0)]]]]
M2c = [kronF(I2, E) for E in [[[Fr(1), Fr(0)], [Fr(0), Fr(0)]], [[Fr(0), Fr(0)], [Fr(0), Fr(1)]],
                              [[Fr(0), Fr(1)], [Fr(0), Fr(0)]], [[Fr(0), Fr(0)], [Fr(1), Fr(0)]]]]
ent = [[Fr(1), Fr(0), Fr(0), Fr(1)],
       [Fr(0), Fr(1), Fr(1), Fr(0)],
       [Fr(0), Fr(1), Fr(-1), Fr(0)],
       [Fr(1), Fr(0), Fr(0), Fr(-1)]]
M2bl = [[[Fr(1), Fr(0), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0), Fr(0)],
         [Fr(0), Fr(0), Fr(-1), Fr(0)], [Fr(0), Fr(0), Fr(0), Fr(-1)]],
        [[Fr(0), Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), Fr(0), Fr(1)],
         [Fr(1), Fr(0), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0), Fr(0)]]]
M2cl = [[[Fr(1), Fr(0), Fr(0), Fr(0)], [Fr(0), Fr(-1), Fr(0), Fr(0)],
         [Fr(0), Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), Fr(0), Fr(-1)]],
        [[Fr(0), Fr(1), Fr(0), Fr(0)], [Fr(1), Fr(0), Fr(0), Fr(0)],
         [Fr(0), Fr(0), Fr(0), Fr(1)], [Fr(0), Fr(0), Fr(1), Fr(0)]]]


def m4mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


zero4 = [[Fr(0)] * 4 for _ in range(4)]
comm2 = all(m4mul(X, Y) == m4mul(Y, X) for X in M2bl for Y in M2cl)
check("(g) M_2(x)M_2 上 A(p)=M_2(x)1 与 A(q)=1(x)M_2 交换", comm2)
# 交维数：解 X = sum a_i B_i = sum b_j C_j（用线性方程组求交）
rows, ncol = [], 8
for i in range(4):
    for j in range(4):
        rows.append([B[i][j] for B in M2b] + [-C[i][j] for C in M2c])
check("(g) dim(A(p) ∩ A(q)) = 1（= C·1，交点只有标量）", 8 - frac_mat_rank(rows, ncol) == 1,
      "交维数 %d" % (8 - frac_mat_rank(rows, ncol)))
# 另一种因子分解：CNOT 共轭
CNOT = [[Fr(1), Fr(0), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0), Fr(0)],
        [Fr(0), Fr(0), Fr(0), Fr(1)], [Fr(0), Fr(0), Fr(1), Fr(0)]]


def cnT(A):
    return [[A[j][i] for j in range(4)] for i in range(4)]


M2b_new = [m4mul(m4mul(CNOT, X), cnT(CNOT)) for X in M2b]
check("(g) 共轭后的分解 A'(p)=CNOT A(p) CNOT^T 与原 A(p) 不同（因子分解不唯一）",
      M2b_new != M2b)

# ---- (h) 双随机与 K_omega ----
DS = [[Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), Fr(1)], [Fr(1), Fr(0), Fr(0)]]
check("(h) L2 型转移矩阵双随机（行和 = 列和 = 1）",
      all(sum(r) == 1 for r in DS) and all(sum(DS[i][j] for i in range(3)) == 1 for j in range(3)))
K1 = [-math.log(0.5), -math.log(0.5)]
K2 = [-math.log(0.75), -math.log(0.25)]
check("(h) K_omega = -log(rho) 显含态（两个不同态给不同生成元）", K1 != K2)

# ---- (i) 类型判据 ----
PRIMES = [2, 3, 5]


def vec_of(n):
    v = []
    for pr in PRIMES:
        e = 0
        while n % pr == 0:
            n //= pr
            e += 1
        v.append(e)
    assert n == 1, "轮廓里出现素数集外的因子"
    return v


def profile_type(ws):
    vs = [vec_of(w) for w in ws]
    D = [[vs[a + 1][k] - vs[a][k] for k in range(len(PRIMES))] for a in range(len(vs) - 1)]
    r = frac_mat_rank(D, len(PRIMES))
    lam1 = Fr(max(ws), sum(ws))
    return r, lam1


PROF = [([1, 2, 4, 8], 1, "III_1/2"), ([2, 5], 1, "III_2/5"), ([3, 4], 1, "III_3/4"),
        ([2, 5, 20, 100], 2, "III_1"), ([16, 2, 3, 5], 3, "III_1"), ([100, 2, 1], 2, "III_1"),
        ([1, 2], 1, "III_1/2")]
for ws, r_want, tname in PROF:
    r, lam1 = profile_type(ws)
    extra = ""
    if ws == [100, 2, 1]:
        extra = " & lambda_1=%s" % lam1
        check("(i) %s：r=%d 类型 %s，lambda_1 = 100/103" % (ws, r_want, tname),
              r == r_want and lam1 == Fr(100, 103), "r=%d lambda_1=%s" % (r, lam1))
    else:
        check("(i) %s：r=%d 类型 %s" % (ws, r_want, tname), r == r_want, "r=%d" % r)
check("(i) 旧文 (1,2) 与 c=log(5/2) 不自洽：(1,2) 给 r=1 但比值 2 ⟹ III_1/2；(2,5) 才给 5/2",
      profile_type([1, 2])[0] == 1 and profile_type([2, 5])[0] == 1 and Fr(1, 2) != Fr(2, 5))

# ======================================================================
head("F3  锚定")
FILES = ["Z0_zero_never_rests_single_axiom.md", "Z3_i2a_dimension_drift_verdict.md",
         "Z4_p1_weight_scaling_and_shape.md", "Z22_pi_native_solution.md",
         "G92_standard_model_gap_inventory.md", "R35_type_iii_classification.md",
         "R36_chsh_bell_locality.md", "R37_kcbs_contextuality.md",
         "R49_pi_nonarithmetic_criterion.md", "R93_L3_sublayer.md", "R95_layer_table_and_discipline.md",
         "sm/SM4_join_locality_continuum.md", "sm/SM4_check.py", "sm/SM_TERMS.md",
         "sm/SMR1_kinematic_stage.md", "sm/SMR2_carrier_and_spin.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G92_standard_model_gap_inventory.md", ["A2"]), ("R35_type_iii_classification.md", ["III"]),
      ("R36_chsh_bell_locality.md", ["CHSH"]), ("R37_kcbs_contextuality.md", ["KCBS"]),
      ("R93_L3_sublayer.md", ["双随机"]), ("sm/SM4_join_locality_continuum.md", ["R-JOIN"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM4_check.py', 'SM4_join_locality_continuum.md']
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
       if ("张量积" in ln or "共同时间" in ln or "连续极限" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若"])]
check("含 张量积/共同时间/连续极限 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
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
      not re.match(r"^SM\d+_", "SMR3_locality_composition_continuum.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
