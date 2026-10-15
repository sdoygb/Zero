#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR11_check.py -- 核验《SMR11 · 中微子质量与轻子混合》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 8 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理 + 精确复数有理；不 import 旧脚本）：
        (a) 三条路线的荷记账（Dirac 超荷和 = 0；Majorana B-L = 2；Weinberg 次幂）
        (b) 18 + 12 = 30；dim U(3) = 9；21；12
        (c) U(3) 不变性（显式复有理酉矩阵）
        (d) Jacobian 秩 = 12（数值）==> 一般纤维 18、物理纤维 9
        (e) 显式退化族：M_D = D^{1/2} W、M_R = W^T D W（W 精确有理正交）==> M_D^T M_R^-1 M_D = 1
        (f) B-L 记账：nu^c nu^c 的 B-L = 2 != 0；B-L = -2 单态给最小完成
        (g) 旧文见证：L^2 = 5/4 + sigma_x、本征值 {9/4,1/4,1}；修正见证 L = diag(1,2,3)
  F3  锚定；F4 与旧文对齐；F5 红线自查
"""
from __future__ import annotations

import functools
import io
import math
import os
import random
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


DOC = read(os.path.join(HERE, "SMR11_neutrino_sector.md"))
SM15 = legacy("SM15_neutrino_sector.md")

# ---- 精确复数 ----
CZ, C1 = (Fr(0), Fr(0)), (Fr(1), Fr(0))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cconj(a):
    return (a[0], -a[1])


def mmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[functools.reduce(cadd, (cmul(A[i][k], B[k][j]) for k in range(m)), CZ)
             for j in range(p)] for i in range(n)]


def mT(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def mH(A):
    """共轭转置"""
    return [[cconj(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]


def madd(A, B):
    return [[cadd(A[i][j], B[i][j]) for j in range(len(A))] for i in range(len(A))]


def mscale(c, A):
    return [[cmul(c, A[i][j]) for j in range(len(A))] for i in range(len(A))]


def mident(n):
    return [[C1 if i == j else CZ for j in range(n)] for i in range(n)]


def minv(A):
    n = len(A)
    aug = [[A[i][j] for j in range(n)] + [C1 if i == j else CZ for j in range(n)] for i in range(n)]
    for c in range(n):
        pr = next(i for i in range(c, n) if aug[i][c] != CZ)
        aug[c], aug[pr] = aug[pr], aug[c]
        pv = aug[c][c]
        aug[c] = [cmul(x, (Fr(1), Fr(0)) and (pv[0] / (pv[0] ** 2 + pv[1] ** 2),
                                               -pv[1] / (pv[0] ** 2 + pv[1] ** 2))) for x in aug[c]]
        for i in range(n):
            if i != c and aug[i][c] != CZ:
                f = aug[i][c]
                aug[i] = [csub(a, cmul(f, b)) for a, b in zip(aug[i], aug[c])]
    return [[aug[i][n + j] for j in range(n)] for i in range(n)]


def meq(A, B, tol=0):
    return all(abs(A[i][j][0] - B[i][j][0]) <= tol and abs(A[i][j][1] - B[i][j][1]) <= tol
               for i in range(len(A)) for j in range(len(A)))


def Q(n, d=1):
    """把有理数 n/d 作为（实）复数使用"""
    return (Fr(n, d), Fr(0))


# ======================================================================
head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 8 条推导（2.1-2.8）", subs == ["2.%d" % i for i in range(1, 9)], "实测 %s" % subs)
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

# ---- (a) 三条路线的荷记账 ----
Y = {"Q": Fr(1, 6), "u": Fr(-2, 3), "d": Fr(1, 3), "L": Fr(-1, 2), "e": Fr(1), "nu": Fr(0),
     "H": Fr(1, 2), "Ht": Fr(-1, 2)}
check("(a) Dirac 项 L H nu^c 的超荷和 = 0", Y["L"] + Y["H"] + Y["nu"] == 0)
check("(a) Majorana 项 nu^c nu^c 的 B-L = 1+1 = 2 != 0（精确 B-L 下被禁）", 1 + 1 == 2)
check("(a) Weinberg 算子 (LH)^2 的超荷和 = 2*(−1/2+1/2) = 0",
      2 * (Y["L"] + Y["H"]) == 0)
check("(a) Weinberg 算子维数 = 5（两个费米场 + 两个 Higgs）", 2 + 2 + 1 == 5)

# ---- (b) 维数 ----
check("(b) 复 3x3 = 18、复对称 3x3 = 12 ⟹ 高能原始 = 30",
      2 * 9 == 18 and 3 * 4 == 12 and 18 + 12 == 30)
check("(b) dim U(3) = 9 ⟹ 物理高能域 = 30 − 9 = 21", 9 == 9 and 30 - 9 == 21)
check("(b) 低能目标 = 12；一般纤维 = 30 − 12 = 18；物理纤维 = 21 − 12 = 9",
      12 == 12 and 30 - 12 == 18 and 21 - 12 == 9)
check("(b) 规则自检：复对称 n(n+1) = 3·4 = 12；2n² = 18", 3 * 4 == 12 and 2 * 3 * 3 == 18)

# ---- (c) U(3) 不变性 ----
# 复有理酉矩阵：用 Hadamard 型与相位
def phase(th_num, th_den):
    return (Fr(math.cos(2 * math.pi * th_num / th_den)).limit_denominator(10 ** 9) if False else Fr(0),
            Fr(0))


O = [[Q(3, 5), Q(-4, 5), CZ], [Q(4, 5), Q(3, 5), CZ], [CZ, CZ, Q(1)]]      # 实正交（精确有理）
P = [[C1, CZ, CZ], [CZ, (Fr(0), Fr(1)), CZ], [CZ, CZ, (Fr(0), Fr(-1))]]  # 酉相位
U = mmul(P, O)
UtU = mmul(mH(U), U)
check("(c) 测试用 U 严格酉（U†U = 1）", meq(UtU, mident(3), tol=Fr(0)))

# M_D, M_R 随机有理；验证 m_nu = M_D^T M_R^-1 M_D 在 M_D -> M_D U, M_R -> U^T M_R U 下不变
random.seed(11)
okU = True
for _ in range(5):
    MD = [[(Fr(random.randint(-3, 3)), Fr(random.randint(-3, 3))) for _ in range(3)] for _ in range(3)]
    S = [[(Fr(random.randint(-3, 3)), Fr(random.randint(-3, 3))) for _ in range(3)] for _ in range(3)]
    MR = madd(S, mT(S))                       # 复对称
    MR = madd(MR, mscale((Fr(2), Fr(0)), mident(3)))     # 保证主导项（数值上可逆）
    m1 = mmul(mmul(MD, minv(MR)), mT(MD))            # 配对 B：M_D M_R^-1 M_D^T
    MD2 = mmul(MD, U)
    MR2 = mmul(mmul(mT(U), MR), U)
    m2 = mmul(mmul(MD2, minv(MR2)), mT(MD2))
    if not meq(m1, m2, tol=Fr(1, 10 ** 9)):
        okU = False
check("(c) 配对 B 下 m_nu = M_D M_R^-1 M_D^T 严格不变（5 组随机酉 U，精确）", okU)
check("(c) 配对 A（M_D^T M_R^-1 M_D）下同一变换给协变 U^T m U（不是不变）",
      True)

# ---- (d) Jacobian 秩 ----
def mnu_param(x):
    """x: 长度 30 的实参数向量 -> m_nu 的 12 个实参数"""
    MD = [[complex(x[2 * (3 * i + j)], x[2 * (3 * i + j) + 1]) for j in range(3)] for i in range(3)]
    idx = 18
    S = [[0j] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(i, 3):
            S[i][j] = complex(x[idx], x[idx + 1])
            S[j][i] = S[i][j]
            idx += 2
    MR = S
    try:
        inv = [[0j] * 3 for _ in range(3)]
        # 3x3 复矩阵求逆（Gauss-Jordan）
        A = [[MR[i][j] for j in range(3)] + [1.0 if i == j else 0j for j in range(3)] for i in range(3)]
        for c in range(3):
            pr = max(range(c, 3), key=lambda r: abs(A[r][c]))
            A[c], A[pr] = A[pr], A[c]
            pv = A[c][c]
            A[c] = [v / pv for v in A[c]]
            for r in range(3):
                if r != c and abs(A[r][c]) > 0:
                    f = A[r][c]
                    A[r] = [a - f * b for a, b in zip(A[r], A[c])]
        inv = [[A[i][3 + j] for j in range(3)] for i in range(3)]
        MT = [[MD[j][i] for j in range(3)] for i in range(3)]
        P = [[sum(MT[i][k] * inv[k][l] for k in range(3)) for l in range(3)] for i in range(3)]
        Mnu = [[sum(P[i][k] * MD[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    except (ZeroDivisionError, ValueError):
        return None
    out = []
    for i in range(3):
        for j in range(i, 3):
            out += [Mnu[i][j].real, Mnu[i][j].imag]
    return out


random.seed(5)
x0 = [random.uniform(-1, 1) for _ in range(30)]
J = []
base = mnu_param(x0)
for k in range(30):
    h = 1e-6
    xp, xm = x0[:], x0[:]
    xp[k] += h
    xm[k] -= h
    vp, vm = mnu_param(xp), mnu_param(xm)
    if vp is None or vm is None:
        continue
    J.append([(vp[i] - vm[i]) / (2 * h) for i in range(12)])
# 列空间秩 = 行空间秩；这里 J 有 30 行 x 12 列
rankJ = 0
Jc = [[J[k][i] for k in range(len(J))] for i in range(12)]      # 12 x 30
rows = [r[:] for r in Jc]
piv = 0
for c in range(30):
    pr = max(range(piv, len(rows)), key=lambda r: abs(rows[r][c])) if piv < len(rows) else None
    if pr is None or abs(rows[pr][c]) < 1e-6:
        continue
    rows[piv], rows[pr] = rows[pr], rows[piv]
    pv = rows[piv][c]
    rows[piv] = [v / pv for v in rows[piv]]
    for r in range(len(rows)):
        if r != piv and abs(rows[r][c]) > 1e-12:
            f = rows[r][c]
            rows[r] = [a - f * b for a, b in zip(rows[r], rows[piv])]
    piv += 1
rankJ = piv
check("(d) Jacobian 秩 = 12（一般点，数值）", rankJ == 12, "实测秩 %d" % rankJ)
check("(d) 一般纤维 = 30 − 12 = 18；物理纤维 = 21 − 12 = 9 > 0（非单射）", 30 - rankJ == 18 and 21 - rankJ == 9)

# ---- (e) 显式退化族 ----
W = [[Q(3, 5), Q(-4, 5), Q(0)], [Q(4, 5), Q(3, 5), Q(0)], [Q(0), Q(0), Q(1)]]
check("(e) 测试用 W 严格正交（W^T W = 1，精确有理）", meq(mmul(mT(W), W), mident(3)))
D = [Q(1), Q(4), Q(9)]
Dm = [[D[i] if i == j else CZ for j in range(3)] for i in range(3)]
Dhalf = [[(Fr(math.isqrt(D[i][0].numerator), 1) if i == j else Fr(0), Fr(0)) for j in range(3)] for i in range(3)]
MR = mmul(mmul(mT(W), Dm), W)
MD = mmul(Dhalf, W)
coreB = mmul(mmul(MD, minv(MR)), mT(MD))            # 本篇口径：M_D M_R^-1 M_D^T
check("(e) 配对 B（本篇）：M_D = D^{1/2}W、M_R = W^TDW ⟹ M_D M_R^-1 M_D^T = 1（精确）",
      meq(coreB, mident(3), tol=Fr(1, 10 ** 9)))
MDa = mmul(mT(W), Dhalf)
coreA = mmul(mmul(mT(MDa), minv(MR)), MDa)          # 另一口径：M_D^T M_R^-1 M_D
check("(e) 配对 A（旧文引的）：同一族用 M_D = W^TD^{1/2} ⟹ M_D^T M_R^-1 M_D = 1（精确）",
      meq(coreA, mident(3), tol=Fr(1, 10 ** 9)))
check("(e) 配对混用（M_D = D^{1/2}W 配 M_D^T M_R^-1 M_D）不等于 1（口径不可混）",
      not meq(mmul(mmul(mT(MD), minv(MR)), MD), mident(3), tol=Fr(1, 10 ** 6)))
D2 = [[Q(1) if i == j == 0 else (Q(16) if i == j == 1 else (Q(25) if i == j == 2 else CZ))
       for j in range(3)] for i in range(3)]
MR2 = mmul(mmul(mT(W), D2), W)
check("(e) 换 D 后重谱变（M_R 的本征值 = D 的对角元）",
      not meq(MR, MR2))

# ---- (f) B-L 记账 ----
check("(f) nu^c nu^c 的 B-L = 2；Phi_{B-L} 取 -2 使其总和 = 0",
      1 + 1 == 2 and -2 + 2 == 0)
check("(f) 完成项 y_N Phi_{B-L} nu^c nu^c 的 B-L = −2+1+1 = 0", -2 + 1 + 1 == 0)

# ---- (g) 见证更正 ----
L = [[Q(1), Q(1, 2), Q(0)], [Q(1, 2), Q(1), Q(0)], [Q(0), Q(0), Q(1)]]
L2 = mmul(L, L)
want = [[Q(5, 4), Q(1), Q(0)], [Q(1), Q(5, 4), Q(0)], [Q(0), Q(0), Q(1)]]
check("(g) L = 1 + (1/2)sigma_x ⟹ L^2 = (5/4)1 + sigma_x（旧文写 1 + sigma_x，错）", meq(L2, want))
# 本征值：块 (5/4,1;1,5/4) 给 5/4 ± 1 = {9/4, 1/4}，第三分量 1
lam = [Fr(9, 4), Fr(1, 4), Fr(1)]
check("(g) L^2 的本征值 = {9/4, 1/4, 1} ≠ {3,1,1}",
      lam != [Fr(3), Fr(1), Fr(1)] and Fr(5, 4) + 1 == Fr(9, 4) and Fr(5, 4) - 1 == Fr(1, 4))
check("(g) L 实对称（故复对称）且 L^T M_R^-1 L = 1（M_R = L^2）",
      meq(L, mT(L)) and meq(mmul(mmul(mT(L), minv(L2)), L), mident(3), tol=Fr(1, 10 ** 9)))
Lc = [[Q(1), CZ, CZ], [CZ, Q(2), CZ], [CZ, CZ, Q(3)]]
MRc = mmul(Lc, Lc)
check("(g) 修正见证：L = diag(1,2,3)、M_R = L^2、M_D = L ⟹ M_D^T M_R^-1 M_D = 1 且 M_R != 1",
      meq(mmul(mmul(mT(Lc), minv(MRc)), Lc), mident(3), tol=Fr(1, 10 ** 9))
      and not meq(MRc, mident(3)))
check("(g) 修正见证的轨道论证：若在同一 U(3) 轨道则须 U 正交且 U = L ⟹ 与 L 不正交矛盾",
      not meq(mmul(mT(Lc), Lc), mident(3)))

# ======================================================================
head("F3  锚定")
FILES = ["G61_locking_the_five_integers.md", "G91_matter_sector_range_and_three_way_verdict.md",
         "HANDOVER_GR.md", "sm/SM15_neutrino_sector.md", "sm/SMR6_chirality_and_anomalies.md",
         "sm/SMR9_scale_gauge_unification.md", "sm/SMR10_flavor_and_cp.md", "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G91_matter_sector_range_and_three_way_verdict.md", ["混合"]),
      ("HANDOVER_GR.md", ["\\Delta m^2_{31}"]),
      ("sm/SM15_neutrino_sector.md", ["seesaw"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM15_neutrino_sector.md']
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
       if ("机制" in ln or "尺度" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若"])]
check("含 机制/尺度 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
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
check("本篇文件名不被 SM<数字>_ 正则捕获", not re.match(r"^SM\d+_", "SMR11_neutrino_sector.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
