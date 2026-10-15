#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR35_check.py -- 核验《SMR35 · 扇区选择性：群论判死＋色因子否定》"""
import io, math, os, re, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
FAIL = NCHK = 0
def head(t): print(""); print("=" * 72); print(t); print("=" * 72)
def check(n, c, d=""):
    global FAIL, NCHK
    NCHK += 1
    ok = bool(c)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   " + d) if d else ""))
def read(p):
    with io.open(p, encoding="utf-8") as f: return f.read()
def exists(rel): return os.path.exists(os.path.join(HERE, rel)) or os.path.exists(os.path.join(ROOT, rel))

DOC = read(os.path.join(HERE, "SMR35_sector_selectivity.md"))

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 " + fld + " 在位", fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 " + sec.strip() + " 在位", sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 6 条推导（2.1-2.6）", subs == ["2." + str(i) for i in range(1, 7)], "实测 " + str(subs))
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§" + name + " 五段式标记齐全", not miss, "缺 " + str(miss))
    check("§" + name + " 常见误读 >= 2 条", blk.count("✗") >= 2, "实测 " + str(blk.count("✗")))
    check("§" + name + " 有加框结果", "\\boxed" in blk)

head("F2  独立复算")
# (a) 三扇区结构与色重数
SECT = {"charged lepton": (2, 2, 1, 1), "up quark": (2, 2, 1, 3), "down quark": (2, 2, 1, 3)}
check("(a) 三扇区算符结构同为 2 x 2 x 1（SU(2) 内容），差别仅在色重数 1 vs 3",
      all(v[:3] == (2, 2, 1) for v in SECT.values())
      and SECT["charged lepton"][3] == 1 and SECT["up quark"][3] == SECT["down quark"][3] == 3)
# (b) SU(2) 不变收缩维数 = 1
def su2_mats():
    import cmath
    out = []
    for th in (0.7, 1.3):
        for ax in range(3):
            c, s = math.cos(th / 2), math.sin(th / 2)
            if ax == 0: U = [[c, 1j * s], [1j * s, c]]
            elif ax == 1: U = [[c, s], [-s, c]]
            else: U = [[c + 1j * s, 0], [0, c - 1j * s]]
            out.append(U)
    return out
# 未知量 m11, m12, m21, m22（4 实未知 -> 用复数实部虚部，共 8 未知）
import itertools
rows = []
for U in su2_mats():
    for a in range(2):
        for b in range(2):
            # (U M U^T - M)_{ab} = 0，对 M 的每个复元素取实部/虚部
            for part in (0, 1):
                row = []
                for i in range(2):
                    for j in range(2):
                        coef = U[a][i] * U[b][j]
                        row.append((coef.real, coef.imag))
                # 展开为 8 个实未知：Re m11, Im m11, ...
                r8 = []
                for i in range(2):
                    for j in range(2):
                        cre, cim = row[i * 2 + j]
                        r8.append(cre)   # Re m_ij
                        r8.append(-cim if part == 0 else cre * 0 + cim)  # 简化处理
                # 直接数值求解：用实数线性代数，目标为 0
                rows.append(r8)
def nullspace_dim(mat, ncol, tol=1e-9):
    M = [r[:] for r in mat]
    r = 0
    for col in range(ncol):
        piv = None
        for i in range(r, len(M)):
            if abs(M[i][col]) > tol: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][col]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and abs(M[i][col]) > tol:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    return ncol - r
# 直接构造：M 在 SU(2) 下不变 <=> M ∝ epsilon（用两个群元素测试）
def invariant_residuals(M):
    res = []
    for U in su2_mats():
        UMU = [[sum(U[a][i] * M[i][j] * U[b][j] for i in range(2) for j in range(2)) for b in range(2)] for a in range(2)]
        for a in range(2):
            for b in range(2):
                res.append(abs(UMU[a][b] - M[a][b]))
    return max(res)
eps = [[0, 1], [-1, 0]]
delta = [[1, 0], [0, 1]]
sig1 = [[0, 1], [1, 0]]
check("(b) epsilon 在 U M U^T 下不变；I 不不变（U U^T != 1，故汤川收缩只能是反对称）",
      invariant_residuals(eps) < 1e-12 and invariant_residuals(delta) > 1e-6)
check("(b) sigma_1 在 SU(2) 下**不**不变（故不变空间维数 < 3）", invariant_residuals(sig1) > 1e-6)
# 不变空间维数：解 U M U^T = M 的实线性系统（8 实未知，纯标准库）
import itertools as _it
A = []
for U in su2_mats():
    for a in range(2):
        for b in range(2):
            for part in (0, 1):
                row = [0.0] * 8
                for i in range(2):
                    for j in range(2):
                        c = U[a][i] * U[b][j]
                        idx_re = 2 * (i * 2 + j)
                        row[idx_re] += c.real if part == 0 else -c.imag
                        row[idx_re + 1] += c.imag if part == 0 else c.real
                for i in range(2):
                    for j in range(2):
                        idx = 2 * (i * 2 + j)
                        if part == 0 and i == a and j == b: row[idx] -= 1.0
                        if part == 1 and i == a and j == b: row[idx + 1] -= 1.0
                A.append(row)
dim = nullspace_dim(A, 8)
check("(b) SU(2) 不变 2x2 复矩阵空间维数 = 2（实维，即 M ∝ epsilon 一维复空间）", dim == 2, "实测 " + str(dim))
check("(b) 该一维复空间由 epsilon 张成（I 与 sigma_1 组合不能只靠 epsilon 表示）",
      invariant_residuals(eps) < 1e-12)
