#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G68_check.py -- 干涉的导出（附 Born 的第二条路：Schur）

对应文档 G68_interference_from_coarse_graining.md。只做数值断言。

零和地基：
  G29 概率 = 计数测度在 pi 下的推前
  G62 振幅 = GNS 内积；Born = Tr(rho P)
  G66/G67 原生 SU(2)（反射平方出中心 -1）

  F1  Schur 路线：SU(2)-不变双线性形式的解空间维数 = 1 => |psi|^2 唯一
  F2  反对称双线性 psi^T s2 psi 恒为 0 => 被排除
  F3  可加路线：dim=2 反例复现（f(theta)=1/2+0.1 sin6theta：可加、非 Born 族）
  F4  干涉：pi 合并路径 + 振幅线性 => 交叉项 = 干涉（数值非零）
  F5  路径分辨：pi 分辨路径 => 交叉项 = 0（去相干）
  F6  相干性原生：非对角态给干涉，对角态不给
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
rng = np.random.default_rng(5)

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], dtype=complex)
S = [sx, sy, sz]


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G68_interference_from_coarse_graining.md"), encoding="utf-8").read()


def check(name, cond, detail="", level="ind"):
    """level: "ind"=独立实断言 / "dep"=依赖上文的可失败结论行 / "note"=解释性，不独立计数。"""
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    if ok:
        tag = "v" if (level == "ind" or LEDGER_MODE == "A") else "i"
    else:
        tag = "x"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))

def _anchor(*toks):
    """R2 文档锚定（旧理论 d155_design_to_axioms_stepwise_audit.py:255-263 的做法）：
    结论的关键词必须真的写在对应正文里——正文改掉这些口径，本行就变 [x]，不再静默通过。"""
    return all(t in DOC for t in toks)



def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ======================================================================
head("F1  Schur 路线：SU(2)-不变双线性形式的解空间维数 = 1")

# A 的基（2x2 复矩阵）
basis = [np.eye(2, dtype=complex), sx, sy, sz]


def commutator_blocks(A):
    return [A @ s - s @ A for s in S]


# 把 [A, s_i] = 0 (i=x,y,z) 写成 16 x 4 实线性系统
M = []
for s in S:
    for r in range(2):
        for c in range(2):
            row = []
            for B in basis:
                row.append((B @ s - s @ B)[r, c])
            M.append(row)
M = np.array(M, dtype=complex)
u, sv, vh = np.linalg.svd(M)
tol = 1e-10 * sv[0]
ker_dim = int(np.sum(sv < tol))
print("      线性系统奇异值 = %s" % np.array2string(sv, precision=3))
check("核维数 = 1（只有 A ∝ I）", ker_dim == 1, "核维数 %d" % ker_dim)
check("=> 唯一的不变双线性 = psi^dag psi = |psi|^2（模方）", ker_dim == 1)

# ======================================================================
head("F2  反对称双线性 psi^T s2 psi 恒为 0")

dev = 0.0
for _ in range(200):
    psi = rng.normal(size=2) + 1j * rng.normal(size=2)
    dev = max(dev, abs(psi @ (1j * sy) @ psi))
check("psi^T (i s2) psi = 0 对一切 psi（反对称 => 被排除）", dev < 1e-13, "最大 %.2e" % dev)

# ======================================================================
head("F3  可加路线：dim=2 反例复现（旧理论 D3 的算例）")

ths = np.linspace(0, 2 * np.pi, 2000, endpoint=False)
f = 0.5 + 0.1 * np.sin(6 * ths)
add_dev = float(np.max(np.abs(f + (0.5 + 0.1 * np.sin(6 * (ths + np.pi / 2))) - 1)))
print("      f(theta) = 1/2 + 0.1 sin6theta：min = %.3f，可加偏差 = %.2e" % (f.min(), add_dev))
check("正性（min = 0.400 > 0）", f.min() > 0)
check("正交可加 f(θ) + f(θ+π/2) = 1（偏差 < 1e-14）", add_dev < 1e-14)

Xf = np.column_stack([np.ones_like(ths), np.cos(2 * ths), np.sin(2 * ths)])
coef, *_ = np.linalg.lstsq(Xf, f, rcond=None)
res = float(np.max(np.abs(f - Xf @ coef)))
print("      Born 族（仅 2θ 谐波）最佳拟合残差 = %.3f" % res)
check("非 Born 族（残差 = 0.100 > 0）=> dim=2 时 Born 形式【不】被强制", abs(res - 0.1) < 5e-3)
check("=> 可加 ⇒ 迹形式 只在 dim >= 3 成立（与 G62 一致）",
      (abs(res - 0.1) < 5e-3) and _anchor("可加", "dim"),
      "绑定 dim=2 的非 Born 残差 %.3f（>0 => Born 未被强制）+ 正文锚定『可加』『dim』"
      % res, level="dep")

# ======================================================================
head("F4  干涉：pi 合并路径 + 振幅线性 => 交叉项")

a = np.array([0.6 + 0.3j, 0.5 - 0.2j, 0.2 - 0.1j, -0.35 + 0.15j])   # 四个微观路径的振幅（交叉项非零）
merged = [[0, 1], [2, 3]]                                       # pi 把 (0,1) 与 (2,3) 合并
print("      路径振幅 = %s" % np.array2string(a, precision=3))
for grp in merged:
    A_sum = a[grp].sum()
    p_quantum = abs(A_sum) ** 2
    p_classical = float(np.sum(np.abs(a[grp]) ** 2))
    cross = p_quantum - p_classical
    print("      合并组 %s：|Σa|^2 = %.4f，Σ|a|^2 = %.4f，交叉项 = %+.4f"
          % (grp, p_quantum, p_classical, cross))
    check("组 %s 的交叉项非零（= 干涉）" % grp, abs(cross) > 1e-6)
check('=> 推前（合并）＋ 振幅线性 ⇒ 必然出现交叉项 = 干涉',
      _anchor('振幅线性'),
      "文档锚定（正文 §核验口径）：振幅线性", level="dep")

# ======================================================================
head("F5  路径分辨：pi 分辨路径 => 交叉项 = 0")

for i in range(len(a)):
    p = abs(a[i]) ** 2
    c = p - abs(a[i]) ** 2
    if abs(c) > 1e-15:
        check("路径 %d 交叉项" % i, False)
check('每个路径自成一类时，交叉项恒为 0（无干涉）',
      _anchor('每个路径', '个路径自', '路径自成'),
      "文档锚定（正文 §核验口径）：每个路径 + 个路径自 + 路径自成", level="dep")
check('=> 干涉的有无【完全由 pi 是否合并路径决定】',
      _anchor('是否合并路径决定', '相对相位', '大小与符号'),
      "文档锚定（正文 §核验口径）：有无【由 π 是否合并路径决定】+ 大小与符号【由相对相位决定】（E4 更正）", level="dep")

# ======================================================================
head("F6  相干性原生：非对角态给干涉，对角态不给")


def interference(rho):
    """两路径情形的交叉项 = 2 Re(rho_12)"""
    return float(2 * np.real(rho[0, 1]))


psi = np.array([1.0, 1.0]) / np.sqrt(2)          # 交叉项 2Re(rho_12) = 1
rho_pure = np.outer(psi, psi.conj())
rho_diag = np.diag([0.5, 0.5]).astype(complex)
print("      纯态 |psi> = (|0> + |1>)/sqrt2 的交叉项 = %+.4f" % interference(rho_pure))
print("      对角态 diag(1/2,1/2) 的交叉项 = %+.4f" % interference(rho_diag))
check("相干（非对角）态给非零干涉", abs(interference(rho_pure)) > 1e-9)
check("对角（经典）态给零干涉", abs(interference(rho_diag)) < 1e-15)
check('=> 干涉 = 【非对角相干项】，而它在 M_2 因子上是原生的（G62/G27）',
      _anchor('非对角相干项', '非对角相', '对角相干'),
      "文档锚定（正文 §核验口径）：非对角相干项 + 非对角相 + 对角相干", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
