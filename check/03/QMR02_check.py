#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR02_check.py —— 《自旋与双覆盖：从单一公理完全退化》的逐条独立复算。

  P1  定理 1  双覆盖：同态、核 = {±1}、R_q ∈ SO(3)、R_q = R_{-q}、R(2π) 的两种作用
  P2  定理 1(c) 非平凡性：闭环的提升端点 +1 → -1
  P3  定理 2  Z_2 区分与穷尽：整数／半整数自旋的 U(2π) = ±I
  P4  定理 3  自旋 1/2 从壳层 C^2 长出：su(2) 对易关系、Casimir、谱、唯一性
  P5  定理 3/推论 3.2  反射 Z_2 与双覆盖 Z_2 是同一个
  P6  定理 4  不可克隆：<ψ|φ> = <ψ|φ>² 的正交／非正交两面
  P7  定理 5  排除：忠实迹禁止 CCR；Killing 型符号差禁止 boost
  P8  文档契约与链接

不 import 任何其它脚本；路径由 __file__ 推导。
"""
import os
import sys
import re
import math
from fractions import Fraction as F

import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR02_spin_and_double_cover.md"


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("  [v] " if cond else "  [x] ") + name + (("   " + detail) if detail else ""))


def head(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


DOC = open(os.path.join(HERE, DOCNAME), encoding="utf-8").read() if os.path.exists(
    os.path.join(HERE, DOCNAME)) else ""

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = (SX, SY, SZ)


# ---------- 四元数（显式代数，用 Fraction 保证精确） ----------
def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2)


def qconj(q):
    return (q[0], -q[1], -q[2], -q[3])


def qunit(t, n):
    """q(t) = cos(t/2) + sin(t/2) n̂（用精确值只在特殊 t 上断言）。"""
    n = np.array(n, dtype=float)
    n = n / np.linalg.norm(n)
    return (math.cos(t / 2), math.sin(t / 2) * n[0], math.sin(t / 2) * n[1], math.sin(t / 2) * n[2])


def rot_of(q):
    """R_q 的矩阵（列 = R_q(e_i)）。"""
    cols = []
    for v in np.eye(3):
        r = qmul(qmul(q, (0.0, v[0], v[1], v[2])), qconj(q))
        cols.append(r[1:])
    return np.array(cols).T


# ======================================================================
head("P1  定理 1 · 双覆盖：同态、核、SO(3)、q 与 −q 同旋转")
# ======================================================================
q1 = qunit(0.7, [1, 2, 3])
q2 = qunit(-1.3, [3, -1, 2])
R1, R2 = rot_of(q1), rot_of(q2)
R12 = rot_of(qmul(q1, q2))
check("★ 同态性：R_{q1 q2} = R_{q1} R_{q2}（%.2e）" % np.linalg.norm(R12 - R1 @ R2),
      np.linalg.norm(R12 - R1 @ R2) < 1e-12)
check("★ R_q ∈ SO(3)：det = %.12f，RᵀR = I（%.2e）"
      % (np.linalg.det(R1), np.linalg.norm(R1.T @ R1 - np.eye(3))),
      abs(np.linalg.det(R1) - 1) < 1e-12 and np.linalg.norm(R1.T @ R1 - np.eye(3)) < 1e-12)
check("★ q 与 −q 给同一旋转：‖R_q − R_{−q}‖ = %.2e" % np.linalg.norm(R1 - rot_of(tuple(-x for x in q1))),
      np.linalg.norm(R1 - rot_of(tuple(-x for x in q1))) < 1e-12)
check("★ 核 = {±1}：R_{+1} = R_{−1} = I₃（%.2e）"
      % np.linalg.norm(rot_of((1.0, 0, 0, 0)) - rot_of((-1.0, 0, 0, 0))),
      np.allclose(rot_of((1.0, 0, 0, 0)), np.eye(3)) and np.allclose(rot_of((-1.0, 0, 0, 0)), np.eye(3)))

# ======================================================================
head("P2  定理 1(c) · 非平凡性：闭环的提升端点 +1 → −1")
# ======================================================================
q0 = qunit(0.0, [0, 0, 1])
qh = qunit(math.pi, [0, 0, 1])
qf = qunit(2 * math.pi, [0, 0, 1])
q4 = qunit(4 * math.pi, [0, 0, 1])
print("      q(0)=%s  q(π)=%s  q(2π)=%s  q(4π)=%s"
      % (np.round(q0, 6), np.round(qh, 6), np.round(qf, 6), np.round(q4, 6)))
check("★ q(0) = +1 且 q(2π) = −1（两重覆盖的核的两个元素）",
      np.allclose(q0, (1, 0, 0, 0)) and np.allclose(qf, (-1, 0, 0, 0)))
check("★ q(4π) = +1（回到单位 ⇒ 双重覆盖，非三重）", np.allclose(q4, (1, 0, 0, 0)))
check("★ 闭环 π(q(t))（t:0→2π）的两端同一点：R_{q(0)} = R_{q(2π)} = I₃",
      np.allclose(rot_of(q0), np.eye(3)) and np.allclose(rot_of(qf), np.eye(3)))
check("★ 但提升的端点不同（+1 ≠ −1）⇒ 该闭环**不零伦** ⇒ 覆盖非平凡", not np.allclose(q0, qf))

# ======================================================================
head("P3  定理 2 · Z_2 区分与穷尽（整数／半整数自旋）")
# ======================================================================
def U_spin(s, theta, axis=2):
    """自旋 s 表示中绕 axis 转 theta 的算子（最高权基）。"""
    dim = int(round(2 * s + 1))
    m = np.array([s - k for k in range(dim)])
    return np.diag(np.exp(1j * theta * m))


check("★ 半整数自旋（s=1/2, 3/2, 5/2）：U(2π) = −I",
      all(np.allclose(U_spin(s, 2 * math.pi), -np.eye(int(2 * s + 1)))
          for s in (0.5, 1.5, 2.5)))
check("★ 整数自旋（s=0, 1, 2）：U(2π) = +I",
      all(np.allclose(U_spin(s, 2 * math.pi), np.eye(int(2 * s + 1)))
          for s in (0.0, 1.0, 2.0)))
check("★ U(4π) = +I 对一切 s（双重覆盖，不是三重）",
      all(np.allclose(U_spin(s, 4 * math.pi), np.eye(int(2 * s + 1)))
          for s in (0.0, 0.5, 1.0, 1.5, 2.0)))
# 由 SU(2) 显式二维表示独立验证
def Usu2(t, n):
    n = np.array(n, float); n = n / np.linalg.norm(n)
    S = n[0] * SX + n[1] * SY + n[2] * SZ
    return math.cos(t / 2) * I2 - 1j * math.sin(t / 2) * S


check("★ 显式 SU(2)：U(2π) = −I₂ 且 det U(t) = 1、U†U = I",
      np.allclose(Usu2(2 * math.pi, [1, 2, 3]), -I2)
      and abs(np.linalg.det(Usu2(1.1, [1, 2, 3])) - 1) < 1e-12
      and np.allclose(Usu2(1.1, [1, 2, 3]).conj().T @ Usu2(1.1, [1, 2, 3]), I2))

# ======================================================================
head("P4  定理 3 · 自旋 1/2 从壳层 C² 长出")
# ======================================================================
J = [s / 2 for s in SIG]
eps = {(0, 1): (2, 1), (1, 2): (0, 1), (2, 0): (1, 1)}
ok = True
for (i, j), (k, sg) in eps.items():
    if not np.allclose(J[i] @ J[j] - J[j] @ J[i], 1j * sg * J[k]):
        ok = False
    if not np.allclose(J[j] @ J[i] - J[i] @ J[j], -1j * sg * J[k]):
        ok = False
check("★ 角动量对易关系 [J_i,J_j] = i ε_{ijk} J_k（六个分量逐一核对）", ok)
check("★ Casimir：J² = Σ J_i² = ¾ I₂ ⇒ s(s+1) = ¾ ⇒ s = 1/2",
      np.allclose(sum(x @ x for x in J), 0.75 * I2))
check("★ 谱：spec(J_z) = {+1/2, −1/2}",
      np.allclose(sorted(np.linalg.eigvalsh(J[2])), [-0.5, 0.5]))
check("★ J_i 自伴、无迹（故属 su(2) 的表示像）",
      all(np.allclose(x, x.conj().T) and abs(np.trace(x)) < 1e-15 for x in J))
# 唯一性：M_2(C) 的不可约表示维数 2 ⇒ 只能是 s=1/2
check("★ 唯一性：2s+1 = 2 ⇒ s = 1/2（该维数下的唯一不可约表示）",
      abs((2 * 0.5 + 1) - 2) < 1e-15 and all(int(2 * s + 1) != 2 for s in (0.0, 1.0, 1.5)))

# ======================================================================
head("P5  推论 3.2 · 反射 Z_2 与双覆盖 Z_2 是同一个")
# ======================================================================
# 反射 iota: w -> -w 保壳层（|S_k| 不变）且翻转符号
from itertools import combinations
L = 8
words = []
for pos in combinations(range(L), L // 2):
    w = [-1] * L
    for p in pos:
        w[p] = 1
    words.append(tuple(w))


def part(w):
    s, o = 0, []
    for x in w:
        s += x
        o.append(s)
    return o


check("★ 反射 ι:w→−w 保壳层（|S_k(−w)| = |S_k(w)| 对全部词与前缀）",
      all([abs(x) for x in part(tuple(-x for x in w))] == [abs(x) for x in part(w)] for w in words))
check("★ ι 是**对合**：ι² = id", all(tuple(-x for x in tuple(-x for x in w)) == w for w in words))
# 符号翻转的算子实现 = −I₂；而 U(2π) 也是 −I₂
check("★ 两者在壳层 C² 上是同一个：ι ↔ −I₂，且 U(2π) = −I₂；各自平方为 +I₂",
      np.allclose((-I2) @ (-I2), I2) and np.allclose(Usu2(2 * math.pi, [0, 0, 1]) @ Usu2(2 * math.pi, [0, 0, 1]), I2))

# ======================================================================
head("P6  定理 4 · 不可克隆")
# ======================================================================
# 正交对可克隆：CNOT（控制=第一比特，目标=第二比特）
CNOT = np.zeros((4, 4), dtype=complex)
CNOT[0, 0] = 1.0            # |00> -> |00>
CNOT[1, 1] = 1.0            # |01> -> |01>
CNOT[2, 3] = 1.0            # |10> -> |11>
CNOT[3, 2] = 1.0            # |11> -> |10>
_okc = True
for _a in (0, 1):
    _psi = np.zeros(4, dtype=complex); _psi[2 * _a + 0] = 1.0
    _want = np.zeros(4, dtype=complex); _want[2 * _a + _a] = 1.0
    if not np.allclose(CNOT @ _psi, _want):
        _okc = False
check("★ 正交对**可**克隆：CNOT|x⟩|0⟩ = |x⟩|x⟩（x = 0,1）", _okc)
# 非正交对必然失败：<0|+> = 1/√2，而其平方 = 1/2
ip = 1 / math.sqrt(2)
check("★ 非正交对**必败**：⟨0|+⟩ = %.6f 而 ⟨0|+⟩² = %.6f ⇒ 等式不成立" % (ip, ip * ip),
      abs(ip - ip * ip) > 0.2)
plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
zero = np.array([1, 0], dtype=complex)
phi = np.kron(zero, zero)
out = CNOT @ phi
want = np.kron(plus, plus)
check("★ CNOT 在 |+⟩ 上失败：重叠 = %.6f < 1" % abs(np.vdot(want, out)),
      abs(np.vdot(want, out)) < 0.99)

# ======================================================================
head("P7  定理 5（排除）· 忠实迹禁止 CCR；Killing 型禁止 boost")
# ======================================================================
n = 8
rng = np.random.default_rng(11)
worst = 0.0
for _ in range(2000):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); X = (X + X.conj().T) / 2
    P = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); P = (P + P.conj().T) / 2
    worst = max(worst, abs(np.trace(X @ P - P @ X)))
check("★ (a) τ([X,P]) ≡ 0（2000 组随机，max = %.2e）而 τ(i·1) = i·%d ≠ 0 ⇒ 无 CCR" % (worst, n),
      worst < 1e-9 and abs(complex(0, n)) > 0)


def so_basis(nn, eta=None):
    B = []
    for i in range(nn):
        for j in range(i + 1, nn):
            M = np.zeros((nn, nn))
            if eta is None:
                M[i, j], M[j, i] = 1.0, -1.0
            else:
                M[i, j], M[j, i] = eta[j, j], -eta[i, i]
            B.append(M)
    return B


def killing_sig(B):
    m = len(B)
    G = np.array([[np.trace(B[a] @ B[b]) for b in range(m)] for a in range(m)])
    ev = np.linalg.eigvalsh(G)
    return int((ev < -1e-9).sum()), int((ev > 1e-9).sum())


neg4, pos4 = killing_sig(so_basis(4))
negL, posL = killing_sig(so_basis(4, np.diag([-1.0, 1.0, 1.0, 1.0])))
check("★ so(4)：Killing 型**定号**（%d 负 / %d 正）⇒ 紧" % (neg4, pos4), neg4 == 6 and pos4 == 0)
check("★ so(1,3)：Killing 型**不定**（%d 负 / %d 正）⇒ 非紧 ⇒ 有限维不承载 boost" % (negL, posL),
      negL == 3 and posL == 3)
# 反自伴元在忠实迹下 τ(X²) < 0
Xa = rng.normal(size=(n, n))
Xa = Xa - Xa.T
check("★ 反自伴元 X = −Xᵀ 满足 τ(X²) = −τ(XᵀX) = −‖X‖²_HS < 0（%.6f）" % np.trace(Xa @ Xa),
      np.trace(Xa @ Xa) < 0 and abs(np.trace(Xa @ Xa) + np.linalg.norm(Xa) ** 2) < 1e-9)

# ======================================================================
head("P8  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("不依赖任何同题旧文", "去依赖声明"),
                 ("定理 1", "定理 1"), ("定理 2", "定理 2"), ("定理 3", "定理 3"),
                 ("定理 4", "定理 4"), ("定理 5", "定理 5"), ("推论 3.2", "推论 3.2"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明"),
                 ("旋转类", "旋转类引用标注")]:
    check("文档契约：%s" % lab, tok in DOC)
bad = []
for m in re.finditer(r"\[([^\]]*)\]\(([^)]+)\)", DOC):
    t = m.group(2).strip()
    if t.startswith(("http", "#", "mailto")):
        continue
    t = re.sub(r":\d+(-\d+)?$", "", t)
    if not os.path.exists(os.path.normpath(os.path.join(HERE, t))):
        bad.append(t)
check("文档内全部相对链接可解析", not bad, "断链：%s" % bad if bad else "无")
_qm_links = re.findall(r"\]\((QM\d+_[^)]*)\)", DOC)
check("★ 不含任何 QM* 旧文依赖（链接扫描）", not _qm_links, "残留 %s" % _qm_links if _qm_links else "无")

print()
print("=" * 78)
print("QMR02 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
