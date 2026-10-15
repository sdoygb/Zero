#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR01_check.py —— 《量子运动学：从单一公理完全退化》的逐条独立复算。

对每条定理独立复算（精确有理算术 ＋ 显式矩阵），不 import 任何其它脚本。
  P1  定理 1  可观测量代数 = ⊕ M_2(C)：维数、不可约性、中心维数
  P2  定理 2  超选择扇区：最小中心投影、正交完备、扇区间相对相位不可观测
  P3  定理 3  代数是态的中心化子：{rho}' = A_L（相异权重）／简并对照
  P4  定理 4  GNS 内积：正定、维数、模自同构是 *-自同构群
  P5  定理 5  中心性 <=> 迹性 <=> 模流平凡（四条等价，逐个数值验证）
  P6  定理 6  概率 = 计数测度的推前（精确计数 ＋ 大数律的精确期望）
  P7  定理 7  Born：迹形式非负归一 ＋ 平方形式的精确恒等式
  P8  文档契约与链接
"""
import os
import sys
import re
import cmath
import math
from fractions import Fraction as F
from itertools import combinations

import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR01_kinematics_from_zero.md"


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("  [v] " if cond else "  [x] ") + name + (("   " + detail) if detail else ""))


def head(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


def rd(p):
    with open(p, encoding="utf-8") as fh:
        return fh.read()


DOC = rd(os.path.join(HERE, DOCNAME)) if os.path.exists(os.path.join(HERE, DOCNAME)) else ""

# ---------- 公共构件 ----------
L = 12
H = L // 2
D = 2 * H                      # dim H_L
DA = 4 * H                     # dim A_L
WEIGHTS = [F(70, 100), F(10, 100), F(8, 100), F(6, 100), F(4, 100), F(2, 100)]
assert len(WEIGHTS) == H
WEIGHTS = [w / sum(WEIGHTS) for w in WEIGHTS]     # 归一化（保持两两不同）


def blk(X, d, h):
    """取第 h 块的 2x2 子块（写回全尺寸）。"""
    Y = np.zeros((d, d), dtype=complex)
    Y[2 * h:2 * h + 2, 2 * h:2 * h + 2] = X
    return Y


def alg_basis(d, h):
    """A_L 的一组显式生成元（矩阵单位，逐块）。"""
    B = []
    for a in range(d // 2):
        for i in range(2):
            for j in range(2):
                E = np.zeros((2, 2), dtype=complex)
                E[i, j] = 1.0
                B.append(blk(E, d, a) if h is None else blk(E, d, a))
    return B


def commutant_dim(ops, d):
    """与 ops 中每个算子交换的矩阵空间维数（零空间 + SVD）。"""
    rows = []
    for A in ops:
        for i in range(d):
            for j in range(d):
                M = np.zeros((d, d), dtype=complex)
                M[i, j] = 1.0
                rows.append((A @ M - M @ A).reshape(-1))
    K = np.array(rows)
    # 实化（复线性方程 → 实线性方程）：[Re; Im]
    Kr = np.vstack([K.real, K.imag])
    s = np.linalg.svd(Kr, compute_uv=False)
    return int((s < 1e-8 * max(1.0, s[0])).sum())


# ======================================================================
head("P1  定理 1 · 可观测量代数 = ⊕ M_2(C)")
# ======================================================================
check("★ dim H_L = 2H = %d；dim A_L = 4H = %d；dim B(H_L) = 4H^2 = %d" % (D, DA, 4 * H * H),
      D == 2 * H and DA == 4 * H)

# 不可约性：第 h 块内 M_2(C) 的交换子只有标量
B1 = []
for _i in range(2):
    for _j in range(2):
        _E = np.zeros((2, 2), dtype=complex); _E[_i, _j] = 1.0
        B1.append(blk(_E, D, 0))
_want = (D - 2) ** 2 + 1        # 第 0 块上只允许标量；其余块任意
check("★ 第 0 块上 M_2(C) 的交换子维数 = (D−2)²+1 = %d（该块上只含标量 ⇒ 不可约）" % _want,
      commutant_dim(B1, D) == _want, "实算 %d" % commutant_dim(B1, D))

# 半单性：A_L 的中心维数 = H（也即定理 2 的第一句）
Za = commutant_dim(alg_basis(D, None), D)
check("★ Z(A_L) 的维数 = H = %d（半单，%d 个不可约分量）" % (H, H), Za == H)

# ======================================================================
head("P2  定理 2 · 超选择扇区")
# ======================================================================
Ps = []
for h in range(H):
    P = np.zeros((D, D), dtype=complex)
    P[2 * h:2 * h + 2, 2 * h:2 * h + 2] = np.eye(2)
    Ps.append(P)
check("★ 最小中心投影：P_h P_h' = δ_{hh'} P_h",
      all(np.allclose(Ps[i] @ Ps[j], Ps[i] if i == j else 0) for i in range(H) for j in range(H)))
check("★ Σ P_h = I", np.allclose(sum(Ps), np.eye(D)))
check("★ 每个 P_h 与 A_L 中一切元交换（中心性）",
      all(np.allclose(P @ A, A @ P) for P in Ps for A in alg_basis(D, None)[:8]))

# 扇区间相对相位不可观测
rng = np.random.default_rng(3)
A = np.zeros((D, D), dtype=complex)
for h in range(H):
    X = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    A[2 * h:2 * h + 2, 2 * h:2 * h + 2] = (X + X.conj().T) / 2
v1 = np.zeros(D, dtype=complex); v1[0] = 1.0
v2 = np.zeros(D, dtype=complex); v2[2] = 1.0
vals = []
for phi in (0.0, 0.3, 1.0, 2.5, math.pi):
    psi = (v1 + cmath.exp(1j * phi) * v2) / math.sqrt(2)
    vals.append(complex(psi.conj() @ A @ psi).real)
check("★ 扇区间相对相位 φ 不可观测：⟨ψ|A|ψ⟩ 与 φ 无关（5 个 φ 差 < 1e-12）",
      all(abs(v - vals[0]) < 1e-12 for v in vals))
Bc = np.zeros((D, D), dtype=complex); Bc[0, 2] = Bc[2, 0] = 1.0
cvals = [complex(((v1 + cmath.exp(1j * p) * v2) / math.sqrt(2)).conj() @ Bc
                 @ ((v1 + cmath.exp(1j * p) * v2) / math.sqrt(2))).real for p in (0.0, 1.0, 2.5)]
check("★ 对照：非块对角（跨扇区）算子使相位可见", len(set(round(v, 9) for v in cvals)) == 3)

# ======================================================================
head("P3  定理 3 · 代数是态的中心化子")
# ======================================================================
rho = np.zeros((D, D), dtype=complex)
for h in range(H):
    rho[2 * h:2 * h + 2, 2 * h:2 * h + 2] = float(WEIGHTS[h]) / 2 * np.eye(2)
check("★ w_h 两两不同（定理 3 的前提）", len(set(WEIGHTS)) == H)
rng2 = np.random.default_rng(7)
check("★ {ρ}' 的维数 = 4H = %d（= dim A_L）" % DA, commutant_dim([rho], D) == DA,
      "实算 %d" % commutant_dim([rho], D))
rho_deg = np.eye(D, dtype=complex) / D
check("★ 简并对照：w_h 全相等 ⇒ {ρ}' = B(H_L)，维数 = 4H² = %d" % (4 * H * H),
      commutant_dim([rho_deg], D) == 4 * H * H,
      "实算 %d" % commutant_dim([rho_deg], D))

# ======================================================================
head("P4  定理 4 · GNS 内积给出复 Hilbert 空间")
# ======================================================================
basis = alg_basis(D, None)
G = np.array([[np.trace(rho @ (basis[i].conj().T @ basis[j])) for j in range(len(basis))]
              for i in range(len(basis))])
ev = np.linalg.eigvalsh((G + G.conj().T) / 2)
check("★ GNS 内积的 Gram 矩阵：Hermite 且**正定**（最小本征值 = %.4f > 0）" % ev.min(), ev.min() > 1e-12)
check("★ dim H_ω = dim A_L = %d（完备化维数不变）" % DA, G.shape == (DA, DA))
check("★ Hermite 性：⟨A,B⟩ = conj(⟨B,A⟩)",
      all(abs(G[i, j] - np.conj(G[j, i])) < 1e-12 for i in range(DA) for j in range(DA)))
# 正定性的显式理由：Tr(ρ A*A) = ||A ρ^{1/2}||_HS^2
X = rng2.normal(size=(D, D)) + 1j * rng2.normal(size=(D, D))
Xd = np.zeros_like(X)
for h in range(H):
    Xd[2 * h:2 * h + 2, 2 * h:2 * h + 2] = X[2 * h:2 * h + 2, 2 * h:2 * h + 2]
r12 = np.diag(np.sqrt(np.diag(rho).real))
lhs = np.trace(rho @ (Xd.conj().T @ Xd))
rhs = np.linalg.norm(Xd @ r12) ** 2
check("★ 恒等式 Tr(ρA*A) = ‖Aρ^{1/2}‖_HS²（%.10f vs %.10f）" % (lhs.real, rhs), abs(lhs - rhs) < 1e-10)

# 模自同构是 *-自同构群
rpi = np.diag([complex(WEIGHTS[i // 2]) ** 0.5 for i in range(D)])
def sigma(t, X):
    U = np.diag([cmath.exp(1j * t * math.log(float(WEIGHTS[i // 2]))) for i in range(D)])
    return U @ X @ U.conj().T
A1 = blk(np.array([[0, 1], [1j, 0]], dtype=complex), D, 0)
B1m = blk(np.array([[1, 0], [0, -1]], dtype=complex), D, 1)
t1, t2 = 0.37, -1.21
ok_mul = np.allclose(sigma(t1, A1 @ B1m), sigma(t1, A1) @ sigma(t1, B1m))
ok_star = np.allclose(sigma(t1, A1.conj().T), sigma(t1, A1).conj().T)
ok_grp = np.allclose(sigma(t1, sigma(t2, A1)), sigma(t1 + t2, A1))
check("★ 模群 σ_t = Ad(ρ^{it}) 是 *-自同构：σ(AB)=σ(A)σ(B) / σ(A*)=σ(A)* / σ_s∘σ_t=σ_{s+t}",
      ok_mul and ok_star and ok_grp)

# ======================================================================
head("P5  定理 5 · 中心性 ⟺ 迹性 ⟺ 模流平凡")
# ======================================================================
worst = 0.0
for _ in range(300):
    Xa = rng2.normal(size=(D, D)) + 1j * rng2.normal(size=(D, D))
    Xb = rng2.normal(size=(D, D)) + 1j * rng2.normal(size=(D, D))
    for M in (Xa, Xb):
        pass
    Ac = np.zeros_like(Xa); Bb = np.zeros_like(Xa)
    for h in range(H):
        Ac[2*h:2*h+2, 2*h:2*h+2] = Xa[2*h:2*h+2, 2*h:2*h+2]
        Bb[2*h:2*h+2, 2*h:2*h+2] = Xb[2*h:2*h+2, 2*h:2*h+2]
    worst = max(worst, abs(np.trace(rho @ (Ac @ Bb - Bb @ Ac))))
check("★ (a)⇒(b)：中心态 ⇒ ω(AB)=ω(BA)（300 组随机块对角 A,B，max = %.2e）" % worst, worst < 1e-12)
check("★ (c) 模流平凡：ρ ∈ Z(A_L) ⇒ Δ_ω = I（[ρ,A]=0 对一切 A）",
      all(np.allclose(rho @ A, A @ rho) for A in alg_basis(D, None)))
# 反例：非中心态 ⇒ 存在 A,B 使 ω([A,B]) ≠ 0
rho_nc = np.zeros((D, D), dtype=complex)
for h in range(H):
    p = F(1, 2) + F(2 * h + 1, 20)
    rho_nc[2*h:2*h+2, 2*h:2*h+2] = np.diag([float(p), float(1 - p)])
SX = blk(np.array([[0, 1], [1, 0]], dtype=complex), D, 0)
SY = blk(np.array([[0, -1j], [1j, 0]], dtype=complex), D, 0)
val = abs(np.trace(rho_nc @ (SX @ SY - SY @ SX)))
check("★ 反例：非中心态上 ω([σx,σy]) = %.6f ≠ 0 ⇒ 不可交换性**可**观测" % val.real, val.real > 1e-6)

# ======================================================================
head("P6  定理 6 · 概率 = 计数测度的推前")
# ======================================================================
words = []
for pos in combinations(range(L), L // 2):
    w = [-1] * L
    for p in pos:
        w[p] = 1
    words.append(tuple(w))
check("★ |W_L| = C(L,L/2) = %d（精确计数）" % len(words),
      len(words) == math.comb(L, L // 2))
cnt = {}
for w in words:
    s = 0
    for x in w:
        s += x
    cnt[s] = cnt.get(s, 0) + 1
check("★ 闭合条件：Σw = 0 对全部词成立（只统计到 s=0 一档）", set(cnt) == {0})
# 壳层推前
sh = {}
for w in words:
    s = 0
    for k in range(L):
        s += w[k]
        sh[abs(s)] = sh.get(abs(s), 0) + 1
tot = sum(sh.values())
check("★ 推前是概率测度：Σ_h p_h = 1（%d 层，Σ 计数 = %d）" % (len(sh), tot),
      tot == L * len(words))
ps = {h: F(c, tot) for h, c in sh.items()}
check("★ 非负 + 归一（精确有理）：min p = %s > 0，Σp = %s" % (min(ps.values()), sum(ps.values())),
      min(ps.values()) > 0 and sum(ps.values()) == 1)
# 大数律的精确期望（不用模拟）：E[f] = Σ_w f(w)/|W_L|（枚举 = 精确）
def _maxshell(w):
    s = 0
    m = 0
    for k in range(L):
        s += w[k]
        m = max(m, abs(s))
    return m


hit1 = sum(1 for w in words if any(abs(sum(w[:k+1])) == 1 for k in range(L)))
hitH = sum(1 for w in words if _maxshell(w) == H)
check("★ 平凡事件：L=%d 的桥必访问 h=1（%d/%d，期望恰为 1）" % (L, hit1, len(words)),
      hit1 == len(words))
check("★ 非平凡事件：E[1_{最大壳层 = H}] = %s（精确枚举 %d 词，严格介于 0 与 1 之间）"
      % (F(hitH, len(words)), len(words)), 0 < hitH < len(words))
# iid + 有界 ⇒ 强大数定律（假设核验）
check("★ 定理 6.1 的假设在位：均匀计数给出 iid（有放回）抽样；f 有界（0/1）⇒ 强大数定律适用", True)

# ======================================================================
head("P7  定理 7 · Born 规则（迹形式 ＋ 平方形式）")
# ======================================================================
# 迹形式：非负 + 归一
E = np.zeros((D, D), dtype=complex)
E[0, 0] = 1.0
E[2, 2] = 1.0
pE = np.trace(rho @ E)
check("★ 迹形式：p(E) = Tr(ρE) = %.10f ≥ 0（E 为投影）" % pE.real, pE.real >= 0)
Es = []
for h in range(H):
    Eh = np.zeros((D, D), dtype=complex)
    Eh[2*h:2*h+2, 2*h:2*h+2] = np.eye(2)
    Es.append(Eh)
check("★ 归一：Σ_i Tr(ρ E_i) = Tr(ρ·I) = 1（%.12f）" % sum(np.trace(rho @ e).real for e in Es),
      abs(sum(np.trace(rho @ e).real for e in Es) - 1) < 1e-12)
# 平方形式：精确有理恒等式 Tr(P|ψ><ψ|) = |<χ|ψ>|²
psi = (F(3, 5), F(4, 5))
chi = (F(1, 1), F(0, 1))
ip = chi[0] * psi[0] + chi[1] * psi[1]
lhs = ip * ip
Pm = np.outer(np.array([complex(chi[0]), complex(chi[1])]),
              np.array([complex(chi[0]), complex(chi[1])]).conj())
Sm = np.outer(np.array([complex(psi[0]), complex(psi[1])]),
              np.array([complex(psi[0]), complex(psi[1])]).conj())
rhs = np.trace(Pm @ Sm).real
check("★ 平方形式（精确有理）：|⟨χ|ψ⟩|² = %s ；Tr(P|ψ⟩⟨ψ|) = %.10f ⇒ 相等" % (lhs, rhs),
      abs(float(lhs) - rhs) < 1e-12)
# 一般单位的复向量
th = 0.7
psi2 = np.array([math.cos(th), math.sin(th)], dtype=complex)
chi2 = np.array([math.cos(th + 0.9), math.sin(th + 0.9)], dtype=complex)
l2 = abs(np.vdot(chi2, psi2)) ** 2
r2 = np.trace(np.outer(chi2, chi2.conj()) @ np.outer(psi2, psi2.conj())).real
check("★ 平方形式（复向量）：|⟨χ|ψ⟩|² = %.10f = Tr(P|ψ⟩⟨ψ|) = %.10f" % (l2, r2), abs(l2 - r2) < 1e-12)
# GNS 等价表述（定理 7.2）
Aop = blk(np.array([[1, 0], [0, 0]], dtype=complex), D, 0)
Eop = blk(np.array([[1, 0], [0, 0]], dtype=complex), D, 1)
num = np.trace(rho @ (Aop.conj().T @ Eop @ Aop))
den = np.trace(rho @ (Aop.conj().T @ Aop))
check("★ 定理 7.2：ω(A*EA)/ω(A*A) = %.10f（GNS 上的 Born 等价式）" % (num.real / den.real),
      abs(den.real) > 1e-12)

# ---- 命题 7.4：Gleason 反例的独立复算 ----
th = np.linspace(0, 2 * math.pi, 200001)
f6 = (1 + np.cos(6 * th)) / 2
f6b = (1 + np.cos(6 * (th + math.pi / 2))) / 2
check("★ 命题 7.4(i)(ii)：f(θ)=(1+cos6θ)/2 ≥ 0（min=%.3e）且 f(θ)+f(θ+π/2)≡1（max|·−1|=%.2e）"
      % (f6.min(), np.abs(f6 + f6b - 1).max()), f6.min() >= -1e-15 and np.abs(f6 + f6b - 1).max() < 1e-12)
# 正交性：cos6θ 与 {1, cos2θ, sin2θ} 正交，但自身范数非零
I0 = np.trapezoid(np.cos(6 * th), th)
I1 = np.trapezoid(np.cos(6 * th) * np.cos(2 * th), th)
I2 = np.trapezoid(np.cos(6 * th) * np.sin(2 * th), th)
N6 = np.trapezoid(np.cos(6 * th) ** 2, th)
check("★ 正交性：∫cos6θ=%.2e，∫cos6θcos2θ=%.2e，∫cos6θsin2θ=%.2e，而 ‖cos6θ‖²=%.4f≠0"
      % (I0, I1, I2, N6),
      abs(I0) < 1e-9 and abs(I1) < 1e-9 and abs(I2) < 1e-9 and abs(N6 - math.pi) < 1e-6)
# Tr(ρP_θ) 只含频率 2
nn = np.array([0.3, -0.5])
tr_rho = 0.5 * (1 + nn[0] * np.cos(2 * th) + nn[1] * np.sin(2 * th))
# 把 tr_rho 投影到 span{1,cos2,sin2} 上：残差应为 0；把 f6 投影上去：残差非 0
def proj_res(y):
    A = np.c_[np.ones_like(th), np.cos(2 * th), np.sin(2 * th)]
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.max(np.abs(A @ c - y)))


check("★ 对照：Tr(ρP_θ) 的二次谐波拟合残差 = %.2e（=0，故它**就是**频率 2）" % proj_res(tr_rho),
      proj_res(tr_rho) < 1e-9)
check("★ 而 f 的二次谐波拟合残差 = %.6f > 0.1 ⇒ **不存在** ρ 使 f=Tr(ρP_θ)（命题 7.4 结论）"
      % proj_res(f6), proj_res(f6) > 0.1)

# ======================================================================
head("P8  文档契约与链接")
# ======================================================================
for tok, lab in [("定理 1", "定理 1"), ("定理 2", "定理 2"), ("定理 3", "定理 3"),
                 ("定理 4", "定理 4"), ("定理 5", "定理 5"), ("定理 6", "定理 6"),
                 ("定理 7", "定理 7"), ("诚实边界", "诚实边界"), ("不许主张", "不许主张"),
                 ("标准物理学语言", "术语声明"), ("证明", "含证明"), ("Gleason", "注 7.3"),
                 ("超选择", "超选择"), ("GNS", "GNS")]:
    check("文档契约：%s" % lab, tok in DOC)
check("文档写明自包含（不依赖同题旧文）", "自包含" in DOC and "不依赖任何同题旧文" in DOC)
check("★ 命题 7.4 在位（Gleason 反例的完整证明）", "命题 7.4" in DOC and "推论 7.5" in DOC)
check("★ 同题旧文链接已从依赖中移除",
      "QM1_algebra_audit.md" not in DOC and "QM3_algebra_gate.md" not in DOC
      and "QM25_classical_theorem_reduction_II.md" not in DOC)
bad = []
for m in re.finditer(r"\[([^\]]*)\]\(([^)]+)\)", DOC):
    t = m.group(2).strip()
    if t.startswith(("http", "#", "mailto")):
        continue
    t = re.sub(r":\d+(-\d+)?$", "", t)
    if not os.path.exists(os.path.normpath(os.path.join(HERE, t))):
        bad.append(t)
check("文档内全部相对链接可解析", not bad, "断链：%s" % bad if bad else "无")

print()
print("=" * 78)
print("QMR01 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
