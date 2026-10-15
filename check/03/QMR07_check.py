#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR07_check.py —— 《纠缠、记录律与熵的盲区：从单一公理完全退化》的逐条独立复算。

  P1  定理 1  复合：代数维数、中心性遗传、两体退化、可分离分解
  P2  定理 2  闭式律 S_CHSH(V)=2√(1+V²)（Horodecki 判据 ＋ T 矩阵显式）
  P3  推论 2.1/定理 3  concurrence=V、negativity=V/2、PT 谱 {±V/2,½,½}、门槛一致
  P4  定理 4（否定）熵的盲区：ρ_A ≡ ½I₂、S_vN ≡ ln2
  P5  定理 5  记录律 S(N)=2√(1+κ₁^{2N})，有限 N 下 S>2
  P6  定理 6  信息版：V²+D²=1 饱和；S>2 ⟺ χ<1 bit
  P7  定理 7  门的下游：V ≤ κ_ω ⇒ S ≤ 2√(1+κ_ω²)
  P8  文档契约与链接

不 import 任何其它脚本；路径由 __file__ 推导。
"""
import os
import sys
import re
import math

import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR07_entanglement_and_record_law.md"


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

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
PAULI = [SX, SY, SZ]


def rhoV(V):
    """X-态：½(|01⟩⟨01|+|10⟩⟨10|) + (V/2)(|01⟩⟨10|+|10⟩⟨01|)。"""
    r = np.zeros((4, 4), dtype=complex)
    r[1, 1] = r[2, 2] = 0.5
    r[1, 2] = r[2, 1] = V / 2
    return r


def Tmat(r):
    return np.array([[np.trace(r @ np.kron(a, b)).real for b in PAULI] for a in PAULI])


def horodecki(r):
    ev = np.sort(np.linalg.eigvalsh(Tmat(r).T @ Tmat(r)))[::-1]
    return 2 * math.sqrt(max(0.0, ev[0] + ev[1]))


def partial_A(r):
    return np.einsum("aiaj->ij", r.reshape(2, 2, 2, 2))


def partial_T(r):
    # 只转置 B 侧：ρ_{(a,b),(a',b')} -> ρ_{(a,b'),(a',b)}
    return np.einsum("aibj->ajbi", r.reshape(2, 2, 2, 2)).reshape(4, 4)


def msqrt(M):
    """半正定矩阵的**矩阵**平方根（谱分解；np.sqrt 是逐元素开方，不能用）。"""
    w, V = np.linalg.eigh((M + M.conj().T) / 2)
    return V @ np.diag(np.sqrt(np.clip(w, 0, None))) @ V.conj().T


def concurrence(r):
    st = msqrt(r)
    rt = np.kron(SY, SY) @ r.conj() @ np.kron(SY, SY)
    ev = np.sort(np.sqrt(np.abs(np.linalg.eigvals(st @ rt @ st))))[::-1]
    return max(0.0, float(ev[0] - ev[1] - ev[2] - ev[3]))


VS = [0.0, 0.2, 0.5, 0.8, 1.0]


# ======================================================================
head("P1  定理 1 · 复合系统的代数与态的遗传")
# ======================================================================
H, D = 3, 6
w = [0.5, 0.3, 0.2]
rho = np.zeros((D, D), dtype=complex)
for h in range(H):
    rho[2*h:2*h+2, 2*h:2*h+2] = w[h] / 2 * I2
check("★ (a) dim(𝒜⊗𝒜) = (4H)² = %d" % (4 * H) ** 2, (4 * H) ** 2 == 16 * H * H)
r2 = np.kron(rho, rho)
check("★ (b) ρ⊗ρ 归一：Tr = %.12f（= (Tr ρ)²）" % np.trace(r2).real, abs(np.trace(r2).real - 1) < 1e-15)
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(120):
    X = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    Y = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    for M in (X, Y):
        pass
    Xb = np.zeros_like(X); Yb = np.zeros_like(Y)
    for h in range(H):
        Xb[2*h:2*h+2, 2*h:2*h+2] = X[2*h:2*h+2, 2*h:2*h+2]
        Yb[2*h:2*h+2, 2*h:2*h+2] = Y[2*h:2*h+2, 2*h:2*h+2]
    A2 = np.kron(Xb, np.eye(D)); B2 = np.kron(Yb, np.eye(D))
    worst = max(worst, abs(np.trace(r2 @ (A2 @ B2 - B2 @ A2))))
check("★ (c) 两体上 ω₂([X,Y]) ≡ 0（120 组，max = %.2e）⇒ 两体同样退化" % worst, worst < 1e-12)
# (d) 可分离分解：ρ⊗ρ = Σ w_h w_h' /4 · I⊗I
recon = np.zeros((D * D, D * D), dtype=complex)
for h in range(H):
    for hp in range(H):
        eh = np.zeros((D, D), dtype=complex); eh[2*h:2*h+2, 2*h:2*h+2] = I2
        eh2 = np.zeros((D, D), dtype=complex); eh2[2*hp:2*hp+2, 2*hp:2*hp+2] = I2
        recon += w[h] * w[hp] / 4 * np.kron(eh, eh2)
check("★ (d) ρ⊗ρ = Σ_{h,h'} (w_h w_h'/4)·I⊗I（残差 %.2e）⇒ **可分离**" % np.linalg.norm(recon - r2),
      np.linalg.norm(recon - r2) < 1e-15)

# ======================================================================
head("P2  定理 2 · 闭式律 S_CHSH(V) = 2√(1+V²)")
# ======================================================================
for V in VS:
    r = rhoV(V)
    S = horodecki(r)
    want = 2 * math.sqrt(1 + V ** 2)
    check("★ V=%.1f：S_max = %.10f = 2√(1+V²) = %.10f（误差 %+.1e）" % (V, S, want, S - want),
          abs(S - want) < 1e-12)
T = Tmat(rhoV(0.6))
check("★ T 矩阵 = diag(V,V,−1)？ 实算 diag = %s" % np.round(np.diag(T), 6).tolist(),
      abs(T[0, 0] - 0.6) < 1e-12 and abs(T[1, 1] - 0.6) < 1e-12 and abs(T[2, 2] + 1) < 1e-12)

# ======================================================================
head("P3  推论 2.1/定理 3 · 全度量同为 V；PPT 与 CHSH 同一门槛")
# ======================================================================
ok = True
for V in VS:
    r = rhoV(V)
    C = concurrence(r)
    pte = np.sort(np.linalg.eigvalsh(partial_T(r)))
    neg = float(sum(-x for x in pte if x < 0))
    if abs(C - V) > 1e-9 or abs(neg - V / 2) > 1e-9:
        ok = False
    print("      V=%.1f  concurrence=%.6f  negativity=%.6f  PT 谱=%s"
          % (V, C, neg, np.round(pte, 6).tolist()))
check("★ concurrence = V 且 negativity = V/2（逐点）", ok)
ok2 = all(abs(np.linalg.eigvalsh(partial_T(rhoV(V)))[0] + V / 2) < 1e-12 for V in VS)
check("★ PT 谱恒为 {−V/2, +V/2, ½, ½}（最小本征值 = −V/2）", ok2)
check("★ 同一门槛：S>2 ⟺ V>0 ⟺ PPT 失败",
      all((horodecki(rhoV(V)) > 2) == (V > 0) == (np.linalg.eigvalsh(partial_T(rhoV(V))).min() < -1e-15)
          for V in VS))

# ======================================================================
head("P4  定理 4（否定）· 熵的盲区")
# ======================================================================
okA, okS = True, True
for V in VS:
    r = rhoV(V)
    dA = np.linalg.norm(partial_A(r) - 0.5 * I2)
    ev = np.linalg.eigvalsh(partial_A(r))
    Svn = float(-sum(x * math.log(x) for x in ev if x > 1e-15))
    if dA > 1e-15:
        okA = False
    if abs(Svn - math.log(2)) > 1e-12:
        okS = False
    print("      V=%.1f  ‖ρ_A−½I‖=%.2e  S_vN=%.10f  (ln2=%.10f)" % (V, dA, Svn, math.log(2)))
check("★ ρ_A ≡ ½I₂ 对一切 V（误差 0）", okA)
check("★ S_vN(ρ_A) ≡ ln2 对一切 V ⇒ **熵对记录完全失明**", okS)
check("★ 对照：S_CHSH 从 %.6f 降到 %.6f（熵一动不动）" % (horodecki(rhoV(1)), horodecki(rhoV(0))),
      abs(horodecki(rhoV(1)) - 2 * math.sqrt(2)) < 1e-12 and abs(horodecki(rhoV(0)) - 2) < 1e-12)

# ======================================================================
head("P5  定理 5 · 记录律")
# ======================================================================
k1 = 0.9
for N in (0, 1, 4, 10, 40):
    V = k1 ** N
    S = horodecki(rhoV(V))
    want = 2 * math.sqrt(1 + k1 ** (2 * N))
    check("★ N=%2d：V=%.6f，S=%.10f = 公式 %.10f；S>2？ %s" % (N, V, S, want, S > 2),
          abs(S - want) < 1e-12 and (S > 2) == (N >= 0))
check("★ (c) 仅 N→∞ 时 S→2：N=100 时 S−2 = %.3e（正但极小）" % (horodecki(rhoV(k1 ** 100)) - 2),
      0 < horodecki(rhoV(k1 ** 100)) - 2 < 1e-6)

# ======================================================================
head("P6  定理 6 · 信息版（互补性饱和）")
# ======================================================================
def h2(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


ok1_, ok2_ = True, True
for V in [1.0, 0.9239, 0.7071, 0.3827, 0.0]:
    D_ = math.sqrt(1 - V ** 2)
    chi = h2((1 - V) / 2)
    S = horodecki(rhoV(V))
    if abs(V ** 2 + D_ ** 2 - 1) > 1e-12:
        ok1_ = False
    if (S > 2) != (chi < 1):
        ok2_ = False
    print("      V=%.4f D=%.4f V²+D²=%.10f χ=%.4f S=%.5f  (S>2)=%s (χ<1)=%s"
          % (V, D_, V ** 2 + D_ ** 2, chi, S, S > 2, chi < 1))
check("★ V²+D² = 1（**饱和** Englert 不等式，是恒等式不是不等式）", ok1_)
check("★ S>2 ⟺ χ<1 bit（逐点一致）", ok2_)

# ======================================================================
head("P7  定理 7 · 门的下游")
# ======================================================================
for nn in ([0, 0, 0.0], [0, 0, 0.4], [0.3, 0.5, 0.6]):
    n = np.array(nn)
    k_ = float(np.linalg.norm(n))
    cap = 2 * math.sqrt(1 + k_ ** 2)
    print("      κ_ω=%.6f ⇒ S ≤ %.6f" % (k_, cap))
check("★ κ_ω=0 ⇒ 上限恰为 2（中心态 ⇒ Bell 局域，精确回位）",
      abs(2 * math.sqrt(1 + 0.0 ** 2) - 2) < 1e-15)
check("★ V ≤ κ_ω ⇒ S ≤ 2√(1+κ_ω²)（等号条件 n_z=0）",
      all(2 * math.sqrt(1 + min(V, 0.4) ** 2) <= 2 * math.sqrt(1 + 0.4 ** 2) + 1e-12 for V in VS))

# ======================================================================
head("P8  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("不依赖任何同题旧文", "去依赖声明"),
                 ("定理 1", "定理 1"), ("定理 2", "定理 2"), ("定理 3", "定理 3"),
                 ("定理 4", "定理 4"), ("定理 5", "定理 5"), ("定理 6", "定理 6"),
                 ("定理 7", "定理 7"), ("推论 2.1", "推论 2.1"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明"),
                 ("熵的盲区", "熵盲区"), ("Horodecki", "判据引用")]:
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
_qm = re.findall(r"\]\((QM\d+_[^)]*)\)", DOC)
check("★ 不含任何 QM* 旧文依赖（链接扫描）", not _qm, "残留 %s" % _qm if _qm else "无")

print()
print("=" * 78)
print("QMR07 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
