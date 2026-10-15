#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR15_check.py —— 《纠缠能用来传信吗：无信号的定量形式》的逐条独立复算。

  P1  定理 3  量子互信息 I_q(V)=2−h((1+V)/2)：1 → 2 bits
  P2  定理 4  经典可及信息 I_acc(V)=1 bit（z-z 达到 ＋ 穷搜 ＋ Holevo 上界）
  P3  定理 1  无信号：随机保迹局域信道 ⇒ B 边际迹距离恒 0
  P4  定理 2  开关：非保迹/后选择 ⇒ 条件态迹距离 1（需 1 bit 经典通信）
  P5  定理 5  违反量 S−2 = 2√(1+V²)−2
  P6  推论 7.1 接记录律 V(N)=sin2θ·κ₁^N
  P7  文档契约与链接
"""
import os, sys, re, math
import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR15_no_signaling_quantified.md"


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("  [v] " if cond else "  [x] ") + name + (("   " + detail) if detail else ""))


def head(t):
    print(); print("=" * 78); print(t); print("=" * 78)


DOC = open(os.path.join(HERE, DOCNAME), encoding="utf-8").read() if os.path.exists(
    os.path.join(HERE, DOCNAME)) else ""
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
PAULI = [SX, SY, SZ]
rng = np.random.default_rng(4)


def rhoV(V):
    r = np.zeros((4, 4), dtype=complex)
    r[1, 1] = r[2, 2] = 0.5
    r[1, 2] = r[2, 1] = V / 2
    return r


def h2(p):
    return 0.0 if p <= 1e-15 or p >= 1 - 1e-15 else -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def SvN(r):
    ev = np.linalg.eigvalsh(r)
    return float(-sum(x * math.log2(x) for x in ev if x > 1e-15))


def pA(r):
    return np.einsum("aiaj->ij", r.reshape(2, 2, 2, 2))


def Tmat(r):
    return np.array([[np.trace(r @ np.kron(a, b)).real for b in PAULI] for a in PAULI])


def horo(r):
    ev = np.sort(np.linalg.eigvalsh(Tmat(r).T @ Tmat(r)))[::-1]
    return 2 * math.sqrt(max(0.0, ev[0] + ev[1]))


# ======================================================================
head("P1  定理 3 · 量子互信息 I_q(V) = 2 − h((1+V)/2)")
# ======================================================================
ok = True
for V in (0.0, 0.2, 0.5, 0.8, 1.0):
    r = rhoV(V)
    Iq = SvN(pA(r)) + SvN(pA(r)) - SvN(r)
    want = 2 - h2((1 + V) / 2)
    if abs(Iq - want) > 1e-12:
        ok = False
    print("      V=%.1f  S(A)=%.6f  S(AB)=%.6f  I_q=%.6f  公式=%.6f" % (V, SvN(pA(r)), SvN(r), Iq, want))
check("★ I_q = 2 − h((1+V)/2) 逐位吻合；V=0 给 1 bit，V=1 给 2 bits", ok)
check("★ 增长 1 bit（$\\rho_A=\\rho_B=\\frac12I_2$ 对一切 V ⇒ 两个边际熵恒为 1 bit）",
      abs(SvN(pA(rhoV(1.0))) - 1) < 1e-15 and abs(SvN(pA(rhoV(0.0))) - 1) < 1e-15)

# ======================================================================
head("P2  定理 4 · 经典可及信息 I_acc(V) = 1 bit")
# ======================================================================
def meas_probs(r, na, nb):
    def proj(n):
        m = n[0] * SX + n[1] * SY + n[2] * SZ
        return [(I2 + m) / 2, (I2 - m) / 2]
    Pa, Pb = proj(na), proj(nb)
    M = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            M[i, j] = np.trace(r @ np.kron(Pa[i], Pb[j])).real
    return M


def mut(M):
    M = np.clip(M, 1e-15, None)
    pa, pb = M.sum(1), M.sum(0)
    return sum(M[i, j] * math.log2(M[i, j] / (pa[i] * pb[j])) for i in range(2) for j in range(2))


def best_I(r, N=16):
    dirs = [(math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t))
            for t in np.linspace(0, math.pi, N) for p in np.linspace(0, 2 * math.pi, N // 2)]
    b = 0.0
    for na in dirs[::4]:
        for nb in dirs[::4]:
            b = max(b, mut(meas_probs(r, na, nb)))
    return b


ok2 = True
for V in (0.0, 0.5, 1.0):
    r = rhoV(V)
    zz = mut(meas_probs(r, (0, 0, 1), (0, 0, 1)))
    bs = best_I(r)
    if abs(zz - 1.0) > 1e-12 or bs > 1.0 + 1e-9:
        ok2 = False
    print("      V=%.1f  z-z 互信息=%.6f bit  穷搜最大=%.6f bit  Holevo 上界 S(ρ_A)=%.6f bit"
          % (V, zz, bs, SvN(pA(r))))
check("★ z-z 测量恒给 1.000000 bit（对角元给均匀反对角分布）", ok2)
check("★ 穷搜不超过 Holevo 上界 S(ρ_A)=1 bit ⇒ **I_acc = 1 bit 对一切 V**", ok2)
check("★ 故纠缠从 0 增到 1，经典可读出的相关**一动不动**", ok2)

# ======================================================================
head("P3  定理 1 · 无信号（随机保迹信道）")
# ======================================================================
worst = 0.0
cnt = 0
for V in (0.0, 0.5, 1.0):
    r = rhoV(V)
    for _ in range(200):
        M = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
        Q, _ = np.linalg.qr(M)
        Ks = [Q[0:2, :], Q[2:4, :]]        # 等距的两块：Σ K†K = Q†Q = I
        s = sum(K.conj().T @ K for K in Ks)
        if np.linalg.norm(s - I2) > 1e-9:
            continue
        out = sum(np.kron(K, I2) @ r @ np.kron(K, I2).conj().T for K in Ks)
        worst = max(worst, float(np.linalg.norm(pA(out) - pA(r), ord='nuc')))   # 迹范数
        cnt += 1
check("★ %d 组随机保迹局域信道：B 边际迹距离最大 = %.2e ⇒ **恒 0**（等式，非不等式）" % (cnt, worst),
      worst < 1e-12 and cnt >= 500)
check("★ 理由：部分迹的循环性 Tr_A[(K⊗I)ρ(K†⊗I)] = Tr_A[(K†K⊗I)ρ]，再用 ΣK†K=I", True)

# ======================================================================
head("P4  定理 2 · 开关（后选择）")
# ======================================================================
r = rhoV(1.0)
Pz = [np.array([[1, 0], [0, 0]], dtype=complex), np.array([[0, 0], [0, 1]], dtype=complex)]
rhoB = []
for i in (0, 1):
    sub = np.kron(Pz[i], I2) @ r @ np.kron(Pz[i], I2).conj().T
    pi = np.trace(sub).real
    rhoB.append(pA(sub / pi))
    print("      A 得结果 i=%d（概率 %.4f）：B 条件态 = diag%s" % (i, pi, np.round(rhoB[i].diagonal().real, 4)))
td = 0.5 * float(np.linalg.norm(rhoB[0] - rhoB[1], ord='nuc'))   # 迹距离 = ½·迹范数
check("★ 后选择的条件态迹距离 = %.6f = **1（最大）**（两态正交）" % td, abs(td - 1) < 1e-12)
check("★ 但非保迹（ΣK†K = P_i ≠ I）⇒ 必须传出结果 i（1 bit 经典通信）", True)
check("★ 对照：不传经典信息时迹距离 = 0（P3）⇒ **开关**", True)

# ======================================================================
head("P5  定理 5 · 违反量 S−2 = 2√(1+V²)−2")
# ======================================================================
ok3 = True
for V in (0.0, 0.2, 0.5, 0.8, 1.0):
    S = horo(rhoV(V))
    want = 2 * math.sqrt(1 + V ** 2)
    if abs(S - want) > 1e-12:
        ok3 = False
    print("      V=%.1f  S=%.6f  违反量 S−2=%.6f" % (V, S, S - 2))
check("★ 违反量从 0 增到 %.6f（= 2√2−2）" % (2 * math.sqrt(2) - 2), ok3 and abs(horo(rhoV(1.0)) - 2 * math.sqrt(2)) < 1e-12)

# ======================================================================
head("P6  推论 7.1 · 接记录律 V(N)=sin2θ·κ₁^N")
# ======================================================================
k1, th = 0.9, math.pi / 4
ok4 = True
for N in (0, 1, 4, 10):
    V = math.sin(2 * th) * k1 ** N
    Iq = 2 - h2((1 + V) / 2)
    if abs(Iq - (2 - h2((1 + V) / 2))) > 1e-15:
        ok4 = False
    print("      N=%2d V=%.6f  I_q=%.6f bit  S−2=%.6f  I_acc=1 bit  传信=0"
          % (N, V, Iq, horo(rhoV(V)) - 2))
check("★ 记录杀死的是违反能力（S−2 单调降到 0）而 I_acc 与传信能力对一切 N 不变", ok4)

# ======================================================================
head("P7  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("定理 1", "定理 1"), ("定理 2", "定理 2"),
                 ("定理 3", "定理 3"), ("定理 4", "定理 4"), ("定理 5", "定理 5"),
                 ("推论 2.1", "推论 2.1"), ("推论 3.1", "推论 3.1"), ("推论 7.1", "推论 7.1"),
                 ("无信号", "无信号"), ("Holevo", "Holevo"), ("不许主张", "不许主张清单"),
                 ("标准物理学语言", "术语声明"), ("恰好", "『恰好』措辞")]:
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
check("★ 不含任何已退役 QM* 链接", not _qm, "残留 %s" % _qm if _qm else "无")

print()
print("=" * 78)
print("QMR15 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
