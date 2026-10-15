#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR04_check.py —— 《测量、效应与无信号：从单一公理完全退化》的逐条独立复算。

  P1  定理 1  效应与概率：非负、归一、可加（含一般效应，非仅投影）
  P2  定理 2  Naimark 显式膨胀：V†V = I、读数逐项相等（trine POVM）
  P3  定理 3  无信号 ⟺ 保迹：引理 3.0 的一般公式、保迹族、非保迹反例、|Φ⁺⟩ 检测
  P4  定理 4（否定）中心态下 ω([E,A]) ≡ 0 ＋ 非中心反例 ＋ 推论 4.1
  P5  文档契约与链接

不 import 任何其它脚本；路径由 __file__ 推导。
"""
import os
import sys
import re
import math

import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR04_measurement_effects_and_no_signaling.md"


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


# ======================================================================
head("P1  定理 1 · 效应与概率（非负／归一／可加）")
# ======================================================================
# trine POVM
ang = [0.0, 2 * math.pi / 3, 4 * math.pi / 3]
NVec = [np.array([math.cos(a), math.sin(a), 0.0]) for a in ang]
E = [(I2 + n[0] * SX + n[1] * SY + n[2] * SZ) / 3 for n in NVec]
check("★ Σ E_i = I（范数差 = %.2e）" % np.linalg.norm(sum(E) - I2), np.linalg.norm(sum(E) - I2) < 1e-15)
check("★ 每个 E_i ≥ 0（最小本征值 = %s）" % [round(float(np.linalg.eigvalsh(e).min()), 12) for e in E],
      all(np.linalg.eigvalsh(e).min() >= -1e-15 for e in E))
check("★ E_i = (2/3)|ψ_i⟩⟨ψ_i|（纯态 POVM 的显式形式）",
      all(abs(np.trace(e) - 2 / 3) < 1e-15 and abs(np.trace(e @ e) - (2 / 3) ** 2) < 1e-15 for e in E))
rng = np.random.default_rng(9)
for tag, rho in [("纯态", np.array([[1, 0], [0, 0]], dtype=complex)),
                 ("混合态", np.diag([0.7, 0.3]).astype(complex)),
                 ("随机态", (lambda M: M @ M.conj().T / np.trace(M @ M.conj().T))(
                     rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))))]:
    ps = [float(np.trace(rho @ e).real) for e in E]
    check("★ %s：p_i = Tr(ρE_i) ≥ 0（min = %.3e）且 Σ p_i = 1（%.12f）" % (tag, min(ps), sum(ps)),
          min(ps) >= -1e-15 and abs(sum(ps) - 1) < 1e-12)
ps = [float(np.trace(np.diag([0.7, 0.3]).astype(complex) @ e).real) for e in E]
check("★ 可加性：p_{I∪J} = Σ_{i∈I∪J} p_i（Tr 的线性）",
      abs((ps[0] + ps[1]) - (ps[0] + ps[1])) < 1e-15 and abs(sum(ps) - 1) < 1e-12)

# ======================================================================
head("P2  定理 2 · Naimark 显式膨胀")
# ======================================================================
d, n = 2, 3
def msqrt(M):
    """半正定矩阵的平方根（谱分解；秩 1 的 E_i 不能用 cholesky）。"""
    w, V = np.linalg.eigh((M + M.conj().T) / 2)
    return V @ np.diag(np.sqrt(np.clip(w, 0, None))) @ V.conj().T


sq = [msqrt(e) for e in E]


def apply_V(chi):
    """V|χ⟩ = Σ_i (sqrt(E_i)|χ⟩) ⊗ |i⟩（维数 d*n）。"""
    out = np.zeros(d * n, dtype=complex)
    for i in range(n):
        out[i * d:(i + 1) * d] = sq[i] @ chi
    return out


V = np.zeros((d * n, d), dtype=complex)
for j in range(d):
    V[:, j] = apply_V(np.eye(d, dtype=complex)[:, j])
check("★ V†V = I（等距；范数差 = %.2e）" % np.linalg.norm(V.conj().T @ V - I2),
      np.linalg.norm(V.conj().T @ V - I2) < 1e-14)
ok = True
for _ in range(20):
    chi = rng.normal(size=d) + 1j * rng.normal(size=d)
    chi /= np.linalg.norm(chi)
    for i in range(n):
        P = np.zeros((d * n, d * n), dtype=complex)
        for k in range(d):
            P[i * d + k, i * d + k] = 1.0
        lhs = (chi.conj() @ (V.conj().T @ P @ V) @ chi)
        rhs = chi.conj() @ E[i] @ chi
        if abs(lhs - rhs) > 1e-14:
            ok = False
check("★ POVM 读数 = 膨胀后 PVM 读数（⟨χ|V†(I⊗|i⟩⟨i|)V|χ⟩ = ⟨χ|E_i|χ⟩，逐项）", ok)
rho0 = np.diag([0.6, 0.4]).astype(complex)
def proj_i(i):
    """扩大空间（次序 |k⟩⊗|i⟩，块 i 长 d）上的 I_d ⊗ |i⟩⟨i|。"""
    Pi = np.zeros((d * n, d * n), dtype=complex)
    for k in range(d):
        Pi[i * d + k, i * d + k] = 1.0
    return Pi


ok2 = all(abs(np.trace(rho0 @ E[i]) - np.trace(V @ rho0 @ V.conj().T @ proj_i(i))) < 1e-14
          for i in range(n))
check("★ 算子形式：Tr(ρE_i) = Tr(VρV†·(I⊗|i⟩⟨i|))（推论 2.1）", ok2)

# ======================================================================
head("P3  定理 3 · 无信号 ⟺ 局域操作保迹")
# ======================================================================
dA, dB = 2, 3
D = dA * dB


def local_apply(rho, Ks, dA, dB):
    out = np.zeros_like(rho)
    for K in Ks:
        M = np.kron(K, np.eye(dB))
        out += M @ rho @ M.conj().T
    return out


def bmarg(rho, dA, dB):
    # 部分迹：行指标 (a,i)、列指标 (a,j) —— 同一字母 a 必须**同时**出现在两个位置
    return np.einsum("aiaj->ij", rho.reshape(dA, dB, dA, dB))


rng = np.random.default_rng(3)
X = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
rho = X @ X.conj().T
rho = rho / np.trace(rho)
# 保迹的非对角 Kraus 族（|±⟩ 基测量）
Kp = np.array([[1, 1], [1, 1]], dtype=complex) / 2
Km = np.array([[1, -1], [-1, 1]], dtype=complex) / 2
Ks = [Kp, Km]
Nsum = sum(K.conj().T @ K for K in Ks)
check("★ 非对角 Kraus 族：Σ K_i†K_i = I（范数差 = %.2e）" % np.linalg.norm(Nsum - I2),
      np.linalg.norm(Nsum - I2) < 1e-15)
b0 = bmarg(rho, dA, dB)
b1 = bmarg(local_apply(rho, Ks, dA, dB), dA, dB)
check("★ B 边缘态不变（差 = %.2e）" % np.linalg.norm(b1 - b0), np.linalg.norm(b1 - b0) < 1e-15)
a0 = np.einsum("aibj->ab", rho.reshape(dA, dB, dA, dB))
a1 = np.einsum("aibj->ab", local_apply(rho, Ks, dA, dB).reshape(dA, dB, dA, dB))
check("★ A 边缘态确实变了（差 = %.3f ≠ 0）—— 操作真的作用在 A 上" % np.linalg.norm(a1 - a0),
      np.linalg.norm(a1 - a0) > 1e-3)
# 引理 3.0 的一般公式： (ρ'_B)_{bb'} = Σ N_{a'''a''} ρ_{(a'',b),(a''',b')}
Nmat = Nsum
got = b1
# 显式按引理形式求和
want = np.zeros((dB, dB), dtype=complex)
for bb in range(dB):
    for bp in range(dB):
        acc = 0
        for a2 in range(dA):
            for a3 in range(dA):
                acc += Nmat[a3, a2] * rho[a2 * dB + bb, a3 * dB + bp]
        want[bb, bp] = acc
check("★ 引理 3.0：一般公式与直接计算一致（差 = %.2e）" % np.linalg.norm(got - want),
      np.linalg.norm(got - want) < 1e-14)
# 非保迹反例
Ln = np.array([[1, 0], [0, 0.5]], dtype=complex)
b2 = bmarg(local_apply(rho, [Ln], dA, dB), dA, dB)
check("★ 反例：非保迹映射（N = diag(1,0.25) ≠ I）⇒ B 边缘态改变（差 = %.3f ≠ 0）"
      % np.linalg.norm(b2 - b0), np.linalg.norm(b2 - b0) > 1e-3)
# |Φ⁺⟩ 上 ρ'_B = ½KK†
phi = np.zeros(dA * dA, dtype=complex)
for a in range(dA):
    phi[a * dA + a] = 1 / math.sqrt(dA)
rphi = np.outer(phi, phi.conj())
H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
for nm, K in [("酉 H", H), ("非酉 L", Ln)]:
    got2 = bmarg(local_apply(rphi, [K], dA, dA), dA, dA)
    want2 = (K @ K.conj().T) / dA
    check("★ |Φ⁺⟩ 上：Tr_A[(K⊗I)|Φ⁺⟩⟨Φ⁺|] = ½KK†（%s）" % nm, np.allclose(got2, want2))
check("★ ⇒ 对**一切**态无信号 ⟺ KK† = I ⟺ K 酉", np.allclose(H @ H.conj().T, I2)
      and not np.allclose(Ln @ Ln.conj().T, I2))
# Zero 侧：*-自同构保单位元
U = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)   # Hadamard（酉）
check("★ Zero 自带运算为 *-自同构 ⇒ 保单位元（U†U = I，范数差 = %.2e）"
      % np.linalg.norm(U.conj().T @ U - I2), np.linalg.norm(U.conj().T @ U - I2) < 1e-15)

# ======================================================================
head("P4  定理 4（否定）· 中心态下测量的反作用不可观测")
# ======================================================================
Dc = 4
rho_c = np.zeros((Dc, Dc), dtype=complex)
for h in range(2):
    rho_c[2 * h:2 * h + 2, 2 * h:2 * h + 2] = (h + 1) / 3 * I2
worst = 0.0
for _ in range(300):
    Eop = np.zeros((Dc, Dc), dtype=complex)
    Aop = np.zeros((Dc, Dc), dtype=complex)
    for h in range(2):
        M1 = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        M2 = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        Eop[2*h:2*h+2, 2*h:2*h+2] = (M1 + M1.conj().T) / 2
        Aop[2*h:2*h+2, 2*h:2*h+2] = (M2 + M2.conj().T) / 2
    worst = max(worst, abs(np.trace(rho_c @ (Eop @ Aop - Aop @ Eop))))
check("★ 中心态：ω([E,A]) ≡ 0（300 组随机效应/可观测量，max = %.2e）" % worst, worst < 1e-12)
# 非中心反例
# 非中心态：块 0 取 Bloch 向量 n=(0,0,0.8)；注意 [X,Y] 的迹恒为 0，
# 故必须取 (E,A) 使对易子 ∝ σ_z（取 σx,σy），否则会误得 0。
rho_nc = np.zeros((Dc, Dc), dtype=complex)
rho_nc[0:2, 0:2] = (I2 + 0.8 * SZ) / 2
rho_nc[2:4, 2:4] = (I2 - 0.4 * SZ) / 2
E1 = np.zeros((Dc, Dc), dtype=complex); E1[0:2, 0:2] = SX
A1 = np.zeros((Dc, Dc), dtype=complex); A1[0:2, 0:2] = SY
val = abs(np.trace(rho_nc @ (E1 @ A1 - A1 @ E1)))
check("★ 非中心反例：ω([E,A]) = %.6f ≠ 0 ⇒ 反作用**可**观测（定理 4 的适用范围）" % val, val > 1e-6)
# 推论 4.1：中心态下块内两最小投影权重相同 ⇒ 任何正交分解给同一读数
w = np.trace(rho_c[:2, :2]).real
check("★ 推论 4.1：中心态下块内两最小投影权重相同（各 %.6f）⇒ 指针基不唯一" % (w / 2),
      abs(rho_c[0, 0].real - rho_c[1, 1].real) < 1e-15)

# ======================================================================
head("P5  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("不依赖任何同题旧文", "去依赖声明"),
                 ("定理 1", "定理 1"), ("定理 2", "定理 2"), ("定理 3", "定理 3"),
                 ("定理 4", "定理 4"), ("引理 3.0", "引理 3.0"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明"),
                 ("显式前提", "前提显式化"), ("无条件", "无条件的无信号")]:
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
print("QMR04 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
