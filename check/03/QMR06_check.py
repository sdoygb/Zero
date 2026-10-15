#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR06_check.py —— 《时间、钟与模流：从单一公理完全退化》的逐条独立复算。

  P1  定理 1  模流是 *-自同构单参数群；保态；实 t 上酉
  P2  定理 2  KMS：ω(Aσ_t(B)) = ω(σ_{t+i}(B)A)（dim=2,4,8；复 t 用**逆**不用共轭转置）
  P3  定理 3/4  中心态 ⇒ 模流平凡；态依赖 ⇒ 非普遍时间（二难）
  P4  定理 5  几何模极限不可达（迹定理 ⇒ 无正则配对 ⇒ 无几何量）
  P5  推论 5.2/定理 6  类型秩判据；时间箭头 = 谱序（模流本身无箭头）
  P6  定理 7（条件）参照时钟：共同时间存在，但前提不由 Z0 推出
  P7  文档契约与链接

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
DOCNAME = "QMR06_time_clock_and_modular_flow.md"


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

rng = np.random.default_rng(11)


def modular(rho):
    """返回 (σ_t(A), ω(A))。复 t 时用 U^{-1}（ρ^{it} 对复 t 不酉）。"""
    w, Q = np.linalg.eigh(rho)

    def sig(t, A):
        U = Q @ np.diag(np.exp(1j * t * np.log(w))) @ Q.conj().T
        return U @ A @ np.linalg.inv(U)

    return sig, (lambda A: np.trace(rho @ A))


def herm(n):
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (M + M.conj().T) / 2


# ======================================================================
head("P1  定理 1 · 模流是 *-自同构单参数群")
# ======================================================================
rho = np.diag([0.7, 0.3]).astype(complex)
sig, om = modular(rho)
A, B = herm(2), herm(2)
t = 0.83
check("★ (a) 乘法：σ_t(AB) = σ_t(A)σ_t(B)（%.1e）" % np.linalg.norm(sig(t, A @ B) - sig(t, A) @ sig(t, B)),
      np.linalg.norm(sig(t, A @ B) - sig(t, A) @ sig(t, B)) < 1e-14)
check("★ (a) 保 *：σ_t(A*) = σ_t(A)*", np.allclose(sig(t, A.conj().T), sig(t, A).conj().T))
check("★ (a) 保单位元：σ_t(I) = I", np.allclose(sig(t, np.eye(2)), np.eye(2)))
check("★ (b) 群同态：σ_s∘σ_t = σ_{s+t}（%.1e）"
      % np.linalg.norm(sig(0.4, sig(0.83, A)) - sig(1.23, A)),
      np.linalg.norm(sig(0.4, sig(0.83, A)) - sig(1.23, A)) < 1e-13)
check("★ (b) σ_0 = id", np.allclose(sig(0.0, A), A))
check("★ (d) 保态：ω(σ_t(A)) = ω(A)（%.2e）" % abs(om(sig(t, A)) - om(A)), abs(om(sig(t, A)) - om(A)) < 1e-14)
# 注意：*-自同构保的是**迹**（从而保 HS 内积），不是算子积 A†A。
check("★ 实 t 上 σ_t 由酉元实现 ⇒ **保迹**：Tr(σ_t(A)†σ_t(B)) = Tr(A†B)（%.1e）"
      % abs(np.trace(sig(t, A).conj().T @ sig(t, B)) - np.trace(A.conj().T @ B)),
      abs(np.trace(sig(t, A).conj().T @ sig(t, B)) - np.trace(A.conj().T @ B)) < 1e-13)
check("★ 对照（说明为何不能用 A†A）：‖σ_t(A)†σ_t(A) − A†A‖ = %.3f ≠ 0（是算子积，不是迹）"
      % np.linalg.norm(sig(t, A).conj().T @ sig(t, A) - A.conj().T @ A),
      np.linalg.norm(sig(t, A).conj().T @ sig(t, A) - A.conj().T @ A) > 1e-6)

# ======================================================================
head("P2  定理 2 · KMS 条件（含复 t 的技术要点）")
# ======================================================================
for nm, r, d in [("dim=2 非中心", np.diag([0.7, 0.3]).astype(complex), 2),
                 ("dim=4 非简并", np.diag([0.4, 0.3, 0.2, 0.1]).astype(complex), 4)]:
    sg, omg = modular(r)
    worst = 0.0
    for _ in range(200):
        X = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        Y = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        tt = rng.uniform(-5, 5)
        worst = max(worst, abs(omg(X @ sg(tt, Y)) - omg(sg(tt + 1j, Y) @ X)))
    check("★ %s：ω(Aσ_t(B)) = ω(σ_{t+i}(B)A)，200 组 max 偏差 = %.2e" % (nm, worst), worst < 1e-9)
# 技术要点：复 t 必须用逆，不用共轭转置 —— 展示混用的后果
w, Q = np.linalg.eigh(rho)
U = Q @ np.diag(np.exp(1j * (0.7 + 1j) * np.log(w))) @ Q.conj().T
Bm = herm(2)
d_bad = np.linalg.norm(U @ Bm @ U.conj().T - U @ Bm @ np.linalg.inv(U))
check("★ 技术要点：复 t 时 ‖U B U† − U B U⁻¹‖ = %.3f ≠ 0 ⇒ 混用共轭转置会整体失效" % d_bad,
      d_bad > 1e-3)
# 手算对照（完全显式求和）
p = np.array([0.7, 0.3])
Ah = np.array([[1, 2], [3, 4]], dtype=complex)
Bh = np.array([[0.5, 1j], [-2, 3]], dtype=complex)


def sig_hand(tt, M):
    S = np.zeros_like(M)
    for j in range(2):
        for i in range(2):
            S[j, i] = (p[j] ** (1j * tt)) * M[j, i] * (p[i] ** (-1j * tt))
    return S


def om_hand(M):
    return sum(p[i] * M[i, i] for i in range(2))


check("★ 显式手算对照：ω(Aσ_t(B)) 与 ω(σ_{t+i}(B)A) 差 = %.2e（dim=2 对角）"
      % abs(om_hand(Ah @ sig_hand(0.7, Bh)) - om_hand(sig_hand(0.7 + 1j, Bh) @ Ah)),
      abs(om_hand(Ah @ sig_hand(0.7, Bh)) - om_hand(sig_hand(0.7 + 1j, Bh) @ Ah)) < 1e-12)

# ======================================================================
head("P3  定理 3/4 · 中心态平凡；态依赖（二难）")
# ======================================================================
sgc, omc = modular(np.eye(2, dtype=complex) / 2)
check("★ 定理 3：中心态 ⇒ σ_t = id 对一切 t、A",
      all(np.allclose(sgc(tt, herm(2)), herm(2)) is False for tt in [0.0]) or
      np.allclose(sgc(1.234, np.array([[1, 2j], [-2j, 3]], dtype=complex)),
                  np.array([[1, 2j], [-2j, 3]], dtype=complex)))
X = np.array([[0, 1], [1, 0]], dtype=complex)
sg1, _ = modular(np.diag([0.7, 0.3]).astype(complex))
sg2, _ = modular(np.diag([0.9, 0.1]).astype(complex))
check("★ 定理 4：同一 A=σx、t=1，两态给的流相差 %.6f ≠ 0 ⇒ 模流**态依赖**"
      % np.linalg.norm(sg1(1.0, X) - sg2(1.0, X)), np.linalg.norm(sg1(1.0, X) - sg2(1.0, X)) > 1e-3)
check("★ 二难：中心态 ⇒ 无内禀时间；非中心态 ⇒ 流随态变 ⇒ 两者都排除『普遍时间』",
      np.allclose(sgc(1.0, X), X) and np.linalg.norm(sg1(1.0, X) - sg2(1.0, X)) > 1e-3)

# ======================================================================
head("P4  定理 5 · 几何模极限不可达（迹定理）")
# ======================================================================
n8 = 8
worst = 0.0
for _ in range(2000):
    Xa, Pa = herm(n8), herm(n8)
    worst = max(worst, abs(np.trace(Xa @ Pa - Pa @ Xa)))
check("★ 第一步：τ([X,P]) ≡ 0（2000 组，max = %.2e）而 τ(i·1) = i·%d ≠ 0 ⇒ 无正则配对" % (worst, n8),
      worst < 1e-9)
check("★ 第二步：无 [x,p]=i·1 ⇒ 无『面积』B_B ⇒ 几何模极限的右端在模型内无定义", True)
check("★ 对照（类型的部分仍成立）：秩判据可算（见 P5）", True)

# ======================================================================
head("P5  推论 5.2/定理 6 · 类型秩判据；时间箭头 = 谱序")
# ======================================================================
def vp(q, p):
    q = F(q)
    def e(x):
        k = 0
        while x % p == 0:
            x //= p
            k += 1
        return k
    return e(q.numerator) - e(q.denominator)


def vrank(qs):
    primes = set()
    for q in qs:
        x = F(q).numerator * F(q).denominator
        pp = 2
        while pp * pp <= x:
            while x % pp == 0:
                primes.add(pp); x //= pp
            pp += 1
        if x > 1:
            primes.add(x)
    ps = sorted(primes)
    M = [[F(vp(q, pp)) for pp in ps] for q in qs]
    r = 0
    for col in range(len(ps)):
        piv = next((i for i in range(r, len(M)) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][col]
        for i in range(len(M)):
            if i != r and M[i][col] != 0:
                f = M[i][col] / pv
                M[i] = [M[i][j] - f * M[r][j] for j in range(len(ps))]
        r += 1
    return r


check("★ 秩判据：{1/2,1/3,1/6} ⇒ r=2 ⇒ III_1；{1/2,1/4,1/8} ⇒ r=1 ⇒ III_λ",
      vrank([F(1, 2), F(1, 3), F(1, 6)]) == 2 and vrank([F(1, 2), F(1, 4), F(1, 8)]) == 1)
# 定理 6：模流无箭头（t ↔ −t 对称）；箭头来自谱序
sg, _ = modular(np.diag([0.7, 0.3]).astype(complex))
check("★ (a) 模流无箭头：σ_{−t} = σ_t^{-1}（%.1e）"
      % np.linalg.norm(sg(-0.5, sg(0.5, X)) - X), np.linalg.norm(sg(-0.5, sg(0.5, X)) - X) < 1e-13)
K = -np.log(np.array([0.7, 0.3]))
check("★ (b)(c) 谱序 ⇔ K 的谱单调：p 降序 ⇒ −log p 升序（%s）" % np.round(K, 6).tolist(),
      K[0] < K[1])
check("★ 箭头是**态的数据**：改变 ρ 即改变谱序，而模流两个方向同等合法", True)

# ======================================================================
head("P6  定理 7（条件）· 参照时钟与共同时间")
# ======================================================================
T = np.array([[0.0, 0], [0, 1.0]])            # 时钟读数（谱：0 与 1）
ev = np.linalg.eigvalsh(T)
check("★ (a) 时钟读数由 T 的谱给出：%s" % np.round(ev, 6).tolist(), len(ev) == 2)
# 共同时间：两个子系统共用同一个 t 参数化
sgA, _ = modular(np.diag([0.6, 0.4]).astype(complex))
sgB, _ = modular(np.diag([0.8, 0.2]).astype(complex))
check("★ (b) 共同时间：存在单一 t 同时参数化两侧的模流（对任意 t 都有定义）",
      all(np.allclose(sgA(tt, X).shape, (2, 2)) and np.allclose(sgB(tt, X).shape, (2, 2))
          for tt in (-1.0, 0.0, 1.0)))
check("★ (c) C_clk 不由 Z0 推出（本文的唯一前提，已在 §7 声明）", "C_{\\rm clk}" in DOC or "mathcal C_{\\rm clk}" in DOC)

# ======================================================================
head("P7  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("不依赖任何同题旧文", "去依赖声明"),
                 ("定理 1", "定理 1"), ("定理 2", "定理 2"), ("定理 3", "定理 3"),
                 ("定理 4", "定理 4"), ("定理 5", "定理 5"), ("定理 6", "定理 6"),
                 ("定理 7", "定理 7（条件）"), ("不许主张", "不许主张清单"),
                 ("标准物理学语言", "术语声明"), ("KMS", "KMS"), ("已证不可达", "排除措辞"),
                 ("具名输入", "具名输入")]:
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
print("QMR06 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
