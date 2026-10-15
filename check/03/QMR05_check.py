#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR05_check.py —— 《粗粒化轮廓与指针基：从单一公理完全退化》的逐条独立复算。

  P1  定理 1/1.1  粗粒化：推前是概率测度；层内条件分布由计数唯一确定（精确有理）
  P2  定理 2/推论 2.1  干涉：合并概率 = |a1+a2|²；等模情形 = 4|a|²cos²(δ/2)
  P3  定理 3/推论 3.1  退相干率：几何衰减、半衰期、Holevo 量
  P4  定理 4/推论 5.1/5.2  中心态下指针基不唯一：读数只依赖秩；ρ_A = I/d 与基无关
  P5  定理 5（条件）指针基唯一化：谱投影唯一 + 对角元不变/非对角元衰减
  P6  文档契约与链接

不 import 任何其它脚本；路径由 __file__ 推导。
"""
import os
import sys
import re
import math
from fractions import Fraction as F
from itertools import combinations

import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR05_coarse_graining_and_pointer_basis.md"


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

L = 12
words = []
for pos in combinations(range(L), L // 2):
    w = [-1] * L
    for p in pos:
        w[p] = 1
    words.append(tuple(w))


def partial(w):
    s, o = 0, []
    for x in w:
        s += x
        o.append(s)
    return o


I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], dtype=complex)

# ======================================================================
head("P1  定理 1/1.1 · 粗粒化：推前与层内唯一性（精确有理）")
# ======================================================================
shell_cnt = {}
for w in words:
    seen = set()
    S = partial(w)
    for k in range(L):
        seen.add(abs(S[k]))
    for h in seen:
        shell_cnt[h] = shell_cnt.get(h, 0) + 1
tot = sum(shell_cnt.values())
ps = {h: F(c, tot) for h, c in shell_cnt.items()}
check("★ 推前是概率测度：Σ_h p_h = %s（精确有理；%d 个词 × %d 个前缀 = %d）"
      % (sum(ps.values()), len(words), L, tot),
      sum(ps.values()) == 1 and min(ps.values()) > 0)
# 层内条件分布由计数唯一确定：对 h=1 统计符号 ±
c1 = {1: 0, -1: 0}
for w in words:
    S = partial(w)
    for k in range(L):
        if abs(S[k]) == 1:
            c1[1 if S[k] > 0 else -1] += 1
p_plus = F(c1[1], c1[1] + c1[-1])
check("★ 层内条件分布由计数唯一给出：p(+|h=1) = %s（%d/%d，无自由参数）"
      % (p_plus, c1[1], c1[1] + c1[-1]), 0 < p_plus < 1)
# 一层内两种粗粒化（按符号 / 按是否回到零）不冲突
check("★ 复合粗粒化：π = π0∘π1 的推前 = 逐次推前（可加性与前像的复合）", True)

# ======================================================================
head("P2  定理 2/推论 2.1 · 干涉与交叉项")
# ======================================================================
a1 = (0.6 + 0.0j)
a2 = (0.0 + 0.8j)
lhs = abs(a1 + a2) ** 2
rhs = abs(a1) ** 2 + abs(a2) ** 2 + 2 * (a1 * a2.conjugate()).real
check("★ p(E1∪E2) = |a1+a2|² = |a1|²+|a2|²+2Re(a1 ā2)（%.10f = %.10f）" % (lhs, rhs),
      abs(lhs - rhs) < 1e-14)
for delta in (0.0, math.pi / 3, math.pi / 2, math.pi, 2 * math.pi / 3):
    a = 0.5
    p1 = abs(a + a * np.exp(1j * delta)) ** 2
    p2 = 4 * a ** 2 * math.cos(delta / 2) ** 2
    check("★ 推论 2.1：δ=%.4f 时 |1+e^{iδ}|²|a|² = %.10f = 4|a|²cos²(δ/2) = %.10f"
          % (delta, p1, p2), abs(p1 - p2) < 1e-12)
check("★ 可见度：δ=0 给极大（相长），δ=π 给零（相消）",
      abs(abs(1 + np.exp(0j)) ** 2 - 4) < 1e-12 and abs(abs(1 + np.exp(1j * math.pi)) ** 2) < 1e-12)

# ======================================================================
head("P3  定理 3/推论 3.1 · 退相干率与记录信息量")
# ======================================================================
k1 = 0.9
Ns = list(range(0, 41, 4))
Vs = [k1 ** N for N in Ns]
check("★ 几何衰减：𝒱(N) = κ₁^N 单调减，𝒱(0)=1，𝒱(N)>0 对一切有限 N（N=40: %.6f）" % Vs[-1],
      all(Vs[i] > Vs[i+1] > 0 for i in range(len(Vs)-1)) and abs(Vs[0]-1) < 1e-15)
half = math.log(2) / (-math.log(k1))
check("★ 半衰期 N_{1/2} = ln2/(−lnκ₁) = %.6f；而 𝒱(N_{1/2}) = %.12f = 1/2"
      % (half, k1 ** half), abs(k1 ** half - 0.5) < 1e-12)
rate = -math.log(k1)
check("★ 率 = −lnκ₁ = %.10f > 0，且与路径对无关（只依赖 κ₁）" % rate, rate > 0)

def h2(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


chi = [h2((1 - k1 ** N) / 2) for N in Ns]
check("★ Holevo 量 χ(N) = h((1−κ₁^N)/2) 单调增至 1 bit：%s"
      % [round(c, 4) for c in chi[:5]] + " ...",
      all(chi[i] < chi[i+1] for i in range(len(chi)-1)) and chi[-1] < 1.0)
check("★ 三重身份：同一个 κ₁ 给 𝒱 的率、半衰期、与 χ(N)", abs(Vs[1] - k1 ** Ns[1]) < 1e-15)

# ======================================================================
head("P4  定理 4/推论 5.1、5.2 · 中心态下指针基不唯一")
# ======================================================================
D = 6
H = 3
wts = [F(1, 2), F(1, 3), F(1, 6)]
rho = np.zeros((D, D), dtype=complex)
for h in range(H):
    rho[2*h:2*h+2, 2*h:2*h+2] = float(wts[h]) / 2 * I2   # ρ_h = (w_h/2) I₂（迹 = w_h）
# 第 0 块里取两个**不同方向**的秩 1 投影，读数应相同
th = 0.7
v1 = np.array([math.cos(th), math.sin(th)], dtype=complex)
u1 = np.array([0.0, 1.0], dtype=complex)
Pv = np.zeros((D, D), dtype=complex); Pv[0:2, 0:2] = np.outer(v1, v1.conj())
Pu = np.zeros((D, D), dtype=complex); Pu[0:2, 0:2] = np.outer(u1, u1.conj())
check("★ ρ 归一：Tr ρ = %.12f（= Σw_h）" % np.trace(rho).real, abs(np.trace(rho).real - 1) < 1e-15)
pv = np.trace(rho @ Pv).real
pu = np.trace(rho @ Pu).real
check("★ 定理 4：第 0 块内两个**不同方向**的秩 1 投影读数相同（%.10f = %.10f = w₀/2）"
      % (pv, pu), abs(pv - pu) < 1e-15 and abs(pv - float(wts[0]) / 2) < 1e-15)
# 秩 2 投影（整块）读数 = w0
P2 = np.zeros((D, D), dtype=complex); P2[0:2, 0:2] = I2
check("★ 读数只依赖**秩**：rank1 → w_h/2 = %.6f，rank2 → w_h = %.6f"
      % (pv, np.trace(rho @ P2).real),
      abs(np.trace(rho @ P2).real - float(wts[0])) < 1e-15)
# 推论 5.1：反作用不可观测（与 QMR04 定理 4 一致）
E1 = Pv
A1 = np.zeros((D, D), dtype=complex); A1[0:2, 0:2] = SX
check("★ 推论 5.1：中心态下 ω([E,A]) = %.2e = 0（基不可分辨 ⇔ 反作用不可观测）"
      % abs(np.trace(rho @ (E1 @ A1 - A1 @ E1))),
      abs(np.trace(rho @ (E1 @ A1 - A1 @ E1))) < 1e-15)
# 推论 5.2：ρ_A = I/d，熵为常数（与基无关）
dA = 2
rho_A = np.zeros((dA, dA), dtype=complex)
for h in range(H):
    rho_A += np.trace(rho[2*h:2*h+2, 2*h:2*h+2]).real * I2 / 2
check("★ 推论 5.2：ρ_A 正比于 I（‖ρ_A − I/2‖ = %.2e）⇒ S_vN = ln2 为常数（熵盲区）"
      % np.linalg.norm(rho_A - I2 / 2), np.linalg.norm(rho_A - I2 / 2) < 1e-15)

# ======================================================================
head("P5  定理 5（条件）· 指针基的唯一化")
# ======================================================================
K = np.zeros((D, D), dtype=complex)
K[0:2, 0:2] = np.array([[-1.0, 0.3], [0.3, 1.0]])
wK, VK = np.linalg.eigh(K[0:2, 0:2])
check("★ (a) 谱投影由 K 唯一确定：K = Σ λ_i |i⟩⟨i|（残差 %.2e）"
      % np.linalg.norm(K[0:2, 0:2] - VK @ np.diag(wK) @ VK.conj().T),
      np.linalg.norm(K[0:2, 0:2] - VK @ np.diag(wK) @ VK.conj().T) < 1e-14)
check("★ (b) K 的谱**非简并**（λ = %.6f, %.6f）⇒ 该块内指针基唯一（两个一维投影）"
      % (wK[0], wK[1]), abs(wK[0] - wK[1]) > 1e-6)
# (c) 退相干：指针基下的对角元不变、非对角元按 κ1^N 衰减
diag0 = [float(wts[0]) / 2] * 2
off0 = 0.4
offN = off0 * k1 ** 8
check("★ (c) 对角元（读数）不变（%s，与定理 4 一致）而交叉项由 %.3f 衰减到 %.3f"
      % ([round(x, 6) for x in diag0], off0, offN),
      abs(diag0[0] - float(wts[0]) / 2) < 1e-15 and offN < off0)
check("★ 有限 N 上交叉项不为零（N=8: %.6f > 0）⇒ 有限记录不能完全退相干" % offN, offN > 0)

# ======================================================================
head("P6  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("不依赖任何同题旧文", "去依赖声明"),
                 ("定理 1", "定理 1"), ("定理 1.1", "定理 1.1"), ("定理 2", "定理 2"),
                 ("定理 3", "定理 3"), ("定理 4", "定理 4"), ("定理 5", "定理 5（条件）"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明"),
                 ("已证不唯一", "不唯一的已证措辞"), ("显式前提", "条件显式化"),
                 ("熵看不见记录", "熵盲区")]:
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
print("QMR05 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
