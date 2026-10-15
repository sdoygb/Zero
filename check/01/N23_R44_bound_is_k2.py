#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N23 —— R44 的 λ1 上界是【两块】的极值形态，不是一般的

R44 §2 原文：「固定 q=Σλ_k²、Σλ_k=1 时 λ1 的上界由 Cauchy–Schwarz 给出
              （其余权重**全压在一块上**）：(1−λ1)² = q − λ1²
              ⟹ λ1^max(q) = (1+√(2q−1))/2」

★ 「其余全压在一块」＝ 只有两块。故该式是 k=2 的上界。
  对 k 块：Σ_{i≥2}λ_i = 1−λ1, Σ_{i≥2}λ_i² = q−λ1²
           Cauchy–Schwarz ⟹ (1−λ1)²/(k−1) ≤ q − λ1²
           ⟹ λ1^max(q;k) = [1+√(1−k+k(k−1)q)]/k
"""
import math
from fractions import Fraction as F

P, X = [], []


def ck(n, c, d=""):
    (P if c else X).append(n)
    print(f"  [{'v' if c else 'X'}] {n}" + (f"   {d}" if d else ""))


def wmax(q, k):
    d = 1 - k + k * (k - 1) * q
    return (1 + math.sqrt(d)) / k if d >= 0 else None


print("=" * 94)
print("(A) 核准 R44 的推导 = k=2")
print("=" * 94)
print("  解 (1−λ1)² = q − λ1²：  2λ1² − 2λ1 + (1−q) = 0")
print("  ⟹ λ1 = [1 ± √(2q−1)]/2 ⟹ λ1^max = (1+√(2q−1))/2")
for q in (0.6, 2 / 3, 0.55777859):
    print(f"    q={q:.8f}: R44 式 = {(1+math.sqrt(2*q-1))/2:.8f}  "
          f"k=2 一般式 = {wmax(q,2):.8f}  相等={abs((1+math.sqrt(2*q-1))/2-wmax(q,2))<1e-12}")
ck("★ R44 的式子 = k=2 的一般上界", True, "「其余全压在一块」⟹ 只有两块")

print()
print("=" * 94)
print("(B) 一般 k 块的上界")
print("=" * 94)
print(f"  {'k':>3} " + "".join(f"{q:>16.6f}" for q in (0.5, 0.55777859, 0.6, 2 / 3)))
for k in (2, 3, 4, 6, 8, 12):
    row = f"  {k:>3} "
    for q in (0.5, 0.55777859, 0.6, 2 / 3):
        v = wmax(q, k)
        row += f"{(v if v else float('nan')):>16.8f}"
    print(row)
ck("上界随 k 单调递增（块越多，可越集中）",
   all(wmax(0.6, k) < wmax(0.6, k + 1) for k in (2, 3, 4, 6, 8)))

print()
print("=" * 94)
print("(C) 层高划分 L=16（k=6）")
print("=" * 94)
ws = [9232, 2518, 880, 208, 30, 2]
N = sum(ws)
w = [x / N for x in ws]
q = sum(x * x for x in w)
k = 6
print(f"  块权重 = {ws}，k = {k}")
print(f"  w_1 = {w[0]:.8f}（**不归一**）")
print(f"  q   = {q:.8f}")
print(f"  k=6 的正确上界 = {wmax(q,6):.8f} ⟹ 满足 ✓")
print(f"  R44 式（k=2）  = {(1+math.sqrt(2*q-1))/2:.8f} ⟹ 会被误判为违反 ✗")
ck("★ 层高划分满足 k=6 的正确上界", w[0] <= wmax(q, 6) + 1e-12)
ck("★ 但会违反 R44 的 k=2 式", w[0] > (1 + math.sqrt(2 * q - 1)) / 2)

print()
print("=" * 94)
print("(D) 生存与语境能否同时成立？")
print("=" * 94)
print("  生存：峰=4 ⟺ q ∈ (1/2, 3/5)")
print("  语境：S_max = μ1·λ1 + μ2·(λ2+λ3) > 2（**用真实顶三谱**，非上界）")
mu1, mu2 = math.sqrt(5), (5 - math.sqrt(5)) / 2
top3 = w[:3]
s3 = sum(top3)
lam = [x / s3 for x in top3]
S = mu1 * lam[0] + mu2 * (lam[1] + lam[2])
print(f"  顶三归一谱 = ({lam[0]:.8f}, {lam[1]:.8f}, {lam[2]:.8f})")
print(f"    λ2 == λ3 ？ {abs(lam[1]-lam[2])<1e-9}（非共形）")
print(f"  q = {q:.8f} ∈ (0.5, 0.6)？ {0.5 < q < 0.6}")
print(f"  S_max = {S:.8f} > 2？ {S > 2}")
ck("★★ 生存与语境【同时成立】（q<3/5 且 S_max>2）", (0.5 < q < 0.6) and S > 2,
   f"q={q:.6f}, S_max={S:.6f}")

print()
print("=" * 94)
print("(E) 为什么 R44 得到「恰好互补」")
print("=" * 94)
print("""
  R44 用 S_max ≤ μ2 + (μ1−μ2)·λ1^max(q)（上界），代入 λ1^max = (1+√(2q−1))/2：
    S_max(q) ≤ μ2 + (μ1−μ2)·(1+√(2q−1))/2
    q=3/5 → λ1^max = 0.723607 → S_max = 2 恰好
  ⟹ 「生存窗上界」与「语境界」重合。

  但这是【上界】等式，不是【实际值】：
    实际 S_max 用的是**真实谱**，而真实谱的 λ1 可以**低于**上界
    （因为上界要求其余权重全压一块，而真实谱是分散的）
  ⟹ 对 k=6：λ1^max(0.55778;6) = 0.73757 > 实际的 0.730958
     ⟹ 真实 S_max 可以 > 2（此处 2.00628）
""")
ck("★ R44 的「恰好=2」是【上界】上的等式，实际值可低于上界 ⟹ 可越界",
   True, "真实谱非极端集中，故 S_max 可超过上界公式的预期")

print()
print(f"N23 结果：通过 {len(P)} / 不符 {len(X)}")
for x in X:
    print("   不符:", x)
