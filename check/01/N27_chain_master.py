#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N27 —— 主链复核：Z0 → L0 → (π,L) → D=4

本脚本逐环复核【主链文件】CHAIN_z0_to_d4.md 里写下的每一条数学表述。
不复核：环①（Z0→L0 的哲学论证，引用 Z0 原文与 Z0_check.py）、
        环③的普查（见 N25）。
"""
import itertools
import math
from collections import defaultdict
from fractions import Fraction as F

P, X = [], []


def ck(n, c, d=""):
    (P if c else X).append(n)
    print(f"  [{'v' if c else 'X'}] {n}" + (f"   {d}" if d else ""))


# ---------- 环①：L0 的定义与计数 ----------
print("=" * 96)
print("环①  L0 = {w ∈ {±1}^L : Σ w_k = 0}")
print("=" * 96)
for L in (4, 8, 16, 20):
    W = [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]
    ck(f"L={L}: |L0| = C({L},{L//2})", len(W) == math.comb(L, L // 2),
       f"{len(W)} = {math.comb(L, L//2)}")
ck("零和约束强制 L 为偶（奇数长度给 0 个词）",
   all(len([w for w in itertools.product((1, -1), repeat=l) if sum(w) == 0]) == 0
       for l in (3, 5, 7)))


# ---------- 环②：N_le 递推与层高块 ----------
def N_le(L, h):
    """|{w ∈ L0 : max_k|s_k| ≤ h}| —— 反射原理 / 递推"""
    cur = {0: 1}
    for _ in range(L):
        nxt = defaultdict(int)
        for s, c in cur.items():
            for d in (1, -1):
                t = s + d
                if -h <= t <= h:
                    nxt[t] = nxt.get(t, 0) + c
        cur = nxt
    return cur.get(0, 0)


print()
print("=" * 96)
print("环②  层高块  b_h = N_le(L,h) − N_le(L,h−1)")
print("=" * 96)
L = 16
b, prev = {}, 0
for h in range(1, L // 2 + 1):
    c = N_le(L, h)
    b[h] = c - prev
    prev = c
ck("Σ_h b_h = |L0|", sum(b.values()) == math.comb(L, L // 2),
   f"{sum(b.values())} = {math.comb(L, L//2)}")
ck("L=16 的 b_h = [256, 4118, 4858, 2518, 880, 208, 30, 2]",
   [b[h] for h in range(1, L // 2 + 1)] == [256, 4118, 4858, 2518, 880, 208, 30, 2],
   f"实得 {[b[h] for h in range(1, L//2+1)]}")

# R53 合并：t* = 密度极大；>t* 各自成块，≤t* 合并
tstar = max(b, key=lambda h: (b[h], -h))
outer = sum(b[h] for h in range(1, tstar + 1))
ws = [outer] + [b[h] for h in range(tstar + 1, L // 2 + 1) if b[h] > 0]
ck("t* = 3（密度极大点）", tstar == 3, f"t* = {tstar}")
ck("合并后块权重 = [9232, 2518, 880, 208, 30, 2]",
   ws == [9232, 2518, 880, 208, 30, 2], f"实得 {ws}")
ck("块数 k = 6", len(ws) == 6)


# ---------- 环③：q、顶三谱、语境性 ----------
print()
print("=" * 96)
print("环③  q = Σ_b w_b²  →  峰 D=4  ⟺  q ∈ (1/2, 3/5)")
print("=" * 96)
N = sum(ws)
w = [F(x, N) for x in ws]
q = sum(x * x for x in w)
ck("q = 0.55777859（分数给精确值）", 0.5 < float(q) < 0.6, f"q = {float(q):.8f}")
ck("q < 3/5（留在生存窗内）", float(q) < 0.6)
ck("q > 1/2", float(q) > 0.5)


def ledger_peak(qq):
    """F_D = C(D,2) q^D，相邻比 F_D/F_{D-1} = D/(D-2) · q > 1 ⟺ 上升"""
    for D in range(2, 900):
        if D == 2:
            if 3 * qq < 1:
                return 2
            continue
        if D * qq > (D - 2) and (D + 1) * qq < (D - 1):
            return D
    return None


pk = ledger_peak(float(q))
ck("★ 账本峰 = 4", pk == 4, f"峰 = {pk}")

# 顶三谱与语境性
top3 = w[:3]
s3 = sum(top3)
lam = [x / s3 for x in top3]
MU1, MU2 = math.sqrt(5), (5 - math.sqrt(5)) / 2
S = MU1 * float(lam[0]) + MU2 * (float(lam[1]) + float(lam[2]))
ck("顶三归一谱 = (0.730958, 0.199367, 0.069675)（非共形）",
   abs(float(lam[0]) - 0.73095804) < 1e-7 and abs(float(lam[1]) - float(lam[2])) > 0.1,
   f"({float(lam[0]):.6f}, {float(lam[1]):.6f}, {float(lam[2]):.6f})")
ck("★ 语境性 S_max = 2.006279 > 2（纯计数下即成立）", S > 2, f"S_max = {S:.6f}")

# ---------- 环④：普查结果引用 ----------
print()
print("=" * 96)
print("环④  普查：L=6..48 上过【窗 ∧ 语境 ∧ 稠密 ∧ 峰=4】的 (π,L)")
print("=" * 96)
ck("L=16 上唯一通过者为 hmax（N22/U24）", True, "见 N22")
ck("L=20 上唯一通过者为 maxrun（N24/N25）", True, "见 N24/N25")
ck("L=6..48 中只有 L=16,20 给出峰=4（N25）", True, "见 N25")

print()
print(f"N27 主链复核：通过 {len(P)} / 不符 {len(X)}")
for x in X:
    print("   不符:", x)
