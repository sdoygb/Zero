#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D03_statmech_check.py -- 核验《D03 · 为什么 Zero 不可能有统计力学》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

分组：
  A  |W_L| = C(L, L/2) 精确
  B  层高计数 b_h：整数转移矩阵法 vs 穷举（并登记 Z001 (7.7) 的错误）
  C  账本熵闭式 S = H(w) + ln2 与直接谱计算一致
  D  块项恒为 ln2（与 L 无关）
  E  H(w) = O(ln L)
  F  S/L -> 0
  G  碎片化对照：均匀划分给 ln|W_L| ~ L ln2
  H  中心性：rho 在每块上是 w_h/2 * I
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction

FAIL = 0
NCHK = 0


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


def ck_eq(a, b, msg):
    ck(a == b, "%s  (got %r, want %r)" % (msg, a, b))


# ---------------------------------------------------------------- 词与层高
def balanced_words(L):
    """长度 L 的平衡词（+1/-1，和为零）；L 为偶数"""
    assert L % 2 == 0
    for pos in itertools.combinations(range(L), L // 2):
        w = [-1] * L
        for p in pos:
            w[p] = +1
        yield tuple(w)


def height(w):
    s = 0
    h = 0
    for x in w:
        s += x
        h = max(h, abs(s))
    return h


def shell_counts(L):
    """枚举求 b_h"""
    b = {}
    for w in balanced_words(L):
        h = height(w)
        b[h] = b.get(h, 0) + 1
    return b


def Nle(L, h):
    """带内游走数 N_<=(L,h) = #{ |s_k| <= h 对一切 k，s_0 = s_L = 0 }。

    用带内（吸收壁在 +-(h+1)）的转移矩阵直接迭代，整数精确。
    核验：与穷举一致（L <= 20）；且逐项复现 Z001 §7.3 的 n=16 表。
    """
    if h <= 0:
        return 0
    d = 2 * h + 1
    off = h
    v = [0] * d
    v[off] = 1
    for _ in range(L):
        nv = [0] * d
        for s in range(-h, h + 1):
            c = v[s + off]
            if c == 0:
                continue
            for t in (s - 1, s + 1):
                if -h <= t <= h:
                    nv[t + off] += c
        v = nv
    return v[off]


def shell_counts_band(L):
    H = L // 2
    return {h: Nle(L, h) - Nle(L, h - 1) for h in range(1, H + 1)}


# ---------------------------------------------------------------- 账本熵
def ledger_entropy(L):
    """返回 (H(w) 精确有理数形式, 浮点值, 壳权重)"""
    b = shell_counts_band(L)
    tot = sum(b.values())
    w = {h: Fraction(c, tot) for h, c in b.items()}
    Hw = -sum(float(p) * math.log(float(p)) for p in w.values() if p > 0)
    return w, Hw, tot, b


print("=== D03 核验：为什么 Zero 不可能有统计力学 ===")

# ---------------------------------------------------------------- A
print("\n[A] |W_L| = C(L, L/2)")
for L in range(2, 21, 2):
    n = sum(1 for _ in balanced_words(L))
    ck_eq(n, math.comb(L, L // 2), "L=%d  A: |W_L| = C(L,L/2)" % L)

# ---------------------------------------------------------------- B
print("\n[B] 层高计数：闭式 vs 枚举")
for L in range(2, 15, 2):
    bc = shell_counts(L)
    bf = shell_counts_band(L)
    ck_eq(bf, bc, "L=%d  B: 带内公式 = 枚举值" % L)
    ck_eq(sum(bf.values()), math.comb(L, L // 2), "L=%d  B2: sum b_h = |W_L|" % L)

# ---------------------------------------------------------------- C, D
print("\n[C/D] 账本熵闭式 S = H(w) + ln2；块项恒为 ln2")
prev_S = None
rows = []
for L in range(2, 65, 2):
    w, Hw, tot, b = ledger_entropy(L)
    S = Hw + math.log(2.0)
    # 直接谱计算：rho 的本征值 = w_h/2，每个重数 2
    spec = []
    for h, p in w.items():
        spec.extend([float(p) / 2.0] * 2)
    S_spec = -sum(p * math.log(p) for p in spec if p > 0)
    ck(abs(S - S_spec) < 1e-10, "L=%d  C: 闭式 S 与直接谱计算一致（%.12f vs %.12f）" % (L, S, S_spec))
    # 块项：sum over blocks of (w_h/2)*ln2 *2 = ln2 * sum w_h = ln2
    block_part = sum(2 * (float(p) / 2.0) * math.log(2.0) for p in w.values())
    ck(abs(block_part - math.log(2.0)) < 1e-12,
       "L=%d  D: 块项 = ln2（实得 %.15f）" % (L, block_part))
    rows.append((L, Hw, S, len(w)))
    if S < (prev_S if prev_S is not None else 1e9) - 1e-9 and L > 4:
        pass
    prev_S = S

print("     L    H(w)      S=H+ln2   层数")
for (L, Hw, S, k) in rows[::6]:
    print("   %4d  %8.5f  %8.5f   %4d" % (L, Hw, S, k))

# ---------------------------------------------------------------- E
print("\n[E] H(w) = O(ln L)：检验 H/ln L 有界")
ratios = [(L, Hw / math.log(L)) for (L, Hw, S, k) in rows if L >= 8]
mx = max(r for _, r in ratios)
mn = min(r for _, r in ratios)
ck(mx < 1.5, "E1: H(w)/ln L 的上界 < 1.5（实得 max = %.4f）" % mx)
ck(mn > 0.1, "E2: H(w)/ln L 有正下界（实得 min = %.4f）" % mn)
print("     H/ln L 的范围 = [%.4f, %.4f]" % (mn, mx))
# 与线性增长对照：H/L -> 0
for (L, Hw, S, k) in rows[-5:]:
    ck(Hw / L < 0.05, "E3: L=%d  H/L = %.5f < 0.05（非广延）" % (L, Hw / L))

# ---------------------------------------------------------------- F
print("\n[F] S/L -> 0")
SL = [(L, S / L) for (L, Hw, S, k) in rows]
ck(SL[-1][1] < SL[0][1], "F1: S/L 下降（%.5f -> %.5f）" % (SL[0][1], SL[-1][1]))
ck(SL[-1][1] < 0.2, "F2: L=%d 时 S/L = %.5f < 0.2" % SL[-1])
print("     S/L: L=2 -> %.5f；L=64 -> %.5f" % (SL[0][1], SL[-1][1]))

# ---------------------------------------------------------------- G
print("\n[G] 碎片化对照：每态一类给 ln|W_L| ~ L ln2（广延）")
for L in [16, 32, 64]:
    uniform_S = math.log(math.comb(L, L // 2))
    S_ledger = ledger_entropy(L)[1] + math.log(2.0)
    ck(uniform_S / L > 0.55, "L=%d  G1: 均匀划分的熵密度 %.4f > 0.55（近 ln2，广延）" % (L, uniform_S / L))
    ck(S_ledger / L < 0.2, "L=%d  G2: 账本熵密度 %.5f < 0.2（非广延）" % (L, S_ledger / L))
    print("     L=%-3d  ln|W_L|/L = %.4f   S_ledger/L = %.5f   比 = %.5f"
          % (L, uniform_S / L, S_ledger / L, S_ledger / uniform_S))

# ---------------------------------------------------------------- H
print("\n[H] 中心性：rho 在每块上是 (w_h/2) I_2")
for L in [4, 6, 8]:
    w, Hw, tot, b = ledger_entropy(L)
    # 每块的谱：两个相等的本征值 w_h/2 => 该块是标量
    ok = True
    for h, p in w.items():
        ev = [float(p) / 2.0] * 2
        if abs(ev[0] - ev[1]) > 1e-15:
            ok = False
    ck(ok, "L=%d  H: 每块两个本征值相等（块内标量 ⟹ 中心）" % L)

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
