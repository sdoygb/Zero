#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z001_hcount_check.py —— Z001 §7 的层高分布闭式 (7.7) 的核验。

核什么：$n=16$ 的层高分布 $b_h$ 由**闭式**算出来，必须与正文表逐位一致。
    b_h = N_≤(n,h) − N_≤(n,h−1)
    N_≤(n,h) = Σ_j [ C(n, n/2+2j(h+1)) − C(n, n/2+(h+1)+2j(h+1)) ]
正文表（Z001 §7）：b_1..b_8 = 256, 4118, 4858, 2518, 880, 208, 30, 2

为什么要有这个脚本：闭式里的 $2j(h+1)$ 与偏移量极易抄错，而**抄错后仍然"像"是对的**——
曾出现一个把 $2j(h+1)$ 写成 $j(h+2)$ 的版本，给出 $1,127,482,484,247,75,13,1$ 而无人察觉。
本脚本把"闭式 ⟺ 表"钉死，谁改动都会立刻红。

用法：python3 _staging/Z001_hcount_check.py    （退出码 0 = 全过）
"""
import math
import sys

FAIL = 0
TABLE = [256, 4118, 4858, 2518, 880, 208, 30, 2]     # Z001 §7 的 n=16 表


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("   " + detail) if detail else ""))
    if not ok:
        FAIL += 1


def N_le(n, h):
    """Z001 (7.7)：层高 ≤ h 的词数。"""
    s = 0
    k = h + 1
    for j in range(-n, n + 1):
        a = n // 2 + 2 * j * k
        b = a + k
        if 0 <= a <= n:
            s += math.comb(n, a)
        if 0 <= b <= n:
            s -= math.comb(n, b)
    return s


print("=" * 68)
print("闭式 (7.7) 复现 Z001 §7 的 n=16 表")
print("=" * 68)
print("   h      N_≤(16,h)      b_h       表值")
got = []
for h in range(1, 9):
    b = N_le(16, h) - N_le(16, h - 1)
    got.append(b)
    print("   %-6d %-14d %-10d %d" % (h, N_le(16, h), b, TABLE[h - 1]))
check("b_1..b_8 与表逐位一致", got == TABLE, "实测 %s" % got)

print("\n" + "=" * 68)
print("反例对照：$2j(h+1)$ 写成 $j(h+2)$ 的版本会给出什么")
print("=" * 68)


def N_le_wrong(n, h):
    s = 0
    k = h + 2
    for j in range(-n, n + 1):
        a = n // 2 + j * k
        b = a + h + 1
        if 0 <= a <= n:
            s += math.comb(n, a)
        if 0 <= b <= n:
            s -= math.comb(n, b)
    return s


wrong = [N_le_wrong(16, h) - N_le_wrong(16, h - 1) for h in range(1, 9)]
print("   错误版本的 b_h = %s" % wrong)
check("错误版本确实与表不符（说明本核验有判别力）", wrong != TABLE, "它给 %s" % wrong)

print("\n" + "=" * 68)
print("旁证：几个 n 的 b_h 之和 = 词数 $\\binom{n}{n/2}$（层高穷尽）")
print("=" * 68)
for n in (4, 8, 12, 16, 20):
    tot = sum(N_le(n, h) - N_le(n, h - 1) for h in range(1, n // 2 + 1))
    want = math.comb(n, n // 2)
    print("   n=%-3d Σb_h=%-8d C(n,n/2)=%-8d" % (n, tot, want))
    check("n=%d 层高穷尽" % n, tot == want)

print("\n汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)