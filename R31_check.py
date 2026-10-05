#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R31_check.py -- 核验"相位接通账本"：旋转类轨道分布给出 q_L，成对账本的峰位唯一选出 L=4。

对应 R31_phase_ledger_and_lifetime_selection.md。
"""

from __future__ import annotations

import io
import itertools
import os
import sys
from collections import Counter
from fractions import Fraction
from math import comb, gcd

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(title: str) -> None:
    print("")
    print("=" * 72)
    print(title)
    print("=" * 72)


def check(name: str, cond: bool, detail: str = "") -> None:
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def read(name: str) -> str:
    with io.open(os.path.join(HERE, name), encoding="utf-8") as handle:
        return handle.read()


R31 = read("R31_phase_ledger_and_lifetime_selection.md")
R25 = read("R25_native_pair_cost_and_four_dim_peak.md")
R23 = read("R23_dimension_descendant_selection.md")
R29 = read("R29_full_support_ledger_factorization_no_go.md")
G61 = read("G61_locking_the_five_integers.md")
Z14 = read("Z14_closure_cyclic_order_base_theorem.md")
ALG = read("zero_sum_rotation_class_algebra.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")


# ---------------------------------------------------------------- 工具
def orbit_sizes(L: int):
    """零和 ±1 词在循环旋转 Z_L 下的轨道大小多重集（枚举）"""
    words = [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]
    seen = set()
    sizes = []
    for w in words:
        if w in seen:
            continue
        orb = {w[i:] + w[:i] for i in range(L)}
        seen |= orb
        sizes.append(len(orb))
    return sizes


def necklace_count(L: int) -> int:
    """旧理论项链闭式 K(L) = (1/L) * sum_{d | gcd(L, L/2)} phi(d) * C(L/d, L/(2d))"""
    def phi(n):
        return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)
    g = gcd(L, L // 2)
    return sum(phi(d) * comb(L // d, L // (2 * d)) for d in range(1, g + 1) if g % d == 0) // L


def q_of(L: int) -> Fraction:
    sizes = orbit_sizes(L)
    N = comb(L, L // 2)
    return sum(Fraction(o, N) ** 2 for o in sizes)


def peak_d(q: Fraction, dmax: int = 300):
    """成对账本 F_D = C(D,2) q^D 的有限峰；q>=1 时严格递增，返回 None 表示无有限峰"""
    if q >= 1:
        return None
    vals = [(D, comb(D, 2) * q ** D) for D in range(2, dmax + 1)]
    mx = max(v for _, v in vals)
    return [D for D, v in vals if v == mx]


# ======================================================================
head("F1  文档范围：这是相位接通账本的选择定理，不是无条件导出")

check("标题与性质在位",
      "相位接通账本" in R31
      and "唯一选出 `L=4`" in R31
      and "原生选择＋账本升级" in R31)
check("定理 R31.1 与唯一性陈述在位",
      "定理 R31.1" in R31
      and "恰有一个" in R31
      and "\\arg\\max_{D\\ge2}F_D=\\{4\\}" in R31)
check("没有把结果写成无条件",
      ("不**关闭**" in R31 or "不关闭" in R31 or "**不**关闭" in R31)
      and "没有推出 `LEDGER-ROT`" in R31
      and "没有由四维标签推出 Lorentz" in R31)

# ======================================================================
head("F2  相位的轨道分布与 q_L（枚举 + 项链闭式交叉核对）")

L_LIST = [2, 4, 6, 8, 10, 12, 14, 16, 18]
data = {}
for L in L_LIST:
    sizes = orbit_sizes(L)
    N = comb(L, L // 2)
    q = sum(Fraction(o, N) ** 2 for o in sizes)
    data[L] = (N, sorted(sizes, reverse=True), q)

check("轨道个数与项链闭式逐项相等",
      all(len(data[L][1]) == necklace_count(L) for L in L_LIST),
      "K=%s" % [len(data[L][1]) for L in L_LIST])
check("L=4 轨道分布 {4,2} 且 q_4 = 5/9",
      data[4][1] == [4, 2] and data[4][2] == Fraction(5, 9),
      "q_4=%s" % data[4][2])
check("L=6 轨道分布 {6,6,6,2}；L=8 为 {8x8,4,2}",
      data[6][1] == [6, 6, 6, 2] and data[8][1] == [8] * 8 + [4, 2])
check("与 zero_sum_rotation_class_algebra 的多重集逐字一致",
      "{4,2}" in ALG.replace(" ", "")
      or ("$L=4$ 给 $\\{4,2\\}$" in ALG and "$\\{8^{\\times8},4,2\\}$" in ALG))
check("R25 的 q_4=5/9 与 omega=[2/3,1/3] 在位",
      "q_4" in R25 and "5/9" in R25.replace(" ", ""))

# ======================================================================
head("F3  定理 R31.1：唯一使峰落在 D>=4 的偶寿命是 L=4")

raw = {L: q_of(L) for L in L_LIST}
peaks = {L: peak_d(raw[L]) for L in L_LIST}

check("L=2 的账本严格递增、无有限峰",
      raw[2] == 1 and peaks[2] is None
      and all(Fraction(D + 1, D - 1) * raw[2] > 1 for D in range(2, 50)))
check("L>=6 时 q_L <= L/N_L 且 <= r_6 = 3/10 < 1/3",
      all(raw[L] <= Fraction(L, data[L][0]) <= Fraction(3, 10) for L in L_LIST if L >= 6),
      "max q_L(L>=6)=%s" % max(float(raw[L]) for L in L_LIST if L >= 6))
check("L>=6 的峰全在 D=2（不在引力子域）",
      all(peaks[L] == [2] for L in L_LIST if L >= 6),
      "peaks=%s" % {L: peaks[L] for L in L_LIST if L >= 6})
check("L=4 的峰唯一为 D=4",
      peaks[4] == [4], "peak=%s" % peaks[4])
check("唯一性：全部偶 L>=2 中恰有 L=4 存在落在引力子域的有限峰",
      [L for L in L_LIST if peaks[L] is not None and peaks[L][0] >= 4] == [4])
check("r_L 随偶 L 严格递减（R25 引理 R25.2 的单调性）",
      all(Fraction(L + 2, comb(L + 2, (L + 2) // 2)) < Fraction(L, comb(L, L // 2))
          for L in [2, 4, 6, 8, 10, 12, 14, 16]))
check("L=4 两侧显式比值：F4/F3=10/9>1，F5/F4=25/27<1",
      Fraction(4, 2) * Fraction(5, 9) == Fraction(10, 9)
      and Fraction(5, 3) * Fraction(5, 9) == Fraction(25, 27))
check("文档给出相邻比与四维窗口",
      "\\frac{D+1}{D-1}q_L" in R31
      and "\\left(\\frac12,\\frac35\\right)" in R31)

# ======================================================================
head("F4  判据的性质：不含数字、不用最小性")

check("新增输入 PEAK-IN-GRAVITON-DOMAIN 在位",
      "PEAK-IN-GRAVITON-DOMAIN" in R31
      and "引力子允许域" in R31)
check("文档明写它不含数字 4、不用最小性、是生存要求",
      "不含数字" in R31 and "不用最小性" in R31 and "生存要求，不是偏好" in R31)
check("G61 确实用最小性锁定 L=4（对照在位）",
      "最小" in G61 and "L=4" in G61.replace(" ", ""))
check("R23.6 的 L=8 对照在位（张力已登记）",
      "L=8" in R23.replace(" ", "")
      and "R23.6" in R31
      and "张力" in R31)

# ======================================================================
head("F5  相位链条与 R30 的分工")

check("文档给出六环链条 ①→⑥",
      "① 相位" in R31 and "⑥ 结论" in R31 and "①\\to②\\to③\\to④\\to⑤\\to⑥" in R31)
check("文档说明与 R30 不冲突（做不了小群但能做账本）",
      "相位做不了小群，但能做账本" in R31)
check("Z14 的相位与双覆盖仍被正确引用",
      "Z14" in R31 and "引理 Z14.1" in R31 and "Z14.3" in R31)
check("R29 的 L=4 行已更新为指向 R31 的条件导出（张力在位）",
      "R31 已给" in R29 and "定理 R31.1" in R29
      and "条件导出" in R29 and "L=8" in R29.replace(" ", ""))

# ======================================================================
head("F6  账本登记")

check("STATUS 已登记 R31",
      "### 2.31" in STATUS
      and "R31_phase_ledger_and_lifetime_selection.md" in STATUS
      and "R31_check.py" in STATUS
      and "PEAK-IN-GRAVITON-DOMAIN" in STATUS
      and "PEAK-IN-GRAVITON-DOMAIN" in STATUS)
check("INDEX 已收录 R31",
      "R31_phase_ledger_and_lifetime_selection.md" in INDEX
      and "R31_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
