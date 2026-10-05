#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R32_check.py -- 核验账本读出的选择（三条路线唯一存活 A）与 R23.6 的 L=8 张力消解。

对应 R32_ledger_readout_selection_and_L8_resolution.md。
"""

from __future__ import annotations

import io
import os
import sys
from fractions import Fraction
from math import comb, e, log

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


R32 = read("R32_ledger_readout_selection_and_L8_resolution.md")
R31 = read("R31_phase_ledger_and_lifetime_selection.md")
R23 = read("R23_dimension_descendant_selection.md")
R27 = read("R27_phase_cochain_quotient_and_dimension_dictionary.md")
G61 = read("G61_locking_the_five_integers.md")
G71 = read("G71_decoherence_from_the_terminal_ledger.md")
G72 = read("G72_kappa1_from_the_ledger.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")


def peak(q: float, dmax: int = 200):
    """成对账本 F_D=C(D,2)q^D 的峰集；q>=1 时返回 None"""
    if q >= 1:
        return None
    vals = [(D, comb(D, 2) * q ** D) for D in range(2, dmax + 1)]
    mx = max(v for _, v in vals)
    return [D for D, v in vals if v == mx]


# ======================================================================
head("F1  文档范围：这是选择定理＋张力消解，不是无条件导出")

check("标题与两条定理在位",
      "账本读出的选择与 `L=8` 张力的消解" in R32
      and "定理 R32.1" in R32
      and "定理 R32.2" in R32)
check("没有写成无条件",
      ("**不**关闭" in R32 or "不关闭" in R32)
      and "没有推出成对账本形式" in R32
      and "没有由四维标签推出 Lorentz" in R32)
check("两条边界在位（L=2 与路线表穷尽性）",
      "`L=2` 的边缘情形" in R32
      and "路线表的穷尽性未证" in R32)

# ======================================================================
head("F2  三条路线的 q 与峰")

check("路线 A 的 q_4=5/9 给峰 {4}",
      peak(5 / 9) == [4])
check("路线 B/C 的 q=1/L 在 L>=4 上峰恒为 {2}",
      all(peak(Fraction(1, L)) == [2] for L in range(4, 21, 2)),
      "L=4..20")
check("路线 D12 的 q=e^-1 给峰 {3}",
      peak(1 / e) == [3], "peak=%s" % peak(1 / e))
check("三条路线中只有 A 落在 D>=4",
      peak(5 / 9) == [4] and peak(Fraction(1, 4)) == [2] and peak(1 / e) == [3])
check("相邻比 (D+1)/(D-1)q 关于 D 严格递减（用于唯一性）",
      all(Fraction(D + 1, D - 1) > Fraction(D + 2, D) for D in range(2, 40)))
check("G71 §4 原文给出时间残类 q=1/L",
      "时间残类 $\\Rightarrow1/L$" in G71 or "时间残类" in G71 and "1/L" in G71)
check("G72 §4 的 e^-1 超越性 no-go 在位",
      "超越" in G72 and "有理数" in G72)
check("文档给出三条路线表",
      "闭类旋转轨道" in R32 and "时间残类" in R32 and "D12" in R32)

# ======================================================================
head("F2.1  命题 R32.3：判据的极小性（可弱化到'峰 != 2'）")

ROUTES = [("A, L=2", 1.0), ("A, L=4", 5 / 9), ("A, L=6", 7 / 25), ("A, L=8", 19 / 175),
          ("B/C, L=4", 0.25), ("B/C, L=6", 1 / 6), ("B/C, L=8", 0.125), ("D12", 1 / e)]


def survivors(smin: int):
    return [name for name, q in ROUTES if peak(q) is not None and peak(q)[0] >= smin]


check("S={D>=4} 唯一存活 (A, L=4)",
      survivors(4) == ["A, L=4"], "survivors=%s" % survivors(4))
check("S={D>=3} 存活 (A,L=4) 与 D12；再用有理纯度 no-go 排除 D12 后唯一",
      survivors(3) == ["A, L=4", "D12"]
      and [n for n in survivors(3) if n != "D12"] == ["A, L=4"])
check("S={D>=2} 不唯一（B/C 存活）",
      len(survivors(2)) > 1 and any("B/C" in n for n in survivors(2)))
check("文档给出有条件弱化形式 PEAK-NOT-GRAVITY-FREE 并限定范围",
      "PEAK-NOT-GRAVITY-FREE" in R32
      and "几何扇区没有动力学" in R32
      and "范围限制" in R32
      and "在固定身份计数下" in R32)
check("R31 已回指该弱化",
      "PEAK-NOT-GRAVITY-FREE" in R31 and "命题 R32.3" in R31)

# ======================================================================
head("F2.2  命题 R32.4：联合唯一性（身份计数 × 路线 × 寿命）")

IDENT = {
    "D": lambda D: D,
    "C(D,2)": lambda D: comb(D, 2),
    "C(D+1,2)": lambda D: comb(D + 1, 2),
}
QCAND = [("A,L=2", 1.0), ("A,L=4", 5 / 9), ("A,L=6", 7 / 25), ("A,L=8", 19 / 175),
         ("B/C,L=4", 0.25), ("B/C,L=6", 1 / 6), ("D12", 1 / e)]


def peak_g(g, q, dmax=200):
    vals = [(D, g(D) * q ** D) for D in range(2, dmax + 1)]
    if q >= 1:
        return None
    mx = max(v for _, v in vals)
    return [D for D, v in vals if v == mx]


joint4 = [(name, qn) for name, g in IDENT.items() for qn, q in QCAND
          if peak_g(g, q) is not None and peak_g(g, q)[0] >= 4]
joint3 = [(name, qn) for name, g in IDENT.items() for qn, q in QCAND
          if peak_g(g, q) is not None and peak_g(g, q)[0] >= 3]

check("S={D>=4} 下联合叉积唯一存活 (C(D,2), A, L=4)",
      joint4 == [("C(D,2)", "A,L=4")], "joint=%s" % joint4)
check("S={D>=3} 下联合叉积不唯一（字典 C(D+1,2) 也存活）",
      len(joint3) > 1 and ("C(D+1,2)", "A,L=4") in joint3,
      "joint=%s" % joint3)
check("字典计数 C(D+1,2) 在 q=5/9 下峰为 D=3（与 R26 一致）",
      peak_g(IDENT["C(D+1,2)"], 5 / 9) == [3])
check("单方向计数 D 在 q=5/9 下峰为 D=2（与 R25 一致）",
      peak_g(IDENT["D"], 5 / 9) == [2])
check("文档给出联合唯一性命题与两条推论",
      "命题 R32.4（联合唯一性）" in R32
      and "PAIR-CARRIER` 也被选出" in R32
      and "联合选择需要判据的完整形式" in R32)

# ======================================================================
head("F3  R23.6 的 L=8 与其三条消解理由")


def r_of(L: int) -> float:
    return (log(L) - 1) / L


rs = {L: r_of(L) for L in range(4, 15, 2)}
check("R23.6 的每步增长率在偶闭圈上唯一极大点为 L=8",
      max(rs, key=rs.get) == 8,
      "r(6)=%.5f r(8)=%.5f r(10)=%.5f" % (rs[6], rs[8], rs[10]))
check("R23.6 的单周期量 P(L)=L/e 严格递增（无有限极大点）",
      all((L / e) < ((L + 2) / e) for L in range(4, 60, 2)))
check("R23 原文含 L=8 与每步增长率 r(L)",
      "L=8" in R23.replace(" ", "")
      and "\\log" in R23
      and "P(L)" in R23.replace(" ", ""))
check("R23.6 的 q_L=e^{-1/L} 不是账本纯度（超越 vs 有理）",
      "e^{-1/L}" in R32.replace(" ", "") or "e^{-1/L}" in R32)
check("R27 §6 的 WIPE-RESET-LEDGER 撤回每步比较量（原文在位）",
      "WIPE-RESET-LEDGER" in R27 and "不再除以串行内部耗时" in R27.replace(" ", ""))
check("文档写明三重理由（非账本 q／已撤回比较量／无极大点）",
      "不是账本纯度" in R32
      and "它的比较量已被撤回" in R32
      and "自陈无极大点" in R32)

# ======================================================================
head("F4  L=2 边缘情形与 G61 的排除（不用最小性）")

check("B/C 在 L=2 给并列峰 {3,4}（边缘情形确实存在）",
      peak(Fraction(1, 2)) == [3, 4], "peak=%s" % peak(Fraction(1, 2)))
check("G61 §1 的 (A)∧(B) 给 L>=4 而不使用最小性",
      "(A)\\wedge(B)" in G61.replace(" ", "")
      and "偶且" in G61.replace(" ", "")
      and "非交换载体" in G61)
check("文档明写该步不用最小性",
      "不用最小性" in R32 and "L\\ge4" in R32)

# ======================================================================
head("F5  账本升级与登记")

check("LEDGER-ROT 升为条件选择",
      "条件选择" in R32 and "LEDGER-ROT" in R32)
check("STATUS 已登记 R32",
      "### 2.32" in STATUS
      and "R32_ledger_readout_selection_and_L8_resolution.md" in STATUS
      and "R32_check.py" in STATUS
      and "定理 R32.1" in STATUS)
check("INDEX 已收录 R32",
      "R32_ledger_readout_selection_and_L8_resolution.md" in INDEX
      and "R32_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
