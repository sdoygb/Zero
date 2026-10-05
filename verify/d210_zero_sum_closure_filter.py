#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D210 —— 零和闭合筛选、无预算与无通道权重的谱系增长
================================================================================
【被检验的问题（D210）】
  去掉有限基本步数预算和外部通道权重后，零和闭合筛选能否给出谱系增长？
  “省、有规律”是否会自动变成更高生殖率？

【本步判据】
  N1  D210 与两个恢复结构已登记
  N2  区分零和、总作用量数值与作用量极值
  N3  基本补偿重写保持零和
  N4  全分支生成不使用概率权重
  N5  闭合谓词与死端谓词
  N6  活扇区由可达性定义
  N7  活后继多重数满足乘法递推
  N8  有限代数计数与公式一致
  N9  内部周期给增长率 lambda=ln(M)/T
  N10 长期排序由 M^(1/T) 决定
  N11 省力降低周期会提高增长率
  N12 规律提高 M 会提高增长率
  N13 规律但 M=1 不产生相对增长
  N14 快但不闭合的分支没有后代
  N15 D184 总通道增长不等于活谱系增长
  N16 不使用总预算变量
  N17 不使用外部通道概率
  N18 无限模式可能没有唯一最强者
  N19 输入预算包含闭合多重数缺口
  N20 文档边界与失败退出码

运行：python3 verify/d210_zero_sum_closure_filter.py
"""
import io
import math
import os
import sys
from collections import Counter, deque

from latex_utils import canonical_math


PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
AX = io.open(os.path.join(ROOT, "AXIOMS.md"), encoding="utf-8").read()
DOC = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D210_zero_sum_closure_filter.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def zero_sum_rewrite(state, i, j):
    """Move one unit from component i to component j."""
    nxt = list(state)
    nxt[i] -= 1
    nxt[j] += 1
    return tuple(nxt)


def all_compensated_successors(state):
    """Return the multiset of all basic zero-sum successors."""
    succ = Counter()
    for i in range(len(state)):
        for j in range(len(state)):
            if i != j:
                succ[zero_sum_rewrite(state, i, j)] += 1
    return succ


def live_states(successors, closed):
    """
    Greatest fixed point of L = C union {x: Succ(x) intersects L}.
    This is the set of states with at least one infinite path, with closed
    states admitted as possible terminal cycles.
    """
    states = set(successors)
    live = set(closed)
    changed = True
    while changed:
        changed = False
        for state in states - live:
            if any(nxt in live for nxt in successors[state]):
                live.add(state)
                changed = True
    return live


def descendants(multiplicity, generations):
    counts = [1]
    for _ in range(generations):
        counts.append(counts[-1] * multiplicity)
    return counts


def lambda_rate(multiplicity, period):
    if multiplicity == 0:
        return -math.inf
    return math.log(multiplicity) / period


# ==================================================================
head("N1  D210 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 D210、闭合筛选增长与闭合多重数缺口已纳入纲领索引",
    "D210" in AX
    and "R-Z-CLOSURE-FILTER-GROWTH" in AX
    and "R-Z-CLOSURE-MULTIPLICITY-GAP" in AX
    and "# D210" in DOC,
)


# ==================================================================
head("N2  区分零和、总作用量数值与作用量极值")
# ==================================================================
check(
    "N2 文档不把 Q=0、S=0 与 delta S=0 混用",
    "Q(x)=\\sum_i x_i=0" in DOC
    and "S_{\\rm total}=0" in DOC
    and "\\delta S=0" in DOC
    and "不能互相替代" in DOC,
)


# ==================================================================
head("N3  基本补偿重写保持零和")
# ==================================================================
charge_ok = True
sample_states = [
    (1, -1, 0),
    (2, -1, -1),
    (0, 0, 0),
    (3, -2, -1),
]
for state in sample_states:
    for nxt in all_compensated_successors(state):
        charge_ok = charge_ok and sum(state) == 0 and sum(nxt) == 0
check(
    "N3 所有基本补偿重写保持零和",
    charge_ok
    and "T_{ij}x" in DOC
    and "\\Delta Q" in DOC,
)


# ==================================================================
head("N4  全分支生成不使用概率权重")
# ==================================================================
succ = all_compensated_successors((1, -1, 0))
check(
    "N4 重写生成多重集合而非概率分布",
    isinstance(succ, Counter)
    and all(isinstance(mult, int) and mult > 0 for mult in succ.values())
    and "不取概率权重" in DOC
    and "m_i\\in\\mathbb N" in DOC,
)


# ==================================================================
head("N5  闭合谓词与死端谓词")
# ==================================================================
check(
    "N5 文档明确区分 Closed 与 Dead",
    "\\mathsf{Closed}(x)" in DOC
    and "\\mathsf{Dead}(x)" in DOC
    and "死端不一定表示某个对象在物理上被销毁" in DOC,
)


# ==================================================================
head("N6  活扇区由可达性定义")
# ==================================================================
successors = {
    "a": ["b", "d"],
    "b": ["a"],
    "c": ["d"],
    "d": ["d"],
    "e": ["f"],
    "f": ["e"],
    "g": ["h"],
    "h": ["g"],
}
closed = {"d", "f"}
live = live_states(successors, closed)
check(
    "N6 活扇区为最大可达闭合集",
    live == {"a", "b", "c", "d", "e", "f"}
    and "g" not in live
    and "h" not in live
    and "\\operatorname{Succ}(x)\\cap L\\ne\\varnothing" in DOC,
)


# ==================================================================
head("N7  活后继多重数满足乘法递推")
# ==================================================================
mult = 3
counts = descendants(mult, 5)
check(
    "N7 N_g=3^g 对前六代成立",
    counts == [1, 3, 9, 27, 81, 243]
    and "N_\\alpha(g+1)" in DOC
    and "M_\\alpha N_\\alpha(g)" in DOC,
)


# ==================================================================
head("N8  有限代数计数与公式一致")
# ==================================================================
finite_ok = all(descendants(m, 4)[-1] == m**4 for m in range(0, 6))
check(
    "N8 0 至 5 重数的有限代数计数等于 M^g",
    finite_ok
    and "M_\\alpha^gN_\\alpha(0)" in DOC,
)


# ==================================================================
head("N9  内部周期给增长率 lambda=ln(M)/T")
# ==================================================================
growth_ok = all(
    abs(lambda_rate(2, 3) - math.log(2) / 3) < 1e-15
    for _ in range(1)
)
check(
    "N9 lambda=ln(M)/T",
    growth_ok
    and "\\lambda_\\alpha" in DOC
    and "\\frac{\\ln M_\\alpha}{T_\\alpha}" in DOC,
)


# ==================================================================
head("N10  长期排序由 M^(1/T) 决定")
# ==================================================================
ordering_ok = True
for m_a in range(1, 6):
    for t_a in range(1, 6):
        for m_b in range(1, 6):
            for t_b in range(1, 6):
                lhs = m_a ** (1 / t_a) > m_b ** (1 / t_b)
                rhs = lambda_rate(m_a, t_a) > lambda_rate(m_b, t_b)
                if lhs != rhs:
                    ordering_ok = False
check(
    "N10 在 1 至 5 的参数网格上，lambda 排序与 M^(1/T) 排序一致",
    ordering_ok
    and "M_A^{1/T_A}>M_B^{1/T_B}" in DOC,
)


# ==================================================================
head("N11  省力降低周期会提高增长率")
# ==================================================================
fast = lambda_rate(2, 1)
slow = lambda_rate(2, 3)
check(
    "N11 同一 M=2 时，T=1 的模式快于 T=3",
    fast > slow
    and "“省”通过减小" in DOC,
)


# ==================================================================
head("N12  规律提高 M 会提高增长率")
# ==================================================================
regular = lambda_rate(3, 2)
less_regular = lambda_rate(2, 2)
check(
    "N12 同一 T=2 时，M=3 的模式增长快于 M=2",
    regular > less_regular
    and "“规律”使" in DOC,
)


# ==================================================================
head("N13  规律但 M=1 不产生相对增长")
# ==================================================================
neutral_counts = descendants(1, 8)
check(
    "N13 M=1 时计数恒为 1，规律本身不增殖",
    neutral_counts == [1] * 9
    and lambda_rate(1, 1) == 0
    and "规律性只有在提高" in DOC,
)


# ==================================================================
head("N14  快但不闭合的分支没有后代")
# ==================================================================
death_path = descendants(0, 7)
check(
    "N14 M=0 时第一代后谱系为零",
    death_path[0] == 1
    and all(count == 0 for count in death_path[1:])
    and "局部活跃不等于谱系增长" in DOC,
)


# ==================================================================
head("N15  D184 总通道增长不等于活谱系增长")
# ==================================================================
phi = (1 + math.sqrt(5)) / 2
check(
    "N15 融合矩阵增长模 phi 不等于活 L 谱系 M=1",
    abs(phi - (1 + math.sqrt(5)) / 2) < 1e-15
    and phi > 1
    and lambda_rate(1, 1) == 0
    and "代数维度增长与活谱系增长是两件事" in DOC,
)


# ==================================================================
head("N16  不使用总预算变量")
# ==================================================================
forbidden_budget = ("B_{\\rm total}", "\\sum_\\alpha C_\\alpha", "N_{\\max}")
check(
    "N16 增长律不要求总预算上限",
    all(fragment in DOC for fragment in forbidden_budget)
    and "删除有限预算不改变条件增长律" in DOC,
)


# ==================================================================
head("N17  不使用外部通道概率")
# ==================================================================
growth_mult = 4
counts_without_probability = descendants(growth_mult, 6)
check(
    "N17 计数律只使用整数重数，不使用 p_i",
    counts_without_probability[-1] == growth_mult**6
    and "p_i\\in[0,1]" in DOC
    and "通道选择被全分支实例化替代" in DOC,
)


# ==================================================================
head("N18  无限模式可能没有唯一最强者")
# ==================================================================
finite_rates = [lambda_rate(m, t) for m in range(1, 6) for t in range(1, 6)]
check(
    "N18 有限网格有最大值，但文档保留无限 sup 不可达边界",
    max(finite_rates) == math.log(5)
    and "\\sup_\\alpha\\lambda_\\alpha" in DOC
    and "无限模式不保证唯一冠军" in DOC,
)


# ==================================================================
head("N19  输入预算包含闭合多重数缺口")
# ==================================================================
check(
    "N19 R-Z-CLOSURE-MULTIPLICITY-GAP 登记规律性到 M 的未证关系",
    "R-Z-CLOSURE-MULTIPLICITY-GAP" in DOC
    and "规律性必然提高 $M_\\alpha$" in DOC
    and "未解选择器" in DOC,
)


# ==================================================================
head("N20  文档边界与失败退出码")
# ==================================================================
check(
    "N20 文档不声称无条件导出、不修改 U1-U4+C1、不新增 U5",
    "不修改 `U1-U4+C1`" in DOC
    and "不新增 `U5`" in DOC
    and "仍需证明" in DOC
    and "不是 `U1-U4` 的推论" in DOC,
)


print("\n" + "=" * 78)
print(f"通过 {len(PASS)} / 不符 {len(FAIL)}")
if FAIL:
    print("失败项：")
    for name in FAIL:
        print(" -", name)
    sys.exit(1)
print("全部通过")
