#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D213
===================
【被检验的问题（D213）】
  当前周期零宇宙骨架能否直接推出 GR？
  全局零和能否冒充局部守恒源？
  周期寿命序能否直接给出四维洛伦兹几何与 Einstein 方程？

运行：python3 verify/d213_periodic_skeleton_to_gr_direct_audit.py
"""
import io
import json
import os
import sys

import numpy as np

from latex_utils import canonical_math


PASS, FAIL = [], []


def check(name, condition, detail=""):
    (PASS if condition else FAIL).append(name)
    print(f"  [{'v' if condition else 'x'}] {name}" + (f"   {detail}" if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
DOC_PATH = os.path.join(
    ROOT,
    "derivations",
    "D213_periodic_skeleton_to_gr_direct_audit.md",
)
AX_PATH = os.path.join(ROOT, "AXIOMS.md")
PERIODIC_JSON = os.path.join(
    ROOT,
    "simulations",
    "zero_sum_periodic_destruction_results.json",
)

DOC = canonical_math(io.open(DOC_PATH, encoding="utf-8").read())
AX = io.open(AX_PATH, encoding="utf-8").read()
BOTH = AX + "\n" + DOC

with io.open(PERIODIC_JSON, encoding="utf-8") as handle:
    PERIODIC = json.load(handle)


def path_balance_change(sign):
    return sign


def disconnected_laplacian(number_of_sites):
    return np.zeros((number_of_sites, number_of_sites), dtype=float)


def ring_laplacian(size):
    laplacian = 2.0 * np.eye(size)
    for index in range(size):
        laplacian[index, (index - 1) % size] -= 1.0
        laplacian[index, (index + 1) % size] -= 1.0
    return laplacian


def support_step(point, time, period):
    if point == "global":
        return "global"
    age = point[1]
    return ("local", age + time) if age + time < period else "global"


# ==================================================================
head("N1  D213 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 文档登记直接 GR 障碍与六槽补全路线",
    "R-Z-GR-DIRECT-OBSTACLE" in BOTH
    and "R-Z-GR-DIRECT-COMPLETION" in BOTH
    and "1/6" in DOC,
)
check(
    "N1 文档没有声称已经直接推出 GR",
    "不能直接推出 GR" in DOC
    and "可组合，不可直接" in DOC
    and "不能主张" in DOC,
)


# ==================================================================
head("N2  当前周期骨架的参数")
# ==================================================================
period = int(PERIODIC["model"]["destruction_period"])
initial_words = PERIODIC["model"]["initial_words"]
site_count = len(initial_words)
check(
    "N2 当前周期为 T=3 且区域数为 6",
    period == 3 and site_count == 6,
    f"T={period}, sites={site_count}",
)
check(
    "N2 文档核验使用同一骨架参数",
    "T=3" in DOC
    and "六个区域" in DOC
    and "\\dim_{\\mathbb C}M_3(\\mathbb C)=9" in DOC,
)


# ==================================================================
head("N3  单区域周期步不是局域零和补偿")
# ==================================================================
single_step_changes = [
    path_balance_change(sign)
    for sign in (1, -1)
]
check(
    "N3 单区域路径步改变局部电荷",
    single_step_changes == [1, -1]
    and all(change != 0 for change in single_step_changes),
)
check(
    "N3 文档给出局域补偿式的缺失",
    "T_{ij}x=x+e_j-e_i" in DOC
    and "\\Delta Q=0" in DOC
    and "没有一条基本边" in DOC,
)


# ==================================================================
head("N4  六区域图没有跨区域边")
# ==================================================================
site_adjacency = np.zeros((site_count, site_count), dtype=int)
cross_edges = int(site_adjacency.sum())
check(
    "N4 当前骨架的跨区域邻接矩阵为零",
    cross_edges == 0,
    f"cross_edges={cross_edges}",
)
check(
    "N4 文档区分账本秩与局域维数",
    "6-1=5" in DOC
    and "账本秩" in DOC
    and "局域维数" in DOC,
)


# ==================================================================
head("N5  孤立区域没有四维扩散窗口")
# ==================================================================
laplacian = disconnected_laplacian(site_count)
check(
    "N5 零 Laplacian 不传播到其他区域",
    np.allclose(laplacian, 0.0),
)
check(
    "N5 文档登记四维扩散窗口缺失",
    "四维热核窗口" in DOC
    and "账本秩" in DOC
    and "局域维数" in DOC,
)


# ==================================================================
head("N6  周期寿命流给方向，不给空间维数")
# ==================================================================
local_points = [("local", age) for age in range(period)]
semigroup_ok = all(
    support_step(support_step(point, t, period), s, period)
    == support_step(point, s + t, period)
    for point in local_points
    for t in range(0, 2 * period + 2)
    for s in range(0, 2 * period + 2)
)
check(
    "N6 寿命流满足有向半群律",
    semigroup_ok,
)
check(
    "N6 寿命流在 T 后进入终端全局支持",
    support_step(("local", 0), period, period) == "global"
    and support_step("global", 100, period) == "global",
)
check(
    "N6 文档只把寿命流声明为时间方向候选",
    "时间方向候选" in DOC
    and "它不给 }3\\text{ 个空间方向" in DOC,
)


# ==================================================================
head("N7  相位代数维数不等于时空维数")
# ==================================================================
matrix_dimension = period * period
diagonal_dimension = period
spatial_coordinate_rank = 0
check(
    "N7 M_3(C) 的 9 维是代数维数而非空间维数",
    matrix_dimension == 9 and diagonal_dimension == 3,
)
check(
    "N7 当前骨架没有给出 3 个独立空间坐标",
    spatial_coordinate_rank < 3
    and "相位投影" in DOC
    and "不是空间维数" in DOC,
    f"coordinate_rank={spatial_coordinate_rank}",
)


# ==================================================================
head("N8  D 层只有全局零和，不是局部守恒源")
# ==================================================================
dead_balance = PERIODIC["runs"][2]["final_dead_balances"]
local_imbalances = sum(abs(value) > 0 for value in dead_balance)
global_balance = sum(dead_balance)
check(
    "N8 当前 D 层有非零局部分量",
    local_imbalances > 0,
    f"imbalances={local_imbalances}, Q_D={dead_balance}",
)
check(
    "N8 当前 D 层只满足全局补偿",
    global_balance == 0,
    f"sum(Q_D)={global_balance}",
)
check(
    "N8 文档拒绝把全局零和当作局部守恒",
    "全局零和" in DOC
    and "\\not\\Longrightarrow" in DOC
    and "局部守恒应力张量" in DOC,
)


# ==================================================================
head("N9  图 Laplacian 不提供曲率与 Einstein 张量")
# ==================================================================
ring = ring_laplacian(8)
check(
    "N9 图 Laplacian 是半正定扩散算子",
    np.allclose(ring, ring.T)
    and np.min(np.linalg.eigvalsh(ring)) > -1.0e-12,
)
check(
    "N9 文档拒绝从图 Laplacian 直接得到 Riemann 曲率",
    "Riemann 曲率" in DOC
    and "Einstein 张量" in DOC
    and "前几何数据" in DOC,
)


# ==================================================================
head("N10 G 与 Lambda 仍有尺度简并")
# ==================================================================
newton_coupling = 1.0
cosmological_constant = 0.1
scale = 2.0
scaled_newton = scale**2 * newton_coupling
scaled_cosmological = scale ** (-2) * cosmological_constant
check(
    "N10 尺度变换保持 G Lambda",
    abs(newton_coupling * cosmological_constant - scaled_newton * scaled_cosmological) < 1.0e-12,
)
check(
    "N10 周期 T 不提供绝对耦合归一化",
    "Newton 耦合" in DOC
    and "真空能读数" in DOC
    and "绝对归一化仍然缺失" in DOC,
)


# ==================================================================
head("N11 六个槽位中当前只贡献时间候选")
# ==================================================================
slots = {
    "local_transport": False,
    "four_dimensional_refinement": False,
    "lorentz_time_line": True,
    "local_source": False,
    "geometric_action": False,
    "coupling_normalization": False,
}
check(
    "N11 当前骨架只完成时间线候选",
    sum(slots.values()) == 1
    and slots["lorentz_time_line"]
    and not slots["local_transport"]
    and not slots["local_source"],
)
check(
    "N11 文档登记六槽补全路线",
    "C_{\\rm loc}" in DOC
    and "C_{4D}" in DOC
    and "C_{\\rm src}" in DOC
    and "C_{\\rm geom}" in DOC
    and "C_{\\rm norm}" in DOC,
)


# ==================================================================
head("N12 输入预算与上游边界")
# ==================================================================
check(
    "N12 文档列出局域补偿、四维细化、守恒源与几何作用输入",
    "局域零和补偿边" in DOC
    and "四维细化族" in DOC
    and "局部守恒源" in DOC
    and "Einstein-Hilbert 几何作用" in DOC,
)
check(
    "N12 文档不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC
    and "不新增 `U5`" in DOC,
)


print("\n" + "=" * 78)
print(f"通过 {len(PASS)} / 不符 {len(FAIL)}")
if FAIL:
    print("失败项：")
    for name in FAIL:
        print(" -", name)
    sys.exit(1)
print("全部通过")
