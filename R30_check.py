#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R30_check.py -- 核验小群—相位路线的对抗审计：空解 no-go、同义改写、循环两读法失败，
以及唯一幸存的桥 LG-FUNCTOR 与依赖账本更正。

对应 R30_little_group_phase_route_audit.md。
"""

from __future__ import annotations

import cmath
import io
import os
import sys

import numpy as np

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


R30 = read("R30_little_group_phase_route_audit.md")
R3 = read("R3_dimension_selection.md")
G11 = read("G11_dimension_as_consistency.md")
G89 = read("G89_dimension_no_go_and_the_balance_condition.md")
G53 = read("G53_connecting_the_rg_axis_to_lattice_refinement.md")
G27 = read("G27_purification_attempt.md")
Z0 = read("Z0_zero_never_rests_single_axiom.md")
Z14 = read("Z14_closure_cyclic_order_base_theorem.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
PROBE = read(os.path.join("simulations", "zero_sum_geometry_probe_results.json"))

# ======================================================================
head("F1  文档范围：本文是撤回稿，不是选维进展")

check("标题与性质为对抗审计＋撤回",
      "小群—相位路线的对抗审计" in R30
      and "撤回" in R30
      and "对抗审计＋no-go＋依赖账本更正" in R30)
check("三条否定结论在位",
      "no-go R30.1" in R30
      and "no-go R30.2" in R30
      and "no-go R30.3" in R30
      and "空解" in R30
      and "同义改写" in R30
      and "两读法全灭" in R30)
check("没有把四维写成导出",
      "没有推出 `D=4`" in R30
      and "四维仍是【条件】" in R30
      and "没有由四维标签推出 Lorentz" in R30)
check("唯一幸存的桥在位",
      "LG-FUNCTOR" in R30
      and "满射到整个" in R30
      and "不许默认这一步" in R30)

# ======================================================================
head("F2  no-go R30.1：全小群读法给空解（ISO(D-2) 非交换）")


def rot(t: float) -> np.ndarray:
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def tr(x: float, y: float) -> np.ndarray:
    return np.array([[1.0, 0.0, x], [0.0, 1.0, y], [0.0, 0.0, 1.0]])


for theta in (np.pi / 2, np.pi / 3, 1.0):
    comm = rot(theta) @ tr(1.0, 0.0) - tr(1.0, 0.0) @ rot(theta)
    mag = float(np.abs(comm).max())
    check("ISO(2) 在 theta=%.3f 非交换" % theta, mag > 1e-6, "max|[R,T]|=%.6f" % mag)

check("字面读法在 D>=4 上存活集为空",
      "ISO(D-2)\\text{ 交换}\\iff D-2\\le1\\iff D\\le3" in R30
      and "存活集}=\\varnothing" in R30)

# ======================================================================
head("F3  no-go R30.2：旋转部分读法是 G89 §4 第 1/2 条的同义改写")

syn_ok = True
detail = []
for D in range(3, 13):
    n = D - 2
    dim_so = n * (n - 1) // 2
    dim_plus = (D - 2) * (D - 3) // 2
    dim_minus = D - 3
    abelian = n <= 2
    minus_le_one = dim_minus <= 1
    if dim_so != dim_plus or abelian != minus_le_one:
        syn_ok = False
    detail.append("D=%d:%s" % (D, "same" if abelian == minus_le_one else "DIFF"))
check("dim SO(D-2) == dim P_D^+ 且交换性 ⟺ dim P_D^- <= 1", syn_ok, " ".join(detail[:4]) + " ...")
check("存活集与 G89 §4 第 1/2 条逐点同真",
      [D for D in range(4, 13) if (D - 2) <= 2] == [4]
      and [D for D in range(4, 13) if (D - 3) <= 1] == [4])
check("文档明写同义改写判定",
      "逐点同真，故不是独立证据" in R30
      and "第 1–10 条那一类" in R30)
check("G11 引理 41 的维数公式被正确引用",
      "dim\\mathcal P_D^{+}=\\frac{(D-2)(D-3)}2=\\dim SO(D-2)" in R30
      and "dim P_D^- ≤ 1" in R30)

# ======================================================================
head("F4  no-go R30.3：循环两读法全灭")


def zp_in_so(n: int, p: int) -> bool:
    A = np.eye(n)
    c, s = np.cos(2 * np.pi / p), np.sin(2 * np.pi / p)
    A[0, 0], A[0, 1], A[1, 0], A[1, 1] = c, -s, s, c
    powered = np.linalg.matrix_power(A, p)
    return bool(np.allclose(powered, np.eye(n), atol=1e-9)
                and not np.allclose(A, np.eye(n), atol=1e-9))


check("Z_p ⊂ SO(n) 对 n=2..10, p=2..8 全部存在（原生循环不选维）",
      all(zp_in_so(n, p) for n in range(2, 11) for p in range(2, 9)))


def u1_divisible(samples=12) -> bool:
    rng = np.random.default_rng(0)
    for _ in range(samples):
        theta = float(rng.uniform(-np.pi, np.pi))
        for n in range(2, 6):
            psi = cmath.exp(1j * theta / n)
            if abs(psi ** n - cmath.exp(1j * theta)) > 1e-9:
                return False
    return True


def z5_not_divisible() -> bool:
    # 有限循环群（加法记号 Z_5）：取 n=|G|=5，则 5x ≡ 0 对一切 x，故只有单位元有 |G| 次根
    for g in range(5):
        roots = [x for x in range(5) if (5 * x) % 5 == g]
        if (len(roots) > 0) != (g == 0):
            return False
    return True


def z_not_divisible() -> bool:
    # 无限循环群 Z（加法）：2x = 1 无整数解
    return not any(2 * x == 1 for x in range(-50, 51))


check("U(1) 可除（每个元素有 n 次根）", u1_divisible())
check("有限循环 Z_5 不可除（单位元外无 |G| 次根）", z5_not_divisible())
check("无限循环 Z 不可除（2x=1 无解）", z_not_divisible())
check("文档给出两读法的失败",
      "存活集}=\\{3,4,5,\\dots,12\\}" in R30
      and "可除" in R30
      and "无解" in R30)

# ======================================================================
head("F5  依赖账本更正：进口是超集，过筛前提不中立")

check("R3 的过筛前提确含「无质量自旋 2」",
      "必须存在传播的无质量自旋 2 引力子" in R3 or "无质量自旋 2" in R3)
check("文档登记过筛前提与判据共用同一外部概念",
      "无质量自旋 2" in R30
      and "同一外部概念的重复使用" in R30
      and "依赖深度没有降低" in R30
      and "\\textbf{不更少}" in R30)
check("依赖归约主张已撤回",
      "前一稿的\"依赖归约\"主张}\\textbf{撤回}" in R30)
check("没有把 (C1) 的进口写成旧路径的严格超集（自我修正）",
      "超集" not in R30
      and "两者共用" in R30
      and "Wigner 分类" in R30)
check("G11 的「无横截平面」原文在位（进口确为外部）",
      "没有横截平面" in G11 and "没有 $O(D-2)$" in G11)
check("G89 §7.1 的「无候选」原文在位",
      "无候选" in G89 and "表示论提升" in G89)

# ======================================================================
head("F6  原生性：Zero 层没有相位，但有有限非交换结构")

check("Z0 §4.3 明写 Zero 层没有相位",
      "Zero 层没有相位" in Z0)
check("文档引用该原文并据此撤回相位桥",
      "Zero 层没有相位" in R30
      and "=\\text{有限 }\\mathbb Z_p" in R30
      and "连续 }U(1)" in R30)
check("G27 的 D_L 非交换原文在位",
      "非交换" in G27 and "二面体" in G27)
check("文档更正「Zero 无原生非交换来源」的表述",
      "仅对}\\textbf{连续}\\text{群成立" in R30
      and "走私的那一步" in R30)
check("Z14 的有限循环与 Z-CONF 边界被准确引用",
      "Z-CONF" in Z14
      and "物理框架确实落在该旋转群上" in R30)

# ======================================================================
head("F7  另两条路线的硬结果")

check("G53 的逆向流发散／无 UV 不动点原文在位",
      "没有 UV 不动点" in G53 and "逆向" in G53)
check("R3 §6 已把谱维数稳定判为排除",
      "排除为当前选维器" in R3
      and "谱维数稳定等于 4" in R3)
check("探针自我判决逐字在位（结果 JSON）",
      "no stable four-dimensional plateau detected" in PROBE
      and "does not fix the number of effective geometric directions" in PROBE)
check("文档总账表登记三条路线的净结果",
      "一条封死、一条有反向硬结果、一条已被排除" in R30
      and "靶值可判性" in R30)

# ======================================================================
head("F8  站得住的边界细化与新增开放项")

check("G11 判定不变，只把缺口定位到零方向",
      "G11\\text{ 的 no-go}\\textbf{ 判定不变}" in R30
      and "零方向" in R30)
check("Z4 终端款证据链的归属张力被登记",
      "Z4-EVIDENCE-SCOPE" in R30
      and "错误归属" in R30
      and "对 }\\Gamma\\text{ 的效力}\\textbf{未证}" in R30)
check("G15 的物质层 no-go 被引用",
      "no-go" in R30 and "G15" in R30)

# ======================================================================
head("F9  账本登记")

check("STATUS 已登记 R30",
      "### 2.30" in STATUS
      and "R30_little_group_phase_route_audit.md" in STATUS
      and "R30_check.py" in STATUS
      and "LG-FUNCTOR" in STATUS
      and "Z4-EVIDENCE-SCOPE" in STATUS)
check("INDEX 已收录 R30",
      "R30_little_group_phase_route_audit.md" in INDEX
      and "R30_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
