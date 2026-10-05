#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D235 —— 年龄到径向映射的重参数化障碍
================================================================================
【被检验的问题（D235）】
  年龄前缀、终端和半流兼容性能否唯一选出年龄到球半径的仿射映射？
  不同单调映射是否给出不同年龄剖面与 likelihood ratio？

【本步判据】
  N1  D235 与恢复结构已登记
  N2  多个 phi_p 都保持零点、终端与单调性
  N3  同一球核拉回后给多个年龄剖面
  N4  全部年龄剖面边界退化并在内部严格下降
  N5  不同年龄剖面的 L1 差异为正
  N6  不同 likelihood ratio 的 L1 差异为正
  N7  p=1 给二次剖面，p=1/2 给线性剖面
  N8  半流共轭保持年龄推进律
  N9  文档登记重参数化无解结果
  N10 文档登记仿射标定缺口
  N11 文档不把 p=1 写成已导出
  N12 上游边界保持
  N13 文档不引用项目外体系名或外部路径

运行：python3 verify/d235_age_radial_reparametrization_no_go.py
"""
import io
import os
import sys

import numpy as np


# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
if os.path.join(ROOT, "verify") not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "verify"))

from latex_utils import canonical_math


PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


AX = io.open(os.path.join(ROOT, "AXIOMS.md"), encoding="utf-8").read()
DOC = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D235_age_radial_reparametrization_no_go.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D235 与恢复结构已登记")
# ==================================================================
check(
    "N1 公理表登记 D235",
    "D235" in AX,
)
check(
    "N1 公理表与文档登记重参数化无解结果",
    "R-Z-AGE-RADIAL-REPARAMETRIZATION-NOGO" in BOTH
    and "重参数化" in DOC,
)
check(
    "N1 公理表与文档登记仿射标定缺口",
    "R-Z-AFFINE-AGE-RADIAL-CALIBRATION-GAP" in BOTH
    and "仿射标定" in DOC,
)


# ==================================================================
head("N2  多个 phi_p 都保持零点、终端与单调性")
# ==================================================================
period = 1.0
radius = 2.0
age = np.linspace(0.0, period, 3001)
p_values = (0.5, 1.0, 2.0, 3.0)
radial_maps = {
    p: radius * (age / period) ** p
    for p in p_values
}

check(
    "N2 所有径向映射保持零点与终端",
    all(
        abs(mapping[0]) < 1e-12
        and abs(mapping[-1] - radius) < 1e-12
        for mapping in radial_maps.values()
    ),
)
check(
    "N2 所有径向映射严格递增",
    all(np.all(np.diff(mapping) > 0.0) for mapping in radial_maps.values()),
)


# ==================================================================
head("N3  同一球核拉回后给多个年龄剖面")
# ==================================================================
age_profiles = {
    p: (radius**2 - radial_maps[p] ** 2) / (2.0 * radius)
    for p in p_values
}
expected_profiles = {
    p: (radius / 2.0) * (1.0 - (age / period) ** (2.0 * p))
    for p in p_values
}

check(
    "N3 拉回公式与直接计算一致",
    all(
        np.max(np.abs(age_profiles[p] - expected_profiles[p])) < 1e-12
        for p in p_values
    ),
)
check(
    "N3 四个年龄剖面不全相同",
    len(
        {
            tuple(np.round(profile, 9))
            for profile in age_profiles.values()
        }
    )
    == len(p_values),
)


# ==================================================================
head("N4  全部年龄剖面边界退化并在内部严格下降")
# ==================================================================
check(
    "N4 所有剖面终端为零且内部为正",
    all(
        abs(profile[-1]) < 1e-12
        and np.all(profile[:-1] > 0.0)
        for profile in age_profiles.values()
    ),
)
check(
    "N4 所有剖面严格下降",
    all(
        np.all(
            -(radius * p / period)
            * (age[1:] / period) ** (2.0 * p - 1.0)
            < 0.0
        )
        for p in p_values
    ),
)


# ==================================================================
head("N5  不同年龄剖面的 L1 差异为正")
# ==================================================================
profile_pairs = ((0.5, 1.0), (1.0, 2.0), (2.0, 3.0))
profile_l1_differences = {
    pair: np.trapezoid(
        np.abs(age_profiles[pair[0]] - age_profiles[pair[1]]),
        age,
    )
    for pair in profile_pairs
}

check(
    "N5 年龄剖面两两 L1 差异为正",
    all(difference > 1e-2 for difference in profile_l1_differences.values()),
    f"L1={ {pair: round(float(value), 6) for pair, value in profile_l1_differences.items()} }",
)


# ==================================================================
head("N6  不同 likelihood ratio 的 L1 差异为正")
# ==================================================================
lambda_0 = -1.0
lambda_1 = 0.0
ratios = {
    p: np.exp((lambda_1 - lambda_0) * profile)
    for p, profile in age_profiles.items()
}
ratio_l1_differences = {
    pair: np.trapezoid(np.abs(ratios[pair[0]] - ratios[pair[1]]), age)
    for pair in profile_pairs
}

check(
    "N6 likelihood ratio 两两 L1 差异为正",
    all(difference > 1e-2 for difference in ratio_l1_differences.values()),
    f"L1={ {pair: round(float(value), 6) for pair, value in ratio_l1_differences.items()} }",
)


# ==================================================================
head("N7  p=1 给二次剖面，p=1/2 给线性剖面")
# ==================================================================
quadratic_reference = (radius / 2.0) * (1.0 - (age / period) ** 2)
linear_reference = (radius / 2.0) * (1.0 - age / period)

check(
    "N7 p=1 精确给二次剖面",
    np.max(np.abs(age_profiles[1.0] - quadratic_reference)) < 1e-12,
)
check(
    "N7 p=1/2 精确给线性剖面",
    np.max(np.abs(age_profiles[0.5] - linear_reference)) < 1e-12,
)


# ==================================================================
head("N8  半流共轭保持年龄推进律")
# ==================================================================
def age_flow(a, s, terminal):
    return np.minimum(a + s, terminal)


sample_ages = np.array([0.0, 0.2, 0.55, 0.9, 1.0])
shift = 0.15
conjugated_flow = {
    p: radial_maps[p][
        np.searchsorted(age, age_flow(sample_ages, shift, period), side="left")
    ]
    for p in p_values
}
direct_flow = {
    p: np.minimum(radial_maps[p][np.searchsorted(age, sample_ages, side="left")] + shift, radius)
    for p in p_values
}

# 共轭半流只保证保序与终端吸收；它不是径向平移。这里核验它确实不等于
# 简单的径向平移，从而暴露仿射映射需要额外假设。
non_affine_p = 2.0
conjugated_sample = conjugated_flow[non_affine_p]
direct_sample = direct_flow[non_affine_p]
check(
    "N8 非仿射共轭半流不是径向平移",
    np.max(np.abs(conjugated_sample - direct_sample)) > 1e-2,
)
check(
    "N8 非仿射共轭半流保持保序",
    np.all(np.diff(conjugated_sample) >= 0.0),
)


# ==================================================================
head("N9  文档登记重参数化无解结果")
# ==================================================================
check(
    "N9 文档登记重参数化无解结果",
    "R-Z-AGE-RADIAL-REPARAMETRIZATION-NOGO" in DOC
    and "\\text{所有 }\\phi_p\\text{ 都保持年龄序" in DOC,
)


# ==================================================================
head("N10  文档登记仿射标定缺口")
# ==================================================================
check(
    "N10 文档登记仿射标定缺口",
    "R-Z-AFFINE-AGE-RADIAL-CALIBRATION-GAP" in DOC
    and "等价于额外规定" in DOC,
)


# ==================================================================
head("N11  文档不把 p=1 写成已导出")
# ==================================================================
check(
    "N11 文档保留 p=1 的附加输入地位",
    "p=1\\text{ 是附加标定，不是现有年龄结构的自然后果。}" in DOC
    and "条件球核 + 仿射年龄标定" in DOC,
)


# ==================================================================
head("N12  上游边界保持")
# ==================================================================
check(
    "N12 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N13  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N13 D235 文档不含项目外体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


print("\n" + "=" * 78)
print(f"通过 {len(PASS)} / 不符 {len(FAIL)}")
print("=" * 78)
if FAIL:
    for name in FAIL:
        print(f"失败：{name}")
    sys.exit(1)
print("全部通过")
