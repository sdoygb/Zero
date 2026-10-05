#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D231 —— 局部模密度剖面缺口
================================================================================
【被检验的问题（D231）】
  区间支持能否决定支持内部的模密度剖面？
  常量、线性退化、二次退化剖面是否给出不同生成元？
  几何 boost 识别还缺什么？

【本步判据】
  N1  D231 与两个恢复结构已登记
  N2  三种剖面共享同一支持
  N3  三种剖面都给出非零局部模生成元
  N4  模流固定三种生成元
  N5  三种生成元谱不同
  N6  常量锐利剖面不满足边界退化
  N7  线性与二次剖面满足边界退化
  N8  同一支持上的剖面差异可测得
  N9  文档登记模密度剖面
  N10 文档登记 boost 剖面缺口
  N11 文档不把剖面选择写成 U1-U4 定理
  N12 文档不声称几何 boost 已识别
  N13 上游边界保持
  N14 文档不引用项目外体系或外部路径

运行：python3 verify/d231_modular_density_profile_gap.py
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
            "D231_modular_density_profile_gap.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D231 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记模密度剖面",
    "D231" in AX
    and "R-Z-MODULAR-DENSITY-PROFILE" in BOTH
    and "K(f)" in DOC,
)
check(
    "N1 公理表与文档登记 boost 剖面缺口",
    "R-Z-BOOST-PROFILE-GAP" in BOTH
    and "f_B" in DOC,
)


# ==================================================================
head("N2  三种剖面共享同一支持")
# ==================================================================
point_count = 2001
x = np.linspace(0.0, 1.0, point_count)
interior = x < 1.0
support_interval = interior

f_constant = np.ones_like(x)
f_linear = 1.0 - x
f_quadratic = 1.0 - x**2

profiles = {
    "constant": f_constant,
    "linear": f_linear,
    "quadratic": f_quadratic,
}

check(
    "N2 三种剖面都在同一支持上非零",
    all(np.all(profile[interior] > 0.0) for profile in profiles.values()),
)
check(
    "N2 三种剖面几乎处处支持区间相同",
    all(
        np.array_equal(profile[interior] > 0.0, support_interval[interior])
        for profile in profiles.values()
    ),
)


# ==================================================================
head("N3  三种剖面都给出非零局部模生成元")
# ==================================================================
r_matrix = 0.75
lambda_0 = np.log(r_matrix)
lambda_1 = np.log(1.0 - r_matrix)
generator_norms = {
    name: np.sqrt(
        np.trapezoid(
            (lambda_0 * profile) ** 2
            + (lambda_1 * profile) ** 2,
            x,
        )
    )
    for name, profile in profiles.items()
}
check(
    "N3 三种剖面生成元范数非零",
    all(norm > 1e-10 for norm in generator_norms.values()),
    f"norms={ {name: round(float(norm), 8) for name, norm in generator_norms.items()} }",
)


# ==================================================================
head("N4  模流固定三种生成元")
# ==================================================================
phase_rate = np.log(r_matrix / (1.0 - r_matrix))
phase = np.exp(1j * phase_rate)
sample_index = int(np.argmin(np.abs(x - 0.5)))
flow_matrix = np.diag([1.0, phase]).astype(complex)


def sample_generator(profile):
    return np.diag(
        [
            lambda_0 * profile[sample_index],
            lambda_1 * profile[sample_index],
        ]
    ).astype(complex)


check(
    "N4 模流固定生成元",
    all(
        np.allclose(
            flow_matrix @ sample_generator(profile) @ flow_matrix.conj().T,
            sample_generator(profile),
        )
        for profile in profiles.values()
    ),
)


# ==================================================================
head("N5  三种生成元谱不同")
# ==================================================================

def generator_levels(profile):
    return np.concatenate([lambda_0 * profile, lambda_1 * profile])


level_sets = {
    name: generator_levels(profile)
    for name, profile in profiles.items()
}

pairwise_differences = {
    "constant_linear": np.max(
        np.abs(np.sort(level_sets["constant"]) - np.sort(level_sets["linear"]))
    ),
    "constant_quadratic": np.max(
        np.abs(np.sort(level_sets["constant"]) - np.sort(level_sets["quadratic"]))
    ),
    "linear_quadratic": np.max(
        np.abs(np.sort(level_sets["linear"]) - np.sort(level_sets["quadratic"]))
    ),
}
check(
    "N5 三种生成元谱两两不同",
    all(difference > 1e-3 for difference in pairwise_differences.values()),
    f"max_diffs={ {name: round(float(value), 6) for name, value in pairwise_differences.items()} }",
)


# ==================================================================
head("N6  常量锐利剖面不满足边界退化")
# ==================================================================
constant_boundary_value = f_constant[-1]
check(
    "N6 常量剖面边界值不为零",
    abs(constant_boundary_value) > 1e-12,
    f"boundary={constant_boundary_value:.3f}",
)


# ==================================================================
head("N7  线性与二次剖面满足边界退化")
# ==================================================================
check(
    "N7 线性与二次剖面右边界为零",
    abs(f_linear[-1]) < 1e-12 and abs(f_quadratic[-1]) < 1e-12,
)
check(
    "N7 退化剖面左边界非零",
    abs(f_linear[0]) > 1e-12 and abs(f_quadratic[0]) > 1e-12,
)


# ==================================================================
head("N8  同一支持上的剖面差异可测得")
# ==================================================================
linear_gap = np.trapezoid(np.abs(f_constant - f_linear), x)
quadratic_gap = np.trapezoid(np.abs(f_constant - f_quadratic), x)
check(
    "N8 剖面 L1 差异为正",
    linear_gap > 1e-3 and quadratic_gap > 1e-3,
    f"L1=({linear_gap:.6f},{quadratic_gap:.6f})",
)


# ==================================================================
head("N9  文档登记模密度剖面")
# ==================================================================
check(
    "N9 文档登记模密度剖面",
    "R-Z-MODULAR-DENSITY-PROFILE" in DOC
    and "K(f)=K_M\\otimes f" in DOC
    and "支持内剖面决定生成元" in DOC,
)


# ==================================================================
head("N10  文档登记 boost 剖面缺口")
# ==================================================================
check(
    "N10 文档登记 boost 剖面缺口",
    "R-Z-BOOST-PROFILE-GAP" in DOC
    and "边界退化" in DOC
    and "几何 boost" in DOC,
)


# ==================================================================
head("N11  文档不把剖面选择写成 U1-U4 定理")
# ==================================================================
check(
    "N11 文档保留选择器地位",
    "U1-U4` 不包含" in DOC
    and "模密度剖面是新的恢复层输入" in DOC,
)


# ==================================================================
head("N12  文档不声称几何 boost 已识别")
# ==================================================================
check(
    "N12 文档保留几何识别缺口",
    "本文没有导出上述任一项" in DOC
    and "缺口已从" in DOC
    and "boost 剖面" in DOC,
)


# ==================================================================
head("N13  上游边界保持")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N14  文档不引用项目外体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N14 D231 文档不含项目外体系名或外部路径",
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
