#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D233 —— 符号年龄对称性无解
================================================================================
【被检验的问题（D233）】
  当前正负延拓重播种能否产生非恒定剖面？
  每个年龄的正负活动计数是否相等？
  非恒定 boost 剖面是否必须增加符号年龄破缺？

【本步判据】
  N1  D233 与两个恢复结构已登记
  N2  每条未闭合路径有符号翻转配对
  N3  每个年龄的正负活动计数相等
  N4  年龄条件符号比例恒为 1/2
  N5  D232 剖面还原为恒定零
  N6  符号年龄协方差为零
  N7  人为偏置 seed 可产生非零剖面
  N8  文档登记符号年龄对称性
  N9  文档登记剖面破缺缺口
  N10 文档不把破缺来源写成已导出
  N11 上游边界保持
  N12 文档不引用项目外体系或外部路径

运行：python3 verify/d233_sign_age_symmetry_no_go_for_profile.py
"""
import io
import itertools
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
            "D233_sign_age_symmetry_no_go_for_profile.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D233 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记符号年龄对称性",
    "D233" in AX
    and "R-Z-SIGN-AGE-SYMMETRY" in BOTH
    and "符号年龄对称性" in DOC,
)
check(
    "N1 公理表与文档登记剖面破缺缺口",
    "R-Z-PROFILE-SYMMETRY-BREAKING-GAP" in BOTH
    and "符号年龄破缺" in DOC,
)


# ==================================================================
head("N2  每条未闭合路径有符号翻转配对")
# ==================================================================
max_steps = 10
paths = list(itertools.product([+1, -1], repeat=max_steps))
path_by_tuple = {path: np.cumsum(path) for path in paths}
pairing_ok = all(
    np.array_equal(path_by_tuple[path], -path_by_tuple[tuple(-step for step in path)])
    for path in paths
)
check(
    "N2 所有路径都有符号翻转配对",
    pairing_ok,
)


# ==================================================================
head("N3  每个年龄的正负活动计数相等")
# ==================================================================
positive_counts = {}
negative_counts = {}
for age in range(1, max_steps + 1):
    positive = 0
    negative = 0
    for path in paths:
        partial_sum = int(np.sum(path[:age]))
        if partial_sum > 0:
            positive += 1
        elif partial_sum < 0:
            negative += 1
    positive_counts[age] = positive
    negative_counts[age] = negative

check(
    "N3 每个年龄正负活动计数相等",
    all(
        positive_counts[age] == negative_counts[age]
        for age in range(1, max_steps + 1)
    ),
    f"counts={[(age, positive_counts[age], negative_counts[age]) for age in range(1, max_steps + 1)]}",
)


# ==================================================================
head("N4  年龄条件符号比例恒为 1/2")
# ==================================================================
conditional_plus = {
    age: positive_counts[age] / (positive_counts[age] + negative_counts[age])
    for age in range(1, max_steps + 1)
}
check(
    "N4 年龄条件正号比例恒为 1/2",
    all(abs(probability - 0.5) < 1e-12 for probability in conditional_plus.values()),
)


# ==================================================================
head("N5  D232 剖面还原为恒定零")
# ==================================================================
r_matrix = 0.75
lambda_0 = np.log(r_matrix)
lambda_1 = np.log(1.0 - r_matrix)
recovered_profiles = {
    age: np.log(
        conditional_plus[age] / (1.0 - conditional_plus[age])
    ) / (lambda_1 - lambda_0)
    for age in range(1, max_steps + 1)
}
check(
    "N5 还原剖面为零",
    all(abs(profile) < 1e-12 for profile in recovered_profiles.values()),
)


# ==================================================================
head("N6  符号年龄协方差为零")
# ==================================================================
age_values = np.array(list(range(1, max_steps + 1)), dtype=float)
imbalance = np.array(
    [
        2.0 * conditional_plus[age] - 1.0
        for age in range(1, max_steps + 1)
    ]
)
covariance = float(
    np.mean((age_values - np.mean(age_values)) * (imbalance - np.mean(imbalance)))
)
check(
    "N6 符号年龄协方差为零",
    abs(covariance) < 1e-12,
    f"cov={covariance:.3e}",
)


# ==================================================================
head("N7  人为偏置 seed 可产生非零剖面")
# ==================================================================
bias_strength = 0.2
biased_plus = {
    age: conditional_plus[age] + bias_strength
    for age in range(1, max_steps + 1)
}
biased_profile = {
    age: np.log(
        biased_plus[age] / (1.0 - biased_plus[age])
    ) / (lambda_1 - lambda_0)
    for age in range(1, max_steps + 1)
}
check(
    "N7 人为偏置给非零但恒定剖面",
    all(abs(profile - next(iter(biased_profile.values()))) < 1e-12
        for profile in biased_profile.values())
    and abs(next(iter(biased_profile.values()))) > 1e-3,
    f"bias_profile={next(iter(biased_profile.values())):.6f}",
)


# ==================================================================
head("N8  文档登记符号年龄对称性")
# ==================================================================
check(
    "N8 文档登记对称性结果",
    "R-Z-SIGN-AGE-SYMMETRY" in DOC
    and "N_+(k)=N_-(k)" in DOC
    and "p_+(k)" in DOC,
)


# ==================================================================
head("N9  文档登记剖面破缺缺口")
# ==================================================================
check(
    "N9 文档登记破缺缺口",
    "R-Z-PROFILE-SYMMETRY-BREAKING-GAP" in DOC
    and "终端吸收对正负支不对称" in DOC
    and "几何时间定向" in DOC,
)


# ==================================================================
head("N10  文档不把破缺来源写成已导出")
# ==================================================================
check(
    "N10 文档保留未解来源",
    "这些都尚未由当前零动力学导出" in DOC
    and "新的缺口是符号年龄相关性的来源" in DOC,
)


# ==================================================================
head("N11  上游边界保持")
# ==================================================================
check(
    "N11 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N12  文档不引用项目外体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N12 D233 文档不含项目外体系名或外部路径",
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
