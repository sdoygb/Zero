#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D241 —— 符号盲配对无解
================================================================================
【被检验的问题（D241）】
  符号盲配对是否能给符号年龄势？
  时间方向是否自动提供正负支耦合差？
  只有哪一种配对差能生成 D238 的线性危险率？

【本步判据】
  N1  D241 与恢复结构已登记
  N2  正负配对势差给出符号年龄势
  N3  符号盲配对给零危险率差
  N4  符号差耦合 K 控制符号比例
  N5  反对称部分在完整方块上抵消
  N6  只保留对称部分 S
  N7  均匀 S 给 k^2 势
  N8  终端中性给 D238 抛物型比例
  N9  时间反演不排除均匀 S
  N10 D220 反射配对给 K=0
  N11 非均匀 S 给不同高阶势
  N12 文档登记与上游边界
  N13 文档不引用项目外体系名或外部路径

运行：python3 verify/d241_sign_blind_pairing_no_go.py
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
            "D241_sign_blind_pairing_no_go.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D241 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D241", "D241" in AX)
check(
    "N1 公理表与文档登记符号盲无解",
    "R-Z-SIGN-BLIND-PAIR-NOGO" in BOTH,
)
check(
    "N1 公理表与文档登记符号差耦合",
    "R-Z-SIGN-DIFFERENTIAL-PAIR-COUPLING" in BOTH,
)
check(
    "N1 公理表与文档登记符号交换破缺缺口",
    "R-Z-PAIR-SIGN-EXCHANGE-BREAKING-GAP" in BOTH,
)


# ==================================================================
head("N2  正负配对势差给出符号年龄势")
# ==================================================================
size = 8
ages = np.arange(0, size + 1)

J_plus = 0.11 * np.ones((size, size))
J_minus = 0.04 * np.ones((size, size))


def block_sum(matrix, k):
    return float(matrix[:k, :k].sum())


phi_plus = np.array([block_sum(J_plus, k) for k in ages])
phi_minus = np.array([block_sum(J_minus, k) for k in ages])
K = J_plus - J_minus
phi_difference = phi_plus - phi_minus
phi_K = np.array([block_sum(K, k) for k in ages])

check(
    "N2 势差等于符号差耦合的势",
    np.allclose(phi_difference, phi_K),
)
check(
    "N2 符号年龄对数比例接势差",
    np.allclose(np.diff(phi_difference), np.diff(phi_K)),
)


# ==================================================================
head("N3  符号盲配对给零危险率差")
# ==================================================================
J_blind = 0.17 * np.ones((size, size))
phi_blind_plus = np.array([block_sum(J_blind, k) for k in ages])
phi_blind_minus = np.array([block_sum(J_blind, k) for k in ages])
blind_log_ratio = phi_blind_plus - phi_blind_minus
blind_hazard_difference = -np.diff(blind_log_ratio)

check(
    "N3 符号盲两支势相等",
    np.allclose(phi_blind_plus, phi_blind_minus),
)
check(
    "N3 符号盲危险率差为零",
    np.allclose(blind_hazard_difference, 0.0),
)


# ==================================================================
head("N4  符号差耦合 K 控制符号比例")
# ==================================================================
terminal_log_ratio = -phi_K[-1]
symbol_log_ratio = phi_K + terminal_log_ratio
symbol_ratio = np.exp(symbol_log_ratio)

check(
    "N4 符号比例只依赖 K",
    np.allclose(
        symbol_log_ratio,
        np.array([block_sum(K, k) for k in ages]) + terminal_log_ratio,
    ),
)
check(
    "N4 终端比例归零",
    abs(symbol_log_ratio[-1]) < 1e-12,
)


# ==================================================================
head("N5  反对称部分在完整方块上抵消")
# ==================================================================
antisymmetric = np.zeros((size, size))
antisymmetric[1, 4] = 0.07
antisymmetric[4, 1] = -0.07
antisymmetric[2, 6] = -0.03
antisymmetric[6, 2] = 0.03
antisymmetric_block_sums = np.array(
    [block_sum(antisymmetric, k) for k in ages]
)

check(
    "N5 反对称矩阵全方块和为零",
    abs(antisymmetric.sum()) < 1e-12,
)
check(
    "N5 所有前缀方块和为零",
    np.allclose(antisymmetric_block_sums, 0.0),
)


# ==================================================================
head("N6  只保留对称部分 S")
# ==================================================================
K_general = K.copy()
K_general[2, 5] += 0.09
K_general[5, 2] -= 0.04
S = 0.5 * (K_general + K_general.T)
A = 0.5 * (K_general - K_general.T)
phi_general = np.array([block_sum(K_general, k) for k in ages])
phi_symmetric = np.array([block_sum(S, k) for k in ages])
phi_antisymmetric = np.array([block_sum(A, k) for k in ages])

check(
    "N6 一般 K 分解为 S+A",
    np.allclose(K_general, S + A),
)
check(
    "N6 符号势只等于对称部分",
    np.allclose(phi_general, phi_symmetric),
)
check(
    "N6 反对称部分势为零",
    np.allclose(phi_antisymmetric, 0.0),
)


# ==================================================================
head("N7  均匀 S 给 k^2 势")
# ==================================================================
delta = 0.052
S_uniform = -delta * np.ones((size, size))
phi_uniform = np.array(
    [block_sum(S_uniform, k) for k in ages]
)

check(
    "N7 均匀负对称耦合给分支势 -delta k^2",
    np.allclose(phi_uniform, -delta * ages**2),
)
check(
    "N7 生成势差给 2k+1",
    np.allclose(
        -np.diff(phi_uniform),
        delta * (2 * ages[:-1] + 1),
    ),
)


# ==================================================================
head("N8  终端中性给 D238 抛物型比例")
# ==================================================================
terminal_age = size
terminal_compensation = -phi_uniform[-1]
log_ratio = phi_uniform + terminal_compensation
ratio = np.exp(log_ratio)
target = np.exp(delta * (terminal_age**2 - ages**2))

check(
    "N8 终端补偿恢复 D238 抛物线",
    np.allclose(ratio, target),
)
check(
    "N8 危险率差线性增长",
    np.allclose(
        -np.diff(log_ratio),
        2.0 * delta * ages[:-1] + delta,
    ),
)


# ==================================================================
head("N9  时间反演不排除均匀 S")
# ==================================================================
reversed_indices = np.arange(size - 1, -1, -1)
S_reversed = S_uniform[np.ix_(reversed_indices, reversed_indices)]

check(
    "N9 均匀 S 在时间反演下不变",
    np.allclose(S_reversed, S_uniform),
)
check(
    "N9 时间反演不产生符号差",
    np.allclose(S_reversed - S_uniform, 0.0),
)


# ==================================================================
head("N10  D220 反射配对给 K=0")
# ==================================================================
J_reflection_plus = 0.13 * np.ones((size, size))
J_reflection_minus = J_reflection_plus.copy()
K_reflection = J_reflection_plus - J_reflection_minus
phi_reflection = np.array(
    [block_sum(K_reflection, k) for k in ages]
)

check(
    "N10 反射配对给 K=0",
    np.allclose(K_reflection, 0.0),
)
check(
    "N10 反射配对给零年龄势",
    np.allclose(phi_reflection, 0.0),
)


# ==================================================================
head("N11  非均匀 S 给不同高阶势")
# ==================================================================
S_nonuniform = S_uniform.copy()
S_nonuniform[5, 5] += 0.06
phi_uniform_small = np.array(
    [block_sum(S_uniform, k) for k in range(1, size + 1)]
)
phi_nonuniform_small = np.array(
    [block_sum(S_nonuniform, k) for k in range(1, size + 1)]
)

check(
    "N11 低阶势保持不变",
    np.allclose(phi_uniform_small[:5], phi_nonuniform_small[:5]),
)
check(
    "N11 高阶级势分叉",
    not np.allclose(phi_uniform_small, phi_nonuniform_small),
)


# ==================================================================
head("N12  文档登记与上游边界")
# ==================================================================
check(
    "N12 文档登记符号盲无解",
    "R-Z-SIGN-BLIND-PAIR-NOGO" in DOC
    and "符号盲配对给零危险率差" in DOC,
)
check(
    "N12 文档登记符号差耦合",
    "R-Z-SIGN-DIFFERENTIAL-PAIR-COUPLING" in DOC
    and "K=J^+-J^-" in DOC,
)
check(
    "N12 文档登记符号交换破缺缺口",
    "R-Z-PAIR-SIGN-EXCHANGE-BREAKING-GAP" in DOC
    and "符号交换破缺来源仍缺" in DOC,
)
check(
    "N12 文档不把新耦合写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "至少要求符号交换破缺" in DOC,
)
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
    "N13 D241 文档不含项目外体系名或外部路径",
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
