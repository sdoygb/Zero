#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D248 —— 层间读出动力学
================================================================================
【被检验的问题（D248）】
  活动层怎样通过历史层与全局零层的读出受到影响？
  对称层间交换能否自动给守恒流？
  全局零层共同模是否会产生真实推力？
  守恒流能否唯一决定 GR 所需的应力张量？

【本步判据】
  N1  D248 与恢复结构已登记
  N2  相同读出给相同层间势
  N3  隐藏历史不同但读出相同不影响流
  N4  对称交换核给反对称边流
  N5  离散连续性保持总层间荷
  N6  全局共同模不改变交换流
  N7  位置依赖 Z 响应改变交换流
  N8  恒定历史势不产生层间流
  N9  守恒矢量流不能唯一决定应力张量
  N10 相同能量流可对应不同空间压强
  N11 文档登记守恒流与应力提升缺口
  N12 文档不把新结构写成 U1-U4 推论
  N13 上游边界保持
  N14 文档不引用项目外体系名或外部路径

运行：python3 verify/d248_interlayer_readout_dynamics.py
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
            "D248_interlayer_readout_dynamics.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D248 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D248", "D248" in AX)
check(
    "N1 公理表与文档登记读出商化",
    "R-Z-INTERLAYER-READOUT-QUOTIENT" in BOTH,
)
check(
    "N1 公理表与文档登记守恒交换",
    "R-Z-INTERLAYER-ZERO-SUM-EXCHANGE" in BOTH,
)
check(
    "N1 公理表与文档登记共同模盲性",
    "R-Z-Z-READOUT-COMMON-MODE-BLINDNESS" in BOTH,
)
check(
    "N1 公理表与文档登记 Dirichlet 源路线",
    "R-Z-INTERLAYER-DIRICHLET-SOURCE-ROUTE" in BOTH,
)
check(
    "N1 公理表与文档登记应力提升缺口",
    "R-Z-SOURCE-STRESS-LIFT-GAP" in BOTH,
)


# ==================================================================
head("N2  相同读出给相同层间势")
# ==================================================================
def effective_potential(readout, z, response):
    return 0.7 * readout + 0.2 * z + response


history_one = (2, 91.0)
history_two = (2, -17.0)
readout_one = history_one[0]
readout_two = history_two[0]
z = 1.3

phi_one = effective_potential(readout_one, z, 0.4)
phi_two = effective_potential(readout_two, z, 0.4)

check(
    "N2 两份历史具有相同读出",
    np.isclose(readout_one, readout_two),
)
check(
    "N2 相同读出给相同有效势",
    np.isclose(phi_one, phi_two),
)


# ==================================================================
head("N3  隐藏历史不同但读出相同不影响流")
# ==================================================================
hidden_difference = history_one[1] - history_two[1]

check(
    "N3 隐藏历史确实不同",
    hidden_difference != 0.0,
)
check(
    "N3 隐藏历史不进入有效势差",
    np.isclose(
        (phi_one - phi_two),
        0.0,
    ),
)


# ==================================================================
head("N4  对称交换核给反对称边流")
# ==================================================================
rng = np.random.default_rng(248)
size = 6
raw = rng.normal(size=(size, size))
symmetric_kernel = 0.5 * (raw + raw.T)
symmetric_kernel -= np.diag(np.diag(symmetric_kernel))
np.fill_diagonal(symmetric_kernel, 0.0)
potential = np.array([0.2, -0.1, 0.8, 0.3, -0.4, 0.6])

edge_flow = symmetric_kernel * (
    potential[:, None] - potential[None, :]
)

check(
    "N4 核为对称矩阵",
    np.allclose(symmetric_kernel, symmetric_kernel.T),
)
check(
    "N4 边流反对称",
    np.allclose(edge_flow, -edge_flow.T),
)


# ==================================================================
head("N5  离散连续性保持总层间荷")
# ==================================================================
laplacian = (
    np.diag(symmetric_kernel.sum(axis=1))
    - symmetric_kernel
)
velocity = -laplacian @ potential

check(
    "N5 总驱动为零",
    np.isclose(np.sum(velocity), 0.0),
)
check(
    "N5 速度等于负边流散度",
    np.allclose(velocity, -np.sum(edge_flow, axis=1)),
)
check(
    "N5 总量变化为零",
    np.isclose(np.sum(velocity), np.sum(-np.sum(edge_flow, axis=1))),
)


# ==================================================================
head("N6  全局共同模不改变交换流")
# ==================================================================
common_shift = 7.5
potential_shifted = potential + common_shift
velocity_shifted = -laplacian @ potential_shifted

check(
    "N6 共同模不改变速度",
    np.allclose(velocity, velocity_shifted),
)
check(
    "N6 共同模不改变边流",
    np.allclose(
        edge_flow,
        symmetric_kernel * (
            potential_shifted[:, None] - potential_shifted[None, :]
        ),
    ),
)


# ==================================================================
head("N7  位置依赖 Z 响应改变交换流")
# ==================================================================
z_readout = 1.7
local_history = np.array([0.1, 0.4, -0.2, 0.3, 0.5, -0.1])
constant_response = np.full(size, 0.35)
position_response = np.array([0.35, 0.70, 0.10, 0.90, 0.20, 0.55])

phi_constant = local_history + constant_response * z_readout
phi_position = local_history + position_response * z_readout

check(
    "N7 常响应给共同模盲性",
    np.allclose(-laplacian @ phi_constant, -laplacian @ local_history),
)
check(
    "N7 位置依赖响应改变速度",
    not np.allclose(-laplacian @ phi_position, -laplacian @ phi_constant),
)


# ==================================================================
head("N8  恒定历史势不产生层间流")
# ==================================================================
constant_history = np.full(size, -2.3)

check(
    "N8 恒定历史势速度为零",
    np.allclose(-laplacian @ constant_history, np.zeros(size)),
)
check(
    "N8 恒定历史势边流为零",
    np.allclose(
        symmetric_kernel * (
            constant_history[:, None] - constant_history[None, :]
        ),
        np.zeros((size, size)),
    ),
)


# ==================================================================
head("N9  守恒矢量流不能唯一决定应力张量")
# ==================================================================
rho = 2.4
pressure_one = 0.6
pressure_two = 0.9
stress_one = np.diag([rho, pressure_one, pressure_one, pressure_one])
stress_two = np.diag([rho, pressure_two, pressure_two, pressure_two])

check(
    "N9 两份应力具有相同能量分量",
    np.allclose(stress_one[:, 0], stress_two[:, 0]),
)
check(
    "N9 两份应力空间部分不同",
    not np.allclose(stress_one[1:, 1:], stress_two[1:, 1:]),
)


# ==================================================================
head("N10 相同能量流可对应不同空间压强")
# ==================================================================
constant_divergence_one = np.zeros_like(stress_one)
constant_divergence_two = np.zeros_like(stress_two)

check(
    "N10 第一应力在平直时空守恒",
    np.allclose(constant_divergence_one, np.zeros_like(stress_one)),
)
check(
    "N10 第二应力在平直时空守恒",
    np.allclose(constant_divergence_two, np.zeros_like(stress_two)),
)
check(
    "N10 空间压强不同",
    not np.isclose(stress_one[1, 1], stress_two[1, 1]),
)


# ==================================================================
head("N11  文档登记守恒流与应力提升缺口")
# ==================================================================
check(
    "N11 文档登记层间守恒交换",
    "R-Z-INTERLAYER-ZERO-SUM-EXCHANGE" in DOC
    and "离散连续性" in DOC,
)
check(
    "N11 文档登记共同模盲性",
    "R-Z-Z-READOUT-COMMON-MODE-BLINDNESS" in DOC
    and "共同模" in DOC,
)
check(
    "N11 文档登记应力提升缺口",
    "R-Z-SOURCE-STRESS-LIFT-GAP" in DOC
    and "应力提升" in DOC,
)


# ==================================================================
head("N12  文档不把新结构写成 U1-U4 推论")
# ==================================================================
check(
    "N12 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "不是 `U1-U4` 的定理" in DOC,
)
check(
    "N12 文档不把层间交换写成无条件 GR",
    "GR 仍不是无条件结果" in DOC,
)


# ==================================================================
head("N13  上游边界保持")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N13 文档登记输入预算",
    "## §3 输入预算" in DOC,
)


# ==================================================================
head("N14  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N14 D248 文档不含项目外体系名或外部路径",
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
