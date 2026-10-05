#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D236 —— 源算子识别的条件嵌入与不唯一性
================================================================================
【被检验的问题（D236）】
  二维矩阵模板 K_M 能否条件嵌入局部源通道？
  源识别需要哪些谱与标度约束？
  为什么 D198 的标量源没有自动给出唯一通道？

【本步判据】
  N1  D236 与恢复结构已登记
  N2  二维源矩阵可以嵌入局域载体
  N3  均匀源通道给出 K_B=2pi I_B W S W^dagger
  N4  谱比例约束 mu0/mu1=lambda0/lambda1
  N5  总体标度约束 2pi I_B mu_i=lambda_i
  N6  不同局域嵌入可以给出相同压缩矩阵
  N7  非中心源扰动改变谱比例
  N8  势能平移改变压缩源通道
  N9  多场混合提供多个嵌入候选
  N10 文档登记条件嵌入与选择缺口
  N11 文档保留接触项缺口
  N12 文档不把源识别写成已导出
  N13 上游边界保持
  N14 文档不引用项目外体系名或外部路径

运行：python3 verify/d236_source_operator_identification_embedding.py
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
            "D236_source_operator_identification_embedding.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D236 与恢复结构已登记")
# ==================================================================
check(
    "N1 公理表登记 D236",
    "D236" in AX,
)
check(
    "N1 公理表与文档登记条件源嵌入",
    "R-Z-TWO-LEVEL-SOURCE-EMBEDDING" in BOTH
    and "条件嵌入" in DOC,
)
check(
    "N1 公理表与文档登记源通道选择缺口",
    "R-Z-SOURCE-EMBEDDING-SELECTION-GAP" in BOTH
    and "源通道" in DOC,
)


# ==================================================================
head("N2  二维源矩阵可以嵌入局域载体")
# ==================================================================
source_eigenvalues = np.array([-1.0, -2.0])
source_matrix = np.diag(source_eigenvalues).astype(complex)
local_dimension = 6

# 取前两个局域基矢为二维源子空间，其余基矢只用于展示非唯一嵌入。
W_one = np.zeros((local_dimension, 2), dtype=complex)
W_one[0, 0] = 1.0
W_one[1, 1] = 1.0
W_two = np.zeros((local_dimension, 2), dtype=complex)
W_two[2, 0] = 1.0
W_two[3, 1] = 1.0

local_source_one = W_one @ source_matrix @ W_one.conj().T
local_source_two = W_two @ source_matrix @ W_two.conj().T

check(
    "N2 两个嵌入都是等距嵌入",
    np.allclose(W_one.conj().T @ W_one, np.eye(2))
    and np.allclose(W_two.conj().T @ W_two, np.eye(2)),
)
check(
    "N2 局域源算子为厄米矩阵",
    np.allclose(local_source_one, local_source_one.conj().T)
    and np.allclose(local_source_two, local_source_two.conj().T),
)


# ==================================================================
head("N3  均匀源通道给出 K_B=2pi I_B W S W^dagger")
# ==================================================================
geometry_integral = 1.0
target_eigenvalues = 2.0 * np.pi * geometry_integral * source_eigenvalues
target_matrix = np.diag(target_eigenvalues).astype(complex)
reconstructed_matrix = 2.0 * np.pi * geometry_integral * source_matrix

check(
    "N3 均匀源通道精确重建目标矩阵",
    np.allclose(reconstructed_matrix, target_matrix),
)
check(
    "N3 重建保留同一二维源子空间",
    np.allclose(
        W_one @ reconstructed_matrix @ W_one.conj().T,
        W_one @ target_matrix @ W_one.conj().T,
    ),
)


# ==================================================================
head("N4  谱比例约束 mu0/mu1=lambda0/lambda1")
# ==================================================================
source_ratio = source_eigenvalues[0] / source_eigenvalues[1]
target_ratio = target_eigenvalues[0] / target_eigenvalues[1]

check(
    "N4 源谱比例与目标谱比例相等",
    abs(source_ratio - target_ratio) < 1e-12,
    f"source_ratio={source_ratio:.12f}",
)


# ==================================================================
head("N5  总体标度约束 2pi I_B mu_i=lambda_i")
# ==================================================================
scale_residuals = (
    2.0 * np.pi * geometry_integral * source_eigenvalues
    - target_eigenvalues
)
check(
    "N5 总体标度约束逐分量满足",
    np.max(np.abs(scale_residuals)) < 1e-12,
)


# ==================================================================
head("N6  不同局域嵌入可以给出相同压缩矩阵")
# ==================================================================
compressed_one = W_one.conj().T @ local_source_one @ W_one
compressed_two = W_two.conj().T @ local_source_two @ W_two

check(
    "N6 两个局域源算子不同",
    not np.allclose(local_source_one, local_source_two),
)
check(
    "N6 两个局域源算子压缩到各自子空间后相同",
    np.allclose(compressed_one, compressed_two),
)


# ==================================================================
head("N7  非中心源扰动改变谱比例")
# ==================================================================
noncentral_perturbation = np.diag([0.15, -0.25]).astype(complex)
perturbed_source = source_matrix + noncentral_perturbation
perturbed_ratio = perturbed_source[0, 0] / perturbed_source[1, 1]

check(
    "N7 非中心扰动改变谱比例",
    abs(perturbed_ratio - source_ratio) > 1e-2,
    f"ratio={perturbed_ratio:.6f}",
)


# ==================================================================
head("N8  势能平移改变压缩源通道")
# ==================================================================
potential_shift = 0.1 * geometry_integral * np.diag([1.0, 1.5])
shifted_potential_source = source_matrix - potential_shift
shifted_ratio = (
    shifted_potential_source[0, 0] / shifted_potential_source[1, 1]
)

check(
    "N8 势能平移改变压缩通道",
    not np.allclose(shifted_potential_source, source_matrix),
)
check(
    "N8 势能平移改变源谱比例",
    abs(shifted_ratio - source_ratio) > 1e-2,
    f"ratio={shifted_ratio:.6f}",
)


# ==================================================================
head("N9  多场混合提供多个嵌入候选")
# ==================================================================
mixing_angle = 0.35
mixing_vector_one = np.array([1.0, 0.0])
mixing_vector_two = np.array([np.cos(mixing_angle), np.sin(mixing_angle)])

check(
    "N9 两个多场混合向量不同",
    not np.allclose(mixing_vector_one, mixing_vector_two),
)
check(
    "N9 两个混合向量都可归一化",
    abs(np.linalg.norm(mixing_vector_one) - 1.0) < 1e-12
    and abs(np.linalg.norm(mixing_vector_two) - 1.0) < 1e-12,
)


# ==================================================================
head("N10  文档登记条件嵌入与选择缺口")
# ==================================================================
check(
    "N10 文档登记条件嵌入与选择缺口",
    "R-Z-TWO-LEVEL-SOURCE-EMBEDDING" in DOC
    and "R-Z-SOURCE-EMBEDDING-SELECTION-GAP" in DOC
    and "条件嵌入" in DOC,
)


# ==================================================================
head("N11  文档保留接触项缺口")
# ==================================================================
check(
    "N11 文档保留接触项缺口",
    "接触项仍不能由该嵌入吸收" in DOC
    and "多通道" in DOC,
)


# ==================================================================
head("N12  文档不把源识别写成已导出")
# ==================================================================
check(
    "N12 文档保留源识别缺口",
    "条件嵌入的存在性不等于源识别的唯一性" in DOC
    and "源通道仍未选择" in DOC,
)


# ==================================================================
head("N13  上游边界保持")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
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
    "N14 D236 文档不含项目外体系名或外部路径",
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
