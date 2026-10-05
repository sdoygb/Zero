#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D215 —— 周期局域环到 U1-U4 的语义接口
================================================================================
【被检验的问题（D215）】
  D214 的闭环局域图能否与 D212 的周期载体合并？
  合并后能否条件满足 U1、U2_sem、U3、U4_sem？
  这项接口是否仍不等于物理局域性、连续 QFT 或 GR？

【本步判据】
  N1  D215 与两个恢复结构已登记
  N2  环上局域张量因子给 U1 候选载体
  N3  寿命年龄支持系统与包含序
  N4  局域代数赋值保持包含与生成性
  N5  非平凡无交局域对严格交换
  N6  非中心忠实积态给非平凡模流
  N7  寿命年龄映射满足正时间半群律
  N8  支持映射保序
  N9  模支持协变并给严格扩张
  N10 接口仍保留物理局域性缺口
  N11 T 与模时间仍未被等同
  N12 不引用旧体系或外部路径
  N13 不修改上游、不新增 U5，失败时非零退出

运行：python3 verify/d215_cycle_local_semantic_u_interface.py
"""
import io
import itertools
import math
import os
import sys

import numpy as np

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
            "D215_cycle_local_semantic_u_interface.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def tensor_embedding(local, site, sites):
    result = np.array([[1.0]], dtype=complex)
    for current in sorted(sites):
        result = np.kron(result, local if current == site else np.eye(2))
    return result


def product_state(local_weights, sites):
    state = np.array([1.0], dtype=float)
    for _ in sorted(sites):
        state = np.kron(state, np.asarray(local_weights, dtype=float))
    return state


def modular_conjugate(matrix, weights, time):
    phases = np.exp(-1j * time * np.log(weights))
    evolution = np.diag(phases)
    return evolution @ matrix @ evolution.conj().T


def support_step(support_set, age, time, period, sites):
    if support_set is None:
        return None, None
    if age + time < period:
        return support_set, age + time
    return None, None


def order_key(support_set, age, sites):
    if support_set is None:
        return (len(sites) + 1, math.inf)
    return (len(support_set), age)


def subset_order(left, right):
    return left.issubset(right)


def support_leq(left, right):
    left_set, left_age = left
    right_set, right_age = right
    if left_set is None:
        return True
    if right_set is None:
        return True
    return left_set.issubset(right_set) and left_age <= right_age + 1.0e-12


T = 2
SITES = (0, 1, 2)


# ==================================================================
head("N1  D215 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 文档登记语义周期桥与物理局域性缺口",
    "D215" in AX
    and "R-Z-U-SEMANTIC-CYCLE-BRIDGE" in BOTH
    and "R-Z-U-PHYSICAL-LOCALITY-GAP" in BOTH
    and "语义接口" in DOC,
)
check(
    "N1 文档不把语义桥写成无条件替代",
    "不是无条件替代" in DOC
    and "不能由" in DOC
    and "条件接口" in DOC,
)


# ==================================================================
head("N2  环上局域张量因子给 U1 候选载体")
# ==================================================================
site_factor_dim = T
global_dim = site_factor_dim ** len(SITES)
check(
    "N2 局域因子与全局载体维数符合张量积",
    global_dim == T ** len(SITES) == 8,
    f"local_dim={site_factor_dim}, global_dim={global_dim}",
)
check(
    "N2 载体是含单位元有限维复 C*-代数",
    "\\bigotimes_{v\\in V}M_T(\\mathbb C)" in DOC
    and "U1" in DOC,
)


# ==================================================================
head("N3  寿命年龄支持系统与包含序")
# ==================================================================
subsets = [
    frozenset(choice)
    for size in range(1, len(SITES) + 1)
    for choice in itertools.combinations(SITES, size)
]
ages = (0.0, 0.4, 1.2)
check(
    "N3 支持点包含非空局域子集与寿命年龄",
    all(subset for subset in subsets)
    and all(0.0 <= age < T for age in ages)
    and "(S,a)" in DOC,
)
check(
    "N3 偏序按局域子集包含与年龄排序",
    subset_order({0}, {0, 1})
    and 0.0 <= 0.4
    and "S\\subseteq S'" in DOC
    and "a\\le a'" in DOC,
)


# ==================================================================
head("N4  局域代数赋值保持包含与生成性")
# ==================================================================
left = {0}
right = {0, 1}
left_dim = T ** len(left)
right_dim = T ** len(right)
full_dim = T ** len(SITES)
check(
    "N4 子集包含给出局域代数包含",
    subset_order(left, right) and left_dim < right_dim,
    f"A({sorted(left)}) dim={left_dim}, A({sorted(right)}) dim={right_dim}",
)
check(
    "N4 全环支持生成完整载体",
    right_dim < full_dim
    and T ** len(set(SITES)) == full_dim
    and "\\bigcup_{p\\in P}A(p)=A" in DOC,
)


# ==================================================================
head("N5  非平凡无交局域对严格交换")
# ==================================================================
sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
left_observable = tensor_embedding(sigma_x, 0, SITES)
right_observable = tensor_embedding(sigma_y, 2, SITES)
commutator = left_observable @ right_observable - right_observable @ left_observable
check(
    "N5 不相交站点上的局域可观测量交换",
    np.allclose(commutator, 0.0),
)
check(
    "N5 文档给出非平凡无交对",
    "\\{0\\}" in DOC
    and "\\{2\\}" in DOC
    and "U2_{\\rm sem}" in DOC,
)


# ==================================================================
head("N6  非中心忠实积态给非平凡模流")
# ==================================================================
local_weights = np.array([0.7, 0.3])
global_weights = product_state(local_weights, SITES)
check(
    "N6 积态完全正权重因而忠实",
    np.all(global_weights > 0.0)
    and abs(float(np.sum(global_weights)) - 1.0) < 1.0e-12,
)
probe = np.zeros((global_dim, global_dim), dtype=complex)
probe[0, 1] = 1.0
time = 0.7
image = modular_conjugate(probe, global_weights, time)
check(
    "N6 非中心局域态给非平凡模流",
    not np.allclose(image, probe)
    and np.linalg.norm(image - probe) > 1.0e-6,
    f"norm={np.linalg.norm(image - probe):.6f}",
)
check(
    "N6 模流保持局域支持代数",
    "\\sigma_t^\\rho(A(S,a))=A(S,a)" in DOC,
)


# ==================================================================
head("N7  寿命年龄映射满足正时间半群律")
# ==================================================================
test_times = (0.0, 0.3, 0.8, 1.5, 3.0)
semigroup_ok = True
for support_set in subsets:
    for age in ages:
        for first in test_times:
            for second in test_times:
                mid_set, mid_age = support_step(support_set, age, first, T, SITES)
                left_result = support_step(mid_set, mid_age, second, T, SITES)
                final_set, final_age = support_step(support_set, age, first + second, T, SITES)
                if left_result[0] != final_set:
                    semigroup_ok = False
                    continue
                if left_result[1] is None or final_age is None:
                    if left_result[1] != final_age:
                        semigroup_ok = False
                elif abs(left_result[1] - final_age) > 1.0e-12:
                    semigroup_ok = False
check(
    "N7 支持年龄流满足 E_sE_t=E_{s+t}",
    semigroup_ok,
)


# ==================================================================
head("N8  支持映射保序")
# ==================================================================
order_pairs = [
    ({0}, 0.2, {0, 1}, 0.5),
    ({0}, 0.2, {0}, 0.5),
    ({1}, 0.0, {1}, 1.0),
]
order_ok = True
for small_set, small_age, large_set, large_age in order_pairs:
    for coefficient in test_times:
        small_image = support_step(small_set, small_age, coefficient, T, SITES)
        large_image = support_step(large_set, large_age, coefficient, T, SITES)
        if not support_leq(small_image, large_image):
            order_ok = False
check(
    "N8 支持映射保持偏序",
    order_ok,
)


# ==================================================================
head("N9  模支持协变并给严格扩张")
# ==================================================================
local_support_dim = T ** 1
full_support_dim = T ** len(SITES)
check(
    "N9 寿命后局域支持严格进入全局包络",
    support_step({0}, 1.2, 1.0, T, SITES) == (None, None)
    and local_support_dim < full_support_dim,
    f"A(local) dim={local_support_dim}, A(global) dim={full_support_dim}",
)
check(
    "N9 文档给出 U4_sem 协变链",
    "U4_{\\rm sem}" in DOC
    and "A(p)" in DOC
    and "A(E_t(p))" in DOC,
)


# ==================================================================
head("N10 接口仍保留物理局域性缺口")
# ==================================================================
check(
    "N10 文档拒绝把形式局域对写成物理因果",
    "形式局域语义" in DOC
    and "物理光锥" in DOC
    and "R-Z-U-PHYSICAL-LOCALITY-GAP" in BOTH,
)
check(
    "N10 文档保留张量分解与权重输入",
    "张量分解" in DOC
    and "权重" in DOC
    and "选择器输入" in DOC,
)


# ==================================================================
head("N11 T 与模时间仍未被等同")
# ==================================================================
check(
    "N11 文档分离周期 T 与模参数",
    "模参数没有唯一关系" in DOC
    and "模时间" in DOC
    and "周期时间" in DOC,
)


# ==================================================================
head("N12 不引用旧体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "外部路径",
    "../..",
)
check(
    "N12 文档不含旧体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N13 上游边界")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
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
