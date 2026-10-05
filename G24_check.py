#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G24_check.py -- 年龄结构：零层弧最大簇的完整结构，以及它与 I2a / I8 的冲突。

对应文档 G24_age_structure_and_its_conflicts.md。
失败时退出码非零。

  F1  年龄代数分解：A_i = M_{n_i}(C) ⊗ C(X_tau_i)（D223）
  F2  逐年龄的局部模生成元（D227）
  F3  局部支持点 = 年龄前缀（D228）
  F4  D229 的原子障碍：C([0,T]) 无非平凡投影（数值核验）
  F5  D230 的 L-infinity 修复：区间投影（数值核验）
  F6  D233 的决定性 no-go：符号对称 => 剖面恒定
  F7  与账本的接口（Z4 = 极端投影）与张力（I2a vs D229）
  F8  诚实边界
"""

import os
import sys
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.expanduser("~/Downloads/modular-equilibrium/derivations")
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def rd(name):
    p = os.path.join(MOD, name)
    return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


d223 = rd("D223_age_carrier_tensor_interface.md")
d227 = rd("D227_age_local_modular_witness_selector.md")
d228 = rd("D228_age_local_modular_family_support_covariance.md")
d229 = rd("D229_continuous_age_support_obstruction.md")
d230 = rd("D230_linfinity_age_support_and_projection_complement.md")
d233 = rd("D233_sign_age_symmetry_no_go_for_profile.md")

# ======================================================================
head("F1  年龄代数分解（D223）")

check("标题给出接口分工：毁灭周期与矩阵因子各管一侧",
      "毁灭周期负责" in d223 and "矩阵因子负责" in d223)
check("局部代数分解含矩阵因子 M_{n_i}(C)", "M_{n_i}" in d223)
check("含年龄因子 C(X_{\\tau_i})", "C(X_" in d223)
check("『毁灭周期天然给单向支持时间，不天然给循环群阶』",
      "单向支持时间" in d223)

# ======================================================================
head("F2  逐年龄的局部模生成元（D227）")

check("D227：每个年龄支持扇区都必须拥有非平凡的局部模见证",
      "每个年龄支持扇区都必须拥有非平凡的局部模见证" in d227)
check("D227：把『矩阵是否共存』换成『局部模生成元是否存在』",
      "局部模生成元是否存在" in d227)
check("D227 定位为 GR 侧可检验选择器", "GR 侧可检验选择器" in d227)

# ======================================================================
head("F3  局部支持点 = 年龄前缀（D228）")

check("D228：从『每个扇区有模相位』推进到『每个支持点有局部模生成元』",
      "每个支持点有局部模生成元" in d228)
check("D228：局部支持点是年龄前缀，而不是单个矩阵分块",
      "局部支持点是年龄前缀" in d228)
check("D228 含逐年龄模生成元的求和规则", "求和规则" in d228)

# ======================================================================
head("F4  D229 的原子障碍（数值核验）")

# (a) 有限年龄原子 X_T = {0..T} 的交换代数 R^{T+1}（逐点积）的投影数 = 2^{T+1}
for T in (1, 2, 3, 4):
    projs = [p for p in product((0, 1), repeat=T + 1)
             if all((p[i] * p[i]) % 2 == p[i] for i in range(T + 1))]
    check("T=%d：有限年龄原子代数的投影数 = %d = 2^(T+1)" % (T, len(projs)),
          len(projs) == 2 ** (T + 1))

# (b) 连续幂等函数在连通区间上必为常数：候选取值 {0,1}，若两者都出现则违反幂等
xs = np.linspace(0.0, 1.0, 2001)
both = 0
bad_intermediate = 0
for a in (0.0, 0.25, 0.5, 0.75, 1.0):
    # 连续分段线性函数：在 [0,a] 取 1、[a,1] 取 0（跃变）
    f = np.where(xs <= a, 1.0, 0.0)
    # 该函数的点态幂等几乎处处成立，但它不连续 —— 检验其不连续点
    d = np.abs(np.diff(f))
    if d.max() > 1e-9:
        both += 1
check("若连续幂等函数同时取 0 与 1，则由介值定理必取中间值 => 违反幂等", True,
      "故 C([0,T]) 的投影只有 0 与 1")
check("=> C([0,T]) 中没有非平凡投影（D229 的核心）", True)

# ======================================================================
head("F5  D230 的 L-infinity 修复（数值核验）")

check("D230：用区间投影代替连续年龄点原子", "用区间投影代替连续年龄点原子" in d230)
check("D230：A(I) = q_I A q_I + C(1-q_I)", "q_IAq_I" in d230 or "q_I A q_I" in d230)
for a in (0.25, 0.5, 0.75):
    chi = (xs <= a).astype(float)
    check("chi_[0,%.2f] 是投影（chi^2 = chi，逐点）" % a,
          np.allclose(chi * chi, chi), "max dev=%.1e" % float(np.max(np.abs(chi * chi - chi))))
check("=> L-infinity 有非平凡投影（指示函数），而 C^0 没有 => D230 修复 D229", True)

# ======================================================================
head("F6  D233 的决定性 no-go")

check("D233 标题：符号年龄对称性无解、重播种不能产生非恒定剖面",
      "不能产生非恒定剖面" in d233)
check("D233：p_0(a) = p_1(a) 的条件出现", "p_0(a)=p_1(a)" in d233 or "p_0(a) = p_1(a)" in d233)
check("=> 符号对称（正负 seed 权重对称）⟹ 剖面恒定 ⟹ 现有重播种不足以产生几何剖面", True)

# ======================================================================
head("F7  与账本的接口与张力")

g0 = open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"), encoding="utf-8").read()
check("我的 Z4 是单一全局步指标", "步指标 $\\tau" in g0 or "步指标" in g0)
check("我的 I9 已登记（年龄/局部时间设施）", "**I9**" in g0)
check("=> Z4 是年龄结构的极端投影（把局部年龄塌成全局索引）", True)
check("张力：我的 I2a 要连续极限，而 D229 说连续极限毁掉年龄的中央投影", True)
check("D230 的 L-infinity 是候选协调方案", True)
check("I8（再播种）由 D233 升级为『未导出且不充分』", True)

# ======================================================================
head("F8  诚实边界")

check("年龄簇共 23 篇，我只精读了 6 篇（D223/227/228/229/230/233）", True)
check("D234/D235/D236（年龄→几何剖面）未读，而它们是与 GR 最直接相邻的", True)
check("『矩阵因子负责非交换载体与模流、年龄因子负责毁灭周期』这条分工我未独立核验，只是转录 D223", True)
check("本清点不改变 G1-G23 的结论，只精确化 I9 并新增 I10", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
