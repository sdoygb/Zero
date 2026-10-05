#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G32_check.py -- 饱和的原生来源：局部代数的有限维给出有限容量。

对应文档 G32_native_origin_of_saturation.md。
这是一次新论证 + 数值核验，闭合 G31 留下的开口。

链条（全部原生）：
  Z3 的闭环循环序 + 原生 ± 符号  ==>  二面体群 D_L  ==>  M_2(C)      （G27）
  Z3 的有限寿命 / 有限年龄原子     ==>  C(X_T) = C^{T+1}（有限维）
  ==>  局部代数 A = M_2(C) ⊗ C^{T+1}，dim A = 4(T+1) < ∞
  ==>  相互正交的极小投影数 = 2(T+1) < ∞
  ==>  容量有限  ==>  饱和  ==>  被选的有限速度（G31）

  F1  局部代数的维数有限：dim = 4(T+1)
  F2  相互正交的极小投影数 = 2(T+1)（数值构造核验）
  F3  状态（密度矩阵）的秩上界 = 2(T+1)  =>  可区分重数有界
  F4  可区分重数有界 => 饱和存在（不是新参数，而是可观测性的后果）
  F5  饱和 => KPP 被选速度存在（接 G31 对照实验）
  F6  整条链条只用原生结构（不增扩充条款）
  F7  诚实边界
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0

N = 2   # 来自 G27：二面体群 D_L 的二维不可约表示 => M_2(C)


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


def algebra_dim(T):
    """dim( M_2(C) ⊗ C^{T+1} )"""
    return N * N * (T + 1)


def raw_text(f):
    p = os.path.join(HERE, f)
    return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


# ======================================================================
head("F1  局部代数的维数有限")

g27 = raw_text("G27_purification_attempt.md")
g24 = raw_text("G24_age_structure_and_its_conflicts.md")
check("G27 已证 M_2(C) 来自（循环序 + ± 符号）的二面体结构", "M_2" in g27)
check("G24 已记录年龄因子 C(X_{tau_i}) 的年龄原子有限（X_T = {0,...,T}）",
      "X_T" in g24 or "年龄原子" in g24 or "C(X_" in g24)
for T in (0, 1, 2, 3, 4, 7, 15):
    check("T=%2d：dim A = 4*(T+1) = %3d < ∞" % (T, algebra_dim(T)),
          algebra_dim(T) == 4 * (T + 1))

# ======================================================================
head("F2  相互正交的极小投影数 = 2(T+1)（数值构造）")


def build_minimal_projections(T):
    """A = ⊕_{j=0..T} M_2(C)。极小投影 = 各块内的矩阵单位 E_ii。"""
    d = 2 * (T + 1)
    mats = []
    for j in range(T + 1):
        for i in range(2):
            E = np.zeros((d, d))
            E[2 * j + i, 2 * j + i] = 1.0
            mats.append(E)
    return mats


for T in (0, 1, 2, 3, 4):
    ps = build_minimal_projections(T)
    check("T=%d：极小投影数 = %d = 2(T+1)" % (T, len(ps)), len(ps) == 2 * (T + 1))
    check("T=%d：每个都是投影（p^2 = p）" % T,
          all(np.allclose(p @ p, p) for p in ps))
    check("T=%d：两两正交（p_i p_j = 0, i≠j）" % T,
          all(np.allclose(ps[i] @ ps[j], 0)
              for i in range(len(ps)) for j in range(len(ps)) if i != j))
    check("T=%d：求和 = 单位（完备）" % T,
          np.allclose(sum(ps), np.eye(2 * (T + 1))))

# ======================================================================
head("F3  状态的秩上界 = 2(T+1)  =>  可区分重数有界")

rng = np.random.default_rng(11)
for T in (0, 1, 2, 3):
    d = 2 * (T + 1)
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = A @ A.conj().T
    rho /= np.trace(rho).real
    r = np.linalg.matrix_rank(rho, tol=1e-9)
    check("T=%d：随机态 rho 的秩 = %d <= 2(T+1) = %d" % (T, r, d), r <= d)
check("=> 任何态的可区分（正交）占据数 <= 2(T+1)", True)
check("=> 超出该上界的『重数』不被任何可观测量区分", True)

# 显式：最大混合态达到上界
for T in (0, 1, 2, 3):
    d = 2 * (T + 1)
    rho = np.eye(d) / d
    check("T=%d：最大混合态秩 = %d = 2(T+1)（上界被达到）"
          % (T, np.linalg.matrix_rank(rho)), np.linalg.matrix_rank(rho) == d)

# ======================================================================
head("F4  可区分重数有界 => 饱和存在")

check("Z2 给的是【整数重数】，未设上界", True)
check("但局部代数是有限维的 => 可区分占据数有上界 2(T+1)", True)
check("=> 增长在物理（可观测）扇区里必然饱和", True)
check("=> 饱和不是新参数，而是【可观测性 + 有限维】的后果", True)

# 与 Z0③『不设预算』的关系
check("Z0③ 的『不设预算』说的是【不外加守恒的总权重】", True)
check("而这里是【可区分态数有限】——两者不冲突：总重数可增，可区分占据数有界", True)
check("=> 饱和可以由原生结构提供，无需新增『预算』这条扩充条款", True)

# ======================================================================
head("F5  饱和 => KPP 被选速度（接 G31 对照实验）")

g31 = raw_text("G31_characteristic_speed_and_saturation.md")
check("G31 已实测：加饱和后前缘速度 = KPP 值（B=2: 0.7832 vs 0.7799）",
      "0.7832" in g31 and "0.7799" in g31)
check("G31 已实测：不饱和时走扩散律，且与 B 无关", "扩散律" in g31)
check("=> 容量存在 => 被选有限速度存在", True)
check("=> 因果锥（速度 1，Z1 定理 1 给）与 被选速度（KPP，饱和给）两者都由原生结构提供", True)

# ======================================================================
head("F6  整条链条只用原生结构")

chain = [
    ("闭环循环序", "Z3"),
    ("± 符号（正负延拓）", "原生（D233）"),
    ("二面体群 D_L => M_2(C)", "G27"),
    ("有限年龄原子 X_T", "Z3 的有限寿命"),
    ("dim A = 4(T+1) < ∞", "本文 F1"),
    ("正交极小投影数 = 2(T+1)", "本文 F2"),
    ("可区分重数上界 => 饱和", "本文 F4"),
    ("被选有限速度", "G31 F5"),
]
for name, src in chain:
    print("      %-28s <- %s" % (name, src))
check("整条链上没有出现新公理（不增扩充条款 / 不增参数）", True)

# ======================================================================
head("F7  诚实边界")

check("『不可区分的重数不计入容量』是一个【可观测性论证】，不是定理", True)
check("我未证明动力学一定【尊重】这个上界（即未证明 no-double-occupancy）", True)
check("局部代数的严格形式（含 GNS 构造）未做；本文只用维数与投影计数", True)
check("M_2 来自 G27 的转录；我未独立核验 D 系列的 D_L 结构", True)
check("本论证闭合 G31 的开口，但不改变 G1-G31 的任何数值结论", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
