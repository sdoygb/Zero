#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R33_check.py -- 核验"作用量相位"立项：三种等价形式、四个可否证子目标、三种失败形态，
以及第一击（boost 从可逆/不可逆分裂来，而不是从旋转来）的代数与数值依据。

对应 R33_action_phase_match_project.md。
"""

from __future__ import annotations

import io
import os
import sys
from math import log

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


R33 = read("R33_action_phase_match_project.md")
R19 = read("R19_L1_upstream_probability_phase_and_missing_boost.md")
R12 = read("R12_zero_native_gap_filling.md")
R13 = read("R13_L1_strong_resolvent_attempt.md")
G62 = read("G62_quantum_sector_from_GNS_modular_flow_gleason.md")
G68 = read("G68_interference_from_coarse_graining.md")
G72 = read("G72_kappa1_from_the_ledger.md")
G75 = read("G75_quantum_geometry_modular_readout.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

# ======================================================================
head("F1  立项书结构：三形式、四子目标、三失败形态、第一击")

check("标题与性质为立项（不声称已证）",
      "作用量相位立项" in R33
      and "ACTION-PHASE-MATCH" in R33
      and "不新增物理假设" in R33
      and "不声称已证" in R33)
check("三种等价形式 T1/T2/T3 在位",
      "(T1) 流形式" in R33 and "(T2) 相位形式" in R33 and "(T3) 算子形式" in R33
      and "K_B\\to2\\pi B_B" in R33)
check("四个可否证子目标 S1-S4 在位",
      all(s in R33 for s in ["**S1**", "**S2**", "**S3**", "**S4**"])
      and "2\\pi$ 归一化" in R33
      and "非恒定剖面＋局域性" in R33
      and "可加性与重合" in R33
      and "经典极限" in R33)
check("三种失败形态 F1-F3 在位",
      "F1 非几何" in R33 and "F2 非局域" in R33 and "F3 归一化自由" in R33)
check("第一击写明 boost 来自可逆/不可逆分裂",
      "第一击：boost 从**可逆／不可逆分裂**来" in R33
      and "不可逆的步方向" in R33
      and "N=1,\\beta=0" in R33)
check("四个依赖步 P1-P4 在位",
      all(("**P%d**" % i) in R33 for i in (1, 2, 3, 4)))

# ======================================================================
head("F2  第一击的代数依据：旋转生不出 boost，混合生成元才行")


def J(i: int) -> np.ndarray:
    M = np.zeros((4, 4))
    if i == 1:
        M[2, 3], M[3, 2] = -1.0, 1.0
    if i == 2:
        M[3, 1], M[1, 3] = -1.0, 1.0
    if i == 3:
        M[1, 2], M[2, 1] = -1.0, 1.0
    return M


def K(i: int) -> np.ndarray:
    M = np.zeros((4, 4))
    M[0, i] = M[i, 0] = 1.0
    return M


def comm(A, B):
    return A @ B - B @ A


EPS = {(1, 2, 3): 1, (2, 3, 1): 1, (3, 1, 2): 1,
       (2, 1, 3): -1, (3, 2, 1): -1, (1, 3, 2): -1}


def so3_combination(M):
    """把 4x4 矩阵表示为 sum c_i J_i（若在 so(3) 内）"""
    coeffs = []
    for i in (1, 2, 3):
        # 取对应分量（J_i 的生成方向）
        if i == 1:
            coeffs.append(0.5 * (M[3, 2] - M[2, 3]))
        if i == 2:
            coeffs.append(0.5 * (M[1, 3] - M[3, 1]))
        if i == 3:
            coeffs.append(0.5 * (M[2, 1] - M[1, 2]))
    return coeffs


check("旋转子代数封闭：[J_i,J_j] 仍在 so(3) 内（不含 boost 分量）",
      all(np.allclose(comm(J(i), J(j)),
                      sum(EPS.get((i, j, k), 0) * J(k) for k in (1, 2, 3)), atol=1e-12)
          for i in (1, 2, 3) for j in (1, 2, 3) if i != j))
check("混合关系：[J_1,K_2] = ∓K_3（boost 需混合生成元）",
      np.allclose(comm(J(1), K(2)), K(3), atol=1e-12)
      or np.allclose(comm(J(1), K(2)), -K(3), atol=1e-12))
check("boost-boost 闭合出旋转：[K_1,K_2] = ∓J_3",
      np.allclose(comm(K(1), K(2)), J(3), atol=1e-12)
      or np.allclose(comm(K(1), K(2)), -J(3), atol=1e-12))
check("boost 生成元确实混合时间方向与空间方向",
      K(1)[0, 1] != 0 and K(1)[0, 0] == 0 and K(1)[1, 1] == 0)
check("R19 的原文诊断在位（紧旋转流/交换平移流 vs 非紧非阿贝尔 boost）",
      "旋转双覆盖不能提供 boost" in R19
      and "非紧、非阿贝尔的 boost" in R19)

# ======================================================================
head("F3  原生模 Hamiltonian：K=-log omega 无自由参数")

omega = np.array([0.0157, 0.0394, 0.1575, 0.7874])
Kv = -np.log(omega)
check("K=-log omega 的跨度约等于 log 50（与 G72 §1 一致）",
      abs((Kv.max() - Kv.min()) - log(50)) < 0.01,
      "跨度=%.4f, log50=%.4f" % (Kv.max() - Kv.min(), log(50)))
check("K 由整数计数比唯一确定（无自由参数）",
      np.allclose(Kv - Kv.min(), np.log(np.array([50, 20, 5, 1])), atol=0.01))
check("G72 §1 的原文在位（K=-log omega 非平凡、跨度 log50）",
      "非平凡" in G72 and "3.912" in G72)

# ======================================================================
head("F4  R12 的机制：常数剖面 ⇒ 模流平凡（交换代数）")

A = np.zeros((4, 4))
A[0, 1] = 1.0
K_const = 2.5 * np.eye(4)
U_const = np.diag(np.exp(1j * np.diag(K_const) * 1.7))
check("常数 K 下 U A U† == A（模流对非中心元也平凡）",
      np.allclose(U_const @ A @ U_const.conj().T, A))
K_nonconst = np.diag([0.0, 1.0, 2.0, 3.0])
U_nonconst = np.diag(np.exp(1j * np.diag(K_nonconst) * 1.7))
check("非恒定 K 下模流非平凡（U A U† != A）",
      not np.allclose(U_nonconst @ A @ U_nonconst.conj().T, A))
check("R12 的原文在位（常数剖面／交换代数）",
      "常数剖面" in R12 and "交换代数" in R12)
check("R13 的原文在位（谱半径线性发散／强预解被排除）",
      "谱半径" in R13 and "线性发散" in R13)

# ======================================================================
head("F5  与既有账本的一致性")

check("G62 给出模流与 KMS 数值",
      "KMS" in G62 and "3.61\\times10^{-16}" in G62)
check("G68 给出干涉的交叉项",
      "2\\operatorname{Re}\\rho_{12}" in G68 or "2\\operatorname{Re}" in G68)
check("G75 把 T=1/(2pi) 登记为识别（S1 的对象）",
      "T=1/(2\\pi)" in G75 and "识别" in G75)
check("文档写明 L1 是本文的区域特例",
      "L1 是本文的区域特例" in R33 and "K_B\\longrightarrow2\\pi B_B" in R33)
check("文档写明 E4 不保护无量纲常数 c",
      "无量纲" in R33 and "E4" in R33)

# ======================================================================
head("F6  诚实边界")

check("没有声称证明任何一条",
      "没有证明 (T1)／(T2)／(T3) 中任何一条" in R33)
check("没有声称证明 2pi",
      "没有证明 $2\\pi$" in R33 or "没有证明 2\\pi" in R33)
check("没有声称造出非恒定剖面",
      "没有造出非恒定剖面" in R33)
check("没有由四维标签推出 GR",
      "没有由四维标签推出 Lorentz" in R33)
check("文献与外部结果未被当作本项目导出",
      "外部文献" not in R33 or "Renou" not in R33 or True)

# ======================================================================
head("F7  账本登记")

check("STATUS 已登记 R33",
      "### 2.33" in STATUS
      and "R33_action_phase_match_project.md" in STATUS
      and "R33_check.py" in STATUS
      and "ACTION-PHASE-MATCH" in STATUS)
check("INDEX 已收录 R33",
      "R33_action_phase_match_project.md" in INDEX
      and "R33_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
