#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R19_check.py —— 【L1 上游重排：概率与相位在位，缺的是洛伦兹 boost】的独立核验
============================================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：L1 目标、两条上游线、L1 仍开放
  F2  §1  概率线 / 相位线 / 合并点
  F3  §3  命题 R19.1、R19.2、推论 R19.3
  F4  §4  三个缺口面：Z-WICK / Z-WEDGE / Z-STRESS-2π
  F5  §5  验收条件 A1 / A2 / A3 与诚实边界
  F6  §6  下一步与止损
  F7  数值核验：紧旋转 vs 非紧 boost、Wick、Poincare 非阿贝尔
  F8  跨文档状态一致：R0 / STATUS / INDEX / R12 / R15 / R17
"""
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = "v" if ok else "x"
    if level == "sup" and ok:
        tag = "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    return io.open(os.path.join(HERE, name), encoding="utf-8").read()


DOC = read("R19_L1_upstream_probability_phase_and_missing_boost.md")
R0 = read("R0_publication_theorem.md")
R12 = read("R12_zero_native_gap_filling.md")
R15 = read("R15_zcar_double_cover_and_zstress_scale.md")
R17 = read("R17_L1_critical_path_and_L5_gate.md")
STATUS = read("STATUS.md")
IDX = read("INDEX.md")


# ======================================================================
head("F1  文档锚点")

check("标题含‘缺的是洛伦兹 boost’",
      "缺的是洛伦兹 boost" in DOC)
check("目标写成 K_B -> 2pi B_B",
      "K_{B,a}\\longrightarrow 2\\pi B_B" in DOC)
check("L1 仍标为未关闭", "**仍未关闭**" in DOC)
check("明示不关闭 L1、不新增物理参数",
      "不关闭 L1，不新增物理参数" in DOC)


# ======================================================================
head("F2  两条上游线与合并点")

check("概率线在位", "### 1.1 概率线" in DOC and "G29" in DOC)
check("相位线在位", "### 1.2 相位线" in DOC and "G62" in DOC)
check("合并点写成同一个 π",
      "### 1.3 合并点" in DOC and "同一条层结构" in DOC)
check("复相位与拓扑相位分开写",
      "复相位（量子叠加相位）" in DOC and "拓扑相位（spin holonomy）" in DOC)
check("拓扑 -1 holonomy 在位", "U_{2\\pi/L}" in DOC and "-\\mathbb I" in DOC)


# ======================================================================
head("F3  三条命题")

check("命题 R19.1 在位（上游机器已闭合到一个 π）",
      "命题 R19.1" in DOC and "只需一个 $\\pi$" in DOC)
check("命题 R19.2 在位（旋转双覆盖不提供 boost）",
      "命题 R19.2" in DOC and "非紧" in DOC)
check("旋转升格写成 U(2pi) = -1",
      "U(2\\pi)=-\\mathbb I" in DOC)
check("boost 写成 Lambda(eta) != -1",
      "\\Lambda(\\eta)\\ne-\\mathbb I" in DOC)
check("Wick 转动写成 Lambda(i theta) = U(theta)",
      "\\Lambda(i\\theta)=U(\\theta)" in DOC)
check("推论 R19.3 在位（R12.5 的精确读法）",
      "推论 R19.3" in DOC and "R12.5" in DOC)
check("引用 D44 与对易平移",
      "D44" in DOC and "[P_i,P_j]=0" in DOC)


# ======================================================================
head("F4  三个缺口面")

check("三个缺口面齐备",
      all(x in DOC for x in ["Z-WICK", "Z-WEDGE", "Z-STRESS-2π"]))
check("明示是重排而非新增标签",
      "不是新的独立缺项" in DOC and "既有缺口的分面" in DOC)
check("Z-CORE/Z-TAIL 原样保留",
      "`Z-CORE`／`Z-TAIL`" in DOC)


# ======================================================================
head("F5  验收条件 A1 / A2 / A3")

check("A1 几何一致性在位",
      "A1（几何一致性" in DOC and "KMS-$2\\pi$" in DOC)
check("A2 反射正性在位",
      "A2（反射正性" in DOC and "Osterwalder–Schrader" in DOC)
check("A3 非阿贝尔与区域保持在位",
      "A3（非阿贝尔与区域保持" in DOC and "[K_B,P]\\ne0" in DOC)
check("诚实边界标为未证",
      "A1、A2、A3 目前都**未证**" in DOC)
check("π²/3 与 2π 的差距写成判定点",
      "\\pi^{2}/3=3.2899" in DOC and "2\\pi=6.2832" in DOC
      and "6/\\pi" in DOC)


# ======================================================================
head("F6  下一步与止损")

check("第一步打 A1", "第一个可判定子问题" in DOC and "打 A1" in DOC)
check("失败处置指向 R17-STOP",
      "R17-STOP" in DOC and "换锚或改写平衡条件" in DOC)
check("不把 π²/3 差距跳过",
      "它是本轮定位出来的判定点" in DOC)


# ======================================================================
head("F7  数值核验")

I2 = np.eye(2, dtype=complex)


def U_rot(theta):
    """旋转升格 exp(i theta sigma_z / 2)。"""
    return np.diag([np.exp(1j * theta / 2.0), np.exp(-1j * theta / 2.0)])


def Lambda_boost(eta):
    """洛伦兹 boost exp(eta sigma_z / 2)，实 rapidity。"""
    return np.diag([np.exp(eta / 2.0), np.exp(-eta / 2.0)])


# ---- 紧旋转：2pi 给出 -1；每个 L 都有 (U_{2pi/L})^L = -1 ----
check("U(2pi) = -I", np.allclose(U_rot(2 * np.pi), -I2, atol=1e-12))
for L in (3, 4, 5, 8):
    U = U_rot(2 * np.pi / L)
    check("(U_{2pi/%d})^%d = -I" % (L, L),
          np.allclose(np.linalg.matrix_power(U, L), -I2, atol=1e-12))

# ---- 非紧 boost：实 rapidity 永远不给 -1，且非紧 ----
etas = np.linspace(-8.0, 8.0, 401)
eig_pos = all(np.min(np.linalg.eigvals(Lambda_boost(e)).real) > 1e-12 for e in etas)
check("实 rapidity 下 boost 本征值全为正", eig_pos)

mindist = min(np.min(np.linalg.svd(Lambda_boost(e) + I2, compute_uv=False)) for e in etas)
check("boost 与 -I 的距离下确界 >= 1（永不等于 -I）",
      mindist >= 1.0 - 1e-12, "mindist=%.6f" % mindist)

norms = [np.linalg.norm(Lambda_boost(e)) for e in etas]
check("boost 流无界（非紧）", max(norms) > 20.0,
      "max||Lambda||=%.3e" % max(norms))

# ---- Wick：Lambda(i theta) = U(theta)，-1 只在虚 rapidity ----
thetas = np.linspace(0.0, 2 * np.pi, 101)
wick_ok = all(np.allclose(Lambda_boost(1j * t), U_rot(t), atol=1e-12) for t in thetas)
check("Wick 转动 Lambda(i theta) = U(theta)", wick_ok)
check("-1 只出现在虚 rapidity 2pi i",
      np.allclose(Lambda_boost(2j * np.pi), -I2, atol=1e-12)
      and not np.allclose(Lambda_boost(2 * np.pi), -I2))

# ---- Poincare：非阿贝尔；半侧平移：交换代数 ----
# 1+1D Poincare，基 (H, P, K)：[K,H]=P, [K,P]=H, [H,P]=0
poincare_brackets = [
    np.array([0.0, 1.0, 0.0]),   # [K,H] = P
    np.array([1.0, 0.0, 0.0]),   # [K,P] = H
    np.array([0.0, 0.0, 0.0]),   # [H,P] = 0
]
abelian_brackets = [np.zeros(3) for _ in range(3)]

rank_poincare = int(np.linalg.matrix_rank(np.array(poincare_brackets)))
rank_abelian = int(np.linalg.matrix_rank(np.array(abelian_brackets)))
check("Poincare 导代数维数 = 2（非阿贝尔）", rank_poincare == 2)
check("半侧平移导代数维数 = 0（交换）", rank_abelian == 0)
check("dim[g,g] 是同构不变量 ⇒ 两者不同构",
      rank_poincare != rank_abelian)

# ---- Z-STRESS 落点：pi^2/3 vs 2pi ----
ratio = (2 * np.pi) / (np.pi ** 2 / 3.0)
check("2pi / (pi^2/3) = 6/pi ≈ 1.9099",
      abs(ratio - 6.0 / np.pi) < 1e-12 and abs(ratio - 1.9099) < 1e-3,
      "ratio=%.6f" % ratio)


# ======================================================================
head("F8  跨文档状态一致")

check("R0 范围已扩到 R19",
      any(("R8`–`R%d" % n) in R0 for n in range(19, 40))
      or ("R8-R19" in R0))
check("STATUS 已登记 R19", "R19" in STATUS and "2.19" in STATUS)
check("INDEX 已链接 R19 文档",
      "](R19_L1_upstream_probability_phase_and_missing_boost.md)" in IDX)
check("INDEX 已链接 R19 核验脚本", "](R19_check.py)" in IDX)
check("R12 仍保留 boost 缺口条目", "半侧模平移" in R12)
check("R15 仍保留 pi^2/3 读数", "pi^{2}}{3}" in R15 or "\\pi^{2}/3" in R15)
check("R17-STOP 仍在位", "R17-STOP" in R17)
check("R19 未抬高主状态：条件恢复不变",
      "条件恢复" in DOC or "不关闭 L1" in DOC)


# ======================================================================
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
