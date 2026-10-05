#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R18_check.py —— 【L5-CERT 第一关：低维通道的目录与归一化判据】的独立核验
=========================================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：L5 目标、V_light、Z-STRESS 必败、L1 仍开放
  F2  §2  R12.2 反例的精确结构：(R18-1)/(R18-2) 与系数 91/8
  F3  §3  命题 R18.1 双方向（破坏方向 / 关闭方向）
  F4  §4  推论 R18.2：Z-STRESS 的约束数上限 m<=2
  F5  §5  推论 R18.3：三条可关路径，不是第四条
  F6  §6  判决与唯一下一步（1+1D 枚举）
  F7  数值核验：BKM 系数、振幅钉死、单约束残留危险方向、秩判据
  F8  跨文档状态一致：R17 / STATUS / R0 / INDEX / R12
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


DOC = read("R18_L5_cert_dimension_normalization_gate.md")
R0 = read("R0_publication_theorem.md")
R12 = read("R12_zero_native_gap_filling.md")
R17 = read("R17_L1_critical_path_and_L5_gate.md")
STATUS = read("STATUS.md")
IDX = read("INDEX.md")


# ======================================================================
head("F1  文档锚点")

check("标题含 L5-CERT 第一关", "L5-CERT 第一关" in DOC)
check("低维通道子空间 V_light 在位",
      "低维通道子空间" in DOC and "V_{\\rm light}" in DOC)
check("目标写成 D(sigma_R || rho) = o(R^d)",
      "D(\\sigma_R\\|\\rho)=o(R^d)" in DOC)
check("Z-STRESS 在 dim V_light >= 3 时必败", "Z-STRESS" in DOC and "必败" in DOC)
check("L1 仍标为未关闭", "**仍未关闭**" in DOC)


# ======================================================================
head("F2  R12.2 反例的精确结构")

check("引用命题 R12.2", "命题 R12.2" in DOC)
check("(R18-1) 展开式在位", "(R18-1)" in DOC and "\\chi_K(X)>0" in DOC)
check("(R18-2) 危险标度在位", "(R18-2)" in DOC and "\\Theta(R^{d})" in DOC)
check("写入 chi_K(x)=91/4 与系数 91/8", "91/4" in DOC and "91/8" in DOC)
check("强调危险来自振幅自由", "振幅的**自由**" in DOC)


# ======================================================================
head("F3  命题 R18.1 双方向")

check("命题 R18.1 在位", "命题 R18.1" in DOC)
check("(i) 破坏方向在位", "(i)（破坏方向）" in DOC)
check("(ii) 关闭方向在位", "(ii)（关闭方向）" in DOC)
check("核条件写成 ker N on V_light",
      "\\ker \\mathcal N|_{V_{\\rm light}}" in DOC)
check("关闭条件写成 N on V_light 单射",
      "\\mathcal N|_{V_{\\rm light}}" in DOC and "单射" in DOC)
check("含 BKM 正定性与边界声明",
      "BKM 二次型在零迹厄米扰动空间上正定" in DOC and "只给充分条件" in DOC)


# ======================================================================
head("F4  推论 R18.2：Z-STRESS 的约束数上限")

check("推论 R18.2 在位", "推论 R18.2" in DOC)
check("约束数写成 m <= 2", "m\\le2" in DOC)
check("三情形表在位：>=3 / =2 / <=1",
      "\\dim V_{\\rm light}\\ge3" in DOC
      and "\\dim V_{\\rm light}=2" in DOC
      and "\\dim V_{\\rm light}\\le1" in DOC)
check(">=3 时 Z-STRESS 不能关闭 L5", "**不能**关闭 L5" in DOC)
check("与 R12 1.3 第 2 条挂钩", "一致 BKM 界" in DOC)


# ======================================================================
head("F5  推论 R18.3：唯一可关路径")

check("推论 R18.3 在位", "推论 R18.3" in DOC)
check("三条路径齐备", all(x in DOC for x in ["DIM-CAT", "SYM-ZERO", "CONTACT"]))
check("明示不是第四条路", "不构成第四条路" in DOC)
check("不把四维判决写进结构判据",
      "它是结构判据，不是四维反例" in DOC)


# ======================================================================
head("F6  判决与唯一下一步")

check("下一步锁定 1+1D 枚举", "1+1D 自由费米子支线" in DOC)
check("枚举四项任务在位",
      all(x in DOC for x in ["规范中性、局域、标量的二次费米型通道",
                             "\\dim V_{\\rm light}",
                             "\\delta V,\\delta T_{00}"]))
check("触发 R17-STOP 的失败处置", "R17-STOP" in DOC and "放弃 Jacobson 2016" in DOC)
check("不把 Z-STRESS 常数当关闭证据",
      "不把 `Z-STRESS` 的常数" in DOC)


# ======================================================================
head("F7  数值核验")

w = np.array([1.0, 2.0, 4.0]) / 7.0


def chi(x):
    return float(np.sum(np.asarray(x) ** 2 / w))


def exact_D(eps, x):
    x = np.asarray(x, dtype=float)
    return float(np.sum((w + eps * x) * np.log(1.0 + eps * x / w)))


# ---- R12.2 的方向：系数 91/8 ----
x_r12 = np.array([1.0, -2.0, 1.0])
check("R12 方向零迹", abs(np.sum(x_r12)) < 1e-12)
check("chi_K(x_R12) = 91/4", abs(chi(x_r12) - 91.0 / 4.0) < 1e-10,
      "chi=%.6f" % chi(x_r12))

d = 4
ratios = []
for R in (0.05, 0.03, 0.02):
    eps = R ** (d / 2.0)          # eps_R = R^{d/2}
    ratio = exact_D(eps, x_r12) / (R ** d)
    ratios.append(ratio)
check("无归一化时 D/R^d -> 91/8 = 11.375",
      all(abs(r - 91.0 / 8.0) < 0.1 for r in ratios),
      "ratios=%s" % ["%.4f" % r for r in ratios])
check("危险标度不为零（L5 失败方向）",
      min(ratios) > 1.0)

# ---- 单条归一化：只去掉一个方向，仍留危险方向 ----
# 零迹空间（sum x = 0）的基
_, _, Vt = np.linalg.svd(np.array([[1.0, 1.0, 1.0]]))
basis = Vt[1:]                      # 2 维：sum(x) = 0
check("零迹空间维数为 2", basis.shape == (2, 3))

nA = np.array([1.0, 1.0, -2.0])     # 单条归一化约束 N_A
rowA = np.array([nA @ b for b in basis])
rankA = int(np.linalg.matrix_rank(rowA.reshape(1, -1)))
check("单条归一化在零迹空间上的秩为 1", rankA == 1)

# kern = span{v1}, v1 = (1,-1,0)
v1 = np.array([1.0, -1.0, 0.0])
check("残留方向满足 N_A(v1)=0", abs(nA @ v1) < 1e-12)
check("残留方向零迹", abs(np.sum(v1)) < 1e-12)
check("chi_K(v1) = 21/2", abs(chi(v1) - 21.0 / 2.0) < 1e-10,
      "chi=%.6f" % chi(v1))

ratios_res = []
for R in (0.05, 0.03, 0.02):
    eps = R ** (d / 2.0)
    ratios_res.append(exact_D(eps, v1) / (R ** d))
check("单约束后仍有 D/R^d -> 21/4 = 5.25 != 0",
      all(abs(r - 21.0 / 4.0) < 0.1 for r in ratios_res),
      "ratios=%s" % ["%.4f" % r for r in ratios_res])
check("单约束确实杀掉了 R12 方向",
      abs(nA @ x_r12) > 1e-9, "N_A(x_R12)=%.3f" % (nA @ x_r12))

# ---- 两条独立归一化：秩 2 -> 单射，无危险方向 ----
nB = np.array([1.0, -1.0, 0.0])
rows = np.array([[nA @ b for b in basis], [nB @ b for b in basis]])
check("两条独立归一化的秩为 2（单射）",
      int(np.linalg.matrix_rank(rows)) == 2)

# ---- 振幅钉死：N 固定后 eps=O(1)，D=O(1)=o(R^d) ----
eps_fixed = 0.01
D_fixed = exact_D(eps_fixed, v1)
check("振幅钉死后 D=O(1)", D_fixed < 1.0, "D=%.6e" % D_fixed)
ratios_fixed = [D_fixed / (R ** d) for R in (0.4, 0.6, 0.8)]
check("固定振幅时 D/R^d 单调趋零",
      all(a < 5e-2 for a in ratios_fixed)
      and ratios_fixed[0] > ratios_fixed[1] > ratios_fixed[2],
      "ratios=%s" % ["%.5f" % r for r in ratios_fixed])


# ======================================================================
head("F8  跨文档状态一致")

check("R17 仍以 L5-CERT 为主攻", "L5-CERT" in R17)
check("R12 仍保留 91/8 反例", "frac{91}{8}" in R12)
check("STATUS 已登记 R18", "R18" in STATUS and "2.18" in STATUS)
check("R0 外部审计范围已扩到 R18 及以后",
      any(("R8`–`R%d" % n) in R0 or ("R8-R%d" % n) in R0
          for n in range(18, 40)))
check("INDEX 已链接 R18 文档",
      "](R18_L5_cert_dimension_normalization_gate.md)" in IDX)
check("INDEX 已链接 R18 核验脚本", "](R18_check.py)" in IDX)
check("R18 未抬高主状态：条件恢复不变",
      "条件恢复" in DOC and "不关闭 L1" in DOC)


# ======================================================================
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
