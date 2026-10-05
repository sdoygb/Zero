#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R15_check.py —— 【Z-CAR 双覆盖升级 + Z-STRESS 归一化常数】的独立核验
====================================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：Z-CAR′、双覆盖、π²/3、J1 未关闭、因子 L 隐患
  F2  R15.1 数值：升格 U(2π/L) 的 L 次幂 = -I（L=6,8,10,12）
  F3  区分：一次移位 U(2π/L) ≠ -I；置换号是另一个对象（L-循环 = -1）
  F4  Z-STRESS：h 对 β-hopping 的最优比例 A_N 稳定 ≈ π²/3，且 ≠ 2π
  F5  交叉引用：R13 探针存在；R14 指向 R15
"""
import cmath
import io
import os
import sys

import mpmath as mp

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


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


def read(name):
    return io.open(os.path.join(HERE, name), encoding="utf-8").read()


DOC = read("R15_zcar_double_cover_and_zstress_scale.md")
R14 = read("R14_L1_from_zero_assembly.md")


# ---------------------------------------------------------------- F1
head("F1  文档锚点")

check("命题 R15.1 在位并标为候选", "R15.1" in DOC and "候选" in DOC)
check("因子 L 隐患被写明", "因子隐患" in DOC or "因子 $L$" in DOC or "因子 L" in DOC)
check("双覆盖机制在位（2:1 覆盖）",
      "双覆盖" in DOC and "2:1" in DOC or "w\\longmapsto w^{2}" in DOC)
check("Z-STRESS 常数 π²/3 在位", "\\pi^{2}}{3}" in DOC or "π²/3" in DOC or "\\pi^{2}/3" in DOC)
check("J1 未关闭被写明", "L1" in DOC and "仍未关闭" in DOC)


# ---------------------------------------------------------------- F2
head("F2  R15.1 数值：升格的 L 次幂 = -I")


def spinor_lift(theta):
    """U(theta)=diag(e^{i theta/2}, e^{-i theta/2})。"""
    return cmath.exp(1j * theta / 2), cmath.exp(-1j * theta / 2)


lift_ok = True
for L in (6, 8, 10, 12):
    theta = 2 * cmath.pi / L
    a, b = spinor_lift(theta)
    aL, bL = a ** L, b ** L
    # (U_theta)^L = diag(e^{i pi}, e^{-i pi}) = -I
    if abs(aL + 1) > 1e-9 or abs(bL + 1) > 1e-9:
        lift_ok = False
check("(U(2π/L))^L = -I（L=6,8,10,12）", lift_ok,
      "每个 L 的两个对角元都 = e^{±iπ} = -1")


# ---------------------------------------------------------------- F3
head("F3  区分升格与置换号")

one_shift_ok = True
for L in (6, 8, 10, 12):
    theta = 2 * cmath.pi / L
    a, _ = spinor_lift(theta)
    if abs(a + 1) < 1e-6:      # 一次移位不应已经是 -1
        one_shift_ok = False
check("一次移位 U(2π/L) ≠ -I（升格不是一次就给 -1）", one_shift_ok)


def perm_sign_cycle(L):
    """L-循环 (1 2 ... L) 的符号 = (-1)^{L-1}。"""
    return (-1) ** (L - 1)


check("对照：L-循环置换号 = (-1)^{L-1} = -1（另一个对象）",
      all(perm_sign_cycle(L) == -1 for L in (6, 8, 10, 12)),
      "它对应‘一次移位’，与双覆盖的‘L 次移位’不是同一量")


# ---------------------------------------------------------------- F4
head("F4  Z-STRESS：h 对 β-hopping 的归一化常数")

mp.mp.dps = 50


def build_h(N):
    C = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            d = i - j
            if d == 0:
                C[i, j] = mp.mpf(1) / 2
            else:
                C[i, j] = mp.sin(mp.pi * d / 2) / (mp.pi * d)
    evals, evecs = mp.eigsy(C)
    mod = [mp.log((1 - evals[k]) / evals[k]) for k in range(N)]
    return evecs * mp.diag(mod) * evecs.T


def profile_scale(N):
    h = build_h(N)
    weights = [abs(h[i, i + 1]) for i in range(N - 1)]
    length = mp.mpf(N)
    beta = [
        (mp.mpf(i) + mp.mpf(1) / 2) * (length - mp.mpf(i) - mp.mpf(1) / 2) / length
        for i in range(N - 1)
    ]
    trim = 1
    w = weights[trim:len(weights) - trim]
    b = beta[trim:len(beta) - trim]
    num = mp.fsum(wi * bi for wi, bi in zip(w, b))
    den = mp.fsum(bi * bi for bi in b)
    return num / den


scales = []
for N in (8, 12, 16, 20):
    scales.append(profile_scale(N))
    print("    N=%2d: A_N = %s" % (N, mp.nstr(scales[-1], 10)))

pi2_3 = mp.pi ** 2 / 3
max_dev = max(abs(s - pi2_3) for s in scales)
check("A_N 稳定，与 π²/3 相合（偏差 < 0.5%）",
      max_dev < mp.mpf("0.005") * pi2_3,
      "max |A_N - π²/3| / (π²/3) = %s" % mp.nstr(max_dev / pi2_3, 4))
check("A_N 明确不等于 2π",
      all(abs(s - 2 * mp.pi) > mp.mpf("2") for s in scales),
      "min A_N = %s, 2π = %s" % (mp.nstr(min(scales), 8), mp.nstr(2 * mp.pi, 8)))


# ---------------------------------------------------------------- F5
head("F5  交叉引用")

check("R13 数值探针在位",
      os.path.exists(os.path.join(HERE, "R13_numeric_probe.py")))
check("R14 §6 已指向 R15（修正 R14.4）",
      "R15" in R14)


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
