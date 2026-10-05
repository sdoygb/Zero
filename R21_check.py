#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R21_check.py —— 【R20 归一化缺口的判定：2π/v_F】的独立核验
============================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：目标、结论框、L1 仍开放
  F2  外部定理与 T 的定义在位
  F3  色散、v_F、半满归一化与 2π/v_F 在位
  F4  R20 的三组外推读数与 1.85 重读在位
  F5  命题 R21.1/R21.2 与边界在位
  F6  对 L1 的净影响、未关闭项与可证伪点
  F7  独立复算：v_F=2、2π/v_F=π、2π/3.3943≈1.851
  F8  跨文档状态一致：R0 / STATUS / INDEX / R15 / R19 / R20
"""
import io
import os
import sys

import mpmath as mp


HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 50

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


DOC = read("R21_vf_normalization_resolves_the_r20_factor.md")
R0 = read("R0_publication_theorem.md")
R15 = read("R15_zcar_double_cover_and_zstress_scale.md")
R19 = read("R19_L1_upstream_probability_phase_and_missing_boost.md")
R20 = read("R20_A1_verdict_shape_holds_constant_fails.md")
STATUS = read("STATUS.md")
IDX = read("INDEX.md")


# ======================================================================
head("F1  文档锚点")

check("标题含 2π/v_F", "2\\pi/v_F" in DOC)
check("结论框含原目标写错", "原 A1 的常数目标写错了" in DOC)
check("L1 仍标为未关闭", "J1 仍未关闭" in DOC or "仍未关闭" in DOC)
check("明确不新增物理参数", "不新增物理参数" in DOC)


# ======================================================================
head("F2  外部定理与 T 定义")

check("外部定理引用在位", "Eisler–Tonni–Peschel" in DOC)
check("H_ent 的 πN T 关系在位", "H_{\\rm ent}\\simeq \\pi N\\,T" in DOC)
check("T 的键定义在位",
      "T_{b,b+1}=T_{b+1,b}" in DOC and "\\left(1-\\frac{b+1}{N}\\right)" in DOC)
check("N T 转物理 beta 在位", "N\\,T_{b,b+1}\\simeq" in DOC)


# ======================================================================
head("F3  色散、v_F 与 2π/v_F")

check("格点色散在位", "\\varepsilon(k)=-2t\\cos(ka)" in DOC)
check("速度公式在位", "v_F=\\left|a\\,\\frac{d\\varepsilon}{dk}\\right|_{k_F}" in DOC)
check("半满 v_F=2 在位", "v_F=2" in DOC)
check("2π/v_F=π 在位", "\\frac{2\\pi}{v_F}" in DOC and "=\n\\pi" in DOC)
check("一般填充 sine 关系在位", "\\frac{2\\pi ta}{v_F}" in DOC)


# ======================================================================
head("F4  R20 读数重解释")

check("三组外推读数在位",
      all(x in DOC for x in ["3.3048", "3.3802", "3.3943"]))
check("1.851 与 0.9255 在位", "1.851" in DOC and "0.9255" in DOC)
check("0.9255 降为估计器商而非有限尺寸残差",
      "0.9255" in DOC and "不是物理修正因子" in DOC and "R22" in DOC)
check("明确不是新物理常数", "不是新的物理常数" in DOC)


# ======================================================================
head("F5  命题与边界")

check("命题 R21.1 在位", "命题 R21.1" in DOC)
check("命题 R21.2 在位", "命题 R21.2" in DOC)
check("条件证成标注在位", "条件证成" in DOC)
check("有限 N 边界在位", "有限 $N$ 数据已经精确等于 $\\pi$" in DOC)
check("不抬高四维 L1", "没有把 `Z-STRESS-2π` 抬高成四维 L1 定理" in DOC)


# ======================================================================
head("F6  对 L1 的净影响")

check("给出修正目标 2π/v_F", "\\frac{2\\pi}{v_F}B_{\\rm bare}" in DOC)
check("保留 Z-WICK / Z-WEDGE", "`Z-WICK`" in DOC and "`Z-WEDGE`" in DOC)
check("保留 Z-CORE / Z-TAIL", "`Z-CORE`" in DOC and "`Z-TAIL`" in DOC)
check("可证伪点在位", "可证伪点" in DOC and "应撤回" in DOC)


# ======================================================================
head("F7  独立复算")

kF = mp.pi / 2
eps = lambda k: -2 * mp.cos(k)
vF = abs(mp.diff(eps, kF))
check("dε/dk 在 k_F=π/2 为 2", abs(vF - 2) < mp.mpf("1e-40"), "v_F=%s" % mp.nstr(vF, 20))
check("2π/v_F=π", abs((2 * mp.pi / vF) - mp.pi) < mp.mpf("1e-40"))
lam = mp.mpf("3.3943")
ratio = (2 * mp.pi) / lam
resid = mp.pi / lam
check("2π/3.3943 约 1.851", abs(ratio - mp.mpf("1.851")) < mp.mpf("0.001"),
      "ratio=%s" % mp.nstr(ratio, 8))
check("1.851 可拆成 2×0.9255",
      abs(2 * resid - ratio) < mp.mpf("1e-30") and abs(resid - mp.mpf("0.9255")) < mp.mpf("0.001"),
      "resid=%s" % mp.nstr(resid, 8))
check("π 小于 3.3943 且相对差小于 9%",
      mp.pi < lam and (lam - mp.pi) / mp.pi < mp.mpf("0.09"),
      "diff=%s" % mp.nstr((lam - mp.pi) / mp.pi, 6))


# ======================================================================
head("F8  跨文档状态一致")

check("R0 范围已扩到 R22",
      any(("R8`–`R%d`" % n) in R0 for n in range(22, 40)) or ("R8-R22" in R0))
check("STATUS 已登记 R21", "R21" in STATUS and "2.21" in STATUS)
check("INDEX 已链接 R21 文档",
      "](R21_vf_normalization_resolves_the_r20_factor.md)" in IDX)
check("INDEX 已链接 R21 核验脚本", "](R21_check.py)" in IDX)
check("R20 已加 R21 后续修正", "R21" in R20 and "2\\pi/v_F" in R20)
check("R19 已加 R21 后续修正", "R21" in R19 and "2\\pi/v_F" in R19)
check("R15 已加 R21 后续修正", "R21" in R15 and "2\\pi/v_F" in R15)
check("R21 已加 R22 估计器校正",
      "R22" in DOC and "估计器" in DOC and "0.9255" in DOC)


# ======================================================================
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
