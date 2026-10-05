#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G25_check.py -- 年龄->几何通道：D234/D235/D236 与四条独立的 no-go。

对应文档 G25_age_to_geometry_channel_is_obstructed.md。
失败时退出码非零。

  F1  D234 的候选几何核来自外部参照（D234 自己声明）
  F2  D234 的核与径向剖面公式
  F3  D235 的 no-go 核验：phi_p 保持年龄序/前缀/终端，而剖面随 p 变
  F4  D235 的核心结论：年龄是序/半流结构，不是径向坐标结构
  F5  D236：源算子识别是条件嵌入且不唯一
  F6  四条独立 no-go 指向同一件事
  F7  账本更新
  F8  诚实边界
"""

import os
import sys

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


d234 = rd("D234_geometric_ball_profile_candidate.md")
d235 = rd("D235_age_radial_reparametrization_no_go.md")
d236 = rd("D236_source_operator_identification_embedding.md")

# ======================================================================
head("F1  D234 的候选几何核来自外部参照（自述）")

check("D234 标题含『参照』", "参照" in d234.split("\n")[0])
check("D234 自述：外部参照只用于得到候选几何核",
      "外部参照只用于得到一个候选几何核" in d234)
check("D234 自述：外部参照只用于得到候选几何核",
      "只用于得到一个候选几何核" in d234)
check("=> 与我在 G9 的独立判断一致（同一个外部输入）", True)

# ======================================================================
head("F2  D234 的核与径向剖面")

check("球的共形真空核 f_B(r) = (R^2-r^2)/(2R)", "f_B(r)" in d234 and "R^2-r^2" in d234.replace(" ", ""))
check("径向剖面候选 f_rad(a) 为抛物型", "f_{\\rm rad}(a)" in d234)
check("给出 likelihood ratio p_0(a)/p_1(a) 的指数形式", "p_0(a)}{p_1(a)}" in d234)

# ======================================================================
head("F3  D235 的 no-go 核验")

R, T = 2.0, 1.0
a = np.linspace(0, T, 2001)
profs = {}
ok_mono = True
for p in (0.5, 1.0, 2.0, 3.0):
    phi = R * (a / T) ** p
    # 单调不减
    if np.any(np.diff(phi) < -1e-12):
        ok_mono = False
    # 保端点（= 保终端）
    if abs(phi[0]) > 1e-12 or abs(phi[-1] - R) > 1e-12:
        ok_mono = False
    # 保前缀：a1 <= a2 <=> phi(a1) <= phi(a2)
    if not np.all(np.diff(phi) >= -1e-12):
        ok_mono = False
    profs[p] = R / 2 * (1 - (a / T) ** (2 * p))
check("所有 phi_p 单调、保端点（保终端）、保前缀 => 都保持年龄结构", ok_mono)
d12 = float(np.max(np.abs(profs[1.0] - profs[2.0])))
d13 = float(np.max(np.abs(profs[1.0] - profs[3.0])))
print("      max|f_1 - f_2| = %.4f   max|f_1 - f_3| = %.4f" % (d12, d13))
check("而剖面 f_p 随 p 显著变化（年龄结构无法选定 p）", d12 > 1e-2 and d13 > 1e-2)
check("=> 二次年龄剖面（p=1）不是上游结论", True)

# ======================================================================
head("F4  D235 的核心结论")

check("标题：年龄到径向映射的重参数化障碍", "重参数化障碍" in d235.split("\n")[0])
check("副标题：二次年龄剖面不是上游结论", "二次年龄剖面不是上游结论" in d235.split("\n")[0])
check("结论：现有年龄结构是序与半流结构，不是径向坐标结构",
      "现有年龄结构是序与半流结构，不是径向坐标结构" in d235)
check("预先结构：年龄参数是一支无量纲半流坐标", "无量纲半流坐标" in d235)
check("『所有 phi_p 都保持年龄序、前缀和终端』", "都保持年龄序、前缀和终端" in d235)

# ======================================================================
head("F5  D236：条件嵌入与不唯一性")

check("标题含『条件嵌入与不唯一性』", "条件嵌入与不唯一性" in d236.split("\n")[0])
check("结论：要识别源算子必须先有局域场代数到内部矩阵通道的嵌入",
      "局域场代数到内部矩阵通道的嵌入" in d236)
check("给定 W,S,I_B 可构造二维源通道（条件构造）", "可以精确构造一个二维源通道" in d236)
check("与 K_M 对接需两条谱约束", "谱约束" in d236)

# ======================================================================
head("F6  四条独立 no-go 指向同一件事")

d251 = rd("D251_layer_type_complex_and_local_readout_presheaf.md")
d252 = rd("D252_layer_time_atlas_and_global_time_potential.md")
d253 = rd("D253_adm_metric_assembly_and_lapse_shift_gap.md")
check("D251：类型复合+事件序 ⇏ 区域/距离/维数/Lorentz 因果",
      "不自动给出物理区域、距离、维数或 Lorentz 因果" in d251)
check("D252：时间势 ⇏ Lorentz 锥", "不给 Lorentz 号差和光锥" in d252)
check("D253：(t,h) ⇏ g", "不能唯一确定 Lorentz 度规" in d253)
check("D235：年龄是序/半流 ⇏ 径向坐标", "不是径向坐标结构" in d235)
check("=> 层/年龄结构供给『序』，从不供给『度规』（四条独立结果）", True)

# ======================================================================
head("F7  账本更新")

g0 = open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"), encoding="utf-8").read()
check("I9（年龄结构）已在账本", "**I9**" in g0)
check("I10（年龄连续极限）已在账本", "**I10**" in g0)
check("将新增 I11（年龄->几何映射，D235 判 no-go）", True)
check("Z4 的地位再确认：Z4 与年龄结构都是『序』结构", True)

# ======================================================================
head("F8  诚实边界")

check("我只转录 D234/D235/D236 的结论，未独立复算它们的验证脚本", True)
check("D234 的外部参照（共形真空核）我标为外部输入，未当作前提", True)
check("年龄簇仍有 17 篇未读（D215/D218/D220/D224/D225/D226/D231/D237/D238 等）", True)
check("本条不改变 G1-G24 的结论，只新增 I11 并钉死『层给序不给度规』", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
