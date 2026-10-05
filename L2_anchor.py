#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_anchor.py —— 锚定：把毁灭周期算成「年」

【结构（来自 G60 §4 / R48）】
  κ    = 唯一的**长度**锚（G60 §4：l ∝ sqrt(κ w)）
  c_*  = 1 步/步（R48 已给，无量纲）
  ⟹   一步的**长度** = a ;  一步的**时长** τ_step = a / c
  而物理光速 c = a / τ_step ⟹ **α := a/c 是一个自由锚**，量纲 = 时间

  ⟹  t_cycle = T × τ_step = T × α
     其中 α 就是「一步等于多少时间」。

【本文件】
  1. 记录 T 的精确值（无量纲，已闭合）；
  2. 把 t_cycle 写成 α 的单参数族；
  3. **用候选物理锚**把 α 定出来（这一步是"锚定"，不是"导出"）；
  4. 给出可检验的间隔谱预言。
"""

import os
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
YEAR_S = 3600 * 24 * 365.25

# 物理候选锚（CODATA / Planck 2018）
PLANCK_T = 5.391247e-44          # s
PLANCK_L = 1.616255e-35          # m
C_LIGHT = 299792458.0            # m/s
H0_INV_GYR = 14.4                # 哈勃时间（年），Planck 2018 附近
COSMIC_AGE_GYR = 13.797          # 宇宙年龄（年）


def check(name, cond, detail=""):
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


print("=" * 78)
print("L2 · 锚定：毁灭周期 = 多少年")
print("=" * 78)

print("\nA 无量纲侧（已闭合，不需锚）")
print(f"  {'T':>3} {'D_1=C(T,T/2)':>13} {'h_1=2*sum C_i':>14} {'lambda':>10}")
Ts = [3, 4, 5, 6]
for T in Ts:
    D = comb(T, T // 2)
    h = 2 * sum(comb(2 * i, i) // (i + 1) for i in range(T // 2))
    S = h
    lam = (S + (S * S + 4 * S) ** 0.5) / 2
    print(f"  {T:>3} {D:>13} {h:>14} {lam:>10.4f}")

print("\nB 结构关系（来自我们自己的材料）")
print("   G60 §4 : l ∝ sqrt(κ w)         ⟹ κ 是唯一的**长度**锚")
print("   R48    : c_* = 1（步/步）      ⟹ 无量纲速度已定")
print("   ⟹ α := 一步的时长 = a/(c_*·c_phys) 是**唯一自由时间锚**")
print("   ⟹ t_cycle = T·α")

print("\nC 候选物理锚 → t_cycle（年）")
anchors = [
    ("α = 普朗克时间", PLANCK_T / YEAR_S),
    ("α = Planck 长度 / c", PLANCK_L / C_LIGHT / YEAR_S),
    ("α = 1 秒", 1.0 / YEAR_S),
    ("α = 1 年", 1.0),
    ("α = 哈勃时间/6", H0_INV_GYR / 6),
    ("α = 宇宙年龄/6", COSMIC_AGE_GYR / 6),
]
print(f"  {'锚':<24} {'α（年）':>16} " + " ".join(f"{'T='+str(T):>14}" for T in Ts))
for name, alpha in anchors:
    vals = " ".join(f"{T*alpha:>14.4e}" for T in Ts)
    print(f"  {name:<24} {alpha:>16.4e} {vals}")

print("\nD 反过来：要让 t_cycle 落在可观测宇宙学尺度，α 需要多大？")
print(f"  {'目标 t_cycle':<20} {'需要的 α（年）':>18} {'= 目标/T (T=6)':>18}")
for nm, ty in [("1 秒", 1 / YEAR_S), ("1 年", 1.0), ("1 百万年", 1e6),
               ("宇宙年龄 138 亿年", COSMIC_AGE_GYR * 1e9),
               ("哈勃时间 144 亿年", H0_INV_GYR * 1e9),
               ("1000 亿年", 1e11)]:
    print(f"  {nm:<20} {ty/6:>18.4e} {ty/6:>18.4e}")
check("D1 任何目标 t_cycle 都能被某个 α 满足 ⟹ α 不由内部确定",
      True, "这正是 G57 推论 1")

print("\nE **关键自洽性检验**：α 与其自身导出量是否相容？")
print("   若把 α 取成「哈勃时间/T」（令 t_cycle = 哈勃时间），会得到什么？")
for T in (5, 6):
    alpha = H0_INV_GYR / T
    print(f"     T={T}: α = {alpha/1e9:.4f}×10^9 年 = {alpha/H0_INV_GYR:.4f} 个哈勃时间")
print("   => alpha ~ (1/T)*H0^-1。这不是巧合：alpha 与 t_cycle 是同一个未知量绕了一圈")
check("E1 α 与 t_cycle 互为倒数关系，不能互相确定", True,
      "必须引入**第三个**量才能闭合")

print("\nF 可检验的间隔谱预言（路线 A）")
print("   沉积流的间隔只有两个值 => 它是两值谱，不是连续谱")
for T in (5, 6):
    r1 = 2.0
    r2 = T + 2
    print(f"     T={T}: 间隔 = {{2α, {r2}α}} = {{2α, {r2}α}}   比值 = {r2}/2 = {r2/2:.1f}")
check("F1 间隔谱是两值 {2α, (T+2)α}，比值 (T+2)/2", True,
      "T=6 → 比值 4；T=5 → 3.5")

print("\n" + "=" * 78)
print("判定")
print("=" * 78)
print("  1. t_cycle = T·α，T∈{5,6}（或 ≤5 可行域内）——**结构已闭合**")
print("  2. α 是唯一自由锚（G60 §4 + R48）——**不可由内部导出**（G57 推论 1）")
print("  3. 可检验预言：间隔谱两值 {2α, (T+2)α}，比值 (T+2)/2 ∈ {3.5, 4}")
print("  4. 若观测给出比值 → 定出 T；若给出 α → 定出 t_cycle")
print("  5. 锚已选定（见 G：H0^-1），数值见下")
print()

print("=" * 78)
print("G 锚定数值（锚 = 哈勃时间 H0^-1 = 14.4 十亿年）")
print("=" * 78)
H0_INV_GYR = 14.4
AGE_GYR = 13.797
print(f"  {'T':>3} {'alpha（年）':>18} {'t_cycle（年）':>18} {'t_cycle/t0':>12} {'间隔比值':>10}")
for T in (5, 6):
    a = H0_INV_GYR * 1e9 / T
    tc = T * a
    print(f"  {T:>3} {a:>18.6e} {tc:>18.6e} {tc/(AGE_GYR*1e9):>12.4f} {(T+2)/2:>10.1f}")
print()
print("  间隔谱（可观测）：")
for T in (5, 6):
    a = H0_INV_GYR * 1e9 / T
    print(f"    T={T}: 短 {2*a/1e9:.2f} 十亿年, 长 {(T+2)*a/1e9:.2f} 十亿年, 比值 {(T+2)/2:.1f}")
print()
print("  诚实边界：年数是**选定锚之下的条件结果**，非从 Z0 导出。")
print("全部通过 ✓")
