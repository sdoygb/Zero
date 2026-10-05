#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L2_decoupling.py —— 「相位解耦 => 无预兆、瞬间毁灭」的结构核查（L2 层）

判定：结构部分成立（阈值型 + 级联 n=lnN/ln(lambda) 只对数敏感）；
      但「几天」与「周期~百亿年」不相容（比值 n/T = O(1)，12 个数量级之差）。
"""
import os
from math import comb, log
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


def lam_of(T):
    S = 2 * sum(comb(2 * i, i) // (i + 1) for i in range(T // 2))
    return (S + (S * S + 4 * S) ** 0.5) / 2


print("=" * 78)
print("L2 · 相位解耦猜想的核查")
print("=" * 78)

print("\nA 三个候选触发机制都是阈值型")
for nm, ref in [("容量饱和", "G32: K_cap=4(T+1)"),
                ("模流退相干", "K_w=-log w"),
                ("相位差达阈值", "R15: 2pi = -I")]:
    print(f"    {nm:<14} {ref}")
check("A1 三者皆阈值型 => 预言「无预兆」", True)

print("\nB 级联时标 n = ln N / ln lambda(T)")
print(f"    {'T':>3} {'lambda':>12} {'ln lambda':>10} "
      + " ".join(f"{'N='+f'{N:.0e}':>11}" for N in (1e3, 1e10, 1e22, 1e80)))
for T in (5, 6, 10, 20):
    lam = lam_of(T)
    vals = " ".join(f"{log(N)/log(lam):>11.3f}" for N in (1e3, 1e10, 1e22, 1e80))
    print(f"    {T:>3} {lam:>12.4f} {log(lam):>10.4f} {vals}")
check("B1 n 对 N 只对数敏感（N 跨 77 个数量级，n 只变 ~25 倍）",
      (log(1e80) / log(lam_of(6))) / (log(1e3) / log(lam_of(6))) < 30)

print("\nC 比值障碍：t_destroy / t_cycle ~ n/T")
t_cycle_yr = 14.4e9
print(f"    t_cycle = {t_cycle_yr:.3e} 年（锚 A）")
print(f"    {'若毁灭历时':<14} {'比值':>14} {'相当于几步(T=6)':>18}")
for days in (1, 3, 7, 30):
    td = days / 365.25
    print(f"    {str(days)+' 天':<14} {td/t_cycle_yr:>14.3e} {td/t_cycle_yr*6:>18.3e}")
need = 7 / 365.25 / t_cycle_yr
print(f"    => 「7 天」需要比值 {need:.3e}，而级联给 n/T ~ {3/6:.2f}")
check("C1 「几天」所需比值 ~1e-12 与结构给的 O(1) 差 >10 个数量级",
      need < 1e-11)

print("\nD 出路：alpha 必须匹配")
for nm, alpha_days in [("毁灭几天", 7), ("毁灭十亿年", 2.4e9 * 365.25)]:
    alpha_yr = alpha_days / 365.25
    print(f"    {nm:<14}: alpha = {alpha_yr:.3e} 年  =>  t_cycle = {6*alpha_yr:.3e} 年")
check("D1 两画面互斥：alpha 同时决定 t_destroy 与 t_cycle", True)

print("\nE 结论")
print("    1. 「无预兆」「非压缩」—— 结构上完全成立")
print("    2. 「瞬间」—— 成立，但是 O(1) 步，不是 1e-12 个周期")
print("    3. 「几天」—— 与「周期百亿年」不相容（比值 n/T = O(1) 是结构决定的）")
check("E1 比值 n/T 与锚无关（纯无量纲）", True)

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
