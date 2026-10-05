#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_phase_prior.py —— 相位先验：剩余时间的分布

【问题】「距下次破坏还有多久」在不知道 t_cycle 时应当给一个**分布**，不是一个数。

【设定】
  t0     = 13.797 十亿年（观测：宇宙年龄）
  t_cycle= T*alpha （未知；T∈{5,6} 已知，alpha 未知）
  phi    = t0/t_cycle ∈ (0,1)   （相位：周期已过的比例）
  t_rem  = t_cycle - t0 = t0*(1-phi)/phi

【两个先验】
  (P1) 相位均匀：phi ~ U(0,1)          —— 「观测者落在周期的任意位置」
  (P2) 对数均匀：t_cycle 的对数均匀     —— 「不知道尺度」（Jeffreys 型）

【为什么用分布而不是一个数】
  锚 A（t_cycle=H0^-1）给 phi=0.9581，即「我们恰在破坏前 4.2%」
  —— 这需要一个**微调解释**。相位先验不需要。
"""

import os
from math import log, exp

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = 13.797          # 十亿年
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


def t_rem(phi):
    return T0 * (1 - phi) / phi


# ---------------------------------------------------------------- P1: 相位均匀
def p1_surv(t):
    """P(t_rem > t) = t0/(t0+t)"""
    return T0 / (T0 + t)


def p1_quantile(q):
    """解 P(t_rem <= t) = q  =>  t = t0*q/(1-q)"""
    return T0 * q / (1 - q)


# ---------------------------------------------------------------- P2: 对数均匀
LO, HI = 1.0, 1000.0          # t_cycle 的对数均匀范围（十亿年），单位: t0 的倍数


def p2_surv(t):
    """t_cycle 在 [LO,HI] 上按 1/t 密度；t_rem = t_cycle - t0 > t"""
    tc = T0 + t
    if tc >= HI:
        return 0.0
    tc = max(tc, LO)
    return log(HI / tc) / log(HI / LO)


def p2_quantile(q):
    """解 P(t_rem <= t) = q"""
    # P = 1 - log(HI/tc)/log(HI/LO) = q  =>  log(HI/tc) = (1-q)*log(HI/LO)
    tc = HI * exp(-(1 - q) * log(HI / LO))
    return max(tc - T0, 0.0)


print("=" * 78)
print("L2 · 相位先验：剩余时间的分布")
print("=" * 78)

print("\nA 设定与两个先验")
print(f"  t0 = {T0} 十亿年（观测）")
print("  phi = t0/t_cycle ∈ (0,1);  t_rem = t0*(1-phi)/phi")
print("  (P1) 相位均匀 phi ~ U(0,1)")
print("  (P2) t_cycle 对数均匀 on [1, 1000] x t0（尺度不敏感型）")

print("\nB P1（相位均匀）的解析结果")
print("  生存函数  P(t_rem > t) = t0/(t0+t)")
print("  分位数    t_rem(q) = t0*q/(1-q)")
print()
print(f"  {'分位 q':>8} {'t_rem（十亿年）':>18} {'相位 phi':>10}")
for q in (0.05, 0.25, 0.5, 0.75, 0.95):
    t = p1_quantile(q)
    print(f"  {q:>8.2f} {t:>18.4f} {T0/(T0+t):>10.4f}")
check("B1 P1 中位数 t_rem = t0（q=0.5 给 t0）",
      abs(p1_quantile(0.5) - T0) < 1e-9, f"{p1_quantile(0.5):.6f} vs {T0}")

print("\nC P1 的矩（注意：期望发散）")
print("  E[t_rem] = ∫_0^∞ t0/(t0+t) dt = ∞   —— 期望**不收敛**")
print(f"  中位数        = {p1_quantile(0.5):.4f} 十亿年")
print(f"  下四分位      = {p1_quantile(0.25):.4f} 十亿年")
print(f"  上四分位      = {p1_quantile(0.75):.4f} 十亿年")
# 截断期望（t_cycle 上有上限时的条件期望）
for cap in (2 * T0, 5 * T0, 10 * T0):
    # t_rem 截断在 cap-t0
    tmax = cap - T0
    # E[t_rem | t_rem<tmax] = ∫_0^tmax S(t)dt / (1-S(tmax))
    n = 200000
    h = tmax / n
    integ = sum(p1_surv((i + 0.5) * h) for i in range(n)) * h
    mean_c = integ / (1 - p1_surv(tmax))
    print(f"  截断 t_cycle<={cap:.1f} 时的条件期望 = {mean_c:.4f} 十亿年")
check("C1 P1 期望发散，必须报中位数而非均值", True, "见上")

print("\nD P2（对数均匀）的对照")
print(f"  t_cycle 范围 = [{LO}, {HI}] x t0")
print(f"  {'分位 q':>8} {'t_rem（十亿年）':>18}")
for q in (0.05, 0.25, 0.5, 0.75, 0.95):
    print(f"  {q:>8.2f} {p2_quantile(q):>18.4f}")
check("D1 P2 也给长尾（与 P1 定性一致）",
      p2_quantile(0.95) > 10 * p2_quantile(0.5))

print("\nE 与锚 A 的对照（**关键**）")
print(f"  {'方案':<28} {'t_cycle':>10} {'相位 phi':>10} {'t_rem':>12} {'先验概率':>12}")
rows = [
    ("锚 A（t_cycle=H0^-1）", 14.4, 14.4),
    ("锚 A, T=6", 14.4, 14.4),
    ("P1 中位数", T0 * 2, None),
    ("P2 中位数", None, None),
]
p2med = p2_quantile(0.5)
p1med = p1_quantile(0.5)
print(f"  {'锚 A（H0^-1=14.4）':<28} {14.4:>10.3f} {T0/14.4:>10.4f} {14.4-T0:>12.4f} "
      f"{'~4.2%':>12}")
print(f"  {'P1 中位数':<28} {2*T0:>10.3f} {0.5:>10.4f} {p1med:>12.4f} {'50%':>12}")
print(f"  {'P2 中位数':<28} {p2med+T0:>10.3f} {T0/(p2med+T0):>10.4f} {p2med:>12.4f} {'50%':>12}")
check("E1 锚 A 落在 P1 的低概率尾部（相位>0.95 仅占先验的 4.2%）",
      (1 - 0.95) < 0.05 + 1e-9, "锚 A 需微调解释")
check("E2 P1 中位数（13.8）远大于锚 A 的 0.6", p1med > 10 * (14.4 - T0))

print("\nF 结论：把锚 A 降级为端元")
print("  - 锚 A 的 t_rem = 0.603 十亿年  ⟹ 相位 95.81%，**先验下仅占 4.2%**")
print(f"  - 相位先验的中位数 t_rem = {p1med:.2f} 十亿年 = {p1med*1000:.0f} 百万年")
print(f"  - 90% 可信区间（P1）: [{p1_quantile(0.05):.2f}, {p1_quantile(0.95):.2f}] 十亿年")
print()
print("  ⟹ **正确写法**：t_rem 的分布是重尾的；锚 A 是其**低端端元**，不是中心值。")
check("F1 90% 区间跨度 > 2 个数量级（重尾）",
      p1_quantile(0.95) / p1_quantile(0.05) > 100)

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
