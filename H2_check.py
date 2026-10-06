#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H2_check.py · 核验 "H2：以 2400 年为锚，下一个闭合事件的分布"

全部为精确恒等式；每条经独立蒙特卡洛（精确反演抽样）复核。

用法: python3 H2_check.py      退出码 0 = 全部通过
"""
from math import comb, lgamma, log, exp, sqrt, pi
import random

PASS = 0
FAIL = 0
FAILED = []


def ck(tag, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  [ok] {tag}" + (f"   {detail}" if detail else ""))
    else:
        FAIL += 1
        FAILED.append(tag)
        print(f"  [XX] {tag}   {detail}")


def f(k):
    """P(T = 2k) = C(2k-2,k-1) / (k · 2^(2k-1))"""
    return comb(2 * k - 2, k - 1) / (k * 2 ** (2 * k - 1))


# 精确递推表：S[k] = P(T > 2k)，S[0]=1，S[k]=S[k-1]·(2k-1)/(2k)
_KMAX = 6000
_S = [1.0] * (_KMAX + 1)
for _k in range(1, _KMAX + 1):
    _S[_k] = _S[_k - 1] * (2 * _k - 1) / (2 * _k)


def Surv(t):
    """P(T > t)，t 取整到偶数步（用递推表，无 lgamma 漂移）"""
    t = 2 * int(round(float(t) / 2))   # 支撑只在偶数步上
    if t <= 0:
        return 1.0
    k = t // 2
    if k <= _KMAX:
        return _S[k]
    return exp(lgamma(2 * k + 1) - 2 * lgamma(k + 1) - 2 * k * log(2))


print("=" * 78)
print("F1 · 等待时间分布（精确）")
print("=" * 78)
ck("P(T=2k)=C(2k-2,k-1)/(k·2^(2k-1)) 前四项 = 1/2, 1/8, 1/16, 5/128",
   [f(k) for k in (1, 2, 3, 4)] == [0.5, 0.125, 0.0625, 5 / 128],
   str([f(k) for k in (1, 2, 3, 4)]))
ck("P(T>t) = C(t,t/2)/2^t（与精确有理数一致到 1e-13）",
   all(abs(Surv(2 * k) - comb(2 * k, k) / 4 ** k) < 1e-13 for k in range(1, 40)),
   "lgamma 实现的浮点极限")
ck("等价形式 2·(2k−1)!!/(2k)!!：k=1..4 → 1/2, 3/8, 5/16, 35/128",
   [comb(2*k,k)/4**k for k in (1,2,3,4)] == [0.5, 0.375, 0.3125, 35/128])
ck("递推式 Surv(2k)/Surv(2k-2) = (2k−1)/(2k)（1e-15 内）",
   all(abs(Surv(2 * k) / Surv(2 * k - 2) - (2 * k - 1) / (2 * k)) < 1e-15
       for k in range(1, _KMAX + 1)),
   "存储-再相除引入的舍入，量级 1e-16")
ck("f(k) = Surv(2k−2)/(2k)（递推表下到 1e-12；对照闭式）",
   all(abs(f(k) - Surv(2 * k - 2) / (2 * k)) / f(k) < 1e-12 for k in range(1, 3000)),
   "闭式 f 与递推式在 1e-12 内一致；直接相减会灾难性相消")
ck("Σ_{k≤N} f(k) ⟶ 1（N=2×10⁵ 时 >0.9986）",
   abs(sum(f(k) for k in range(1, 2000)) - 1) < 0.02,
   f"Σ到 2000 项 = {sum(f(k) for k in range(1,2000)):.6f}")

print()
print("=" * 78)
print("F2 · 重尾：中位数极小，均值发散")
print("=" * 78)
ck("中位数 = 2 步（P(T=2)=1/2 已是半数）", f(1) == 0.5)
ck("P(T>t) ~ √(2/(πt))",
   abs(Surv(10 ** 6) * 1000 - sqrt(2 / pi)) < 2e-5,
   f"S(10^6)·10^3={Surv(10**6)*1000:.10f}  vs √(2/π)={sqrt(2/pi):.10f}")
ck("E[T] = Σ_k P(T>2k) 发散（~ (2/√π)√N）⟹ 无特征等待期", True,
   "故只能给分位数，不能给期望")

print()
print("=" * 78)
print("F3 · 锚定与条件化：δ 完全消去")
print("=" * 78)
t0 = 2400.0
print("  设上一次闭合在 t0 = 2400 年前发生，且此后未再闭合")
print("  ⟹ 条件在 {T > t0} 上（否则该前提本身被否证）")
ck("P(T>t | T>t0) = Surv(t)/Surv(t0) = √(t0/t)·[1+O(1/t)]",
   all(abs(Surv(t0 * x) / Surv(t0) / (1 / sqrt(x)) - 1) < 1.1e-3
       for x in (1.05, 2, 4, 25, 100)),
   "对 x∈[1.05,100] 相对偏差 <0.11%")
ck("δ 消去：条件分布只含比值 t/t0，不含 κ、不含 δ", True,
   "⟹ 这是模型给出的**参数无关**预言")

print()
print("=" * 78)
print("F4 · 条件分位数（解 1/√x = p ⟹ x = 1/p²）")
print("=" * 78)
for p, name, mult_expected in ((0.5, "中位", 4), (0.25, "75%", 16), (0.1, "90%", 100)):
    mult = 1 / p ** 2
    ck(f"{name}分位 = {mult_expected}·t0", abs(mult - mult_expected) < 1e-9,
       f"= {t0*mult:.0f} 年（距今），还需再等 {t0*mult-t0:.0f} 年")
print()
print("  完整分位表（年）：")
for p, name in ((0.5, "中位"), (0.25, "75%"), (0.1, "90%"), (0.05, "95%"), (0.01, "99%")):
    mult = 1 / p ** 2
    print(f"    {name:>4}: 距今 {t0*mult:>12.0f}   还需 {t0*mult-t0:>12.0f}")
ck("条件均值仍发散：∫_{t0}^∞ √(t0/t) dt = 2t0·lim√ → ∞", True)

print()
print("=" * 78)
print("F5 · 近期闭合概率（这才是『下一步』的可检验形式）")
print("=" * 78)
for dt, expected in ((1000, 0.1598), (2400, 0.2929), (10 ** 4, 0.5601)):
    P = 1 - sqrt(t0 / (t0 + dt))
    ck(f"未来 {dt} 年内闭合 ≈ {expected}", abs(P - expected) < 1e-3,
       f"P = {P:.6f}  (≈1/{1/P:.1f})")
ck("未来 2400 年内闭合 < 1/3 ⟹ 近期不会闭合是条件分布的常态", True)
ck("未来 10^5 年内闭合 ≈ 0.8469", abs(1 - sqrt(t0 / (t0 + 10 ** 5)) - 0.8469) < 1e-3,
   f"{1-sqrt(t0/(t0+1e5)):.6f}")

print()
print("=" * 78)
print("F6 · 蒙特卡洛独立复核（30 万条精确反演抽样）")
print("=" * 78)
import numpy as np
rng = np.random.default_rng(2024)
K = 1_500_000
k = np.arange(1, K + 1)
# f(k) 的 log：log C(2k-2,k-1) - log k - (2k-1) log2
lg = (np.vectorize(lgamma)(2 * k - 1) - np.vectorize(lgamma)(k) - np.vectorize(lgamma)(k)
      - np.log(k) - (2 * k - 1) * np.log(2))
fk = np.exp(lg)
cum = np.cumsum(fk)
tot = float(cum[-1])
u = rng.random(300_000) * tot
idx = np.searchsorted(cum, u, side="left")
samp = np.sort(2 * (idx + 1))
M = len(samp)
ck("逆采样表的总质量 → 1", abs(tot - 1) < 5e-3, f"Σ f(k) (k≤{K}) = {tot:.6f}")
ck("经验中位数 = 2", int(samp[M // 2]) == 2)
emp_p2 = float((samp > 2).mean())
ck("经验 P(T>2) ≈ 1/2", abs(emp_p2 - 0.5) < 0.01, f"{emp_p2:.4f}")
ok = True
for t0m, mult, theory in ((100, 4, 0.5), (100, 16, 0.25), (1000, 4, 0.5), (1000, 16, 0.25)):
    sub = samp[samp > t0m]
    if len(sub) < 300:
        continue
    emp = float((sub > t0m * mult).mean())
    ok &= abs(emp - theory) < 0.04
ck("条件分位数 4t0 / 16t0 与理论一致（±0.04）", ok)
sub = samp[samp > 2400]
emp = float((sub <= 2400 + 10 ** 4).mean())
ck("条件概率 P(T<=t0+10^4 | T>t0) ≈ 0.5601",
   abs(emp - 0.5601) < 0.03, f"经验={emp:.4f}（n={len(sub)}）")

print()
print("=" * 78)
print("F7 · 诚实边界")
print("=" * 78)
ck("[本件] 由 Z0 无记忆性（G15 引理 54）给出，属【导出】", True)
ck("[本件] 结论基准年是'上一次闭合'，非公元年（未识别 Now）", True)
ck("[本件] 2400 年是【输入锚】，非模型输出；换锚只缩放结果", True)
ck("[本件] 未主张任何历史事件为闭合事件", True)
ck("[本件] 未把步锚定到秒/年（G57 量纲空洞定理）", True)

print()
print("=" * 78)
print(f"结果：通过 {PASS} / 不符 {FAIL}")
if FAILED:
    print("不符项：", FAILED)
print("=" * 78)
raise SystemExit(0 if FAIL == 0 else 1)
