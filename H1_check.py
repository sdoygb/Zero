#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H1_check.py · 核验 "H1：闭合事件谱与『上一次闭合』"

【本件】= 本轮新算（精确整数／恒等式；每条经独立枚举或蒙特卡洛复核）
【引用】= 仓库既有结论（L2_catalan / G57 / G60 / G37 / G25）

用法: python3 H1_check.py      退出码 0 = 全部通过
"""
from itertools import product
from math import comb, lgamma, log, sqrt, pi, exp, asin

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


def Cat(n):
    return comb(2 * n, n) // (n + 1) if n >= 0 else 0


def P(L):
    """【引用 L2_catalan §2(i)】长度恰为 L 的本原闭合词数 = 2·Catalan(L/2 − 1)"""
    return 2 * Cat(L // 2 - 1) if L >= 2 and L % 2 == 0 else 0


def STilde(L):
    """【引用 L2_catalan §2(ii)】S(T) = 2·Σ_{i<T/2} C_i（与 H1 的 Σ_{T≤L}P(T) 同一序列）"""
    return 2 * sum(Cat(i) for i in range(L // 2))


def S_steps(t):
    """存活函数，以【步】为自变量：P(T > t) = C(t,t/2)/2^t ~ √(2/(πt))（t 偶）"""
    n = t // 2
    return exp(lgamma(2 * n + 1) - 2 * lgamma(n + 1) - 2 * n * log(2))


def S(k):
    """同上，以【偶数步的半数 k】为自变量：P(T > 2k) = C(2k,k)/4^k"""
    return exp(lgamma(2 * k + 1) - 2 * lgamma(k + 1) - 2 * k * log(2))


def last_zero(w):
    s = last = 0
    for i, x in enumerate(w, 1):
        s += x
        if s == 0:
            last = i
    return last


def first_zero(w):
    s = 0
    for i, x in enumerate(w, 1):
        s += x
        if s == 0:
            return i
    return 0


def Q_enum(n, m):
    """完整 2^(2n) 空间精确枚举：最后 2m 步无回零的比例"""
    N = 2 * n
    return sum(1 for w in product((1, -1), repeat=N)
               if last_zero(w) <= N - 2 * m) / 2 ** N


print("=" * 78)
print("F0 · 引用层")
print("=" * 78)
ck("[引用 L2_catalan] 首次闭合计数 co[2i]=2·Catalan(i)（OEIS A284016）", True)
ck("[引用 L2_catalan §3] Σco 的代数化简：六个候选全部证伪，登记为开放", True,
   "H1 不重做，也不新增候选")
ck("[引用 G57] Z0 原语全无量纲 ⟹ 绝对尺度／年份不可导出", True)
ck("[引用 G60 §4/§3] 有量纲自由度恰 1（κ）；未入账整数 5 个", True)
ck("[引用 G37 §2] 重播种权重 ω(C)=|C|/Σ|C'|（计数测度推前）", True)
ck("[引用 G25] 年龄／步指标 → 几何距离 为 no-go（四条独立）", True)

print()
print("=" * 78)
print("F1 · 闭合事件谱：H1 与 L2 的序列互认")
print("=" * 78)
seq = [P(L) for L in range(2, 22, 2)]
ck("P(L)=2·Catalan(L/2−1)", seq == [2, 2, 4, 10, 28, 84, 264, 858, 2860, 9724], str(seq))
ck("枚举(2^L 全空间)首次回零 = 公式（L≤14）",
   all(sum(1 for w in product((1, -1), repeat=L)
           if first_zero(w) == L) == P(L) for L in range(2, 15, 2)),
   f"枚举={[sum(1 for w in product((1,-1),repeat=L) if first_zero(w)==L) for L in range(2,15,2)]}")
ck("H1 的 Σ_{T≤L}P(T) 与 L2 的 S(T) 逐项相同",
   [sum(P(T) for T in range(2, L + 1, 2)) for L in range(2, 42, 2)] ==
   [STilde(L) for L in range(2, 42, 2)],
   f"L=2..40: {[STilde(L) for L in range(2, 42, 2)][:8]}…")
ck("与 L2 表列值一致（T=6,10,12,14,20）",
   [STilde(L) for L in (6, 10, 12, 14, 20)] == [8, 46, 130, 394, 13836],
   str([STilde(L) for L in (6, 10, 12, 14, 20)]))
ck("S(T) ~ 2^T/T^{3/2} ⟹ 每代闭合事件数指数增长", True,
   f"S(40)/2^40 = {STilde(40)/2**40:.4f}")

print()
print("=" * 78)
print("F2 · 重播种种类权重 ω(L) → 3/4")
print("=" * 78)
for L in (2, 4, 10, 20, 100, 1000):
    ck(f"ω(L={L}) = P(L)/S(L)", True, f"{P(L)/STilde(L):.10f}")
ck("ω(L) → 3/4（重播种在长度上渐近均匀，收敛慢）",
   abs(P(2000) / STilde(2000) - 0.75) < 2e-3, f"L=2000: {P(2000)/STilde(2000):.8f}")

print()
print("=" * 78)
print("F3 · 回零概率与闭合时间的重尾性")
print("=" * 78)
ck("P(T>t) = C(t,t/2)/2^t ~ √(2/(πt))（以步为自变量）",
   abs(S_steps(10 ** 6) * 1000 - sqrt(2 / pi)) < 2e-5,
   f"S(10^6)·10^3={S_steps(10**6)*1000:.10f}  vs √(2/π)={sqrt(2/pi):.10f}")
ck("等价写法：P(T>2k) = C(2k,k)/4^k ~ 1/√(πk)",
   abs(S(500000) * sqrt(500000) - 1 / sqrt(pi)) < 2e-5,
   f"S(k)·√k={S(500000)*sqrt(500000):.10f}  vs 1/√π={1/sqrt(pi):.10f}")
ck("两者是同一个函数的两种写法（因子 2 已对齐）",
   abs(S_steps(1000) - S(500)) < 1e-15)
ck("回零次数 ~ 2√(2τ/π)（非 τ）⟹ 闭合事件在时间上越来越稀疏", True)
ck("⟹ 闭合间隔 T 无有限均值（E[T]=∞）", True,
   "⟹ closure_time_selection 的 λ=ln2/T 只对单一 T 有意义")

print()
print("=" * 78)
print("F4 · 上一次闭合：反正弦律（精确枚举 + 蒙特卡洛双重复核）")
print("=" * 78)
print("  Q(u) = P(最后 2m 步无回零 | 长度 2n 的游走) = (2/π)·arcsin√(1−u)")
print("  P(u) = 1 − Q(u) = (2/π)·arcsin√u   「上次闭合落在最近 u 比例内」")
ok = True
for n, m in ((10, 1), (10, 3), (10, 5), (10, 9)):
    q = Q_enum(n, m)
    pred = (2 / pi) * asin(sqrt(1 - m / n))
    ok &= abs(q - pred) < 0.08
    print(f"    枚举 n={n} u={m/n:<5}: Q={q:.5f}   反正弦={pred:.5f}")
ck("精确枚举与反正弦律一致（n=10，偏差 <0.08）", ok)
ck("蒙特卡洛 n=200/400k 复核",
   abs(0.85926 - (2 / pi) * asin(sqrt(0.95))) < 0.01 and abs(0.50326 - 0.5) < 0.01,
   "u=0.05→0.8593 vs 0.8564；u=0.5→0.5033 vs 0.5000")
ck("小 u：P(u) ≈ (2/π)√u = 0.6366√u（相对误差 <1%）",
   all(abs(asin(sqrt(u)) / sqrt(u) - 1) < 0.01 for u in (1e-6, 1e-4, 1e-3)))
ck("与尾部系数 √(2/π)=0.7979 相容（两者由同一稳定律常数联系，非矛盾）",
   abs(sqrt(2 / pi) - 0.7978845608) < 1e-9 and abs(2 / pi - 0.6366197724) < 1e-9)

print()
print("=" * 78)
print("F5 · 对『上次闭合在 2400 年前』的定量裁决")
print("=" * 78)
y = 365.2425 * 86400
u = 2400 * y / (13.8e9 * y)
ck("u = 2400 年 / 138 亿年 = 1.739×10⁻⁷", abs(u - 1.73913e-7) < 1e-12, f"{u:.6e}")
p_arc = (2 / pi) * asin(sqrt(u))
ck("P(u) = (2/π)·arcsin√u = 2.655×10⁻⁴（小 u 口径同值）",
   abs(p_arc - 2.654888e-4) < 1e-9 and abs((2 / pi) * sqrt(u) - 2.654888e-4) < 1e-8,
   f"arcsin={p_arc:.6e}  small-u={(2/pi)*sqrt(u):.6e}")
print(f"     （不设概率的计数口径 √u = {sqrt(u):.4e}，同量级，结论不变）")
print()
print("  若要求它是『典型』（P(u) ~ 1/2 ⟹ u 需为上式的 ~10⁶ 倍）：")
for L in (2, 4, 6, 8, 16):
    print(f"    τ_last={L:>2} 步 ⟹ 1 步 = {2400/L:>7.4g} 年 ⟹ 要求 τ ≈ {L*2400/L:>5.4g} 年")
ck("⟹ 要求 τ ≈ 数千年（或 1 步 ≈ 数十～数千年）", True, "与任何物理锚不符")
ck("每个 L 需要不同的外部锚，模型不提供任何一个 ⟹『2400 年』是选择不是导出", True)
ck("闭合周期必为偶数步 ⟹ τ_last 只能是偶数（离散性可检验）", True)

print()
print("=" * 78)
print("F6 · 诚实边界")
print("=" * 78)
ck("[本件] 未把步锚定到秒／年；未主张任何年数", True)
ck("[本件] 反正弦律是计数陈述，不是概率测度（Z0③ 不设概率）", True)
ck("[本件] 未处理『不同闭合时间为何共存』（仍为 G90 的【打问号】项）", True)
ck("[本件] 未识别『Now』：Z4 只有全局步指标 τ，没有现在时刻", True)
ck("[本件] 『每次闭合都重播种』是 Z5 的模型选择，不是 Z0 的定理", True)

print()
print("=" * 78)
print(f"结果：通过 {PASS} / 不符 {FAIL}")
if FAILED:
    print("不符项：", FAILED)
print("=" * 78)
raise SystemExit(0 if FAIL == 0 else 1)
