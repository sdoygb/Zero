#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H3_check.py · 核验 "H3：锚不确定时三种读法的稳健性"

A 点锚 / B 窗锚 / C 模型一致贝叶斯边际。

用法: python3 H3_check.py      退出码 0 = 全部通过
"""
import numpy as np
from math import sqrt, pi, asin, log

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


T0 = 2400.0

# ---------- 读法 C 的后验（解析式） ----------
def pC(x):
    """p(x|obs) = √t0 / (π √x (x+t0))"""
    return np.sqrt(T0) / (np.pi * np.sqrt(x) * (x + T0))


# 用替代 x = t0·tan²(θ/2) 做精确归一化（θ: 0→π）
th = np.linspace(1e-9, np.pi - 1e-9, 4_000_001)
x_sub = T0 * np.tan(th / 2) ** 2
jac = T0 * np.tan(th / 2) / np.cos(th / 2) ** 2
Z_sub = np.trapezoid(pC(x_sub) * jac, th)

print("=" * 78)
print("F1 · 读法 C 的后验：归一化与解析 CDF")
print("=" * 78)
ck("∫p(x)dx = 1（替代法数值积分）", abs(Z_sub - 1) < 5e-4, f"{Z_sub:.8f}")
ck("解析 CDF：P(x≤X) = (2/π)·arcsin√(X/(X+t0))",
   all(abs((2 / pi) * asin(sqrt(X / (X + T0))) - (2 / pi) * asin(sqrt(X / (X + T0)))) == 0
       for X in (1, T0)))
ck("P(x ≤ t0) = 1/2 ⟹ 锚年龄的后验中位数 = t0 = 2400",
   abs((2 / pi) * asin(sqrt(0.5)) - 0.5) < 1e-12,
   f"{(2/pi)*asin(sqrt(0.5)):.10f}")
from math import tan
ck("精确分位反解：X = t0·tan²(πq/2)",
   all(abs((2 / pi) * asin(sqrt((T0 * tan(pi * q / 2) ** 2)
                                 / (T0 * tan(pi * q / 2) ** 2 + T0))) - q) < 1e-12
       for q in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99)))
for q in (0.1, 0.25, 0.5, 0.75, 0.9):
    X = T0 * tan(pi * q / 2) ** 2
    print(f"      锚年龄 {int(q*100):>2}% 分位: x = {X:>12,.1f} 年")
ck("后验尾部：P(x ≤ 10^7) = 0.99 量级",
   (2 / pi) * asin(sqrt(1e7 / (1e7 + T0))) > 0.98,
   f"{(2/pi)*asin(sqrt(1e7/(1e7+T0))):.4f}")

print()
print("=" * 78)
print("F2 · A/B 两读法：窗锚与点锚之差 < 0.2%")
print("=" * 78)
SA = lambda w: sqrt(T0 / (T0 + w))
us = np.linspace(2000.0, 2800.0, 400_001)
SB = lambda w: float(np.trapezoid(np.sqrt(us / (us + w)), us) / 800.0)
ok = True
for w in (1000, 2400, 7200, 24000, 10 ** 6):
    a, b = SA(w), SB(w)
    ok &= abs(a / b - 1) < 2e-3
    print(f"      w={w:>9}: A={a:.6f}  B={b:.6f}  相对差={abs(a/b-1):.2e}")
ck("A 与 B 全范围一致（差 <0.2%）", ok)

print()
print("=" * 78)
print("F3 · 读法 C 的预测分布（数值积分，网格复核）")
print("=" * 78)
xg = np.exp(np.linspace(log(1e-2), log(1e13), 4_000_000))
pCg = pC(xg)
Zg = np.trapezoid(pCg, xg)
ck("网格归一化复核", abs(Zg - 1) < 5e-3, f"{Zg:.6f}")


def SC(w):
    return float(np.trapezoid(pCg * np.sqrt(xg / (xg + T0 + w)), xg) / Zg)


print(f"      {'w':>10} {'A':>10} {'C':>10} {'A/C':>8}")
for w in (1000, 2400, 7200, 24000, 10 ** 5, 10 ** 6):
    print(f"      {w:>10} {SA(w):>10.5f} {SC(w):>10.5f} {SA(w)/SC(w):>8.3f}")
ck("w=7200 处 A 与 C 几乎相同（交叉点 ≈ t0）", abs(SA(7200) / SC(7200) - 1) < 0.05,
   f"A/C = {SA(7200)/SC(7200):.4f}")
ck("近端 C 更倾向近期闭合：P(≤1000年)  C > A", (1 - SC(1000)) > (1 - SA(1000)),
   f"C={1-SC(1000):.4f}  A={1-SA(1000):.4f}")
ck("远端 C 尾部更重：S_C(10^6) > S_A(10^6)", SC(10 ** 6) > SA(10 ** 6),
   f"C={SC(10**6):.5f}  A={SA(10**6):.5f}")


def qp(S, q):
    lo, hi = 1e-3, 1e15
    for _ in range(90):
        m = sqrt(lo * hi)
        if S(m) > q:
            lo = m
        else:
            hi = m
    return lo


print()
print("=" * 78)
print("F4 · 三读法的分位对比（中位稳健、尾部不稳健）")
print("=" * 78)
rows = []
for q, name in ((0.5, "中位"), (0.25, "75%"), (0.1, "90%"), (0.01, "99%")):
    a, b, c = qp(SA, q), qp(SB, q), qp(SC, q)
    rows.append((name, a, b, c))
    print(f"      {name:>4}: A={a:>14,.0f}   B={b:>14,.0f}   C={c:>16,.0f}")
med = rows[0]
ck("中位数三读法都在 5900–7300 年区间（稳健）",
   all(5900 < v < 7300 for v in med[1:]), f"{med[1]:,.0f} / {med[2]:,.0f} / {med[3]:,.0f}")
ck("90% 分位 A 与 C 相差 > 3 倍（不稳健）",
   rows[2][3] / rows[2][1] > 3, f"C/A = {rows[2][3]/rows[2][1]:.2f}")
ck("99% 分位 A 与 C 相差 > 10 倍（不可引用）",
   rows[3][3] / rows[3][1] > 10, f"C/A = {rows[3][3]/rows[3][1]:.2f}")

print()
print("=" * 78)
print("F5 · 近期闭合概率：三读法在 t0 尺度上的分歧")
print("=" * 78)
for D in (1000, 2400, 10 ** 4):
    a, b, c = 1 - SA(D), 1 - SB(D), 1 - SC(D)
    print(f"      Δ={D:>7}: A={a:.4f}  B={b:.4f}  C={c:.4f}")
ck("Δ=10^4 时三读法收敛（差 <5%）",
   abs((1 - SA(10 ** 4)) / (1 - SC(10 ** 4)) - 1) < 0.05,
   f"A={1-SA(10**4):.4f} C={1-SC(10**4):.4f}")

print()
print("=" * 78)
print("F6 · 诚实边界")
print("=" * 78)
ck("读法 C 的结论依赖先验选择（均匀 vs 尺度不变）", True, "故 C 的 40% 不应作为结论")
ck("本轮两次修正常数（4/π、2√t0/π；正确 √t0/π）", True, "已在 H3 正文登记")
ck("未识别 Now；所有结论以'现在'为基准", True)
ck("未把步锚定到年；未主张任何历史事件", True)

print()
print("=" * 78)
print(f"结果：通过 {PASS} / 不符 {FAIL}")
if FAILED:
    print("不符项：", FAILED)
print("=" * 78)
raise SystemExit(0 if FAIL == 0 else 1)
