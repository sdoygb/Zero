# -*- coding: utf-8 -*-
# Z086 核验：Σw=0 ⟹ 2 自由度（模 b + 辐角 θ = 同一复数 ŝ₁）；载体 = U(1)
# 口径：质量 MeV，同 Z082 §1.3；DFT 同 Z082 §2。
import cmath, math

# 实测质量（MeV）
m = (0.51099895, 105.6583755, 1776.86)
s = tuple(math.sqrt(x) for x in m)

# 平均与偏离
sbar = sum(s) / 3.0
w = tuple(x / sbar - 1.0 for x in s)

# --- 断言 1：Σw = 0（精确到浮点）---
assert abs(sum(w)) < 1e-12, sum(w)

# --- 断言 2：w_2 = -(w_0 + w_1)（第三个是前两个的负和）---
assert abs(w[2] + (w[0] + w[1])) < 1e-12

# --- 断言 3：三项互异（不全在 1 维，秩 = 2）---
assert len(set(round(x, 6) for x in w)) == 3

# --- DFT ---
N = 3
s0 = sum(s) / 3.0
s1 = sum(s[i] * cmath.exp(-2j * math.pi * 1 * i / 3) for i in range(3)) / 3.0
s2 = sum(s[i] * cmath.exp(-2j * math.pi * 2 * i / 3) for i in range(3)) / 3.0

# --- 断言 4：ŝ₂ = conj(ŝ₁)（2 自由度：1 个实数 ŝ₀ + 1 个复数 ŝ₁）---
assert abs(s2 - s1.conjugate()) < 1e-9

# --- 断言 5：arg ŝ₁（辐角 θ）---
theta = cmath.phase(s1)
assert abs(theta - 2.316624733883) < 1e-9

# --- 断言 6：|ŝ₁|/|ŝ₀| = b/2，b 为模 ---
ratio = abs(s1) / abs(s0)
assert abs(ratio - 0.707100252592) < 1e-9
b = 2.0 * ratio
assert abs(b - 1.414200505185) < 1e-9

# --- 断言 7：b 与 θ 是同一复数 ŝ₁ 的模与辐角 ---
assert abs(abs(s1) - ratio * abs(s0)) < 1e-9
assert abs(cmath.phase(s1) - theta) < 1e-12

# --- 断言 8：θ 与 2π/3 + 2/9 的偏差（开放/输入侧）---
twopiw3 = 2 * math.pi / 3
delta = theta - twopiw3 - 2.0 / 9.0
assert abs(delta - 7.409e-6) < 1e-8

print("m  =", m)
print("s  =", [round(x, 6) for x in s])
print("<s> =", round(sbar, 9))
print("w  =", [round(x, 6) for x in w])
print("sum w      =", round(sum(w), 12))
print("w0+w1      =", round(w[0] + w[1], 12))
print("-w2        =", round(-w[2], 12))
print("ŝ0 =", round(s0.real, 9))
print("ŝ1 =", round(s1.real, 9), "+", round(s1.imag, 9), "i")
print("arg ŝ1     =", theta, "rad")
print("2π/3       =", round(twopiw3, 9))
print("2/9        =", round(2.0 / 9.0, 9))
print("δθ (arg-2π/3-2/9) =", round(delta, 12))
print("|ŝ1|/|ŝ0|  =", round(ratio, 12))
print("b          =", round(b, 12))
print("自由度：3 个偏离 + 1 个约束(Σw=0) = 2；ŝ₁ 提供模+辐角 = 2 参数")
print("ALL ASSERTIONS PASSED")
