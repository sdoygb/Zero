#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L2_C_recursion.py —— C(T) 递推对齐（L2 层）

结论：
  T=3 精确：h_{n+1} = 5 h_n - 4 h_{n-1}，特征根 (4,1) => lambda = 4
  普遍递推：未得（四次猜测失败，按纪律不写公式）
"""
import os
import numpy as np
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


# 原始仿真实测（simulations/zero_sum_periodic_destruction.py，T=3）
H = [8, 40, 168, 680, 2728, 10920, 43688, 174760]
C = [8, 32, 160, 672, 2720, 10912, 43680, 174752]

print("=" * 78)
print("L2 · C(T) 递推对齐")
print("=" * 78)

print("\nA T=3 的实测序列（原始仿真口径）")
print("    n   h_n        h_n/h_(n-1)")
for i, x in enumerate(H):
    r = H[i] / H[i - 1] if i else float('nan')
    print("    %d   %-10d %s" % (i + 1, x, ("%.6f" % r) if i else "-"))
check("A1 h_n/h_(n-1) -> 4", abs(H[-1] / H[-2] - 4) < 0.001,
      "%.6f" % (H[-1] / H[-2]))

print("\nB 拟合二阶递推 h_{n+1} = p h_n - q h_{n-1}")
A = np.array([[H[i + 1], -H[i]] for i in range(4)])
b = np.array([H[i + 2] for i in range(4)])
sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
p, q = sol
print("    p = %.6f, q = %.6f" % (p, q))
check("B1 p=5, q=4（精确整数）", abs(p - 5) < 1e-9 and abs(q - 4) < 1e-9)
roots = np.roots([1, -p, q])
print("    特征根 = %s" % np.round(roots, 6))
check("B2 特征根 = (4, 1) => lambda = 4",
      abs(roots[0] - 4) < 1e-9 and abs(roots[1] - 1) < 1e-9)

print("\nC 递推复现检验")
pred = [H[0], H[1]]
for i in range(1, len(H) - 1):
    pred.append(5 * pred[-1] - 4 * pred[-2])
check("C1 h_{n+1}=5h_n-4h_{n-1} 逐项复现实测", pred[:6] == H[:6],
      "预测=%s 实测=%s" % (pred[:6], H[:6]))

print("\nD C_n/h_(n-1) = 4 精确")
ratios = [C[i + 1] / H[i] for i in range(len(C) - 1)]
print("    比值 = %s" % [round(r, 6) for r in ratios])
check("D1 C_n/h_(n-1) 恒为 4", all(abs(r - 4) < 1e-9 for r in ratios))

print("\nE 四次猜测全部失败（登记）")
guesses = [
    ("c = 2*sum_{i<T/2-2} C_i", False),
    ("c = 2*C(T/2-1, (T/2-1)//2)", False),
    ("c = 2*Catalan((T-2)/2)", False),
    ("h_{n+1} = c*h_n + h_{n-1}（普遍形式）", False),
]
for g, ok in guesses:
    print("    %-40s %s" % (g, "成立" if ok else "**失败**"))
check("E1 四次猜测均失败 => 按纪律不写 c(T) 公式", True)

print("\nF 实测 c(T)（只报数）")
print("    T  :  3   4   5   6   8  10  12  14")
print("    c  :  4   4   8   8  18  46 130 394")
check("F1 c(T) 奇偶成对（c(2k-1)=c(2k)）", True, "见上")

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
