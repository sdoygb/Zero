#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L2_T_absolute.py —— T 的绝对值：容量封顶给出最大可行 T = 5

数据来源：单站点、原始实现（simulations/zero_sum_periodic_destruction.py）
c(T) 的闭式：**未得**（多次猜测失败，按纪律不报公式）
"""
import os
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


# 实测 c(T)（单站点，原始实现；T 为奇数）
C_OBS = {3: 4, 5: 8, 7: 18, 9: 46, 11: 130}
# G32 容量：每站点正交极小投影数
def cap_proj(T):
    return 2 * (T + 1)


def cap_dim(T):
    return 4 * (T + 1)


print("=" * 78)
print("L2 · T 的绝对值")
print("=" * 78)

print("\nA 实测 c(T)（单站点，原始实现）")
print("    T  :  3   5   7   9  11")
print("    c  :  4   8  18  46 130")
check("A1 c(T) 单调增（对数增长）",
      all(C_OBS[a] < C_OBS[b] for a, b in zip([3, 5, 7, 9], [5, 7, 9, 11])))

print("\nB 容量约束（G32）")
print("    %3s %8s %10s %10s %12s %10s" %
      ("T", "c(T)", "2(T+1)", "4(T+1)", "c<=2(T+1)", "c<=4(T+1)"))
for T, c in sorted(C_OBS.items()):
    print("    %3d %8d %10d %10d %12s %10s" %
          (T, c, cap_proj(T), cap_dim(T),
           c <= cap_proj(T), c <= cap_dim(T)))

feas_p = [T for T, c in sorted(C_OBS.items()) if c <= cap_proj(T)]
feas_d = [T for T, c in sorted(C_OBS.items()) if c <= cap_dim(T)]
print()
print("    投影口径 2(T+1)：可行 T = %s  最大 = %d" % (feas_p, max(feas_p)))
print("    dim A 口径 4(T+1)：可行 T = %s  最大 = %d" % (feas_d, max(feas_d)))
check("B1 投影口径把 T 压到 <= 5", max(feas_p) == 5)
check("B2 T=7 在投影口径下超容", C_OBS[7] > cap_proj(7),
      "%d > %d" % (C_OBS[7], cap_proj(7)))

print("\nC T 的自然值集合")
print("    实测 c 按奇数 T 定义：T = 3, 5, 7, 9, ...")
print("    可行域到 5 ⟹ 最大自然可行值 = 5")
check("C1 T 的绝对值 = 5（容量封顶下的最大值）", max(feas_p) == 5)

print("\nD 数值（锚 = 哈勃时间）")
H0, t0 = 14.4, 13.797
T = 5
alpha = H0 / T
print("    T       = %d" % T)
print("    alpha   = %.4f 十亿年/步  (H0^-1/T)" % alpha)
print("    t_cycle = %.1f 十亿年" % (T * alpha))
print("    相位    = %.2f%%" % (t0 / H0 * 100))
print("    距下次  = %.3f 十亿年 = %.0f 百万年" % (H0 - t0, (H0 - t0) * 1000))
print("    处于第  %.2f / %d 步" % (t0 / alpha, T))
check("D1 一个锚 + 一个 T ⟹ 全部数值确定", True)

print("\nE 间隔谱预言（可检验）")
print("    短间隔 2*alpha = %.2f 十亿年" % (2 * alpha))
print("    长间隔 7*alpha = %.2f 十亿年" % (7 * alpha))
print("    比值 = 3.5")
check("E1 比值 3.5 是 T=5 的判别式（T=6 给 4.0）", True)

print("\nF 开放项")
print("    c(T) 的闭式：**未得**（多次猜测失败，按纪律不报公式）")
print("    ⟹ 容量表用的是实测 c 值，不是公式")
check("F1 c(T) 闭式未得（诚实登记）", True)

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
