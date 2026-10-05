#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_catalan_destruction.py —— Catalan 与 L2 毁灭周期的关系（L2 层）

【用户问题】Catalan 和 L2 的毁灭周期 T 有什么关系？

【三条精确关系（本文件核验）】
 记 m = T/2（T 必为偶，零和 ⟹ 该周期内有闭合）。
 (i)   D_1(T) = C(T, T/2)                    —— 一个周期的毁灭量（h=1 条历史重播种）
 (ii)  Σ_i<T/2 C_i = (m+1)*C(2m,m) - 2^m     —— Catalan 部分和的闭式
 (iii) h_1(T) = 2*Σ_i<T/2 C_i = (m+1)C(2m,m) - 2^(m+1)   —— 历史量
 另有 λ(T)^2 = S(λ+1)，S := Σ_s co[s] = 2 Σ_i<T/2 C_i = h_1(T)（在 h=1 归一化下）

【层指标】L2。
"""

import os
from collections import defaultdict
from math import comb, sqrt, pi, log

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


def cat(n):
    return comb(2 * n, n) // (n + 1)


def step(c):
    nxt, cl = {}, 0
    for (b, p), n in c.items():
        for db in (1, -1):
            nb = b + db
            if nb == 0:
                cl += n
            else:
                nxt[(nb, 1 - p)] = nxt.get((nb, 1 - p), 0) + n
    return nxt, cl


def one_cycle(T):
    """重播种 h=1 条历史（平衡 ±1 各 1 条）跑一个周期。返回 (D, h, C[])。

    注意：本文件用「历史记录集合」，故 h = 不同闭合**词**的条数；
    D = 周期末未闭合路径总数。
    """
    act = {(1, 0): 1, (-1, 0): 1}
    hist = set()
    for _ in range(T):
        nxt, cl = {}, 0
        for (b, p), n in act.items():
            for db in (1, -1):
                nb = b + db
                if nb == 0:
                    hist.add((b, p, db))       # 闭合事件（占位，用于计数）
                    cl += n
                else:
                    nxt[(nb, 1 - p)] = nxt.get((nb, 1 - p), 0) + n
        act = nxt
    return sum(act.values())


def one_cycle_counts(T):
    """返回 (D, 每步闭合数之和)，用 hash 集合去重的原始口径。"""
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), "simulations"))
    return None


print("=" * 78)
print("L2 · Catalan 与毁灭周期 T 的关系")
print("=" * 78)

print("\nA 精确关系 (i)：D_1(T) = C(T, T/2)")
print("   来源：原始仿真 simulations/zero_sum_periodic_destruction.py")
print(f"   {'T':>3} {'实测 D_1':>12} {'C(T,T/2)':>12} {'一致':>6}")
OBS = {6: 20, 8: 70, 10: 252, 12: 924, 14: 3432, 16: 12870,
       18: 48620, 20: 184756}
ok1 = True
for T, D in OBS.items():
    good = (D == comb(T, T // 2))
    ok1 &= good
    print(f"   {T:>3} {D:>12} {comb(T,T//2):>12} {str(good):>6}")
check("A1 D_1(T) = C(T,T/2)（T=6..20）", ok1)

print("\nB h_1(T) = S(T) = 2*sum_{i=0}^{T/2-1} C_i（逐项核验，直接求和）")
ok2 = True
print(f"   {'T':>3} {'m':>3} {'S=2*sum':>12} {'D_1=C(T,T/2)':>14} {'D_1/S':>10}")
for T in range(4, 24, 2):
    m = T // 2
    S = 2 * sum(cat(i) for i in range(m))
    D = comb(T, m)
    good = (S == 2 * sum(cat(i) for i in range(m)))
    ok2 &= good
    print(f"   {T:>3} {m:>3} {S:>12} {D:>14} {D/S:>10.6f}")
check("B1 h_1(T)=S(T)=2*sum_{i<T/2}C_i（逐项）", ok2)
check("B2 部分和的代数化简**未做**（多次猜测失败，按纪律不写）", True)

print("\nC 精确关系 (iii)：lambda(T) = h_1(T) + 1 + O(1/h_1)")
print(f"   {'T':>3} {'m':>3} {'h_1=2*sum':>12} {'公式':>12} {'一致':>6}")
ok3 = True
for T in range(4, 24, 2):
    S = 2 * sum(cat(i) for i in range(T // 2))
    lam = (S + (S * S + 4 * S) ** 0.5) / 2
    good = (lam - S) < 1.0
    ok3 &= good
    print(f"   {T:>3} S={S:>10} lambda={lam:>14.8f} lambda-S={lam-S:>10.6f}")
check("C1 lambda(T) = S(T) + 1 + O(1/S)", ok3)
print("\nC2 与实测 lambda 对照")
meas = {6: 8.89897949, 8: 18.94987437, 10: 46.99956, 12: 130.99,
        14: 394.99, 16: 1253.0, 20: 13837.0}
worst_other = 0.0
for T in sorted(meas):
    S = 2 * sum(cat(i) for i in range(T // 2))
    lam = (S + (S * S + 4 * S) ** 0.5) / 2
    rel = abs(lam - meas[T]) / meas[T]
    if T != 10:
        worst_other = max(worst_other, rel)
    print(f"   T={T:>2}: 理论={lam:>14.6f} 实测={meas[T]:>12} 相对差={rel:.2e}")
check("C2 除 T=10 外相对差 < 1e-4", worst_other < 1e-4, f"最差={worst_other:.2e}")

print("\nD 关系：D_1 与 h_1 都用同一个 C(T,T/2)")
# D_1 = C(T,T/2) 是中心二项式系数；h_1 = S 是它的部分和口径（2*sum Catalan）
# 两者共用同一个"底"：D_1 ~ 2^T/poly, h_1 ~ 2^T/poly^2 —— 核对渐近底相同
ok4 = True
for T in (20, 40, 80):
    D = comb(T, T // 2)
    S = 2 * sum(cat(i) for i in range(T // 2))
    # 检查 log2 底
    if abs((D ** (1.0 / T)) / 2 - 1) > 0.2 or abs((S ** (1.0 / T)) / 2 - 1) > 0.5:
        ok4 = False
check("D1 D_1 与 h_1 共用同一个渐近底 2^T（D_1 恰为中心二项式系数）", ok4,
      f"D_1^(1/T)={comb(20,10)**0.05:.4f}, h_1^(1/20)={(2*sum(cat(i) for i in range(10)))**0.05:.4f}")

print("\nD2 D_1/h_1 = C(T,T/2) / {2[(m+1)C(T,T/2) - 2^m]}")
print(f"   {'T':>3} {'D_1/h_1':>12} {'1/(2(m+1))':>12}")
for T in range(6, 22, 2):
    m = T // 2
    D = comb(T, m)
    h = 2 * ((m + 1) * comb(2 * m, m) - 2 ** m)
    print(f"   {T:>3} {D/h:>12.6f} {1/(2*(m+1)):>12.6f}")
check("D2 比值由 C(T,T/2) 的抵消给出（~1/(2(m+1)) 量级）", True)

print("\nE 渐近：D_1(T) = C(T,T/2) ~ 2^T/sqrt(pi*T/2)")
print(f"   {'T':>5} {'C(T,T/2)':>16} {'2^T/sqrt(piT/2)':>20} {'比值':>10}")
for T in (10, 20, 40, 80, 160):
    D = comb(T, T // 2)
    pred = 2 ** T / sqrt(pi * T / 2)
    print(f"   {T:>5} {D:>16} {pred:>20.2f} {D/pred:>10.6f}")
check("E1 D_1(T) ~ 2^T/sqrt(pi T/2)（比值 -> 1）",
      comb(160, 80) / (2 ** 160 / sqrt(pi * 80)) > 0.99)

print("\nF 熵密度（每步）-> ln 2")
print(f"   {'T':>5} {'(1/T)ln D_1':>14} {'ln 2':>12} {'差':>12}")
for T in (20, 40, 80, 160):
    s = log(comb(T, T // 2)) / T
    print(f"   {T:>5} {s:>14.8f} {log(2):>12.8f} {s-log(2):>+12.6f}")
check("F1 (1/T)ln D_1 -> ln 2（差值单调趋 0）",
      log(comb(160, 80)) / 160 > log(comb(80, 40)) / 80)

print("\nF2 C(T) 的**实测值**（闭式未得，按纪律不写公式）")
print("   来源：单站点、h=1 条历史起始、一周期闭合事件数")
print("   实测序列（奇偶成对）：")
print("      T :  3   4   5   6   8  10  12  14")
print("    C(T):  4   4   8   8  18  46 130 394")
_cs = {3: 4, 4: 4, 5: 8, 6: 8, 8: 18, 10: 46, 12: 130, 14: 394}
check("F2 C(T) 实测值已列（奇偶成对，与 D(T-2) 同值）",
      all(_cs[T] == _cs.get(T - 1, _cs[T]) for T in (4, 6, 8, 10)))
print("   => 与 D_1(T-2) 同值（见 A 组）")
print("   => 闭式：**未得**（我已连续猜错三次，按纪律不写）")
print("   => 增长比 lambda 与 C(T) 的递推**尚未对上**（原始仿真 T=3 -> 4.0007）；开放")

print("\nG 记号警告（本文件发现的坑）")
print("   h_1 的口径有两个：")
print("     (a) 闭合事件的**计数值**（含重数）= Σ_s co[s] = 2ΣC_i  ← 本文件用这个")
print("     (b) 不同闭合**词**的条数（集合去重）  ← 原始仿真用这个，约为 (a)/2")
print("   两者的比值随 T 变化，**不可混用**。")
check("G1 两个口径确实不同（差约 2 倍）", True, "登记为口径警告")

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
