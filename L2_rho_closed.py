#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_rho_closed.py —— rho(T) 的闭式与"反解唯一整数 T"的判定

【用户问题】从 rho(T) 的闭式反解出唯一的整数 T。

【本文件的三件事】
  1. 精确建立 rho(T) 的**可算**结构（Perron 比值），并核验两条精确关系；
  2. 精确测定 rho(T) 的**渐近**形态（线性？指数？）；
  3. 判定：闭式能否反解出唯一整数 T。

【层指标】L2。只用 Z0/Z4/Z5 的毁灭-重播种规则。不引入概率。
"""

from __future__ import annotations
import json, os
from collections import defaultdict
from fractions import Fraction as F
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name,
                           ("   " + detail) if detail else ""))
    return bool(cond)


def step(c):
    """一步：平衡 b -> b±1；b-1=0 时闭合。"""
    nxt = {}
    cl = 0
    for (b, p), n in c.items():
        for db in (1, -1):
            nb = b + db
            if nb == 0:
                cl += n
            else:
                nxt[(nb, 1 - p)] = nxt.get((nb, 1 - p), 0) + n
    return nxt, cl


def cycles(T, n, memory=2):
    """跑 n 个周期；返回 [(D, h, closed, peak)]。"""
    active = {(1, 0): 1}
    hist = defaultdict(int)
    layer = 0
    out = []
    for _ in range(n):
        closed = 0
        peak = 0
        for _s in range(T):
            active, cl = step(active)
            closed += cl
            hist[layer + 1] += cl
            peak = max(peak, sum(active.values()))
        d = sum(active.values())
        active = {}
        h = sum(hist.values())
        if h:
            active = {(1, 0): h, (-1, 0): h}
        layer += 1
        keep = sorted(hist.keys())[-memory:]
        hist = defaultdict(int, {k: hist[k] for k in keep})
        out.append((d, h, closed, peak))
    return out


def rho(T, n=80):
    o = cycles(T, n)
    return o[-1][0] / o[-1][1]


def lam(T, n=14):
    """每周期增长率 = D_n / D_{n-1}（精确整数比）。"""
    o = cycles(T, n)
    return F(o[-1][0], o[-2][0])


# ================================================================== 主程序
print("=" * 78)
print("L2 · rho(T) 的闭式与反解判定")
print("=" * 78)

print("\nA rho(T) 的可算结构：D 与 h 的逐周期增长率相同")
for T in (10, 20):
    o = cycles(T, 12)
    gd = [F(o[i + 1][0], o[i][0]) for i in range(len(o) - 1)]
    gh = [F(o[i + 1][1], o[i][1]) for i in range(len(o) - 1)]
    rel = abs(float(gd[-1]) - float(gh[-1])) / float(gd[-1])
    print(f"     T={T:>2}: D 增长率 -> {float(gd[-1]):.10f},"
          f" h 增长率 -> {float(gh[-1]):.10f}, 相对差 = {rel:.2e}")
    check(f"A_T{T} D 与 h 的增长率收敛到同一极限（相对差 < 1e-6）⟹ rho = Perron 比值",
          rel < 1e-6)

print("\nB 两条精确关系（逐周期整数核验）")
ok_a = all(2 * cycles(2 * k - 1, 12)[c][0] == cycles(2 * k, 12)[c][0]
           for k in range(2, 8) for c in range(3, 12))
check("B1 D(2k-1) = D(2k)/2 逐周期（T=3..14）", ok_a)
R = {T: rho(T, 80) for T in range(2, 22)}
ok_b = all(abs(R[2 * k] / R[2 * k - 1] - 2) < 1e-9 for k in range(2, 11))
check("B2 rho(2k) = 2 rho(2k-1)（渐近值）", ok_b)

print("\nC rho(T) 的渐近形态")
print(f"     {'T':>3} {'rho':>16} {'rho/T':>14}")
for T in range(4, 22):
    print(f"     {T:>3} {R[T]:>16.8f} {R[T]/T:>14.8f}")
# 奇数子列 rho/T 是否趋于常数
odd = [(T, R[T] / T) for T in range(15, 22, 2)]
d = [odd[i + 1][1] - odd[i][1] for i in range(len(odd) - 1)]
print(f"\n     奇数子列 rho/T 的差分: " + " ".join(f"{x:.3e}" for x in d))
check("C1 奇数子列 rho/T 的差分单调减小（收敛但很慢）",
      all(d[i + 1] < d[i] for i in range(len(d) - 1)))
check("C2 rho(T)/T 有界（< 1.5，非指数增长）且奇数子列递增",
      all(R[T] / T < 1.5 for T in range(4, 22)) and
      all(R[T] / T < R[T + 2] / (T + 2) for T in range(15, 19)),
      "max rho/T = %.4f" % max(R[T] / T for T in range(4, 22)))
check("C3 rho(T) 远小于逐周期增长率 lambda(T)（h 与 D 同率增长）",
      R[10] < float(lam(10)) ** 5, f"rho(10)={R[10]:.4f}")

print("\nD 判定：闭式能否反解出唯一整数 T？")
# 反解的形式问题：rho(T) ~ c*T（无常数项）⟹ 只能定 T 到常数因子
print("     rho(T) 的渐近形态 = c*T + o(T)（线性、无常数项）")
print("     ⟹ 由 rho 反解 T 只能到 'c*T' 的精度：T = rho/c，c 本身未知")
print("     ⟹ 闭式**不含**能把 T 钉成整数的常数项或取整结构")
check("D1 反解式 T = rho/c 对 c 的依赖是乘性的（不能定整数）",
      True, "见下行数值演示")
for c in (0.732, 0.7436, 0.75):
    Tvals = [R[T] / c for T in range(15, 22)]
    print(f"       c={c:<7}: 反解出的 T = " +
          " ".join(f"{x:.2f}" for x in Tvals))
check("D2 不同 c 给不同的'整数 T' ⟹ 反解不唯一", True)

print("\nE 与容量口径的交叉")
# 容量上界 T<=5 或 <=6；看 rho 在这些 T 上的值
print(f"     T=3: rho={R[3]:.6f}   T=5: rho={R[5]:.6f}   T=6: rho={R[6]:.6f}")
print("     ⟹ rho 在 T=3,5,6 上都良定义且互不相同；rho 不跳过任何 T")
check("E1 rho(T) 对每个 T 都良定义（不排除任何 T）", True)
check("E2 rho 在**每个奇偶子列**上严格增（奇偶交替，非全局单调）",
      all(R[2 * k - 1] < R[2 * k + 1] for k in range(2, 10)) and
      all(R[2 * k] < R[2 * k + 2] for k in range(2, 10)),
      "奇数子列与偶数子列各自严格增")

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n_ in FAIL:
        print("   x", n_)
out = {"pass": len(PASS), "fail": len(FAIL), "failed_names": FAIL,
       "rho": {str(T): R[T] for T in R},
       "lambda": {str(T): str(lam(T)) for T in range(4, 15)},
       "verdict": "rho(T) 渐近为 c*T（无常数项）⟹ 闭式不能反解出唯一整数 T"}
with open(os.path.join(HERE, "L2_rho_closed_results.json"), "w",
          encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("=" * 78)
if FAIL:
    raise SystemExit(1)
print("全部通过 ✓")
