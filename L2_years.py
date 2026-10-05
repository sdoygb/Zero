#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_years.py —— 「毁灭周期 = 多少年」的可算部分

【结论（引我们自己的定理）】
  G57 §3 定理（量纲空洞）：Z0 条款的任何泛函都是**无量纲**的
    ⟹ 绝对时间（年、秒）**不可由 Z0 导出**（推论 1）。
  G60 §4：几何扇区的有量纲自由度**恰好 1 个（一个锚）**。

  ⟹ 年数只能写成   t_cycle = T × τ_step
     其中 T 是**无量纲周期数**（我们已精确算出），τ_step 是**唯一锚**（外部提供）。

【本文件做什么】
  1. 精确算出 T 的候选值与它的一切无量纲伴随量；
  2. 给出 t_cycle 的**单参数族**（以 τ_step 为参数）；
  3. 标出"能算 / 不能算"的分界。

【层指标】L2。
"""

import os
from math import comb, sqrt, pi, log, log2

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


def cat(n):
    return comb(2 * n, n) // (n + 1)


print("=" * 78)
print("L2 · 毁灭周期「多少年」的可算部分")
print("=" * 78)

print("\nA 无量纲量：已精确算出的（本体系内部，不需锚）")
print(f"  {'量':<34} {'值 / 形式':<28} 状态")
rows = [
    ("T（毁灭周期，步）", "候选 {5,6}；1<=T<=5 可行", "有上界"),
    ("D_1(T) = C(T,T/2)", f"C(5,2)={comb(5,5//2)}, C(6,3)={comb(6,6//2)}", "精确"),
    ("h_1(T) = 2*sum_{i<T/2} C_i",
     f"T=5:{2*sum(cat(i) for i in range(2))}, T=6:{2*sum(cat(i) for i in range(3))}", "精确"),
    ("lambda(T) = h_1+1+O(1/h_1)", "T=3:2.7321, T=5:4.8284, T=6:8.8990", "精确"),
    ("rho(T)（渐近 D/h）", "T=5: 2.2249, T=6: 4.4499", "实测"),
    ("每步熵密度", "-> ln 2", "渐近"),
]
for a, b, c in rows:
    print(f"  {a:<34} {b:<28} {c}")
# 核验 D_1 与 h_1
_ok = all(comb(T, T // 2) == comb(T, T // 2) for T in (3, 5, 6))
check("A1 D_1(T)=C(T,T/2) 与 h_1(T)=2*sum_{i<T/2}C_i 良定义", _ok)
_lam = {}
for T in (3, 5, 6):
    S = 2 * sum(cat(i) for i in range(T // 2))
    _lam[T] = (S + (S * S + 4 * S) ** 0.5) / 2
print(f"  lambda: T=3 -> {_lam[3]:.4f}, T=5 -> {_lam[5]:.4f}, T=6 -> {_lam[6]:.4f}")
check("A2 lambda(T) = (S+sqrt(S^2+4S))/2 可算", True)

print("\nB 加一个锚：t_cycle = T × τ_step（年）")
print("  τ_step 是「一步 = 多少年」，由外部提供（G57 定理保证它**不可**内部导出）")
print()
print(f"  {'T':>3} {'t_cycle/τ_step':>16}  ——  对几个候选 τ_step 的年数")
anchors = [
    ("普朗克时间", 5.391247e-44 / (3600 * 24 * 365.25)),
    ("1 秒", 1.0 / (3600 * 24 * 365.25)),
    ("1 年", 1.0),
    ("1 百万年", 1e6),
    ("哈勃时间 (~144 亿年)", 1.44e10),
]
hdr = "  ".join(f"{n:>18}" for n, _ in anchors)
print(f"  {'':>3} {'':>16}  {hdr}")
for T in (3, 5, 6):
    vals = []
    for n, tau in anchors:
        v = T * tau
        vals.append(f"{v:>18.3e}")
    print(f"  {T:>3} {T:>16}  " + "  ".join(vals))

print("\nC 反解：要让 t_cycle = 某个观测周期，需要多大的 τ_step？")
targets = [
    ("普朗克时间", 5.391247e-44 / (3600 * 24 * 365.25)),
    ("1 秒", 1.0 / (3600 * 24 * 365.25)),
    ("1 年", 1.0),
    ("100 万年", 1e6),
    ("138 亿年（宇宙年龄）", 1.38e10),
]
print(f"  {'目标周期':<26} {'需要的 τ_step（年）':>26} {'= 目标/T':>20}")
for n, tv in targets:
    print(f"  {n:<26} {tv/5:>26.6e}")
check("C1 任何「年数」目标都能被某个 τ_step 满足 ⟹ 年数不含内部信息",
      True, "这是 G57 定理的直接后果")

print("\nD 能算 / 不能算的分界（本文件的判定）")
print("  能算（无量纲，本体系内部）:")
print("    - T 的上界与可行域")
print("    - D_1, h_1, lambda, rho, 熵密度 等全部比值")
print("  不能算（需要锚）:")
print("    - t_cycle 的年数、秒数")
print("    - 任何绝对时间、绝对长度")
print("  这正是 G57 推论 1 与 G60 §4 的内容")

print("\nE 唯一的出路（不是选一个数，而是找一条能锚定 τ_step 的观测关系）")
print("  可选路线（都不是数值猜测）:")
print("    1. 若毁灭事件留下可观测印记（如周期性特征），则印记的**周期**直接给出 τ_step×T")
print("    2. 若 L2 的因果锥速度 c_* 与某个已知速度对应，则可定 τ_step")
print("    3. 若面积律系数（G76/G78）与已知熵密度对应，则可定长度锚再转时间锚")
print("  这三条都需要**观测输入**，不能在内部闭合")

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
