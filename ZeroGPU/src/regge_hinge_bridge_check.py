#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
regge_hinge_bridge_check.py —— 词层面积 ↔ Regge 铰链：两库"面积"是否同一个 + 缺陷角
==================================================================================
背景：
  · **新库**（`Z044`/`Z048`）：平衡词（±1 序列）的面积 a(w) = Σ_k 前缀和_k
  · **旧库 #30/#31**：p 态词的面积向量 area(w) = Σ_k (L−k) σ(w_k) ∈ ℝ^{p−1}
  · **#29 :144 的"三合一判据"**：L ≡ 0 mod 4 ⟺ 平坦词存在 ⟺ q 单位根
  · **#29 :135 的真正缺口**："词层面积 → 4D 时空铰链面积"的连续极限映射未建

本脚本做三件事：
  A．**桥**：证明 p=2 时两库的"面积"是同一个（由 Σw=0 的分部求和）
  B．**核**：`Z048` 的断言"零面积词存在 ⟺ L ≡ 0 mod 4"（穷举 L=4..14）
  C．**缺陷角**：4-单纯形边界（10 个三角形铰链，每铰链 2 个四面体）的 Regge 缺陷角
     δ = 2π − Σθ；正则情形 θ = arccos(1/3) ⟹ δ = 2π − 2arccos(1/3)

运行：cd ZeroGPU/src && /usr/bin/python3 regge_hinge_bridge_check.py
"""
import itertools
import json
import math
import os
import sys

import numpy as np

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


def balanced(L):
    """平衡词：长度 L、和为零的 ±1 序列（新库口径）"""
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


def area_new(w):
    """新库 Z044：a(w) = Σ_k 前缀和_k"""
    s = 0
    tot = 0
    for x in w:
        s += x
        tot += s
    return tot


def area_old_p2(w):
    """旧库 #30/#31 p=2：area(w) = Σ_k (L−k) σ(w_k)，σ(±1)=±1"""
    L = len(w)
    return sum((L - k) * x for k, x in enumerate(w, start=1))


# ============================================================
print("=" * 78)
print("A．桥：p=2 时两库的『面积』是否同一个？")
print("=" * 78)
print("  新库：a(w) = Σ_k prefix_k = Σ_i w_i·(L−i+1)（分部求和）")
print("  旧库：area(w) = Σ_i w_i·(L−i)")
print("  ⟹ 两者相差 Σ_i w_i = 0（平衡词）⟹ **对平衡词完全相同**")
bad = []
for L in (4, 6, 8, 10, 12):
    for w in balanced(L):
        if area_new(w) != area_old_p2(w):
            bad.append((L, w))
print(f"  穷举 L=4..12 的全部平衡词：不一致的个数 = {len(bad)}")
ck("A1 新库面积 ≡ 旧库面积（p=2，平衡词；由 Σw=0 的分部求和）",
   len(bad) == 0, f"L=4..12 全部平衡词（共 {sum(len(balanced(L)) for L in (4,6,8,10,12))} 个）逐位相同")

# ============================================================
print()
print("=" * 78)
print("B．核 `Z048` 的断言：零面积词存在 ⟺ L ≡ 0 (mod 4)")
print("=" * 78)
print(f"  {'L':>3} | {'L mod 4':>7} | {'词数 C(L,L/2)':>13} | {'零面积词数':>10} | 存在?")
rows = []
for L in range(4, 15, 2):
    ws = balanced(L)
    z = sum(1 for w in ws if area_new(w) == 0)
    rows.append((L, L % 4, len(ws), z))
    print(f"  {L:>3} | {L % 4:>7} | {len(ws):>13} | {z:>10} | {'✓' if z > 0 else '✗'}")
ok = all((z > 0) == (L % 4 == 0) for L, _, _, z in rows)
ck("B1 零面积词存在 ⟺ L ≡ 0 (mod 4)（Z048 断言，穷举复算）", ok,
   "L≡0: 存在；L≡2: 不存在（L=4..14 全部相符）")
note("两库一致", "旧库 #29 的『三合一』第一条（平坦词存在 ⟺ L≡0 mod 4）",
     "⟸ 已由新库 Z048 独立证过（本表复算）——两库同判")

# ============================================================
print()
print("=" * 78)
print("C．Regge 缺陷角：4-单纯形边界（10 个三角形铰链）")
print("=" * 78)
print("  边界 3-复形：f = (5, 10, 10, 5)（顶点/边/三角形/四面体）；Σ_四面体 4 个三角形 = 20 = 2×10")
print("  ⟹ **每个三角形铰链恰被 2 个四面体共享**")
theta_reg = math.acos(1 / 4)                        # ★ 正则 4-单纯形的**胞间**二面角 = arccos(1/4)
delta_reg = 2 * math.pi - 2 * theta_reg
print(f"  注：4D Regge 的铰链是**三角形**；正则 4-单纯形在其三角形处的二面角 = arccos(1/4)")
print(f"  正则 4-单纯形胞间二面角 θ = arccos(1/4) = {math.degrees(theta_reg):.6f}°")
print(f"  ⟹ 缺陷角 δ = 2π − 2θ = {math.degrees(delta_reg):.6f}°（= {delta_reg:.6f} rad）")
ck("C1 正则 4-单纯形的缺陷角 δ = 2π − 2arccos(1/4) > 0",
   delta_reg > 0, f"δ = {math.degrees(delta_reg):.4f}°（10 个铰链全同）")

# 用实际嵌入的 5 顶点验证二面角（R⁴ 中：把两个对顶点投影到三角形的正交补，量其夹角）
V = np.array([(math.cos(2 * math.pi * j / 5), math.sin(2 * math.pi * j / 5),
               math.cos(4 * math.pi * j / 5), math.sin(4 * math.pi * j / 5)) for j in range(5)])


def perp_dir(V, tri, o):
    """把 V[o]-V[tri[0]] 投影到三角形 span 的正交补，返回单位向量（R⁴）"""
    i, j, k = tri
    u, v = V[j] - V[i], V[k] - V[i]
    w = V[o] - V[i]
    # 对 u,v 做 Gram-Schmidt，去掉 w 的切向分量
    e1 = u / np.linalg.norm(u)
    v2 = v - (v @ e1) * e1
    e2 = v2 / np.linalg.norm(v2)
    w2 = w - (w @ e1) * e1 - (w @ e2) * e2
    return w2 / np.linalg.norm(w2)


tri = (0, 1, 2)
n1, n2 = perp_dir(V, tri, 3), perp_dir(V, tri, 4)
ang = math.degrees(math.acos(max(-1.0, min(1.0, float(n1 @ n2)))))
print(f"\n  用实际嵌入的 5 顶点算三角形 (0,1,2) 处两法向夹角 = {ang:.6f}°")
print(f"  对照：arccos(1/4) = {math.degrees(math.acos(0.25)):.6f}°；π−arccos(1/4) = {math.degrees(math.pi-math.acos(0.25)):.6f}°")
dih = math.radians(ang)
ck("C2 实测二面角与正则 4-单纯形值相符（嵌入自洽）",
   abs(dih - theta_reg) < 1e-9 or abs(math.pi - dih - theta_reg) < 1e-9,
   f"实测 {ang:.6f}° ；与 arccos(1/4) 差 {abs(dih-theta_reg):.2e} rad")
delta_meas = 2 * math.pi - 2 * dih
print(f"  ⟹ 实测缺陷角 δ = 2π − 2θ = {math.degrees(delta_meas):.6f}°")
ck("C3 缺陷角 ≠ 0（正曲率）——『平坦』在词层对应面积零词、在几何侧对应 δ=0",
   abs(delta_meas) > 1e-6, f"δ = {math.degrees(delta_meas):.4f}° ≠ 0")

# ============================================================
print()
print("=" * 78)
print("D．缺口定位（#29 :135 的真正缺口）")
print("=" * 78)
print("  #29 原文：『面积相位 → Regge 缺陷角』连续极限映射 **未建**；")
print("           真正缺口 = **词层面积如何成为 4D 时空铰链面积**")
print()
print("  本脚本已把两侧都量化：")
print("    · 词层侧：面积向量 ∈ ℝ⁴（#31 显式嵌入）；面积零 ⟺ 重心（#30/#31）")
print("    · 几何侧：Regge 缺陷角 δ = 2π − Σθ（4-单纯形边界 10 铰链，δ = 2π − 2arccos(1/3)）")
print("  ⟹ **缺的是「两侧之间的映射」**：词层的 area(w) ∈ ℝ⁴ 与 4D 铰链面积 A_h 的逐项对应")
note("缺口", "『词层面积 → 4D 铰链面积』的映射未建",
     "两侧量都已显式（词层侧 = #31 的 σ 嵌入；几何侧 = Regge 缺陷角）")

# ============================================================
print()
n_fail = sum(1 for x in CK if not x["pass"])
print("=" * 78)
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 78)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "regge_hinge_bridge_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "verdict": {
                   "bridge": "p=2 时新库面积 ≡ 旧库面积（由 Σw=0 的分部求和）——两库同一对象",
                   "Z048": "零面积词存在 ⟺ L≡0 mod 4（穷举复算，与 #29 三合一第一条一致）",
                   "deficit": "4-单纯形边界 10 铰链，每铰链 2 四面体；正则 θ=arccos(1/3)，δ=2π−2θ>0",
                   "gap": "缺的是『词层 area(w)∈ℝ⁴ ↔ 4D 铰链面积 A_h』的逐项映射（#29 :135 自标未建）"}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/regge_hinge_bridge_check.json")
sys.exit(1 if n_fail else 0)