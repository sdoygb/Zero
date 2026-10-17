#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
curvature_from_defects_check.py —— 《曲率从哪来》独立复算（规矩 3）
==================================================================
复算对象：旧库 `Zero-261010/unify/_src_db/` 的两篇 Zero-DB 攻线文档
    #30 `连续极限映射_词层单纯形复形与Regge铰链面积.md`
    #31 `显式等距嵌入_Zeta5到R4_正4单纯形与词层面积向量.md`
断言（待独立复算，不引用其数字）：
    QQ：面积零词 = 重心词（五态 L=10 恰 120 = 5!·1；三态 L=12 给 192 = 6·T₄）
    RR：面积向量范数 = 单纯形格点范数（三态 = 六角格范数 x²−xy+y²）
    嵌入：σ(ζ^j) = (cos 2πj/5, sin 2πj/5, cos 4πj/5, sin 4πj/5) ∈ ℝ⁴，5 顶点两两距离² = 5 恒等
    L=5：120 词面积向量范数² 恒等于 25
    连续极限"谱密度→0.5"：密度 = (Λ+1)/(2Λ) = (n²+1)/(2n²)  ← 本脚本判定它是否为重言式

运行：cd ZeroGPU/src && /usr/bin/python3 curvature_from_defects_check.py
"""
import itertools
import json
import os
import sys
from math import cos, sin, pi, factorial

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


def sigma5(j):
    """ℤ[ζ₅] → ℝ⁴：σ(ζ^j) = (ζ^j, ζ^{2j}) 的实虚部"""
    t = 2 * pi * j / 5
    return np.array([cos(t), sin(t), cos(2 * t), sin(2 * t)])


def sigma3(j):
    """ℤ[ζ₃] → ℝ²：σ(ζ^j) 的实虚部（ℤ[ζ₃] 的六角格）"""
    t = 2 * pi * j / 3
    return np.array([cos(t), sin(t)])


def area_of(word, sig):
    """area(w) = Σ_k (L−k) σ(w_k)，k 从 1 起（L = len(word)）"""
    L = len(word)
    out = np.zeros_like(sig(0))
    for k, s in enumerate(word, start=1):
        out = out + (L - k) * sig(s)
    return out


# ============================================================
print("=" * 78)
print("A．嵌入（篇 #31）：σ(ζ^j) ∈ ℝ⁴，5 顶点两两距离² 与重心")
print("=" * 78)
V = np.array([sigma5(j) for j in range(5)])
print(f"  顶点 |σ(ζ^j)|² = {[round(float(v@v), 10) for v in V]}")
d2 = [float((V[i] - V[j]) @ (V[i] - V[j])) for i in range(5) for j in range(i + 1, 5)]
print(f"  两两距离² = {[round(x, 10) for x in d2]}")
ck("A1 5 顶点两两距离² 恒等于 5（正 4-单纯形，边长 √5）",
   all(abs(x - 5) < 1e-12 for x in d2), f"min={min(d2):.12f} max={max(d2):.12f}")
s = V.sum(axis=0)
ck("A2 重心 Σσ(ζ^j) = 0", float(np.abs(s).max()) < 1e-15, f"|Σ| = {float(np.abs(s).max()):.2e}")

# 仿射无关性（面积=0 ⟹ 系数全等的前提）
M = np.vstack([V.T, np.ones(5)])
ck("A3 5 顶点在 ℝ⁴ 中仿射无关（rank = 5）⟹ 面积=0 ⟺ 系数全等",
   np.linalg.matrix_rank(M) == 5, f"rank = {np.linalg.matrix_rank(M)}")

# ============================================================
print()
print("=" * 78)
print("B．L=5（五态，每符号恰一次）：120 词，面积范数² 恒等 25？")
print("=" * 78)
norms = []
for w in itertools.permutations(range(5)):
    a = area_of(w, sigma5)
    norms.append(float(a @ a))
print(f"  词数 = {len(norms)}（= 5! = {factorial(5)}）")
print(f"  范数² 取值集合 = {sorted(set(round(x, 9) for x in norms))}")
ck("B1 L=5 的 120 词面积向量范数² 恒等于 25（全部落半径 5 球面）",
   len(set(round(x, 9) for x in norms)) == 1 and abs(norms[0] - 25) < 1e-9,
   f"唯一值 = {round(norms[0], 9)}")

# ============================================================
print()
print("=" * 78)
print("C．L=10（五态，每符号恰两次）：总词数 / 面积零词数")
print("=" * 78)
words = []
for pos in itertools.combinations(range(10), 2):          # 符号 0 的两个位置
    rest = [p for p in range(10) if p not in pos]
    for pos1 in itertools.combinations(rest, 2):          # 符号 1
        r2 = [p for p in rest if p not in pos1]
        for pos2 in itertools.combinations(r2, 2):        # 符号 2
            r3 = [p for p in r2 if p not in pos2]
            for pos3 in itertools.combinations(r3, 2):    # 符号 3
                r4 = [p for p in r3 if p not in pos3]
                occ = [pos, pos1, pos2, pos3, tuple(r4)]
                w = [0] * 10
                for sym, ps in enumerate(occ):
                    w[ps[0]] = sym
                    w[ps[1]] = sym
                words.append(tuple(w))
print(f"  总词数 = {len(words)}（理论 10!/(2!⁵) = {factorial(10)//32}）")
ck("C1 总词数 = 10!/(2!⁵) = 113400",
   len(words) == factorial(10) // 32, f"{len(words)}")

zero = [w for w in words if float(np.abs(area_of(w, sigma5)).max()) < 1e-12]
print(f"  面积零词数 = {len(zero)}")
ck("C2 面积零词恰 120 个（= 5!·1）", len(zero) == 120, f"{len(zero)}")
note("缺陷词数", len(words) - len(zero), "＝总词数 − 重心词数（篇 #30 报 113,280）")

# 显式刻画：每个符号的两次出现位置之和 = 11？
# ★ 修正 1：位置用 1 基（0 基时应为 9）
pred = [w for w in words
        if all(sum(p + 1 for p, s in enumerate(w) if s == sym) == 11 for sym in range(5))]
print(f"  按『每符号位置和 = 11（关于中点对称）』预测的个数 = {len(pred)}")
ck("C3 面积零 ⟺ 每符号两次出现位置互补（p₁+p₂ = L+1 = 11）",
   len(pred) == 120 and set(zero) == set(pred),
   f"预测 {len(pred)} 个；与实测集合{'相同' if set(zero)==set(pred) else '不同'}")

# ============================================================
print()
print("=" * 78)
print("D．三态（p=3，ℤ[ζ₃] → ℝ²）：面积向量是否落六角格（范数 x²−xy+y²）？")
print("=" * 78)
print("  注：坐标取 σ(ζ^j) = (cos 2πj/3, sin 2πj/3)；改用六角格基 (1,0),(−1/2,√3/2) 时")
print("      面积向量坐标应为整数 ⟹ 范数 = x²−xy+y²")
B = np.array([[1.0, 0.0], [-0.5, np.sqrt(3) / 2]])          # 六角格基
Binv = np.linalg.inv(B)
words3 = [w for w in itertools.product(range(3), repeat=6) if all(w.count(s) == 2 for s in range(3))]
hexnorms, ok_int = set(), True
for w in words3:
    a = area_of(w, sigma3)
    xy = np.linalg.solve(B.T, a)     # ★ 修正 2：基向量按行放 ⟹ 解 Bᵀc = a（原写 Binv@a）
    if max(abs(xy[0] - round(xy[0])), abs(xy[1] - round(xy[1]))) > 1e-9:
        ok_int = False
    m, n_ = round(xy[0]), round(xy[1])
    hexnorms.add(m * m - m * n_ + n_ * n_)
print(f"  三态 L=6 词数 = {len(words3)}（= 6!/(2!³) = {factorial(6)//8}）")
print(f"  面积向量在六角格基下的坐标全为整数？ {ok_int}")
print(f"  六角格范数 x²−xy+y² 的取值（前 12 个）= {sorted(hexnorms)[:12]}")
ck("D1 三态面积向量落六角格（坐标为整数 ⟹ 范数 = x²−xy+y²，∈{0,1,3,4,7,9,…}）",
   ok_int and all(v >= 0 for v in hexnorms), f"取值 {sorted(hexnorms)[:8]}")

# ============================================================
print()
print("=" * 78)
print("E．连续极限『谱密度→0.5』是什么？（判定是否重言）")
print("=" * 78)
print("  篇 #30 的表：谱点数 = Λ+1 = n²+1；范围 = 2Λ = 2n²；密度 = 点数/范围")
for n in (4, 6, 8, 10):
    a = (n * n + 1) / (2 * n * n)
    print(f"    n = {n:>2}: 点数 = {n*n+1:>3}，范围 = {2*n*n:>3}，密度 = {a:.4f}")
print()
print("  ⟹ 这是**算术恒等式**：任何 Λ+1 个点落在长度 2Λ 的区间里，密度都 → 1/2。")
print("     它**不依赖词层结构**（换任何『点数=Λ+1、范围=2Λ』的谱都成立）")
note("判定", "『谱密度→0.5』= 计数恒等式 (Λ+1)/(2Λ)→1/2；不构成连续极限的证据",
     "篇 #30 自己把 Donsker 型严格收敛标为【未证】（§6）——本判定与其自标一致")

# ============================================================
print()
n_fail = sum(1 for x in CK if not x["pass"])
print("=" * 78)
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 78)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "curvature_from_defects_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "verdict": {
                   "A_embedding": "σ(ζ^j)∈ℝ⁴ 的 5 顶点两两距离²=5、重心为零、仿射无关（rank 5）",
                   "B_L5": "L=5 的 120 词面积范数² 恒为 25",
                   "C_L10": "总词数 113400；面积零词恰 120；且 ⟺ 每符号位置和 = L+1（关于中点对称）",
                   "D_p3": "三态面积向量落六角格，范数 = x²−xy+y²",
                   "E_density": "『谱密度→0.5』= 计数恒等式 (Λ+1)/(2Λ)→1/2，不是连续极限的证据"}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/curvature_from_defects_check.json")
sys.exit(1 if n_fail else 0)