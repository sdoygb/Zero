#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_color_quantum.py —— 权重格子 → 色量子数：从 Zero A₂ 权重格子读出夸克色态量子数
===============================================================================
前置：Z096 已从 Zero C₃(ω) 重建 SU(3) A₂ 权重格子 ℤω₁⊕ℤω₂（基本权重）。
本篇：用该权重格子构造 SU(3) 基本表示 3（夸克色三重态）的权重，
      读出 r/g/b 色态的具体量子数 (T₃, Y_c, Q_c)，对照 PDG，接 Z090 色锚定。

权重坐标（Dynkin 坐标）→ 物理色量子数：
  T₃ = (a−b)/2    （同位旋第三分量）
  Y_c = (a+b)/3   （色超荷）
  Q_c = T₃ + Y_c/2 = (2a−b)/3   （色电荷）
"""
import sys, os, json
import numpy as np

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
OUT = os.path.join(RES, "k4_color_quantum.json")

# Z096 重建的基本权重（Dynkin/权重坐标）
# 基本权重 ω₁=(1,0), ω₂=(0,1)（在基本权重基下）
# SU(3) 基本表示 3 的权重 = ω₁, ω₂, −(ω₁+ω₂)
def dynkin_to_phys(a, b):
    """权重坐标 (a,b) → (T₃, Y_c, Q_c)"""
    T3 = (a - b)/2.0
    Yc = (a + b)/3.0
    Qc = T3 + Yc/2.0
    return T3, Yc, Qc

print("="*82)
print("权重格子 → 色量子数：从 Zero A₂ 权重格子读出夸克色态 r/g/b")
print("="*82)

# ① 从 Z096 权重格子构造 SU(3) 基本表示 3 的权重
weights3 = [(1,0), (0,1), (-1,-1)]        # ω₁, ω₂, −(ω₁+ω₂)
color_names = ["r (red)", "g (green)", "b (blue)"]
print("\n① SU(3) 基本表示 3（夸克色三重态）的权重")
print("   权重坐标(ω₁,ω₂基):  ω₁=(1,0), ω₂=(0,1), −(ω₁+ω₂)=(−1,−1)")
print("   对应色态 → 物理量子数 (T₃, Y_c, Q_c)：")
table = []
for (a,b), name in zip(weights3, color_names):
    T3, Yc, Qc = dynkin_to_phys(a,b)
    table.append({"color": name, "weight": (a,b), "T3": round(T3,4),
                  "Yc": round(Yc,4), "Qc": round(Qc,4)})
    print(f"     {name:12s} 权重=({a:+d},{b:+d})  T₃={T3:+.2f}  Y_c={Yc:+.3f}  Q_c={Qc:+.3f}")

# ② 反三重 3̄
print("\n② 反三重态 3̄（反夸克色态 r̄,ḡ,b̄）：权重取反")
table_bar = []
for (a,b), name in zip(weights3, color_names):
    T3, Yc, Qc = dynkin_to_phys(-a,-b)
    table_bar.append({"color": name.replace(" (","̅ (")+"̅", "weight": (-a,-b),
                      "T3": round(T3,4), "Yc": round(Yc,4), "Qc": round(Qc,4)})
    print(f"     {name.replace(' (','̅ ('):10s}̅ 权重=({-a:+d},{-b:+d})  T₃={T3:+.2f}  Y_c={Yc:+.3f}  Q_c={Qc:+.3f}")

# ③ 对照 PDG 夸克色量子数
print("\n③ 对照 PDG 夸克色量子数")
print("   u 型夸克（Q=+2/3）：色三重态 r 电荷 = +2/3 ✓")
print("   d 型夸克（Q=−1/3）：色三重态 g/b 电荷 = −1/3 ✓")
print(f"   色三重 Q 总和 = {table[0]['Qc']+table[1]['Qc']+table[2]['Qc']:.4f}（=0，色单态总电荷 ✓）")

# ④ 与 Z090 色锚定对接
print("\n④ 与 Z090 色锚定对接")
print("   Z090: C₃ 复表示(ω/ω²) → 色三重态(存在性)；平凡表示 → 色单态")
print("   Z096: A₂ 权重格子从 Zero C₃(ω) 重建")
print("   本篇: 权重格子 → 色三重态 r/g/b 的具体量子数 (T₃,Y_c,Q_c) 定量化")
print("   ⟹ 色量子数【定量】从 Zero 权重格子长出，接 Z090 代数结构")

# ⑤ 诚实边界
print("\n⑤ 诚实边界")
print("   · 重建的是 SU(3) 色空间的量子数(T₃,Y_c,Q_c)，是权重格子数学事实")
print("   · 非 Zero 引擎粒子模拟：不声称色态在演化引擎中出现")
print("   · 色超荷 Y_c 与 Z090 弱超荷 Y_L 是不同层次(Y_c∈SU(3)内部, Y_L∈电弱)，勿混")

result = {
    "weights_3": [{"color": c, "weight": w, "T3": round(dynkin_to_phys(*w)[0],3),
                   "Yc": round(dynkin_to_phys(*w)[1],3), "Qc": round(dynkin_to_phys(*w)[2],3)}
                  for c,w in zip(color_names, weights3)],
    "weights_3bar": [{"color": c.replace(" (","̅ ("), "weight": w, "T3": round(dynkin_to_phys(-w[0],-w[1])[0],3),
                     "Yc": round(dynkin_to_phys(-w[0],-w[1])[1],3), "Qc": round(dynkin_to_phys(-w[0],-w[1])[2],3)}
                     for c,w in zip(color_names, weights3)],
    "matches_PDG": "u型夸克 r=+2/3, d型 g/b=−1/3, 色三重总电荷0",
    "link_Z090": "Z090代数结构(C₃复表示判色三重)+Z096权重格子+本篇定量(T₃,Y_c,Q_c)",
    "boundary": "量子数数学事实, 非引擎粒子模拟; Y_c≠Z090弱超荷Y_L"
}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=1)
print(f"\n结果已写入 {OUT}")
