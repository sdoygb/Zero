#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_flavor_gauge.py —— 补最后 2 项：味量子数 + 规范玻色子种类
================================================================
按 Z094 信息层检验的剩余 2 项"部分"，尝试从 Zero 结构锚定，诚实标记能/缺。

① 味量子数（奇异数S/粲数C/底数B'/顶数T）：
   结构依据：3 代 = Z₃/S₃ 循环（Z093）。上型/下型夸克（电荷区分，Z090）。
   锚定：代1=u/d(味0)，代2=c/s(C=+1/S=-1)，代3=t/b(T=+1/B'=-1)。
   诚实：Zero 结构给"3 代循环"，具体味名/味数值是 PDG 映射（约定），
   味守恒性质（弱破缺）是动力学层，不在信息层。

② 规范玻色子种类（γ/W±/Z/胶子/Higgs）：
   Zero 3 个玻色子：自旋0/1、电荷0/±1、色单态、有质量（3 频率槽）。
   匹配：W±(Q=±1,色单,有质量)、Z(Q=0,色单,有质量)、H(自旋0,色单,有质量)。
   缺口：γ(无质量) 不符(Zero玻色子有质量)；胶子(色八重) 不符(Zero玻色子色单态)。
   诚实：3 玻色子可锚定电弱 W±/Z/H 候选；γ 与胶子是 Zero 结构缺口。
"""
import json, os

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
OUT = os.path.join(RES, "k4_flavor_gauge.json")

# ---------- ① 味量子数 ----------
# 3 代 → 具体味（上型/下型 × 代）
flavor_map = [
    # (代, 上型味, 下型味)
    (1, "u", "d"), (2, "c", "s"), (3, "t", "b"),
]
# PDG 味量子数：(S, C, B', T)  下划线 B'=底数（区别于重子数 B）
flavor_PDG = {
    "u": (0, 0, 0, 0), "d": (0, 0, 0, 0),
    "c": (0, 1, 0, 0), "s": (-1, 0, 0, 0),
    "t": (0, 0, 0, 1), "b": (0, 0, -1, 0),
}
# 上型/下型夸克质量（3 频率槽，代1,2,3）
m = {1: 0.413, 2: 0.707, 3: 1.000}

print("="*76)
print("① 味量子数锚定（3 代 Z₃ 循环 → 味标记）")
print("="*76)
flavor_rows = []
for gen, up, down in flavor_map:
    for name, fl in (("上型", up), ("下型", down)):
        S, C, Bd, T = flavor_PDG[fl]
        flavor_rows.append({"代": gen, "味": fl, "上/下": name, "质量": m[gen],
                            "S": S, "C": C, "B'": Bd, "T": T})
        print(f"  代{gen} {name}夸克 {fl:<2} m={m[gen]:.3f}  "
              f"S={S:+d} C={C:+d} B'={Bd:+d} T={T:+d}   (PDG)")

print(f"  → 3 代循环给味结构；味名/数值为 PDG 映射（约定），非内生守恒")

# ---------- ② 规范玻色子种类 ----------
print("\n" + "="*76)
print("② 规范玻色子种类锚定（3 个色单态有质量玻色子）")
print("="*76)
# PDG 规范玻色子属性：(自旋, 电荷, 色, 质量)
gauge_PDG = {
    "γ":    dict(自旋=1, Q=0.0,  色="单态", 质量="0"),
    "W±":   dict(自旋=1, Q="±1", 色="单态", 质量=">0"),
    "Z":    dict(自旋=1, Q=0.0,  色="单态", 质量=">0"),
    "胶子": dict(自旋=1, Q=0.0,  色="八重", 质量="0"),
    "H":    dict(自旋=0, Q=0.0,  色="单态", 质量=">0"),
}
# Zero 玻色子：色单态、有质量、自旋0/1、电荷0/±1
zero_boson = dict(自旋=[0,1], Q=[0.0, "±1"], 色="单态", 质量=">0")

def match(b):
    """Zero 玻色子能否承载该规范玻色子。"""
    if b["色"] != zero_boson["色"]:        # 色八重 → 缺（Zero 色单态）
        return "✗ 缺(色八重)"
    if b["质量"] == "0":                    # 无质量 γ/胶子 → 缺（Zero 有质量）
        return "✗ 缺(需无质量)"
    if b["自旋"] in zero_boson["自旋"] and b["Q"] in zero_boson["Q"]:
        return "✓ 可锚定"
    return "◐ 部分"

gauge_rows = []
for name, b in gauge_PDG.items():
    s = match(b)
    gauge_rows.append({"玻色子": name, **b, "Zero承载": s})
    print(f"  {name:<3} 自旋{b['自旋']} Q={str(b['Q']):<3} 色{b['色']} 质量{b['质量']:<3}"
          f" → Zero 承载: {s}")

print(f"  → 3 玻色子可锚定 W±/Z/H 候选；γ(无质量)、胶子(色八重) 为结构缺口")

# ---------- 汇总 ----------
print("\n" + "="*76)
print("信息层覆盖更新：10/12 → 味量子数可补 + 规范玻色子部分")
print("="*76)
covered = 10 + 1   # + 味量子数
print(f"  明确覆盖: {covered}/12  (电荷/色/T₃/Y/手性/自旋/三代/B/L/C/味)")
print(f"  部分: 规范玻色子种类 (W±/Z/H 可, γ/胶子缺)")
print(f"  结论: 信息层到 11/12 明确; 规范玻色子种类是真实结构缺口")
print(f"    - γ: Zero 玻色子有质量(3频率槽), 无质量光子需零质量槽")
print(f"    - 胶子: Zero 玻色子色单态, 色八重需额外 SU(3) 结构")

result = {
    "flavor_anchor": flavor_rows,
    "gauge_anchor": gauge_rows,
    "coverage": {"明确": 11, "部分": 1, "总数": 12},
    "gauge_gap": ["γ(需零质量)", "胶子(需色八重)"],
    "conclusion": "味量子数可补(3代→味); 规范玻色子种类部分(W±/Z/H 可, γ/胶子缺)"
}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=1)
print(f"\n结果已写入 {OUT}")
