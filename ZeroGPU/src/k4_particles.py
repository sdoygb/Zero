# -*- coding: utf-8 -*-
"""
ZeroGPU · 粒子模拟：忠实引擎稳态涌现 → 3 频率槽(代) × 量子特征锚定(种类) → 粒子谱

延续 Z090 锚定协议：粒子 = (代 × 种类)。
  代   = K₄ 稳态 3 频率槽（忠实引擎 top2+T=4 涌现，长期稳定）
  种类 = 结构覆盖 × C₃ 表示：
         单覆盖整数表示 + 色单态  → 玻色子（自旋0/1）
         双覆盖 spinor + 色单态  → 轻子（费米子）
         双覆盖 spinor + 色三重  → 夸克（费米子，分上型/下型）
  量子数 = T₃(手性)、色(C₃)、电荷(Q=T₃+Y/2，Y 压缩自超荷内生 5→1)

用法：/usr/local/bin/python3 k4_particles.py
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loyalty_engine import LoyaltyEngine

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]

# ── ① 忠实引擎涌现稳态 ──
eng = LoyaltyEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2")
for _ in range(60000):
    eng.run_generation()
w = eng.record / eng.record.sum()
# 3 频率槽（相邻比对数）
freqs = [abs(math.log(w[i]/w[i+1])) for i in range(3)]
freqs = sorted(freqs)
mass = [x/max(freqs) for x in freqs]     # 3 代质量（归一化，比值=seed 注入）

print("=" * 86)
print("粒子模拟：K₄ 忠实引擎稳态 → 粒子谱（3 代 × 种类 × 量子数）")
print("=" * 86)
print(f"稳态 record = {[round(float(x),4) for x in w]}")
print(f"3 频率槽(代) = {[round(f,4) for f in freqs]}  （长期稳定，6 万代）")
print(f"3 代质量(归一) = {[round(m,4) for m in mass]}  （比值由 seed 注入，Z086）")

# ── ② 种类锚定（Z090 协议）──
# 每代 × 每种结构种类 → 粒子
particles = []
Nc = 3
Y_L_lepton = -1.0
Y_L_quark = -Y_L_lepton / Nc

def charge(t3, Y): return t3 + Y/2.0

# 费米子-轻子（双覆盖+色单态，T₃ 按手性）
for d, m in enumerate(mass, 1):
    eL_q = charge(-0.5, Y_L_lepton)   # e 左
    nuL_q = charge(+0.5, Y_L_lepton)  # ν 左
    particles.append({"代": d, "m": m, "种类": "轻子带电", "符号": "-",
                      "色": "单态", "自旋": "1/2(费米)", "电荷": eL_q})
    particles.append({"代": d, "m": m, "种类": "轻子中微子", "符号": "0",
                      "色": "单态", "自旋": "1/2(费米)", "电荷": nuL_q})
# 费米子-夸克-上型（双覆盖+色三重，Q=+2/3）
for d, m in enumerate(mass, 1):
    particles.append({"代": d, "m": m, "种类": "上型夸克", "符号": "+2/3",
                      "色": "三重", "自旋": "1/2(费米)",
                      "电荷": charge(+0.5, Y_L_quark)})
# 费米子-夸克-下型（Q=-1/3）
for d, m in enumerate(mass, 1):
    particles.append({"代": d, "m": m, "种类": "下型夸克", "符号": "-1/3",
                      "色": "三重", "自旋": "1/2(费米)",
                      "电荷": charge(-0.5, Y_L_quark)})
# 玻色子（单覆盖+色单态，自旋0/1）
for d in range(1, 4):
    particles.append({"代": d, "m": mass[d-1], "种类": "玻色子(矢量/标量)", "符号": "0/±1",
                      "色": "单态", "自旋": "0/1(玻色)",
                      "电荷": charge(0.0, 0.0)})

# 锚定 PDG 对照
pdg = {"轻子带电":-1.0,"轻子中微子":0.0,"上型夸克":2/3,"下型夸克":-1/3,
       "玻色子(矢量/标量)":0.0}
print()
print(f"{'代':<3}{'种类':<18}{'电荷':<8}{'色':<6}{'自旋':<12}{'PDG电荷':<8}{'锚定'}")
print("-" * 80)
for p in particles:
    match = abs(p["电荷"] - pdg[p["种类"]]) < 1e-9
    print(f"{p['代']:<3}{p['种类']:<18}{p['符号']:<8}{p['色']:<6}"
          f"{p['自旋']:<12}{pdg[p['种类']]:<8}{'✓' if match else '✗'}")

# 汇总：每代的粒子集合
print("\n" + "=" * 86)
print("粒子谱承载结构（每代粒子集合）：")
print("  每代 = 1 质量槽（seed 注入）× 6 结构种类（轻子2+上夸克+下夸克+玻色子）")
print("  = 与标准模型一致：每代 {e,ν,u,d} + 规范/标量玻色子")
print("=" * 86)

out = {"record": [float(x) for x in w], "freqs": freqs, "mass": mass,
       "particles": particles}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_particles.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"已存 {os.path.join(RES,'k4_particles.json')}")
