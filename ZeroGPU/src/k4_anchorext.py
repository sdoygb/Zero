#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_anchorext.py —— 锚定扩展：把 ν_R 纳入粒子谱（每代 8 场）
================================================================
Z091 锚定把"轻子中微子"作一类(电荷0)，未区分 ν_L/ν_R。
本篇扩展锚定协议：按【手性×色】展开每代费米子为 8 场，读出手性/超荷，
把右手中微子 ν_R（T₃=0, Y=0, Q=0）纳入粒子谱。

超荷体系（Z090）：
  轻子多重态 Y_L=−1（输入单位）; 夸克多重态 Y_L=+1/3（N_c=3 导出）
  Y_R = Y_L + 2·T₃(多重态)（内生关系）
  电荷 Q = T₃ + Y/2（右手 T₃=0）
每代费米子 8 场：ν_L, e_L, e_R, ν_R, u_L, u_R, d_L, d_R（含 ★新粒子 ν_R）
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loyalty_engine import LoyaltyEngine

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]

# ── ① 忠实引擎涌现稳态（Z091）──
eng = LoyaltyEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2")
for _ in range(60000):
    eng.run_generation()
w = eng.record / eng.record.sum()
freqs = sorted([abs(math.log(w[i]/w[i+1])) for i in range(3)])
mass = [x/max(freqs) for x in freqs]

print("="*86)
print("锚定扩展：把 ν_R 纳入粒子谱（每代 8 场，含右手惰性中微子）")
print("="*86)
print(f"稳态 record = {[round(float(x),4) for x in w]}")
print(f"3 频率槽(代) = {[round(f,4) for f in freqs]}")

# ── ② 超荷体系 + 手性展开（Z090 + 锚定扩展）──
Nc = 3
Y_L_lepton = -1.0
Y_L_quark = -Y_L_lepton / Nc   # +1/3（N_c=3 导出）

def Y_R(Y_L, T3_mult): return Y_L + 2*T3_mult   # 右手超荷（多重态 T₃）
def charge(t3, Y): return t3 + Y/2.0

# 轻子（色单态，费米）：多重态 (ν:+1/2, e:−1/2)，右手 T₃=0
leptons = [
    ("ν_L", +0.5, Y_L_lepton),           # 左手：T₃=+1/2, Y=−1, Q=0
    ("e_L", -0.5, Y_L_lepton),           # 左手：T₃=−1/2, Y=−1, Q=−1
    ("e_R", -0.5, Y_R(Y_L_lepton,-0.5)), # 右手：Y_R=−2, Q=−1
    ("ν_R", +0.5, Y_R(Y_L_lepton,+0.5)), # 右手：Y_R=0, Q=0 ★新粒子
]
# 夸克（色三重，费米）：多重态 (u:+1/2, d:−1/2)，右手 T₃=0
quarks = [
    ("u_L", +0.5, Y_L_quark),
    ("d_L", -0.5, Y_L_quark),
    ("u_R", +0.5, Y_R(Y_L_quark,+0.5)),
    ("d_R", -0.5, Y_R(Y_L_quark,-0.5)),
]

def fermi_field(name, T3m, Y, col):
    T3 = T3m if "L" in name else 0.0      # 右手 T₃=0
    return {"name": name, "色": col, "T3": T3, "Y": Y, "Q": charge(T3, Y)}

# 每代 8 场
particles = []
for d, m in enumerate(mass, 1):
    for nm, t3m, y in leptons:
        f = fermi_field(nm, t3m, y, "单态"); f["代"]=d; f["m"]=m; particles.append(f)
    for nm, t3m, y in quarks:
        f = fermi_field(nm, t3m, y, "三重"); f["代"]=d; f["m"]=m; particles.append(f)
# 玻色子（非世代，色单态）；W± 是 SU(2) 三重态 T₃=±1, Y=0
bosons = [("γ", 0.0, 0.0, 1), ("Z", 0.0, 0.0, 1),
          ("W+", 1.0, 0.0, 1), ("W−", -1.0, 0.0, 1), ("Higgs", 0.0, 0.0, 0)]
for nm, t3, y, sp in bosons:
    q = charge(t3, y)
    particles.append({"name": nm, "色":"单态", "T3":t3, "Y":y, "Q":q,
                      "代":"—", "自旋":f"{sp}玻色"})

# ── ③ 对照 PDG + ν_R 标注 ──
sm_Q = {"ν_L":0.0,"e_L":-1.0,"e_R":-1.0,"u_L":2/3,"u_R":2/3,"d_L":-1/3,"d_R":-1/3,
        "γ":0.0,"Z":0.0,"W+":1.0,"W−":-1.0,"Higgs":0.0}   # 不含 ν_R（SM 无）
print("\n每代费米子 8 场（含 ★ν_R 新粒子）：")
print(f"{'代':<3}{'场':<6}{'色':<5}{'T₃':<7}{'Y':<8}{'Q':<7}{'SM'}")
print("-"*52)
for p in particles:
    if "代" in p and p["代"] != "—":
        is_sm = p["name"] in sm_Q and abs(p["Q"]-sm_Q[p["name"]])<1e-9
        tag = "✓" if is_sm else "★新粒子(ν_R)" if p["name"]=="ν_R" else "?"
        print(f"{p['代']:<3}{p['name']:<6}{p['色']:<5}{p['T3']:<7.1f}{p['Y']:<8.2f}{p['Q']:<7.2f}{tag}")
# 玻色子
print("\n玻色子：")
for p in particles:
    if p["代"]=="—":
        print(f"   {p['name']:<7} 色单态 Q={p['Q']:.1f} {p['自旋']}")

# ── ④ 汇总 ──
nuR_count = sum(1 for p in particles if p["name"]=="ν_R")
print("\n" + "="*86)
print("锚定扩展结果：")
print(f"  每代费米子 = 8 场（ν_L,e_L,e_R,★ν_R,u_L,u_R,d_L,d_R）×3代 = 24 场")
print(f"  ν_R 已纳入粒子谱（每代 1 个 ×3）—— SM 无, 新粒子候选")
print(f"  超荷/手性由锚定读出：ν_R(T₃=0,Y=0,Q=0) vs e_R(T₃=0,Y=−2,Q=−1) 区分")
print("  玻色子：γ,Z,W±,Higgs（色单态）；胶子色八重另载（Z095）")
print("="*86)

out = {"record": [float(x) for x in w], "freqs": freqs,
       "particles": [{"name":p["name"],"代":p["代"],"色":p["色"],"T3":p["T3"],
                       "Y":round(p["Y"],3),"Q":round(p["Q"],3),
                       "自旋":p.get("自旋","1/2费米")} for p in particles],
       "nuR_included": True, "nuR_count": nuR_count,
       "conclusion": "锚定扩展: 每代8场(含ν_R), 手性×色区分读出, ν_R纳入粒子谱"}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_anchorext.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"已存 {os.path.join(RES,'k4_anchorext.json')}")
