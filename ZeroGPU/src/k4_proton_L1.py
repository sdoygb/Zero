#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_L1.py —— 从 L0-L1-L2 层间动力学(L1历史记录层=合力) 推导质子质量 0.35% 偏移
======================================================================================
用户方向: 不从力的数值算, 而从三层闭合动力学推导; L1历史记录层与"合力"关系大.
Z088 基础: 六边闭环; 三色荷共享 L1 历史层 → 记录归一化(Z088 §3.1 重数加权划分) = 结合力.
目标: 质子的有效合成从 ΣM=3(L2独立) 受 L1 记录共享修正, 内生给出 3→2.9895 (Δ=0.010532).
"""
from __future__ import annotations
import os, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

me, m_p = 0.51099895, 938.272088
theta_e = 2*np.pi/3 + 2/9
A = np.sqrt(me)/(1+np.sqrt(2)*np.cos(theta_e)); A2 = A**2
S0 = m_p/A2      # 2.989468
D  = 3 - S0      # 0.010532  (合成修正)
drel = D/3       # 0.0035107 (相对修正 0.351%)
P_("=" * 96)
P_("质子 0.35% 修正 · L1 历史记录层(合力) 机制推导 · k4_proton_L1.py")
P_("=" * 96)
P_("目标: 合成 3→%.6f (Δ=%.6f, 相对 %.6f%%)" % (S0, D, drel*100))

# ---- Z088 内生量 (L1 记录层) ----
two9 = 2/9; d0 = np.pi/12 - two9
q_hist = 5/9        # Z088 稳态集中度
w2, w1 = 2/3, 1/3   # Z088 二分支 2:1
four9 = 4/9
log2 = np.log(2)    # Z087 内禀频率
L1q = {"q_hist=5/9": q_hist, "1−q=4/9": four9, "2/3": w2, "1/3": w1,
       "2/9=(2/3)(1/3)": two9, "1/9": 1/9, "5/9−1/2=1/18": 5/9-1/2,
       "q²=25/81": q_hist**2, "(4/9)²": four9**2, "log2": log2,
       "δ₀": d0, "π/12": np.pi/12}

P_("\n" + "="*96)
P_("【A. 从 Z088 L1 内生量找 Δ=0.010532 (合成修正)】")
hits = []
for nm, v in L1q.items():
    for k in [2,3,4,5,6,8,9,12,16,24,25,27,32,36,48,64,81,100]:
        for name, val in [(f"{nm}/{k}", v/k), (f"{nm}×{k}", v*k),
                          (f"{nm}·(1/3)", v*(1/3)), (f"{nm}·(2/9)", v*two9)]:
            rel = abs(val-D)/D*100 if val>0 else 1e4
            hits.append((rel, name, val))
# 组合
for a in range(1,5):
    for b in range(1,5):
        v = a*w2*b*w1/ (10**2)
        hits.append((abs(v-D)/D*100 if v>0 else 1e4, f"{a}·(2/3)×{b}·(1/3)/100", v))
hits.sort(key=lambda x:x[0])
P_("  Top12 最接近 Δ=0.010532:")
for rel, nm, v in hits[:12]:
    P_("    %-22s = %.7f  (偏差 %.2f%%)" % (nm, v, rel))

P_("\n" + "="*96)
P_("【B. 机制候选: L1 记录归一化 → 合成修正】")
P_("  候选1: 合成 S' = ΣM × (1−R), 需 R=%.6f (%.4f%%)" % (drel, drel*100))
P_("     R 候选: q_hist 相关 1−q=%.6f; (1−q)/q=%.4f; (1/3)²=%.4f; (1/9)=%.4f" % (four9, four9/q_hist, 1/9, 1/9))
for nm, v in [("1−q_hist=4/9", four9), ("1/9", 1/9), ("(1/3)²", 1/9), ("(2/3)(1/9)", w2/9), ("(4/9)/q", four9/q_hist)]:
    rel = abs(v-drel)/drel*100
    P_("      R≈%-10s = %.6f  (偏差 %.1f%%)" % (nm, v, rel))

P_("\n" + "="*96)
P_("【C. 三色荷共享 L1: 记录重数归一化的直接机制】")
# 三色荷在 L1 的记录: 若按 Z088 稳态二分支 2:1, 三份记录可能 (2/3,1/3,?) 
# 机制: 三独立色荷 L2 合成 ΣM=3; 共享 L1 后, 记录归一化使合成按重数加权重排
# 试: 三色荷 L1 记录 ω=(2/3,1/3,0)归一化 或 (ω1,ω2,1-ω1-ω2)
# 质子的 M 合成 = 3×Σ(ω_i·M_i)/Σω? 或 3×(重数加权M平均)
for w3 in [0, 1/3, 1/6]:
    w = np.array([w2, w1, w3]); w = w/w.sum()
    P_("    L1 记录 ω=%s (归一化) → 重数加权系数 (Σω_i²=%.4f)" % (np.round(w,4), np.sum(w**2)))

out = {"delta": D, "s0": S0, "drel": drel, "best_simple": hits[0],
       "note": "从 Z088 L1 内生量(q_hist=5/9, 2:1) 构建, 未见精确内生 0.0105"}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_L1.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_L1.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_L1.json / .txt")
