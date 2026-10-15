# -*- coding: utf-8 -*-
"""
ZeroGPU · 超荷 Y 内生来源推导

目标：5 个标准模型超荷值（轻子 -1/-2、夸克 +1/3,+4/3,-2/3）中，多少是真正
独立的？其余能否从 Zero 结构（手性 T₃、色重数 N_c=3）内生导出？

候选结构关系（待验证）：
  ① Y_R = Y_L + 2T₃^L   （右手超荷 = 左手超荷 + 弱同位旋，手性内生）
  ② Y_L(夸克) = -Y_L(轻子)/N_c   （夸克超荷 = 轻子超荷/色重数，色内生）
  ③ Q = T₃ + Y/2        （Gell-Mann-Nishijima 自洽）

若①②③全成立 ⟹ Y 表 5 值中仅 1 个独立（轻子 Y_L=-1，即 U(1) 归一化）。

用法：/usr/local/bin/python3 k4_hypercharge.py
"""
from __future__ import annotations

Nc = 3  # K₄ 色重数（内生）

# ── 假设：唯一独立输入 = 轻子左手超荷（U(1) 归一化，Z086）──
Y_L_lepton = -1.0

# ── 关系②：夸克左手超荷 ──
Y_L_quark = -Y_L_lepton / Nc          # = +1/3

# ── 弱同位旋（手性内生）：
# 左手双重态 (ν,e)L: T3=+1/2,-1/2； (u,d)L: +1/2,-1/2
T3 = {"nuL": +0.5, "eL": -0.5, "uL": +0.5, "dL": -0.5}

# ── 关系①：右手超荷 ──
def Y_R(Y_L, t3): return Y_L + 2 * t3

print("=" * 80)
print("超荷 Y 内生来源：仅 1 个独立输入（轻子 Y_L=-1）")
print("=" * 80)
print(f"独立输入：Y_L(轻子) = {Y_L_lepton}（U(1) 归一化，Z086）")
print(f"关系②：Y_L(夸克) = -Y_L(轻子)/N_c = -({Y_L_lepton})/{Nc} = {Y_L_quark}")
print()
print(f"{'粒子':<8}{'T₃(手性)':<12}{'Y_L':<10}{'Y_R(导出)':<14}{'Q_左':<10}{'Q_右':<10}{'PDG'}")
print("-" * 80)
pdg = {"e":-1.0,"nu":0.0,"u":+2.0/3,"d":-1.0/3}
allok = True
for name, t3, yl in [("nu", +0.5, Y_L_lepton), ("e", -0.5, Y_L_lepton),
                     ("u", +0.5, Y_L_quark),  ("d", -0.5, Y_L_quark)]:
    yr = Y_R(yl, t3)
    qL = t3 + yl/2
    qR = 0.0 + yr/2        # 右手 T3=0
    ok = abs(qL - pdg[name]) < 1e-9 and abs(qR - pdg[name]) < 1e-9
    allok &= ok
    print(f"{name:<8}{t3:<12}{yl:<10}{yr:<14}{qL:<10}{qR:<10}{pdg[name]}{' ✓' if ok else ' ✗'}")

print("-" * 80)
print(f"右手超荷（{['-2,-2','-2','+4/3','-2/3']}）全部由 Y_R=Y_L+2T₃ 导出：{'✓' if allok else '✗'}")
print()
print("=" * 80)
print("结论：")
print("  ① Y_R = Y_L + 2T₃ —— 右手超荷由左手+同位旋内生（不独立）")
print("  ② Y_L(夸克) = -Y_L(轻子)/N_c —— 夸克超荷由轻子+色重数内生")
print("  ③ Q = T₃ + Y/2 自洽 ⟹ 电荷全锚定")
print("  ⟹ 超荷 5 值中仅 <轻子 Y_L=-1> 1 个独立（U(1) 归一化）")
print("  ⟹ Zero 内生：手性(T₃)、色重数(3)、1/3 量子化、结构关系①②")
print("  ⟹ 残留唯一输入：轻子超荷单位 -1（Z086 U(1)，需 GUT 归一化才完全内生）")
print("=" * 80)
