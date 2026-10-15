# -*- coding: utf-8 -*-
"""
ZeroGPU · K₄ 三频率槽 → 三代质量比：Koide 与观测检验

背景：K₄（Γ 原生来源）无反馈多频率精确并存 = 3 频率槽（Z089），频率由 seed
相邻比决定（Z086 U(1) 具名输入）。本脚本检验：这 3 频率槽在合理的质量映射下，
能否再现 Koide 公式（Q=2/3）与三代观测质量比。

质量映射假设（标口径，非导出）：
  (A) m_i = f_i       质量=频率（Z087 粒子质量=频率）
  (B) m_i = r_i       质量=seed 相邻比
  (C) m_i = e^{f_i} = r_i （=B）

Koide 公式：Q = (Σm_i)/(Σ√m_i)² = 2/3（对带电轻子精确）。
Koide 族参数解：m_i ∝ (1+√2 cos(θ+2πi/3))²，对任意 θ 满足 Q=2/3。

观测（带电轻子，PDG）：m_e=0.511, m_μ=105.66, m_τ=1776.86 MeV。
三代夸克质量比也检验（u,c,t / d,s,b）。

用法：/usr/local/bin/python3 k4_koide.py
"""
from __future__ import annotations
import math

# ── Koide 公式 ──
def koide_q(m):
    s = sum(m); sr = sum(math.sqrt(x) for x in m)
    return s / (sr * sr)

# ── Koide 族：m_i = (1+√2 cos(θ+2πi/3))² ──
def koide_family(theta):
    sq2 = math.sqrt(2.0)
    return [(1 + sq2 * math.cos(theta + 2 * math.pi * i / 3)) ** 2 for i in range(3)]

print("=" * 82)
print("一、Koide 族验证：任意 θ 是否恒 Q=2/3")
print("=" * 82)
qs = [koide_q(koide_family(t)) for t in [0, 0.3, 1.0, 2.0, math.pi / 4]]
print("  θ∈{0,0.3,1,2,π/4} 的 Q =", [round(q, 6) for q in qs], "（应全=2/3）")

# ── 观测 ──
print("\n" + "=" * 82)
print("二、观测三代质量比（PDG, MeV）")
print("=" * 82)
obs = {
    "轻子(e,μ,τ)":  (0.51099895, 105.6583755, 1776.86),
    "上型夸克(u,c,t)": (2.16, 1.27e3, 1.73e5),
    "下型夸克(d,s,b)": (4.67, 93.0, 4.18e3),
}
for name, m in obs.items():
    print(f"  {name:<14} 质量={m}  Q={koide_q(m):.6f}  "
          f"（Koide 目标 0.666667）")

# ── K₄ 三频率槽：几种 seed 的频率 / 相邻比 ──
print("\n" + "=" * 82)
print("三、K₄ 三频率槽 → 质量映射 → Koide Q")
print("=" * 82)
k4_seeds = {
    "trip [8,4,1.5,1]": (8.0, 4.0, 1.5, 1.0),
    "quad [16,8,3,1]":  (16.0, 8.0, 3.0, 1.0),
    "geom [8,4,2,1]":   (8.0, 4.0, 2.0, 1.0),
    "slow [9,5,2,1]":   (9.0, 5.0, 2.0, 1.0),
    "rand [10,4,1.7,1]":(10.0, 4.0, 1.7, 1.0),
}
for name, (s0, s1, s2, s3) in k4_seeds.items():
    # 频率槽 = 3 相邻比的对数（单调 seed 假定）
    r = (s0 / s1, s1 / s2, s2 / s3)
    f = tuple(math.log(x) for x in r)
    qA = koide_q(f)                      # 假设 A：m=f
    qB = koide_q(r)                      # 假设 B：m=r（=e^f）
    print(f"  [{name:<14}] 相邻比={[round(x,3) for x in r]} "
          f"频率={[round(x,3) for x in f]}")
    print(f"       Q(质量=频率) = {qA:.6f}   Q(质量=相邻比) = {qB:.6f}")

# ── 反求：K₄ seed 相邻比需满足什么才能 Koide=2/3 ──
print("\n" + "=" * 82)
print("四、反求：Koide 成立（Q=2/3）的三质量比族（Koide 参数解）")
print("=" * 82)
print("  Koide 族 m_i ∝ (1+√2 cos(θ+2πi/3))² 恒满足 Q=2/3。其质量比随 θ：")
for t in [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]:
    fam = koide_family(t)
    mm = [x / min(fam) for x in fam]
    print(f"    θ={t:.1f}: 质量比(小:中:大)={[round(x,3) for x in mm]}  "
          f"Q={koide_q(fam):.6f}")

# ── K₄ 频率槽能否落入 Koide 族（质量=相邻比 r，且 r 取自 Koide 族）──
print("\n" + "=" * 82)
print("五、若 seed 相邻比取自 Koide 族，则 K₄ 三频率槽自动满足 Q=2/3")
print("=" * 82)
print("  即：存在 seed（s 单调）使 3 相邻比 = Koide 族质量比 ⟹ K₄ 三频率槽再现 Koide。")
print("  seed 从质量比重建：s=[1, r1, r1*r2, r1*r2*r3]（归一化），r_i = 相邻比。")
for t in [0.4, 0.6, 0.8, 1.0]:
    fam = koide_family(t)
    r = tuple(fam[0] / fam[1], ) if False else tuple(fam[i] / fam[(i + 1) % 3] for i in range(3))
    # 只取有序（小→中→大）的三代，相邻比按排序后取
    srt = sorted(fam)
    r_sorted = tuple(srt[i] / srt[i + 1] for i in range(2))
    r3 = (r_sorted[0], r_sorted[1], 1.0)
    seed = [1.0, r3[0], r3[0] * r3[1], r3[0] * r3[1] * r3[2]]
    print(f"    θ={t:.1f}: 三代质量={[round(x,3) for x in srt]}  "
          f"相邻比={[round(x,3) for x in r3]}  → seed(末两个)=[...,{round(seed[2],3)},{round(seed[3],3)}]")

print("\n" + "=" * 82)
print("六、K₄ 三频率槽能否直接匹配轻子观测质量比")
print("=" * 82)
me, mmu, mtau = obs["轻子(e,μ,τ)"]
ratio_obs = sorted((me / mme if False else 1.0, mmu / me, mtau / me))
print(f"  轻子质量比(归一化 e=1)：1 : {mmu/me:.3f} : {mtau/me:.1f}  Q={koide_q((me,mmu,mtau)):.6f}")
# K₄ 频率槽直接取轻子质量比作为 seed 相邻比，检验自洽
seed_l = [1.0, 1.0, mmu / me, mtau / me]  # 单调
r_l = (seed_l[1] / seed_l[0], seed_l[2] / seed_l[1], seed_l[3] / seed_l[2])
print(f"  若 seed 单调取轻子质量 e,μ,τ：相邻比={[round(x,4) for x in r_l]}")
print(f"    → 质量=相邻比映射下 Q={koide_q(r_l):.6f}（应≈2/3 若自洽）")
