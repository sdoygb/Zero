# -*- coding: utf-8 -*-
"""
ZeroGPU · 粒子模拟（数量扩展）：种类之外，加"数量/丰度"

延续 k4_particles.py（种类谱：3 代 × 种类 × 量子数），补充第四维——数量。

【数量的来源（机制声明）】
  引擎里唯一携带"产生多少"的量 = 稳态占据权重 w = record/Σrecord。
  4 个 K₄ 态槽 → 3 个相邻频率槽（=3 代质量）：
    代 d 的频率槽来自 |log(w[d-1]/w[d])|，故代 d 的"人口"取源顶点权重 w[d-1]。
  每个粒子的相对数量 = 代人口 × 简并度：
    轻子(带电/中微子) ×1；夸克(上/下) ×Nc=3（色三重）；玻色子 ×1。

【诚实边界】"数量=稳态占据权重×简并度"是【一种内生选择】。备选：
  (a) 热平衡丰度 n ∝ exp(-m/T)×简并度 —— 需引入温度 T（外部/锚定），轻代更多；
  (b) Z129 非对角流源/汇净收支 —— 源顶点净产、汇顶点净收，给出"生产/吸收"方向；
  三者对"哪一代最多"的排序一致（第 1 代最多、第 3 代最少），但具体比例不同。
  本文默认取 (0) 占据权重，零新参数、完全内生；换机制需重跑。

用法：/usr/local/bin/python3 k4_particles_quantity.py
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loyalty_engine import LoyaltyEngine

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
Nc = 3                       # 色三重简并
DEG = {"轻子带电":1, "轻子中微子":1, "上型夸克":Nc, "下型夸克":Nc, "玻色子(矢量/标量)":1}

# ── ① 忠实引擎稳态（与 k4_particles.py 一致）──
eng = LoyaltyEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2")
for _ in range(60000):
    eng.run_generation()
w = eng.record / eng.record.sum()
freqs = sorted([abs(math.log(w[i]/w[i+1])) for i in range(3)])
mass = [x/max(freqs) for x in freqs]

# 代人口 = 源顶点权重（频率槽 d 的源 = w[d-1]）
gen_pop = np.array([w[0], w[1], w[2]])          # 未归一，相对即可
gen_pop_n = gen_pop / gen_pop.sum()             # 归一（占三代总人口比例）
print("=" * 88)
print("粒子模拟（数量扩展）：种类 + 数量/丰度（数量=稳态占据权重×简并度）")
print("=" * 88)
print(f"稳态 record 占据权重 w = {[round(float(x),4) for x in w]}")
print(f"3 代质量(归一)      = {[round(m,4) for m in mass]}  （代1最轻→代3最重）")
print(f"3 代人口(源顶点)    = {[round(float(x),4) for x in gen_pop]}  "
      f"（归一 {[round(float(x),4) for x in gen_pop_n]}）→ 轻代最多")
print(f"简并度: 轻子×1, 夸克×{Nc}(色), 玻色子×1")

# ── ② 组装带数量的粒子谱 ──
particles = []
Y_L_lepton = -1.0; Y_L_quark = -Y_L_lepton / Nc
def charge(t3, Y): return t3 + Y/2.0

def add(d, mval, kind, symbol, color, spin, q):
    deg = DEG[kind]
    n = gen_pop[d-1] * deg                 # 相对数量
    particles.append({"代": d, "m": mval, "种类": kind, "符号": symbol,
                      "色": color, "自旋": spin, "电荷": q,
                      "简并度": deg, "数量": round(float(n), 6)})

for d, m in enumerate(mass, 1):
    add(d, m, "轻子带电", "-", "单态", "1/2(费米)", charge(-0.5, Y_L_lepton))
    add(d, m, "轻子中微子", "0", "单态", "1/2(费米)", charge(+0.5, Y_L_lepton))
    add(d, m, "上型夸克", "+2/3", "三重", "1/2(费米)", charge(+0.5, Y_L_quark))
    add(d, m, "下型夸克", "-1/3", "三重", "1/2(费米)", charge(-0.5, Y_L_quark))
    add(d, m, "玻色子(矢量/标量)", "0/±1", "单态", "0/1(玻色)", charge(0.0, 0.0))

# ── ③ 输出带数量的谱 ──
tot = sum(p["数量"] for p in particles)
print(f"\n{'代':<3}{'种类':<18}{'电荷':<7}{'色':<5}{'简并':<5}{'数量(相对)':<12}{'数量%':<9}{'质量'}")
print("-" * 88)
for p in particles:
    print(f"{p['代']:<3}{p['种类']:<18}{p['符号']:<7}{p['色']:<5}"
          f"{p['简并度']:<5}{p['数量']:<12.4f}{100*p['数量']/tot:<9.2f}{p['m']:.4f}")

# 每代小计
print("\n每代小计（含色/自旋简并）：")
for d in range(1, 4):
    sub = sum(p["数量"] for p in particles if p["代"] == d)
    print(f"  第 {d} 代: 数量 = {sub:.4f}  ({100*sub/tot:.1f}% of total)")
print(f"  总相对数量 = {tot:.4f}（归一化: 100%）")

# 丰度排序
print("\n丰度排序（相对数量降序，前 5）：")
for p in sorted(particles, key=lambda x: -x["数量"])[:5]:
    print(f"  {p['种类']}(代{p['代']})  数量={p['数量']:.4f}")

print("\n" + "=" * 88)
print("机制说明：数量 = 稳态占据权重(源顶点) × 简并度；零新参数、完全内生。")
print("备选机制（需重跑）：(a)热平衡 exp(-m/T)；(b)Z129 非对角流源/汇净收支。")
print("三者对'第1代最多/第3代最少'的排序一致，比例不同。")
print("=" * 88)

# ── ④ 落盘 ──
out = {"record": [float(x) for x in w],
       "gen_mass": [float(m) for m in mass],
       "gen_pop": [float(x) for x in gen_pop],
       "gen_pop_norm": [float(x) for x in gen_pop_n],
       "degeneracy": DEG,
       "quantity_mechanism": "稳态占据权重(源顶点)×简并度（内生，零参数）",
       "quantity_alternatives": ["热平衡 exp(-m/T)×简并(需T)",
                                  "Z129 非对角流源/汇净收支"],
       "particles": particles,
       "total_relative": float(tot)}
os.makedirs(RES, exist_ok=True)
fp = os.path.join(RES, "k4_particles_quantity.json")
with open(fp, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"已存 {fp}")
