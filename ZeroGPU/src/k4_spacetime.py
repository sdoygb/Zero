#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_spacetime.py —— 时空动态生成：一代一代迭代忠实地从 L2 演化长出势（模流）
============================================================================
依据 Zero 文章库：
  Z021: 维数 D = 通道数−1; 通道几何=正则单纯形(D+1顶点)
  Z033: 离散泊松 Lφ=−x 在连通图上精确可解唯一(零和约束) → 势 φ 可算
        （轨道/位置/时间被否证: 缺位置算符与共同时间 → 只生成"势"，非经典坐标时空）
  Z088: 物理世界在 L2, 时空出现在 L2（滑动窗/格点/模流）, 闭合回灌=重数累积
设计（忠实）：
  每代: L2 演化(忠实引擎) → 荷分布 x(record 归一化零和)
        → 解 Lφ=−x 得势 φ（时空形状）
        → 随代迭代, φ 逐代重算 = 时空势的动态生成（模流）
  时空量: 势分布 φ, 梯度(力场 g=∇φ), 二阶差分(曲率代理 κ≈−Δφ)
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loyalty_engine import LoyaltyEngine

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]

def laplacian(adj):
    n = len(adj)
    L = np.zeros((n,n))
    for i, nb in enumerate(adj):
        L[i,i] += len(nb)
        for j in nb:
            L[i,j] -= 1
    return L

def solve_potential(L, x):
    """解 Lφ = −x（模常数），零和约束。用伪逆（零空间投影）。"""
    n = L.shape[0]
    # 投影到零和子空间（Σx=0, Σφ=0）
    x0 = x - x.mean()
    Lpinv = np.linalg.pinv(L)          # 图 Laplacian 伪逆（零空间=常数向量）
    phi = -Lpinv @ x0
    phi -= phi.mean()                  # 模常数
    return phi

print("="*84)
print("时空动态生成：一代一代从 L2 演化长出势 φ（模流）")
print("="*84)

n = len(K4)
L = laplacian(K4)
print(f"时空骨架图 = K₄（{n} 顶点）")
print(f"通道结构 → 维数 D = 通道数−1（Z021）")
print(f"图 Laplacian L = diag(deg)−A，Z033 离散泊松 Lφ=−x 可解唯一(零和)")

# 忠实引擎，每代采当前闭合分布（归一化零和）→ 势
eng = LoyaltyEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2")
gens = 300
pot_traj = []
for g in range(gens):
    eng.run_generation()
    if (g+1) % 50 == 0 or g < 5:
        # 当前代闭合分布（动态），归一化零和 → 势量级可控
        x = eng.gen_closures.copy()
        x = x / x.sum() - 1.0/len(x)          # 归一化 + 零和
        phi = solve_potential(L, x)
        # 时空量
        grad = np.abs(np.diff(np.sort(phi)))        # 势梯度（力场代理）
        curva = -L @ phi                            # 二阶差分（曲率代理 ≈ −Δφ）
        pot_traj.append({"gen": g+1, "x": [round(v,4) for v in x],
                         "phi": [round(v,4) for v in phi],
                         "grad_norm": round(float(np.linalg.norm(grad)),5),
                         "curv": [round(v,4) for v in curva]})
        print(f"  代{g+1:>3}: 势φ={[round(v,4) for v in phi]}  "
              f"梯度范数={np.linalg.norm(grad):.4f} 曲率={[round(v,3) for v in curva]}")

# 势是否收敛（模流是否到稳态）
final_x = eng.gen_closures / eng.gen_closures.sum() - 1.0/len(eng.gen_closures)
final_phi = solve_potential(L, final_x)
print("\n势的收敛（模流稳定性）:")
print(f"  最终势 φ = {[round(v,3) for v in final_phi]}")
print(f"  势方差 = {np.var(final_phi):.4f}")

# 时空量的物理对照（诚实）
print("\n" + "="*84)
print("时空动态生成结果")
print("="*84)
print(f"  每代: 荷分布 x(record零和) → 解 Lφ=−x → 势 φ（时空形状）随代迭代 = 模流")
print(f"  时空量: 势φ（引力势）, 梯度g=∇φ（力场）, 曲率κ=−Δφ（弯曲代理）")
print(f"  模流: 势随代演化（动态生成）; 若收敛→时空稳态（不动点）")
print("  【诚实边界】生成图上的势+通道结构(Z033允许), 非位置+时间的经典时空(Z033否证)")
print("  Z033: 算得出势, 算不出轨道——缺位置算符与共同时间")

out = {"config": {"graph": "K4", "n_verts": n, "gens": gens},
       "spacetime_definition": {"dim": "D=通道数−1(Z021)", "potential": "Lφ=−x(Z033)"},
       "potential_traj": pot_traj,
       "final_phi": [round(v,3) for v in final_phi],
       "conclusion": "每代从L2演化解出势φ随代迭代=时空动态生成(模流)",
       "boundary": "图上的势+通道结构, 非经典时空坐标(Z033否证轨道/位置/时间)"}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_spacetime.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"\n已存 {os.path.join(RES,'k4_spacetime.json')}")
