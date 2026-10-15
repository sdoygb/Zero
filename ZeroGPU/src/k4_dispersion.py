#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_dispersion.py —— 验证色散关系：模流时空分解, 找 Zero 内生的传播模式(波)
====================================================================================
用 Zero 结构代替经典位置/时间/相位后, "波" = 势在图上随演化步传播的时空模式:
  φ(v, gen) ~ Σ A_m cos(ω_m·gen − k_m·v)
  时间频率 ω = 沿代(gen) FFT
  空间波数 k = 图 Laplacian 特征值(内生空间频率)
验证: 模流是否存在 ω−k 色散关联的传播模式(势扰动在空间模式上以确定频率振荡/传播)
"""
from __future__ import annotations
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from k4_inhomogeneity_scan import InhomogEngine, laplacian, solve_potential
from loyalty_engine import path_adj

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def run_dispersion(adj, p, seed, label, N=3000, lam=1.0, feedback=0.5):
    n = len(adj); L = laplacian(adj)
    eng = InhomogEngine(adj, seed=seed, T=4, retention="top2", p=p, beta=0.2, feedback=feedback)
    phis = []
    for _ in range(N):
        eng.run_generation()
        m = eng.mass_dist(eng.record)
        xs = m - 1.0/n
        phis.append(solve_potential(L, xs))
    phis = np.array(phis)                       # (N, n)
    # 空间分解: Laplacian 特征分解
    w, Vv = np.linalg.eigh(L)                    # 特征值 w(升序), 特征向量列 Vv
    coeffs = phis @ Vv                           # (N, n): 每列=某空间模式随代演化
    print(f"\n{'='*80}\n{label}（n={n}顶点, p={p}）: 时空分解")
    rows = []
    for m in range(n):
        k = w[m]
        sig = coeffs[:, m]
        sig0 = sig - sig.mean()
        spec = np.abs(np.fft.rfft(sig0))
        freqs = np.fft.rfftfreq(N, 1)
        kk = np.argsort(spec)[::-1]
        # 主导频率(排除~0漂移)
        peakf = freqs[kk[0]]; peakamp = spec[kk[0]]
        # 排除常数模式漂移: 用次主导正频率
        for i in kk[:5]:
            if freqs[i] > 0.005:   # 排除慢漂移低频
                peakf = freqs[i]; peakamp = spec[i]; break
        amp_norm = np.linalg.norm(spec)
        energy_share = (peakamp**2)/(np.sum(spec**2)+1e-12)
        rows.append({"mode": m, "k": round(float(k),5), "omega": round(float(peakf),4),
                     "amp": round(float(peakamp),5), "energy_share": round(float(energy_share),4)})
        print(f"  模式{m}: 空间波数k={k:.4f}  主导频率ω={peakf:.4f}  幅度={peakamp:.4f}  占比={energy_share:.3f}")
    return w, Vv, phis, rows

# ---- K₄ ----
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
w4, V4, ph4, rows4 = run_dispersion(K4, p=1.0, seed=[8.0,4.0,1.5,1.0], label="K₄")

# 色散关系判断
print("\nK₄ 色散关系 (ω vs k):")
for r in rows4:
    print(f"  k={r['k']:.4f} → ω={r['omega']:.4f}")
k_set = set(round(r["k"],3) for r in rows4)
omega_set = set(round(r["omega"],3) for r in rows4)
if len(k_set) <= 1:
    print("  ⚠ K₄ 谱简并(特征值{0,4,4,4}), 空间频率单一 → 需更大图验证色散")
else:
    print("  K₄ 有多个空间频率, 检查 ω−k 关联")

# ---- 更大图: path8 ----
print("\n" + "="*80)
print("换 path8 图（8 顶点, 谱丰富）验证色散关系")
adj8 = path_adj(8); n8 = len(adj8)
seed8 = [float(2**(n8-1-v)) for v in range(n8)]
w8, V8, ph8, rows8 = run_dispersion(adj8, p=1.0, seed=seed8, label="path8")

print("\npath8 色散关系 (ω vs k):")
disp = []
for r in rows8:
    if r["k"] > 1e-6:   # 排除常数模式
        disp.append((round(r["k"],3), round(r["omega"],3)))
for k, om in disp:
    print(f"  k={k} → ω={om}")

# 色散关联: 相关系数(k vs ω, 排除漂移)
if len(disp) > 2:
    ks = np.array([d[0] for d in disp]); oms = np.array([d[1] for d in disp])
    corr = np.corrcoef(ks, oms)[0,1]
    print(f"  色散 ω vs k 相关系数 = {corr:.4f}")
    if abs(corr) > 0.5:
        print("  → 存在 ω−k 色散关联: 势扰动在空间模式上有确定频率 = Zero 内生传播模式(波)候选")
    else:
        print("  → 无强 ω−k 色散关联: 传播模式不明显, 记录为收敛缺口")

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_dispersion.json"), "w", encoding="utf-8") as f:
    json.dump({"K4": rows4, "path8": rows8}, f, ensure_ascii=False, indent=1)
print(f"\n已存 {os.path.join(RES,'k4_dispersion.json')}")
