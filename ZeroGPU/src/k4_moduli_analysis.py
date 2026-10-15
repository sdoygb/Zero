#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_moduli_analysis.py —— 深挖时空动态模流的性质：漂移 vs 振荡、频谱、周期结构
====================================================================================
对 Z101 的持续动态模流（λ=1.0+反馈）做时间序列分析：
  1) 漂移: 每顶点势的线性趋势斜率 → 是收敛(斜率→0) 还是持续漂移
  2) 频谱: 去趋势后 FFT → 是否有主导频率(特征时间尺度), 对比粒子频率槽
  3) 周期: 毁灭周期 T=4 是否引入周期结构; 自相关
  4) 顶点差异: 哪个顶点动/静
"""
from __future__ import annotations
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from k4_particle_spacetime import ParticleSpacetimeEngine, laplacian, solve_potential, freq_mass

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
L4 = laplacian(K4); n = 4

print("="*88)
print("深挖时空动态模流性质：漂移 vs 振荡、频谱、周期结构")
print("="*88)

eng = ParticleSpacetimeEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2",
                              lam=1.0, beta=0.2, feedback=0.5)
N = 3000
phis = []
for g in range(N):
    eng.run_generation()
    m = freq_mass(eng.record, n)
    x = eng.gen_closures / eng.gen_closures.sum() - 1.0/n
    xs = 1.0*(m/m.sum() if m.sum()>0 else np.ones(n)/n); xs -= xs.mean()
    phis.append(solve_potential(L4, xs))
phis = np.array(phis)          # (N, 4)
t = np.arange(N)

# 1) 漂移: 每顶点线性趋势斜率 + 末段漂移率
print("\n【1】漂移分析（每顶点势的线性趋势）")
slopes = []
for v in range(n):
    sl = np.polyfit(t, phis[:, v], 1)[0]
    slopes.append(sl)
    # 末段(后1000代)漂移总量
    drift_tail = abs(phis[-1, v] - phis[-1000, v])
    print(f"  顶点{v}: 线性斜率={sl:+.6f}/代  末1000代漂移={drift_tail:.5f}  "
          f"势末值={phis[-1,v]:+.4f}")
slopes = np.array(slopes)
print(f"  → 全顶点斜率均值为{slopes.mean():+.6f}, 若|斜率|持续>0 → 持续漂移(未收敛)")

# 2) 去趋势 + 频谱
print("\n【2】频谱分析（去趋势后 FFT，找主导频率）")
detrend = phis - t[:, None]*slopes[None, :]
for v in range(n):
    sig = detrend[:, v] - detrend[:, v].mean()
    spec = np.abs(np.fft.rfft(sig))
    freqs = np.fft.rfftfreq(N, 1)            # 频率单位: 每代
    k = np.argsort(spec)[::-1]
    peakf = [round(freqs[i],4) for i in k[:4] if freqs[i] > 0]
    peakamp = [round(spec[i],5) for i in k[:4] if freqs[i] > 0]
    # 能量集中度: 主导频率占比
    energy = np.sum(spec**2)
    top_share = spec[k[0]]**2/energy if energy > 0 else 0
    print(f"  顶点{v}: 主导频率={peakf}  幅度={peakamp}  主导占比={top_share:.3f}")

# 3) 周期结构: 自相关 + 毁灭周期 T=4
print("\n【3】周期结构（自相关 + 毁灭周期 T=4）")
sig = detrend[:, 0]; sig -= sig.mean()
ac = np.correlate(sig, sig, mode='full')[N-1:] / (np.var(sig)*N + 1e-12)
for lag in [1, 2, 3, 4, 5, 8, 12, 16, 32, 64]:
    print(f"  lag={lag}: 自相关={ac[lag]:+.3f}")

# 4) 收敛性: 末段波动是否衰减(趋势外推)
print("\n【4】收敛性（末段趋势 vs 持续波动）")
seg_slope = np.polyfit(t[-1000:], detrend[-1000:, 0], 1)[0]
print(f"  末1000代去趋势残差斜率={seg_slope:+.6f} (→0 则振荡稳定; ≠0 则仍有慢趋势)")
print(f"  全程势动态总量={np.sum(np.linalg.norm(np.diff(phis,axis=0),axis=1)):.4f}")

# 粒子频率槽对比
w = np.sort(eng.record/eng.record.sum())[::-1]
freq_slots = [abs(np.log(w[i]/(w[i+1]+1e-12))) for i in range(min(3,n-1))]
print(f"\n粒子频率槽(Z089)={[round(f,4) for f in freq_slots]}")

os.makedirs(RES, exist_ok=True)
out = {"config": {"graph": "K4", "lam": 1.0, "feedback": 0.5, "N": N},
       "slopes": [round(s,6) for s in slopes],
       "particle_freq_slots": [round(f,4) for f in freq_slots]}
with open(os.path.join(RES, "k4_moduli_analysis.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"\n已存 {os.path.join(RES,'k4_moduli_analysis.json')}")
