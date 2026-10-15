# -*- coding: utf-8 -*-
"""Z131 核验：Zero 原生演化律 = 模流 σ_t(A)=ρ^{it}Aρ^{-it}，非中心态以 log(p/q) 振荡。
A: 2:1 非中心态 → 实测角频率 = log2
B: K4 稳态记录 = seed 归一化
C: 记录相邻权重比 log = K4 三频率槽
D: K4 稳态记录的模流，各相邻对实测角频率 = log 相邻比
独立复算：纯 numpy，不 import 项目脚本、不写文件。FFT 循环频率 ×2π 换算角频率。
"""
import numpy as np

def modflow_freq(rho, a, b, T=400.0, N=8000):
    """模流 σ_t(A)=ρ^{it}Aρ^{-it} 在非中心态 ρ 上，离对角元(a,b)的实测角频率(FFT)。"""
    ev, vec = np.linalg.eigh(rho)
    t = np.linspace(0, T, N)
    def el(s):
        Rp = (vec * np.exp(1j * s * np.log(ev))).dot(vec.T)      # ρ^{it}（本征分解，复数指数）
        Rm = (vec * np.exp(-1j * s * np.log(ev))).dot(vec.T)     # ρ^{-it}
        A = np.zeros((len(rho), len(rho))); A[a, b] = 1.0
        At = Rp @ A @ Rm
        return np.real(At[a, b])                                 # 离对角矩阵元(迹对离对角恒0)
    sig = np.array([el(s) for s in t])
    sp = np.fft.rfft(sig - np.mean(sig)); fr = np.fft.rfftfreq(N, t[1] - t[0])
    k = np.argmax(np.abs(sp)[1:]) + 1
    return 2 * np.pi * abs(fr[k])                                # 角频率

print("=" * 84)
print("Z131 核验：模流在非中心态上以角频率 log(p/q) 振荡；K4 频率槽即此产物")
print("=" * 84)
ok_all = True

# A: 2:1 非中心态
rho = np.diag([2/3, 1/3]); fA = modflow_freq(rho, 0, 1); theoA = np.log(2/3/(1/3))
mA = abs(fA - theoA) < 0.02; ok_all &= mA
print(f"[A] 2:1非中心态: 实测角频率={fA:.4f} vs log2={theoA:.4f}  {'✓' if mA else '✗'}")

# B: 稳态记录 = seed 归一化
seed = np.array([8.0, 4.0, 1.5, 1.0]); rec = np.array([0.5517, 0.2759, 0.1034, 0.069])
ns = seed / seed.sum(); mB = np.allclose(rec, ns, atol=2e-3); ok_all &= mB
print(f"[B] seed归一{np.round(ns,4)} vs 记录{rec}  {'✓' if mB else '✗'}")

# C: 记录相邻比 log = K4 三频率槽
slots = np.array([0.4055, 0.6931, 0.9808])
adjl = np.sort(np.log(rec[:3] / rec[1:]))
mC = np.allclose(adjl, slots, atol=5e-3); ok_all &= mC
print(f"[C] 相邻比log {np.round(adjl,4)} vs K4槽{slots}  {'✓' if mC else '✗'}")

# D: K4 稳态记录模流，相邻对
print("[D] K4稳态记录 ρ=diag(rec) 的模流，各相邻对实测角频率:")
rho4 = np.diag(rec)
for a, b in [(0, 1), (1, 2), (2, 3)]:
    f = modflow_freq(rho4, a, b); theo = np.log(rec[a] / rec[b]); m = abs(f - theo) < 0.02
    ok_all &= m
    print(f"    对({a},{b}): 实测={f:.4f} vs log({rec[a]:.4f}/{rec[b]:.4f})={theo:.4f}  {'✓' if m else '✗'}")

print("=" * 84)
print(f"【判定】{'✓ 演化律成立: 非中心态模流以 log(权重比) 振荡, K4 频率槽即此原生律产物' if ok_all else '✗ 不成立'}")
print("=" * 84)
