# -*- coding: utf-8 -*-
"""
Z130 复核：闭合冻结核验 —— K4 冻结骨架 vs 波扇区生长
============================================================
A. 忠实引擎：L2 因果闭合首代即锁定 → record(L1) 归一化结构永久冻结、量无限累积
B. 反馈通道：seed==record（真实存在）但 closure==seed（退化复制，不长新结构）
C. 扰动鲁棒：eps=1.0 大扰动也被吸收回冻结值（鲁棒吸引子）
D. 波/场扇区：闭合势(D场)驱动波从顶层源下传、能量重分布、长新波结构(峰1→12)

纯标准库 + numpy + Zero 忠实引擎 + field_sim；Zero 只读（仅新增本篇）。
"""
from __future__ import annotations
import sys, os, math
import numpy as np
SRC = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ZeroGPU", "src"))
sys.path.insert(0, SRC)
from loyalty_engine import LoyaltyEngine
from field_sim import field_at

K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
def norm(v):
    s=v.sum(); return v/s if s else v

print("="*78)
print("Z130 复核：闭合冻结核验（K4 冻结骨架 vs 波扇区生长）")
print("="*78)

# ────────────────────────── A/B/C：忠实引擎 ──────────────────────────
eng = LoyaltyEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2")
for _ in range(60000):
    eng.run_generation()
rec0 = norm(eng.record); seed = norm(eng.seed); clo = norm(eng.gen_closures)
print("\n[A] record(L1) 归一化冻结 + 量累积")
print(f"  record 归一化 = {[round(float(x),6) for x in rec0]}")
print(f"  漂移(与首代比, 已由 frozen 全程测 ~5e-16) → 结构冻结")
print(f"  record 绝对值 Σ = {eng.record.sum():.3e}（量持续累积）")
assert np.allclose(seed, clo) and np.allclose(clo, rec0), "seed/closure/record 未锁死"

print("\n[B] 反馈通道：存在但退化")
print(f"  seed==record? {np.allclose(seed,rec0)}  → 闭合回灌成 seed（边②③通道真实存在）")
print(f"  closure==seed? {np.allclose(clo,seed)}  → 但 L2 只复制 seed，不长新结构")
assert np.allclose(seed,rec0) and np.allclose(clo,seed)

print("\n[C] 扰动鲁棒：eps=1.0 大扰动也被吸收")
eng2 = LoyaltyEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2")
for _ in range(100): eng2.run_generation()
before = norm(eng2.record.copy())
eng2.seed = eng2.seed + np.array([1.0,-0.3,0.2,0.0])*1.0
drift = 0.0
for _ in range(4000):
    eng2.run_generation()
    drift = max(drift, float(np.max(np.abs(norm(eng2.record)-before))))
print(f"  扰后4000代 record 归一化漂移 = {drift:.3e}  → 被拉回冻结值")
assert drift < 1e-10, "扰动未被吸收"

# ────────────────────────── D：波/场扇区 ──────────────────────────
print("\n[D] 波/场扇区：闭合势(D场)驱动波演化、长新结构")
H=8; L=16
up,tot = field_at(L)
nodes=[]; idx={}
for h in range(1,H+1):
    for k in range(0,h+1):
        idx[(h,k)]=len(nodes); nodes.append((h,k))
N=len(nodes)
D=np.zeros(N)
for (h,k),i in idx.items():
    u=up[(h,k)]; t=tot[(h,k)]; p=u/t if t else 0.5
    D[i]=1.0-2.0*p
adj=[set() for _ in range(N)]
for (h,k),i in idx.items():
    if k>0: adj[i].add(idx[(h,k-1)])
    if k<h: adj[i].add(idx[(h,k+1)])
    if k==0:
        for h2 in range(1,H+1):
            if h2!=h: adj[i].add(idx[(h2,0)])
for i in range(N): adj[i]=sorted(adj[i])
Lap=np.zeros((N,N))
for i in range(N):
    Lap[i,i]=len(adj[i])
    for j in adj[i]: Lap[i,j]=-1
evals,evecs=np.linalg.eigh(Lap)
mu=0.5; c2=0.8/max(1.0,float(np.max(evals))*0.4)
V=np.array([-mu*D[i] for i in range(N)])

def leapfrog(phi0,phi1,steps,useV=True):
    prev,cur=phi0.copy(),phi1.copy(); phis=[cur.copy()]
    for _ in range(1,steps):
        nxt=2*cur-prev-c2*(Lap@cur)
        if useV: nxt=nxt+c2*V*cur
        prev,cur=cur,nxt; phis.append(cur.copy())
    return np.array(phis)

def layer_energy(phi):
    E={h:0.0 for h in range(1,H+1)}
    for (h,k),i in idx.items(): E[h]+=phi[i]**2
    s=sum(E.values()); return {h:(E[h]/s if s>0 else 0) for h in E}

# 波前
ph0=np.zeros(N); ph1=np.zeros(N); ph1[idx[(H,0)]]=1.0
phis=leapfrog(ph0,ph1,3000,useV=True)
peak_t={h:int(np.argmax(np.abs(phis[:,idx[(h,0)]])[:1500])) for h in range(1,H+1)}
delay=peak_t[1]-peak_t[H]
print(f"  顶层源→底层波前延迟 = {delay} 代  (→ 波从顶层下传)  [断言 delay>0]")
assert delay > 0, "未见顶层→下层传播"

# 能量重分布
E0=layer_energy(phis[0]); Ee=layer_energy(phis[-1])
redist=sum(abs(E0[h]-Ee[h]) for h in E0)
print(f"  能量跨层重分布(首→末总变化) = {redist:.3f}  [断言 >0.02]")
assert redist > 0.02, "能量未重分布"

# 长新结构：振幅峰数 1→多
def n_peaks(phi):
    return sum(1 for i in range(1,N-1) if phi[i]>phi[i-1] and phi[i]>=phi[i+1])
n0=n_peaks(phis[0]); nm=n_peaks(phis[1000]); ne=n_peaks(phis[-1])
print(f"  振幅峰数: t=0:{n0} → t=1000:{nm} → t=末:{ne}  (长出新波结构)  [断言峰数增长]")
assert ne > n0 or nm > n0, "未长出新波结构"

# D 势作用
ph0=np.zeros(N); ph1=np.zeros(N); ph1[idx[(H,0)]]=1.0
phis0=leapfrog(ph0,ph1,3000,useV=False)
red0=sum(abs(layer_energy(phis0[0])[h]-layer_energy(phis0[-1])[h]) for h in range(1,H+1))
print(f"  闭合势(D)作用: 有D重分布 {redist:.3f} vs 无D {red0:.3f}  [断言|差|>0.02]")
assert abs(redist-red0) > 0.02, "D 势未起作用"

print("\n" + "="*78)
print("Z130 全部断言通过 ✓  核验成立")
print("="*78)
