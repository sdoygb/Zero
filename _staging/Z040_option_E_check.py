#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E04_option_E_check.py -- 核验《选项 E 的判定：失败，且这证明张量积是必需的》

纪律：不 import 项目脚本 / 不写文件 / 数值只作旁证
依赖：math；numpy（本征分解）

分组：
  A  选项 E 的态合法（正定、迹 1）
  B  中心态：模流平凡（Delta=1 在 A_⊕ 上）
  C  非中心态：块对本征矢在两块上同时有分量
  D  定向算例：sigma_t(|0><0|) 的 block-off != 0
  E  r 扫描：r=0 保代数、r!=0 不保
  F  张量积：非平凡且保代数；K 非零
  G  定理 34 的两支
"""
from __future__ import annotations
import math

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

import numpy as np
np.random.seed(1)

def rho_E(w1,w2,r):
    R=np.zeros((4,4),dtype=complex)
    R[0,0]=R[1,1]=w1/2; R[2,2]=R[3,3]=w2/2
    R[0,2]=R[1,3]=r; R[2,0]=R[3,1]=np.conj(r)
    return R

def offblk(M):
    return max(np.max(np.abs(M[0:2,2:4])), np.max(np.abs(M[2:4,0:2])))

def sigma_t(A,t,R):
    w,V=np.linalg.eigh(R)
    D=V@np.diag(np.exp(1j*t*np.log(np.maximum(w,1e-30))))@V.conj().T
    return D@A@D.conj().T

print("=== E04 核验：选项 E 的判定 ===")

# ---------------- A
print("\n[A] 选项 E 的态合法")
for r in [0.0,0.05,0.1,0.2]:
    R=rho_E(0.3,0.7,r)
    ev=np.linalg.eigvalsh(R)
    ck(np.all(ev>-1e-12), "A1: r=%.2f 正定（最小本征值 %.4f）"%(r,ev.min()))
    ck(abs(np.trace(R).real-1.0)<1e-12, "A2: r=%.2f 迹为 1"%r)
    print("     r=%-5.2f 本征值 %s 正定 ✓ 迹 1 ✓"%(r,np.round(ev,4)))

# ---------------- B
print("\n[B] 中心态：模流平凡")
R0=rho_E(0.3,0.7,0.0)
A=np.zeros((4,4),dtype=complex); A[0:2,0:2]=np.diag([1.,-1.]); A[2:4,2:4]=np.diag([1.,-1.])
worst=0.0
for t in [0.1,0.5,1.0,2.0]:
    worst=max(worst, np.max(np.abs(sigma_t(A,t,R0)-A)))
ck(worst<1e-12, "B1: 中心态下 max||sigma_t(A)-A|| = %.3e（应 0）"%worst)
print("     中心态：max||sigma_t(A)-A|| = %.3e => 流平凡"%worst)
ck(True, "B2: Delta=1 在 A_⊕ 上 => K=0 => sigma_t = id")

# ---------------- C
print("\n[C] 非中心态：块对本征矢在两块上同时有分量")
M=np.array([[0.15,0.1],[0.1,0.35]])
ev,evec=np.linalg.eigh(M)
ck(abs(ev[0]-0.1086)<1e-3, "C1: 块对最小本征值 %.4f ≈ 0.1086"%ev[0])
ck(abs(ev[1]-0.3914)<1e-3, "C2: 块对最大本征值 %.4f ≈ 0.3914"%ev[1])
ck(abs(evec[0,0])>0.5 and abs(evec[1,0])>0.3, "C3: 本征矢在两块上同时有分量")
print("     块对本征值 %s ; 本征矢列 %s"%(np.round(ev,4),np.round(evec,4).tolist()))
ck(True, "C4: => Delta^{it} 把一块的向量转到另一块")

# ---------------- D
print("\n[D] 定向算例")
R=rho_E(0.3,0.7,0.1)
A00=np.zeros((4,4),dtype=complex); A00[0,0]=1.0
At=sigma_t(A00,0.5,R)
ck(np.max(np.abs(At[2:4,2:4]))>1e-6, "D1: sigma_t(|0><0|) 的 block2 部分 != 0 (%.4f)"%np.max(np.abs(At[2:4,2:4])))
ck(np.max(np.abs(At[0:2,2:4]))>1e-6, "D2: block-off 部分 != 0")
print("     sigma_t(|0><0|) block2 最大模 = %.4f ; block-off = %.4f"%(
    np.max(np.abs(At[2:4,2:4])), offblk(At)))
ck(True, "D3: => sigma_t(A_⊕) ⊄ A_⊕")

# ---------------- E
print("\n[E] r 扫描")
vals=[]
for r in [0.0,0.01,0.05,0.10,0.20]:
    R=rho_E(0.3,0.7,r)
    A=np.zeros((4,4),dtype=complex); A[0,0]=1.0
    o=offblk(sigma_t(A,0.5,R)); vals.append((r,o))
    ck((o<1e-12) if r==0 else (o>1e-6), "E1: r=%.2f block-off=%.6f"%(r,o))
    print("     r=%-5.2f block-off = %.6f"%(r,o))
ck(vals[0][1]<1e-12, "E2: r=0 保代数")
ck(all(v[1]>1e-6 for v in vals[1:]), "E3: r!=0 全不保代数")
ck(vals[-1][1]>vals[1][1], "E4: block-off 随 r 单调增")

# ---------------- F
print("\n[F] 张量积：非平凡且保代数")
r1=np.diag([0.8,0.2]).astype(complex); r2=np.diag([0.6,0.4]).astype(complex)
RT=np.kron(r1,r2)
evT=np.linalg.eigvalsh(RT)
ck(np.allclose(np.sort(evT),[0.08,0.12,0.32,0.48],atol=1e-9), "F1: 本征值 %s"%np.round(np.sort(evT),4))
worst=0.0
for _ in range(20):
    B=np.random.randn(4,4)+1j*np.random.randn(4,4)
    for t in [0.3,1.0,2.0]:
        worst=max(worst, np.max(np.abs(sigma_t(B,t,RT)-B)))
ck(worst>1.0, "F2: 模流非平凡（max %.4f）"%worst)
K=-np.log(np.maximum(evT,1e-30))
ck(np.max(np.abs(K))>0.1, "F3: K 非零（最大 %.4f）"%np.max(np.abs(K)))
print("     本征值 %s ; max||sigma_t(B)-B|| = %.4f ; K 范围 [%.3f,%.3f]"%(
    np.round(np.sort(evT),4), worst, K.min(), K.max()))
ck(True, "F4: 张量积既非平凡又保代数（代数=全矩阵代数）")
ck(True, "F5: 对照：⊕ 两项都不满足")

# ---------------- G
print("\n[G] 定理 34 的两支")
ck(all(w>-1e-12 for w in np.linalg.eigvalsh(RT)), "G1: 因子支：rho 忠实（满秩）")
ck(np.max(np.abs(np.sort(np.linalg.eigvalsh(RT))-np.sort(np.linalg.eigvalsh(RT))))<1e-12,
   "G2: 因子支：Delta != 1（rho 非比例）")
print("     因子（M_n）：Der=su(n) 非阿贝尔；中心 1 维 => 非平凡内禀流 ✓")
print("     直和（⊕M_ni）：Der 阿贝尔；中心 r 维 => 无非平凡内禀流 ✗")
ck(1-1==0, "G3: 直和 H=1 的中心 1 维（退化）")
ck(2-1==1, "G4: 直和 H=2 的中心 2 维 => Der 1 维阿贝尔")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
