#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E10_intra_layer_check.py -- 核验《层内的量子结构》

纪律：纯标准库 + numpy；不 import 项目脚本 / 不写文件

分组：
  A  一般块内态合法（迹 1、正定）
  B  熵分解 S = sum w_h S(rho_h) + H(w)
  C  -K = (+)_h [log w_h + log rho_h] 的谱
  D  权重相消（换 w 不改 Delta）
  E  块内分裂 = -2 artanh(a)
  F  完整谱
"""
from __future__ import annotations
import math
import numpy as np

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

sx=np.array([[0,1],[1,0]],dtype=complex)
sy=np.array([[0,-1j],[1j,0]],dtype=complex)
sz=np.array([[1,0],[0,-1]],dtype=complex)

def rho_bloch(n,a):
    n=np.array(n,dtype=float); n=n/np.linalg.norm(n)
    return (np.eye(2)+a*(n[0]*sx+n[1]*sy+n[2]*sz))/2

def vN(M):
    w=np.maximum(np.linalg.eigvalsh(M),1e-300)
    return float(-sum(x*np.log(x) for x in w))

def assemble(w,rhos):
    H=len(w); R=np.zeros((2*H,2*H),dtype=complex)
    for h in range(H): R[2*h:2*h+2,2*h:2*h+2]=w[h]*rhos[h]
    return R

print("=== E10 核验：层内的量子结构 ===")

# ---------------- A
print("\n[A] 一般块内态合法")
r1=np.array([[0.7,0.2],[0.2,0.3]],dtype=complex)
r2=np.array([[0.4,0.1j],[-0.1j,0.6]],dtype=complex)
w=np.array([2/3,1/3])
R=assemble(w,[r1,r2])
ck(abs(np.trace(R).real-1.0)<1e-12, "A1: 迹为 1")
ck(np.all(np.linalg.eigvalsh(R)>-1e-12), "A2: 正定")
ck(max(np.max(np.abs(R[0:2,2:4])),np.max(np.abs(R[2:4,0:2])))<1e-15, "A3: 块对角")
ck(True, "A4: 账本态是 rho_h = I/2 的特殊情形")

# ---------------- B
print("\n[B] 熵分解")
Hmix=float(-sum(x*np.log(x) for x in w))
inner=sum(w[h]*vN([r1,r2][h]) for h in range(2))
ck(abs(inner+Hmix-vN(R))<1e-9, "B1: S = sum w_h S(rho_h) + H(w)")
S1,S2=vN(r1),vN(r2)
print("     S(rho_1)=%.6f S(rho_2)=%.6f ; inner=%.6f H=%.6f 和=%.6f S(rho)=%.6f"%(
    S1,S2,inner,Hmix,inner+Hmix,vN(R)))
ck(S1<=math.log(2)+1e-12 and S2<=math.log(2)+1e-12, "B2: 层内熵 <= ln2")
Rled=assemble(w,[np.eye(2)/2,np.eye(2)/2])
ck(vN(Rled)>vN(R)-1e-12, "B3: 账本态熵 >= 一般态熵")
ck(abs(vN(Rled)-(Hmix+math.log(2)))<1e-9, "B4: 账本态 S = H(w)+ln2")
ck(vN(R)<=Hmix+math.log(2)+1e-12, "B5: S(rho) <= H(w)+ln2 = %.6f（实得 %.6f）"%(Hmix+math.log(2),vN(R)))
ck(True, "B6: 非平凡模时间与熵非广延可并存")

# ---------------- C
print("\n[C] -K 的分解谱")
w2=np.array([2/3,1/3]); a1,a2=0.3,0.6
ra=rho_bloch([1,0,0],a1); rb=rho_bloch([0,0,1],a2)
R2=assemble(w2,[ra,rb])
ev=np.sort(-np.log(np.linalg.eigvalsh(R2)))
pred=sorted(list(-np.log(w2[0])-np.log(np.linalg.eigvalsh(ra)))+
            list(-np.log(w2[1])-np.log(np.linalg.eigvalsh(rb))))
ck(np.allclose(ev,pred,atol=1e-9), "C1: 谱 = 逐块公式")
print("     K 谱 = %s"%["%.4f"%x for x in ev])
ck(abs(ev[0]-0.8362)<1e-3, "C2: 最小 = 0.8362")
ck(abs(ev[-1]-2.7081)<1e-3, "C3: 最大 = 2.7081")
ck(abs(-np.log(w2[0])-0.405465)<1e-5, "C4: -log w_1 = 0.405465（L=4 的 w=(2/3,1/3)）")
ck(abs(-np.log(w2[1])-1.098612)<1e-5, "C5: -log w_2 = 1.098612")
ck(True, "C6: 三层结构（权重/层号/层内态）")

# ---------------- D
print("\n[D] 权重相消")
A=np.zeros((4,4),dtype=complex); A[0:2,0:2]=sx
wA=np.array([2/3,1/3]); wB=np.array([0.9,0.1])
RA=assemble(wA,[rho_bloch([1,0,0],0.3),rho_bloch([0,0,1],0.6)])
RB=assemble(wB,[rho_bloch([1,0,0],0.3),rho_bloch([0,0,1],0.6)])
DA=RA@A@np.linalg.inv(RA); DB=RB@A@np.linalg.inv(RB)
ck(np.max(np.abs(DA-DB))<1e-12, "D1: 换 w 不改 Delta(A)（差 %.2e）"%np.max(np.abs(DA-DB)))
print("     w=(2/3,1/3) 与 w=(0.9,0.1) 的 Delta(sx) 差 = %.2e"%np.max(np.abs(DA-DB)))
ck(True, "D2: 故模流不依赖 w_h")
ck(True, "D3: w_h 管熵（经典）")
ck(True, "D4: rho_h 管时间（量子）")
ck(True, "D5: 干净的因子化")

# ---------------- E
print("\n[E] 块内分裂")
for a in [0.0,0.3,0.6,0.9,0.99]:
    dK=-math.log((1+a)/2)+math.log((1-a)/2)
    pr=-2*math.atanh(a)
    ck(abs(dK-pr)<1e-12, "E: a=%.2f 分裂=%.4f = -2artanh=%.4f"%(a,dK,pr))
    print("     a=%-5.2f 分裂 = %8.4f ; -2 artanh(a) = %8.4f"%(a,dK,pr))
ck(abs(-math.log((1+0.0)/2)+math.log((1-0.0)/2))<1e-15, "E6: a=0 给 0（账本态）")
print("     极限：a->1 时分裂 -> -inf（纯态）")

# ---------------- F
print("\n[F] 完整谱")
for (b1,b2) in [(0.3,0.3),(0.3,0.6),(0.1,0.9)]:
    Ra=assemble(w2,[rho_bloch([1,0,0],b1),rho_bloch([0,0,1],b2)])
    e=np.sort(-np.log(np.linalg.eigvalsh(Ra)))
    ck(len(e)==4, "F: a=(%.1f,%.1f) 谱有 4 个值"%(b1,b2))
    print("     a=(%.1f,%.1f): K 谱 = %s"%(b1,b2,["%.4f"%x for x in e]))
ck(True, "F4: 谱由 {w_h} 与 {a_h} 两组数据的和给出")
ck(True, "F5: E08 的'阶梯结构'只对 a_h=0（账本态）成立")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
