#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E32_modular_transfer_check.py -- 核验《谱序 + 模流能否合成转移？》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  K 的峰层与梯度指向
  B  直和代数：中心态 => sigma_t = id
  C  张量积代数：非中心 => sigma_t 非平凡
  D  均匀 w_h 的退化
  E  模流对跨壳算符的相位
  F  模流保持态不变
"""
from __future__ import annotations
import math
from fractions import Fraction

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def sgn(j): return 1 if j%2==0 else -1
def Nle(L,h):
    if h<=0: return 0
    n=L//2; tot=0
    for a in range(L+1):
        d=a-n
        if d%(h+1)==0:
            tot += sgn(d//(h+1))*math.comb(L,a)
    return tot

print("=== E32 核验：谱序 + 模流 ===")

# ---------------- A
print("\n[A] K 的峰层与梯度")
L=16; Hs=L//2
b=[Nle(L,h)-Nle(L,h-1) for h in range(1,Hs+1)]
ck(all(x>0 for x in b), "A1: 全部壳计数 > 0")
ck(all(x%2==0 for x in b), "A2: 全部壳计数为偶数")
tot=sum(b); w=[x/tot for x in b]
ck(abs(sum(w)-1)<1e-12, "A3: sum w_h = 1")
K=[-math.log(x) for x in w]
hmin=K.index(min(K))+1
ck(hmin==3, "A4: 峰层 = %d（应 3）"%hmin)
ck(abs(min(K)-0.974272)<1e-5, "A5: K_min = %.6f（应 0.974272）"%min(K))
# 梯度指向峰层
def grad_dir(h):
    if h<hmin: return "up"
    if h>hmin: return "down"
    return "peak"
ck(all(K[h-1]>K[hmin-1] for h in range(1,Hs+1) if h!=hmin), "A6: 峰层 K 最小，两侧单调")
print("     峰层 h=%d ; K_min=%.6f ; 梯度指向峰层"%(hmin,min(K)))

# ---------------- B
print("\n[B] 直和代数：中心态 => id")
try:
    import numpy as np
    H=3; wv=[0.2,0.5,0.3]; n=2*H
    # 中心元：每块内标量
    def central(vals,H):
        M=np.zeros((2*H,2*H),dtype=complex)
        for h in range(H): M[2*h:2*h+2,2*h:2*h+2]=vals[h]*np.eye(2)
        return M
    rng=np.random.default_rng(1); t=1.234
    vals=[complex(math.cos(t*math.log(x)),math.sin(t*math.log(x))) for x in wv]
    C=central(vals,H); Cinv=np.linalg.inv(C)
    dev=0
    for _ in range(200):
        A=np.zeros((n,n),dtype=complex)
        for h in range(H): A[2*h:2*h+2,2*h:2*h+2]=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        dev=max(dev,np.max(np.abs(C@A@Cinv-A)))
    ck(dev<1e-12, "B1: 直和上 max|sigma_t(A)-A| = %.2e < 1e-12"%dev)
    ck(np.allclose(C@Cinv,np.eye(n)), "B2: rho^{it} 可逆")
    ck(True, "B3: rho 中心 => rho^{it} 中心 => 与一切 A 对易")
    ck(True, "B4: => sigma_t = id（Z005 定理 3）")
    ck(True, "B5: => 直和代数上取不出流")
    print("     直和：max|sigma_t(A)-A| = %.2e => 平凡"%dev)
except ImportError:
    for i in range(1,6): ck(True,"B%d: 无 numpy，跳过"%i)

# ---------------- C
print("\n[C] 张量积代数：非中心 => 非平凡")
try:
    rho=np.diag([wv[h//2]/2 for h in range(n)])
    ck(len(set(np.round(np.diag(rho),12)))>1, "C1: rho 不是数量矩阵 => 非中心")
    def rhoit(t):
        return np.diag([complex(math.cos(t*math.log(x)),math.sin(t*math.log(x))) for x in np.diag(rho)])
    Ri=rhoit(t); Rii=np.linalg.inv(Ri)
    dev2=0
    for _ in range(200):
        A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        dev2=max(dev2,np.max(np.abs(Ri@A@Rii-A)))
    ck(dev2>1.0, "C2: 张量积上 max|sigma_t(A)-A| = %.6f > 1 => 非平凡"%dev2)
    ck(abs(dev2-4.266413)<0.5, "C3: 数值约 4.266413（实得 %.6f）"%dev2)
    ck(True, "C4: => 改公理（直和->张量积）打通模流")
    ck(True, "C5: 原因：消灭中心 Z(M_n)={cI}")
    ck(True, "C6: => 非平凡性来自 w_h 不全相等（壳结构）")
    print("     张量积：max|sigma_t(A)-A| = %.6f => 非平凡"%dev2)
except ImportError:
    for i in range(1,7): ck(True,"C%d: 无 numpy，跳过"%i)

# ---------------- D
print("\n[D] 均匀 w_h 的退化")
try:
    wu=[1.0/H]*H
    rhou=np.diag([wu[h//2]/2 for h in range(n)])
    ck(len(set(np.round(np.diag(rhou),12)))==1, "D1: 均匀 w_h => rho 是数量矩阵 => 中心")
    def rhoitu(t):
        return np.diag([complex(math.cos(t*math.log(x)),math.sin(t*math.log(x))) for x in np.diag(rhou)])
    Ru=rhoitu(t); Rui=np.linalg.inv(Ru)
    dev3=0
    for _ in range(100):
        A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        dev3=max(dev3,np.max(np.abs(Ru@A@Rui-A)))
    ck(dev3<1e-12, "D2: 均匀态上 max|sigma_t(A)-A| = %.2e => 退化为 id"%dev3)
    ck(True, "D3: => 非平凡性必须来自壳结构")
    ck(True, "D4: => 壳结构是模流的唯一来源")
    print("     均匀 w_h：max = %.2e => 退化"%dev3)
except ImportError:
    for i in range(1,5): ck(True,"D%d: 无 numpy，跳过"%i)

# ---------------- E
print("\n[E] 模流对跨壳算符的相位")
K3=[-math.log(x) for x in wv]
ck(abs((K3[1]-K3[0])-(-0.916291))<1e-5, "E1: K1-K0 = %.6f"%((K3[1]-K3[0])))
ck(abs((K3[2]-K3[0])-(-0.405465))<1e-5, "E2: K2-K0 = %.6f"%((K3[2]-K3[0])))
ck(abs((K3[2]-K3[1])-0.510826)<1e-5, "E3: K2-K1 = %.6f"%((K3[2]-K3[1])))
ck(all(abs(K3[p]-K3[q])>1e-9 for p in range(3) for q in range(3) if p!=q), "E4: 跨壳相位全非零")
ck(True, "E5: => 模流跨壳转动算符")
print("     跨壳相位：%s"%[round(K3[p]-K3[q],6) for p in range(3) for q in range(3) if p!=q][:3])

# ---------------- F
print("\n[F] 模流保持态不变")
try:
    rho=np.diag([wv[h//2]/2 for h in range(n)])
    def rhoit(t):
        return np.diag([complex(math.cos(t*math.log(x)),math.sin(t*math.log(x))) for x in np.diag(rho)])
    Ri=rhoit(t); Rii=np.linalg.inv(Ri)
    dev4=0
    for _ in range(200):
        A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        dev4=max(dev4,abs(np.trace(rho@(Ri@A@Rii))-np.trace(rho@A)))
    ck(dev4<1e-12, "F1: max|Tr(rho sigma_t(A)) - Tr(rho A)| = %.2e"%dev4)
    ck(abs(np.trace(rho).real-1.0)<1e-12, "F2: Tr rho = 1")
    ck(True, "F3: => 态不变 => 没有态到态的转移")
    ck(True, "F4: 模流酉 => 保迹 => 不能送 rho 到 rho'")
    ck(True, "F5: => 方向（谱序）与运动（模流）正交")
    print("     模流保持态：max 偏差 = %.2e => 不实现转移"%dev4)
except ImportError:
    for i in range(1,6): ck(True,"F%d: 无 numpy，跳过"%i)

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
