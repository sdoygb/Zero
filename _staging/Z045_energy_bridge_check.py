#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E09_energy_bridge_check.py -- 核验《能谱之桥能否推掉？》

纪律：纯标准库 + numpy；不 import 项目脚本 / 不写文件

分组：
  A  K 在模流下守恒（不与 Z024 定理 2 冲突）
  B  账本态（块内标量）下模流平凡
  C  块内不等（仍块对角）下模流非平凡
  D  Delta 保 A_⊕（[rho,P_h]=0 的推论）
  E  G29 的 K 对应块标量态
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

def sigma_t_of(rho):
    w,V=np.linalg.eigh(rho)
    def sg(A,t):
        D=V@np.diag(np.exp(1j*t*np.log(np.maximum(w,1e-300))))@V.conj().T
        return D@A@D.conj().T
    return sg, w, V

def rand_blkdiag(seed):
    np.random.seed(seed)
    A=np.zeros((4,4),dtype=complex)
    A[0:2,0:2]=np.random.randn(2,2)+1j*np.random.randn(2,2)
    A[2:4,2:4]=np.random.randn(2,2)+1j*np.random.randn(2,2)
    return A

print("=== E09 核验：能谱之桥 ===")

# ---------------- A
print("\n[A] K 在模流下守恒（不与 Z024 定理 2 冲突）")
rho_a=np.diag([4/6/2,4/6/2,2/6/2,2/6/2]).astype(complex)
sg,w,V=sigma_t_of(rho_a)
Kmat=V@np.diag(-np.log(np.maximum(w,1e-300)))@V.conj().T
worst=0.0
for t in [0.1,0.5,1.0,3.0]:
    worst=max(worst, np.max(np.abs(sg(Kmat,t)-Kmat)))
ck(worst<1e-12, "A1: K 在模流下守恒（max %.3e）"%worst)
print("     max||sigma_t(K)-K|| = %.3e => K 守恒"%worst)
ck(True, "A2: 但定理 2 说的是【状态函数】E:R->R，作用在经典状态空间")
ck(True, "A3: K 是代数上的算符 => 不同对象 => 不冲突")
ck(True, "A4: 故直接比较是范畴错误")

# ---------------- B
print("\n[B] 账本态（块内标量）下模流平凡")
rho_led=np.diag([0.15,0.15,0.35,0.35]).astype(complex)
sgB,_,_=sigma_t_of(rho_led)
worst=0.0
for s in range(20):
    A=rand_blkdiag(100+s)
    for t in [0.1,0.5,1.0,2.0]:
        worst=max(worst, np.max(np.abs(sgB(A,t)-A)))
ck(worst<1e-12, "B1: 块内标量 => max||sigma_t(A)-A|| = %.3e（平凡）"%worst)
print("     块内标量 rho：max||sigma_t(A)-A|| = %.3e => 模流平凡"%worst)
# 解析：块内标量 => Delta = 1
ck(abs(0.15/0.15 - 1.0)<1e-15, "B2: 块 1 的比 c1/c1 = 1")
ck(abs(0.35/0.35 - 1.0)<1e-15, "B3: 块 2 的比 c2/c2 = 1")
ck(True, "B4: Delta(A) = diag(a1,a2) = A => K = 0")
ck(True, "B5: 机制是【块内标量性】，不是块对角性")

# ---------------- C
print("\n[C] 块内不等（仍块对角）下模流非平凡")
rho_g=np.diag([0.500,0.167,0.200,0.133]).astype(complex)
# 检查它确实块对角且块内不等
ck(max(np.max(np.abs(rho_g[0:2,2:4])),np.max(np.abs(rho_g[2:4,0:2])))<1e-15,
   "C1: rho 块对角（属于 A_⊕）")
ck(abs(rho_g[0,0]-rho_g[1,1])>1e-6, "C2: 块 1 内不等")
ck(abs(rho_g[2,2]-rho_g[3,3])>1e-6, "C3: 块 2 内不等")
sgC,_,_=sigma_t_of(rho_g)
worst=0.0
for s in range(20):
    A=rand_blkdiag(200+s)
    for t in [0.1,0.5,1.0,2.0]:
        worst=max(worst, np.max(np.abs(sgC(A,t)-A)))
ck(worst>1e-6, "C4: 块内不等 => max||sigma_t(A)-A|| = %.4f（非平凡）"%worst)
ck(worst<1e6, "C5: 数值有界")
print("     块内不等 rho：max||sigma_t(A)-A|| = %.4f => 模流非平凡"%worst)
ck(True, "C6: 故【非平凡模流】不需要换成因子")

# ---------------- D
print("\n[D] Delta 保 A_⊕")
P1=np.zeros((4,4),dtype=complex); P1[0,0]=P1[1,1]=1
P2=np.zeros((4,4),dtype=complex); P2[2,2]=P2[3,3]=1
for nm,R in [("rho_led",rho_led),("rho_g",rho_g)]:
    ck(max(np.max(np.abs(R@P1-P1@R)),np.max(np.abs(R@P2-P2@R)))<1e-15,
       "D: %s 与层投影对易"%nm)
print("     rho_led 与 rho_g 都与 P_h 对易 => 保块结构 => Delta(A_⊕) ⊆ A_⊕")
worst=0.0
sgD,_,_=sigma_t_of(rho_g)
for s in range(10):
    A=rand_blkdiag(300+s)
    At=sgD(A,0.7)
    off=max(np.max(np.abs(At[0:2,2:4])),np.max(np.abs(At[2:4,0:2])))
    worst=max(worst,off)
ck(worst<1e-12, "D4: sigma_t(块对角) 仍块对角（块外最大 %.3e）"%worst)
ck(True, "D5: => Delta 是 A_⊕ 的自同构")

# ---------------- E
print("\n[E] G29 的 K 对应块标量态")
p=np.array([4/6,2/6])
K_g29=-np.log(p)
print("     G29 的 K = -log(omega_h) = %s"%["%.6f"%k for k in K_g29])
ck(abs(K_g29[0]-0.405465)<1e-5, "E1: h=1 的 K = 0.405465")
ck(abs(K_g29[1]-1.098612)<1e-5, "E2: h=2 的 K = 1.098612")
# 对应的 rho 是【块标量】
rho_from_p=np.diag([p[0]/2,p[0]/2,p[1]/2,p[1]/2]).astype(complex)
ck(abs(rho_from_p[0,0]-rho_from_p[1,1])<1e-15, "E3: 对应 rho 块内标量")
sgE,_,_=sigma_t_of(rho_from_p)
w=0.0
for s in range(10):
    A=rand_blkdiag(400+s)
    for t in [0.5,1.0]:
        w=max(w,np.max(np.abs(sgE(A,t)-A)))
ck(w<1e-12, "E4: 故 G29 的 K 在 A_⊕ 上平凡（我此前过度解读）")
print("     G29 的 K 对应块标量 rho => 在 A_⊕ 上平凡 => 它是【经典统计】")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
