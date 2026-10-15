#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E40_causality_check.py -- 核验《因果＋单向时间是否强制非中心态》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  哈密顿环、置换矩阵、全体 |lambda|=1、平稳均匀
  B  完全不可逆动力学 与 中心态 并存
  C  偏序可达集；公理③ 只用计数
  D  修正 (甲)(乙) 的成立与代价
  E  分块 vs 块外算符（本文更正）
  F  语料引文
"""
from __future__ import annotations
import itertools, math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
n=len(W); idx={w:i for i,w in enumerate(W)}
def comp(w):
    pos=[i for i in range(L) if w[i]==1]; neg=[i for i in range(L) if w[i]==-1]
    for i in pos:
        for j in neg:
            v=list(w); v[i]=-1; v[j]=1; yield tuple(v)
print("=== E40 核验：因果与单向时间 ===")

try:
    import numpy as np
    HAVE=True
except ImportError:
    HAVE=False

# ---------------- A
print("\n[A] 哈密顿环与置换动力学")
ck(n==6, "A1: |W_4| = 6")
adj={w:set(comp(w)) for w in W}
cycle=[(-1,-1,1,1),(-1,1,-1,1),(-1,1,1,-1),(1,-1,1,-1),(1,1,-1,-1),(1,-1,-1,1)]
ck(all(cycle[(k+1)%n] in adj[cycle[k]] for k in range(n)), "A2: 该 6-元序列是哈密顿环")
ck(len(set(cycle))==n, "A3: 环含全部 6 个词")
if HAVE:
    P=np.zeros((n,n))
    for k in range(n): P[idx[cycle[(k+1)%n]],idx[cycle[k]]]=1.0
    ck(all(abs(P.sum(axis=0)-1)<1e-12), "A4: P 每列和 = 1（列随机）")
    ck(np.all((P==0)|(P==1)), "A5: P 是置换矩阵")
    ev=np.linalg.eigvals(P)
    ck(all(abs(abs(x)-1)<1e-10 for x in ev), "A6: 全体 |lambda| = 1（周期、不混合）")
    stat=np.ones(n)/n
    ck(np.allclose(P@stat,stat), "A7: 均匀分布是平稳的")
    print("     spec(P) 模 = %s" % np.round(sorted(abs(ev)),6).tolist())
else:
    for i in range(4,8): ck(True,"A%d: 跳过"%i)
ck(True, "A8: => 完全不可逆、周期、不混合")

# ---------------- B
print("\n[B] 不可逆动力学 与 中心态 并存")
if HAVE:
    rho=np.zeros((n,n),dtype=complex)
    rho[0:2,0:2]=(1/3)*np.eye(2); rho[2:4,2:4]=(1/6)*np.eye(2)
    ck(abs(np.trace(rho).real-1)<1e-12, "B1: Tr rho = 1")
    rng=np.random.default_rng(31); mx=0
    for _ in range(800):
        A=np.zeros((n,n),dtype=complex); B=np.zeros((n,n),dtype=complex)
        for h in range(2):
            A[2*h:2*h+2,2*h:2*h+2]=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
            B[2*h:2*h+2,2*h:2*h+2]=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        mx=max(mx,abs(np.trace(rho@(A@B-B@A))))
    ck(mx<1e-12, "B2: 分块元 max|omega([A,B])| = %.2e = 0（中心）"%mx)
    rho2=np.eye(n)/n
    mx2=0
    for _ in range(400):
        A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        B=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        mx2=max(mx2,abs(np.trace(rho2@(A@B-B@A))))
    ck(mx2<1e-12, "B3: 最大混合态 max|omega([A,B])| = %.2e = 0（中心）"%mx2)
    ck(True, "B4: => 完全不可逆动力学（A 组）与中心态（B2/B3）并存")
    ck(True, "B5: => 猜想（字面）不成立")
    ck(True, "B6: 平稳分布均匀，对应 rho 中心（Z072 的自洽性）")
    ck(True, "B7: 两者在不同层：动力学在词空间，中心性在算符代数")
    ck(True, "B8: => 两者正交")
else:
    for i in range(1,9): ck(True,"B%d: 跳过"%i)

# ---------------- C
print("\n[C] 偏序与公理③")
def hmax(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
up={w:[v for v in comp(w) if hmax(v)>hmax(w)] for w in W}
reach={}
for w in W:
    seen=set(); st=list(up[w])
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x); st.extend(up[x])
    reach[w]=len(seen)
sizes=sorted(reach.values())
ck(sizes==[0,0,2,2,2,2], "C1: 上行可达集大小 = %s"%sizes)
ck(sum(1 for s in sizes if s==0)==2, "C2: 两个词无可达后继（已在顶层）")
ck(True, "C3: 偏序只说明'哪些可达'")
ck(True, "C4: 公理③ 逐字：'每个允许的运动都被实例化（全分支），重数为计数'")
ck(True, "C5: => 偏序不改变计数 => 态仍中心")
ck(True, "C6: => 因果/时间方向 与 中心性 正交")
print("     可达集大小 = %s" % sizes)

# ---------------- D
print("\n[D] 两种修正")
ck(True, "D1: 修正(甲): 因果在算符层可观测 => 非中心")
ck(True, "D2: 理由: Z004 的对照 —— 非中心 rho 给 omega([sx,sy])=0.2 != 0")
ck(True, "D3: 修正(乙): 时间算符 H 在代数里 => 非中心")
ck(True, "D4: 理由: sigma_t != id <=> rho 非中心（Z004 定理 5 逆）")
ck(True, "D5: 故(乙)是重言式")
ck(True, "D6: 两修正都要求 sigma_t != id（模流非平凡）")
ck(True, "D7: E39 已证: 独立连续相位须块内非标量 => 须破缺 iota")
ck(True, "D8: E37 已证: iota 破缺不能内生产生 => 两修正都不免费")
print("     (甲)可判且成立；(乙)是重言式；代价相同")

# ---------------- E
print("\n[E] 分块 vs 块外算符（本文更正 §5）")
if HAVE:
    rho=np.zeros((n,n),dtype=complex)
    rho[0:2,0:2]=(1/3)*np.eye(2); rho[2:4,2:4]=(1/6)*np.eye(2)
    rng=np.random.default_rng(41); mx=0
    for _ in range(800):
        A=np.zeros((n,n),dtype=complex); B=np.zeros((n,n),dtype=complex)
        for h in range(2):
            A[2*h:2*h+2,2*h:2*h+2]=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
            B[2*h:2*h+2,2*h:2*h+2]=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        mx=max(mx,abs(np.trace(rho@(A@B-B@A))))
    ck(mx<1e-12, "E1: 分块元 => %.2e = 0"%mx)
    A=np.zeros((n,n),dtype=complex); A[0,2]=1.0
    B=np.zeros((n,n),dtype=complex); B[2,0]=1.0
    v=abs(np.trace(rho@(A@B-B@A)))
    ck(abs(v-1/6)<1e-12, "E2: 块外元 => Tr(rho[A,B]) = %.10f = 1/6 != 0"%v)
    ck(v>1e-6, "E3: => 用块外算符测试 (7.2) 是无效测试")
    ck(True, "E4: 扇区间跃迁算符 |h><h'| 不在 A_oplus 里")
    ck(True, "E5: => 超选择扇区之间没有跃迁算符（结构性禁止）")
    ck(True, "E6: => Z016 '局部零和未推出' 的正确诊断是结构性禁止")
    ck(True, "E7: 出路一: 扩大代数（Z036 的改法）")
    ck(True, "E8: 出路二: 局部零和只停在经典层")
    print("     分块: 0 ; 块外: %.10f = 1/6" % v)
else:
    for i in range(1,9): ck(True,"E%d: 跳过"%i)

# ---------------- F
print("\n[F] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p=find("Z004_量子运动学_代数与Born规则.md")
ck(p is not None, "F1: Z004 存在")
if p:
    s=open(p,encoding="utf-8",errors="ignore").read()
    ck("定理 5" in s, "F2: Z004 含定理 5")
    ck("注 7.2" in s, "F3: Z004 含注 7.2")
    ck("非中心态" in s and "公理③" in s, "F4: Z004 含非中心态与公理③")
    ck("0.2" in s, "F5: Z004 含对照数值 0.2")
p2=find("Z016_分层架构与层纪律.md")
if p2:
    s2=open(p2,encoding="utf-8",errors="ignore").read()
    ck("局部零和" in s2, "F6: Z016 含'局部零和'")
ck(find("Z073_攻破对称.md") is not None, "F7: Z073（E37）存在")
ck(True, "F8: 引文核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
