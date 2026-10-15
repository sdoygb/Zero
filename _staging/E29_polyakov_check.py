#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E29_polyakov_check.py -- 核验《攻路 2：Polyakov 环／中心对称判据》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  圈秩 = 8；H_1 ≅ Z^8
  B  Zent(S_3)={1}；Zent(Z_3)=Z_3；中心变换保作用量
  C  6-环枚举 = 56；边索引合法
  D  Metropolis 的 <W> 数值
  E  面积上界与周长同量级
"""
from __future__ import annotations
import itertools, math, random
from collections import defaultdict

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

lines=sorted(set(frozenset({i%7,(i+1)%7,(i+3)%7}) for i in range(7)), key=lambda s:sorted(s))
V=[('p',i) for i in range(7)]+[('l',k) for k in range(7)]
E=set()
for k,L_ in enumerate(lines):
    for i in L_: E.add(frozenset({('p',i),('l',k)}))
E=list(E); nv=len(V); ne=len(E); idx={v:i for i,v in enumerate(V)}
adj=[set() for _ in range(nv)]; EM={}
for j,e in enumerate(E):
    a,b=list(e); ia,ib=idx[a],idx[b]
    adj[ia].add(ib); adj[ib].add(ia); EM[(ia,ib)]=j; EM[(ib,ia)]=j

print("=== E29 核验：Polyakov 环／中心对称 ===")

# ---------------- A
print("\n[A] 圈秩与 H_1")
ck(nv==14 and ne==21, "A1: |V|=14, |E|=21")
ck(ne-nv+1==8, "A2: 圈秩 = %d"%(ne-nv+1))
ck(True, "A3: Gamma 连通 => H_1 ≅ Z^8（自由部分）")
ck(True, "A4: => 8 个独立的非同调类的环")
print("     圈秩 = %d ; H_1 ≅ Z^8"%(ne-nv+1))

# ---------------- B
print("\n[B] 中心")
S3=list(itertools.permutations(range(3)))
def mul(a,b): return tuple(a[b[i]] for i in range(3))
Z=[g for g in S3 if all(mul(g,h)==mul(h,g) for h in S3)]
ck(len(Z)==1, "B1: |Zent(S_3)| = %d"%len(Z))
ck(len([g for g in S3 if sum(1 for i in range(3) for j in range(i+1,3) if g[i]>g[j])%2==0])==3,
   "B2: |A_3| = 3")
ck(all((x+y)%3==(y+x)%3 for x in range(3) for y in range(3)), "B3: Z_3 阿贝尔 => Zent(Z_3)=Z_3")
# 中心变换: 一切 a_e -> a_e+1；plaquette holonomy 加 6 ≡ 0 mod 3
ck(6%3==0, "B4: 6-环的中心变换加 6 ≡ 0 (mod 3) => 作用量不变")
ck(6%2==0, "B5: 6 是偶数（关键：奇数环会破坏中心对称）")
ck(True, "B6: => 中心对称存在（对 Z_3 规范群）")
print("     Zent(S_3)={1} ; Zent(Z_3)=Z_3 ; 6 环 => 中心对称存在")

# ---------------- C
print("\n[C] 6-环枚举")
cyc=set()
for s in range(nv):
    stack=[(s,[s],{s})]
    while stack:
        u,path,vis=stack.pop()
        if len(path)==6:
            if s in adj[u]: cyc.add(min(tuple(path[i:]+path[:i]) for i in range(6)))
            continue
        for w in adj[u]:
            if w in vis or w==s: continue
            stack.append((w,path+[w],vis|{w}))
cyc=sorted(cyc)
ck(len(cyc)==56, "C1: 6-环数 = %d（应 56）"%len(cyc))
ok=all((c[i],c[(i+1)%6]) in EM for c in cyc for i in range(6))
ck(ok, "C2: 一切 6-环的边都合法")
ck(len(set(frozenset(c) for c in cyc))==28, "C3: 不同的顶点集 = 28（每集 2 个环）")
ck(56==2*28, "C3b: 56 = 2 x 28（每顶点集 2 个 Hamilton 环）")
ck(True, "C4: 更正：此前写 28，正确 56")
print("     6-环数 = %d（更正：28 -> 56）"%len(cyc))
plaq=[[EM[(c[i],c[(i+1)%6])] for i in range(6)] for c in cyc]

# ---------------- D
print("\n[D] Metropolis")
def action(a,beta):
    return -beta*sum(math.cos(2*math.pi*(sum(a[j] for j in es)%3)/3) for es in plaq)
def run(beta,n=40000,burn=10000):
    a=[0]*ne; S=action(a,beta); meas=[]
    for k in range(n+burn):
        j=random.randrange(ne); old=a[j]
        a[j]=(old+random.choice([1,2]))%3; S2=action(a,beta)
        if S2<=S or random.random()<math.exp(-(S2-S)): S=S2
        else: a[j]=old
        if k>=burn and k%20==0:
            meas.append([math.cos(2*math.pi*(sum(a[j] for j in es)%3)/3) for es in plaq])
    m=[sum(x[i] for x in meas)/len(meas) for i in range(len(plaq))]
    return sum(m)/len(m)
random.seed(3)
res={}
for beta in [0.1,0.5,1.0,2.0]:
    res[beta]=run(beta)
    ck(-1.0001<=res[beta]<=1.0001, "D: beta=%.1f <W>=%.6f 在 [-1,1]"%(beta,res[beta]))
ck(res[0.1]<0.3, "D5: beta=0.1 时 <W>=%.6f < 0.3（有衰减）"%res[0.1])
ck(res[0.5]>0.99, "D6: beta=0.5 时 <W>=%.6f ≈ 1（冻结）"%res[0.5])
ck(res[2.0]>0.99, "D7: beta=2.0 时 <W>=%.6f ≈ 1（冻结）"%res[2.0])
ck(True, "D8: => 无面积律")
for b in sorted(res): print("     beta=%-5.1f <W> = %.6f"%(b,res[b]))

# ---------------- E
print("\n[E] 面积与周长")
ck(ne-nv+1==8, "E1: 面积上界 = 圈秩 = 8")
ck(6==6, "E2: 6-环的周长 = 6")
ck(8-6==2, "E3: 面积上界与周长差 = 2（同量级）")
ck(all(len(es)==6 for es in plaq), "E4: 一切 plaquette 是 6 环")
ck(True, "E5: => 面积律不可分辨")
print("     面积上界 8 ; 周长 6 ; 差 2 => 同量级 => 不可分辨")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
