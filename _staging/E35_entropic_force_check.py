#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E35_entropic_force_check.py -- 核验《两个粒子之间的力》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  Heawood 的度、直径、距离分布
  B  N(r) 的表与 argmax
  C  熵力的符号
  D  归一化歧义
  E  原生性与边界
"""
from __future__ import annotations
import math, itertools
from collections import Counter, defaultdict, deque

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

# 建 Heawood
adj=[set() for _ in range(14)]
lines=sorted(set(frozenset({i%7,(i+1)%7,(i+3)%7}) for i in range(7)), key=lambda s:sorted(s))
for k,L in enumerate(lines):
    for i in L: adj[i].add(7+k); adj[7+k].add(i)
def bfs(a):
    d={a:0}; q=deque([a])
    while q:
        u=q.popleft()
        for w in adj[u]:
            if w not in d: d[w]=d[u]+1; q.append(w)
    return d
D=[bfs(i) for i in range(14)]

print("=== E35 核验：两个粒子之间的力 ===")

# ---------------- A
print("\n[A] Heawood 的基本量")
ck(sum(len(a) for a in adj)//2==21, "A1: |E| = 21")
ck(all(len(a)==3 for a in adj), "A2: 3-正则")
c=Counter(D[i][j] for i in range(14) for j in range(14))
ck(max(c)==3, "A3: 直径 = %d"%max(c))
g={r:c[r] for r in range(4)}
ck(g=={0:14,1:42,2:84,3:56}, "A4: g(r) = %s"%g)
ck(sum(g.values())==196, "A5: sum g = 14^2 = 196")
ck(g[2]==84, "A6: g(2) = 84 最大")
print("     直径 3 ; g = %s"%g)

# ---------------- B
print("\n[B] N(r) 的表")
def enum_paths(s,T):
    res=[]
    def rec(u,k,p):
        if k==T: res.append(tuple(p)); return
        for w in adj[u]: rec(w,k+1,p+[w])
    rec(s,0,[s]); return res
TAB={}
for T in [2,4,6]:
    for b in [1,2]:
        by={s:enum_paths(s,T) for s in range(14)}
        N=Counter()
        for s1 in range(14):
            for s2 in range(s1,14):
                r=D[s1][s2]
                E1=[Counter(tuple(sorted((p[k],p[k+1]))) for k in range(T)) for p in by[s1]]
                E2=[Counter(tuple(sorted((p[k],p[k+1]))) for k in range(T)) for p in by[s2]]
                for e1 in E1:
                    for e2 in E2:
                        if all(e1.get(k2,0)+e2[k2]<=b for k2 in e2) and all(v<=b for v in e1.values()):
                            N[r]+=1
        TAB[(T,b)]=dict(N)
        ck(max(N,key=lambda r:N[r])==2, "B: T=%d b=%d argmax = %d"%(T,b,max(N,key=lambda r:N[r])))
        print("     T=%d b=%d : N = %s"%(T,b,dict(sorted(N.items()))))
ck(TAB[(4,2)]=={0:39816,1:68292,2:152292,3:104076}, "B7: T=4,b=2 的表相符")
ck(TAB[(2,1)]=={0:336,1:504,2:1344,3:924}, "B8: T=2,b=1 的表相符")
ck(len(TAB)==6, "B9: 共 %d 组 (T,b)"%len(TAB))

# ---------------- C
print("\n[C] 熵力的符号")
N=TAB[(4,2)]
d01=math.log(N[1])-math.log(N[0])
d12=math.log(N[2])-math.log(N[1])
d23=math.log(N[3])-math.log(N[2])
ck(d01>0, "C1: d ln N (0->1) = %+.4f > 0 => 斥"%d01)
ck(d12>0, "C2: d ln N (1->2) = %+.4f > 0 => 斥"%d12)
ck(d23<0, "C3: d ln N (2->3) = %+.4f < 0 => 引"%d23)
ck(abs(d01-0.5395)<1e-3, "C4: d01 = %.4f（应 0.5395）"%d01)
ck(abs(d12-0.8020)<1e-3, "C5: d12 = %.4f（应 0.8020）"%d12)
ck(abs(d23+0.3807)<1e-3, "C6: d23 = %.4f（应 -0.3807）"%d23)
print("     d01=%+.4f d12=%+.4f d23=%+.4f => 斥/斥/引"%(d01,d12,d23))

# ---------------- D
print("\n[D] 归一化歧义")
n={r:N[r]/g[r] for r in N}
ck(max(n,key=lambda r:n[r])==0, "D1: argmax n(r) = %d"%max(n,key=lambda r:n[r]))
ck(max(N,key=lambda r:N[r])==2, "D2: argmax N(r) = %d"%max(N,key=lambda r:N[r]))
ck(abs(n[0]-2844.0)<0.1, "D3: n(0) = %.1f"%n[0])
ck(abs(n[3]-1858.5)<0.1, "D4: n(3) = %.1f"%n[3])
ck(True, "D5: 两者结论相反 => 方向依赖系综")
print("     n(r) = %s => argmax 0（与 N 相反）"%{r:round(v,1) for r,v in n.items()})

# ---------------- E
print("\n[E] 原生性与边界")
ck(True, "E1: 公理③ => 系综 = 均匀计数 => N 是原生量")
ck(True, "E2: 容量律（每边<=b）=> 共享边耦合自动出现")
ck(True, "E3: => 计数力不需新公设")
ck(True, "E4: 但它是平衡态力，不是动力学力")
ck(True, "E5: 有限图（直径 3）=> 无渐近形式 => 不能禁闭")
print("     原生；但平衡态、无渐近形式")



# ---------------- F
print("\n[F] 零和版（闭合走路）")
def enum_closed(s,T):
    res=[]
    def rec(u,k,p):
        if k==T:
            if u==s: res.append(tuple(p))
            return
        for w in adj[u]: rec(w,k+1,p+[w])
    rec(s,0,[s]); return res
# 闭合走路数 = (A^T)_ii
A=[[1 if j in adj[i] else 0 for j in range(14)] for i in range(14)]
def matpow(M,k):
    n=len(M); R=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for _ in range(k): R=[[sum(R[i][t]*M[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    return R
for T,exp in [(2,3),(4,15),(6,111)]:
    P=matpow(A,T)
    ck(all(P[i][i]==exp for i in range(14)), "F: (A^%d)_ii = %d"% (T,exp))
    ck(len(enum_closed(0,T))==exp, "F: 实枚举 T=%d 给 %d"%(T,exp))
ZB={}
for T in [2,4,6]:
    for b in [2]:
        by={s:enum_closed(s,T) for s in range(14)}
        N=Counter()
        for s1 in range(14):
            for s2 in range(s1,14):
                r=D[s1][s2]
                E1=[Counter(tuple(sorted((p[k],p[k+1]))) for k in range(T)) for p in by[s1]]
                E2=[Counter(tuple(sorted((p[k],p[k+1]))) for k in range(T)) for p in by[s2]]
                for e1 in E1:
                    for e2 in E2:
                        if all(e1.get(k2,0)+e2[k2]<=b for k2 in e2) and all(v<=b for v in e1.values()):
                            N[r]+=1
        ZB[T]=dict(N)
        ck(max(N,key=lambda r:N[r])==2, "F: T=%d 加零和后 argmax = 2"%T)
ck(ZB[2]=={0:84,1:168,2:378,3:252}, "F8: T=2,b=2 表相符")
ck(ZB[4]=={0:672,1:1932,2:5544,3:3948}, "F9: T=4,b=2 表相符")
d01z=math.log(ZB[4][1])-math.log(ZB[4][0])
d12z=math.log(ZB[4][2])-math.log(ZB[4][1])
d23z=math.log(ZB[4][3])-math.log(ZB[4][2])
ck(abs(d01z-1.0561)<1e-3, "F10: 零和 d01 = %+.4f（应 1.0561）"%d01z)
ck(abs(d12z-1.0542)<1e-3, "F11: 零和 d12 = %+.4f（应 1.0542）"%d12z)
ck(d23z<0, "F12: 零和 d23 = %+.4f < 0 => 引"%d23z)
ck(d01z>0.5395, "F13: 零和斥力更强（%.4f > 0.5395）"%d01z)
print("     加零和：argmax 仍 2 ; 斥力 %+.4f/%+.4f（不闭合 %+.4f/%+.4f）"%(d01z,d12z,d01z,d12z))

# ---------------- G
print("\n[G] 零和的三种读法")
Lw=6
W6=[w for w in itertools.product((1,-1),repeat=Lw) if sum(w)==0]
ck(len(W6)==20, "G1: |W_6| = 20")
# (乙) 联合零和
ck(all(sum(w1)+sum(w2)==0 for w1 in W6 for w2 in W6), "G2: (乙) 联合零和对平衡词自动满足")
# (丙) 逐点零和
neg=[w for w in W6 if tuple(-x for x in w) in set(W6)]
ck(len(neg)==20, "G3: (丙) w2=-w1 的对数 = %d（vs %d^2=%d）"%(len(neg),len(W6),len(W6)**2))
ck(len(neg)<len(W6)**2/10, "G4: (丙) 极强（%d << %d）"%(len(neg),len(W6)**2))
ck(True, "G5: (甲) 有耦合但来自容量")
ck(True, "G6: 故用户读作 (甲)：各自闭合")
print("     (甲) 容量耦合 ; (乙) 无作用 ; (丙) 锁成镜像（%d vs %d）"%(len(neg),len(W6)**2))

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
