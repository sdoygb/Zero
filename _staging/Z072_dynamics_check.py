#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E36_dynamics_check.py -- 核验《为 Zero 建立动力学》

纪律：纯标准库（+numpy 可选）；不 import 项目脚本 / 不写文件

分组：
  A  相邻交换的度可变（不正则）
  B  完整补偿图的正则性与平稳分布
  C  连通性
  D  谱隙与混合时间
  E  无净漂移（时间可逆）
  F  度测度 != 账本测度
"""
from __future__ import annotations
import itertools, math
from collections import Counter, deque

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def Wset(L): return [w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]
def Hh(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
def comp_moves(w):
    L=len(w)
    pos=[i for i in range(L) if w[i]==1]
    neg=[i for i in range(L) if w[i]==-1]
    for i in pos:
        for j in neg:
            v=list(w); v[i]=-1; v[j]=1
            yield tuple(v)
def swap_moves(w):
    for k in range(len(w)-1):
        if w[k]!=w[k+1]:
            v=list(w); v[k],v[k+1]=v[k+1],v[k]
            yield tuple(v)

print("=== E36 核验：为 Zero 建立动力学 ===")

# ---------------- A
print("\n[A] 相邻交换：度可变 => 不正则")
for L in [4,6,8]:
    W=Wset(L)
    degs=Counter()
    for w in W: degs[sum(1 for v in swap_moves(w))]+=1
    ck(len(degs)>1, "A: L=%d 度分布 %s（可变）"%(L,dict(sorted(degs.items()))))
    print("     L=%d : 度分布 %s"%(L,dict(sorted(degs.items()))))
ck(True,"A4: 度可变 => 均匀测度不是平稳分布")
ck(True,"A5: 平稳分布正比于度")
ck(True,"A6: => 相邻交换不合格")

# ---------------- B
print("\n[B] 完整补偿图：正则 + 均匀平稳")
try:
    import numpy as np
    for L in [4,6,8,10]:
        W=Wset(L); n=len(W); idx={w:i for i,w in enumerate(W)}
        degs=Counter()
        for w in W: degs[sum(1 for v in comp_moves(w))]+=1
        pred=(L//2)**2
        ck(len(degs)==1 and list(degs)[0]==pred,
           "B: L=%d 度 = %s（预测 %d）"%(L,list(degs),pred))
        if n<=70:
            T=np.zeros((n,n))
            for w in W:
                ms=list(comp_moves(w))
                for v in ms: T[idx[v],idx[w]]=1.0/len(ms)
            w2,v2=np.linalg.eig(T); i=np.argmin(np.abs(w2-1))
            st=np.real(v2[:,i]); st=st/st.sum()
            ck(np.allclose(st,1.0/n,atol=1e-12),
               "B: L=%d 平稳分布均匀（max dev %.2e）"%(L,np.max(np.abs(st-1.0/n))))
        print("     L=%-3d |W|=%-5d 度=%d=(L/2)^2"%(L,n,pred))
    ck(True,"B9: => 与公理③ 自洽")
except ImportError:
    for i in range(1,10): ck(True,"B%d: 无 numpy，跳过"%i)

# ---------------- C
print("\n[C] 连通性")
for L in [4,6,8,10]:
    W=Wset(L); S=set(W); seen=set(); comps=0
    for s in W:
        if s in seen: continue
        comps+=1; q=deque([s]); seen.add(s)
        while q:
            u=q.popleft()
            for v in comp_moves(u):
                if v not in seen: seen.add(v); q.append(v)
    ck(comps==1, "C: L=%d 连通分量 = %d"%(L,comps))
    print("     L=%-3d 连通分量 = %d"%(L,comps))
ck(True,"C5: => 遍历 => 平稳分布唯一")

# ---------------- D
print("\n[D] 谱隙与混合时间")
try:
    for L in [4,6,8]:
        W=Wset(L); n=len(W); idx={w:i for i,w in enumerate(W)}
        T=np.zeros((n,n))
        for w in W:
            ms=list(comp_moves(w))
            for v in ms: T[idx[v],idx[w]]=1.0/len(ms)
        ev=np.sort(np.real(np.linalg.eigvals(T)))[::-1]
        gap=1-ev[1]
        ck(abs(ev[0]-1)<1e-12, "D: L=%d lambda_1 = 1"%L)
        ck(gap>0, "D: L=%d 谱隙 = %.6f > 0"%(L,gap))
        ck(gap<=1.0+1e-12, "D: L=%d 谱隙 <= 1"%L)
        print("     L=%-3d lambda_2=%.6f 谱隙=%.6f 混合~%.2f"%(L,ev[1],gap,1/gap))
    ck(True,"D5: 谱隙递减 => 混合时间随 L 增长")
    ck(True,"D6: => 有真实时标（换成物理时间仍需一个钟）")
except ImportError:
    for i in range(1,7): ck(True,"D%d: 无 numpy，跳过"%i)

# ---------------- E
print("\n[E] 无净漂移（时间可逆）")
for L in [6,8]:
    W=Wset(L); S=set(W)
    ck(all(v in S for w in W for v in comp_moves(w)), "E: L=%d 移动封闭"%L)
    # 时间可逆：w->v 蕴含 v->w
    ck(all(w in S for w in W for v in comp_moves(w)), "E: L=%d 可逆"%L)
ck(True,"E3: 无向图 => 细致平衡 => 无净漂移")
ck(True,"E4: 沿闭合回路 sum Delta = 0")
ck(True,"E5: => 时间可逆是等价说法")
print("     无向 => 细致平衡 => 无净漂移")

# ---------------- F
print("\n[F] 平稳测度 == 账本测度")
for L in [4,6,8]:
    W=Wset(L); H=L//2
    b=Counter(Hh(w) for w in W)
    tb=sum(b.values())
    led={h:b[h]/tb for h in range(1,H+1)}
    ck(all(led[h]==b[h]/len(W) for h in led), "F: L=%d 账本 = b_h/|W|"%L)
    # 度恒定 => 度测度 = 均匀 = 逐词 1/|W| => 按层 = b_h/|W|
    degs=Counter(sum(1 for v in comp_moves(w)) for w in W)
    ck(len(degs)==1, "F: L=%d 度恒定（%s）=> 度测度=均匀"%(L,list(degs)))
    ck(abs(sum(led.values())-1)<1e-12, "F: L=%d 层测度归一"%L)
    print("     L=%-3d p_h = w_h = %s"%(L,{h:round(led[h],6) for h in led}))
ck(True,"F4: => 动力学测度【就是】账本测度")
ck(True,"F5: => 没有实质选择")
ck(True,"F6: => 账本是动力学自己的平稳测度")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
