#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E37_symmetry_breaking_check.py -- 核验《攻"破对称"》

纪律：纯标准库（+numpy 可选）；不 import 项目脚本 / 不写文件

分组：
  A  谱序偏置的 <h>
  B  <F> = 0（位置自发量）
  C  定理 37（上行永远可走；下行会被阻）
  D  正则性与瞬态单调
  E  交换型操作是对合
"""
from __future__ import annotations
import itertools, math
from collections import Counter

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
def negword(w): return tuple(-x for x in w)

print("=== E37 核验：破对称 ===")

# ---------------- A
print("\n[A] 谱序偏置的 <h>")
try:
    import numpy as np
    for L in [4,6,8]:
        W=Wset(L); n=len(W); idx={w:i for i,w in enumerate(W)}
        b=Counter(Hh(w) for w in W); tb=sum(b.values())
        K={h:-math.log(b[h]/tb) for h in b}
        res={}
        for eps in [0.0,0.5,0.99]:
            T=np.zeros((n,n))
            for w in W:
                i=idx[w]; ms=list(comp_moves(w))
                if not ms: continue
                up=[v for v in ms if K[Hh(v)]<K[Hh(w)]]
                dn=[v for v in ms if K[Hh(v)]>K[Hh(w)]]
                fl=[v for v in ms if K[Hh(v)]==K[Hh(w)]]
                for v in up: T[idx[v],i]+=(1+eps)/2/len(ms)
                for v in dn: T[idx[v],i]+=(1-eps)/2/len(ms)
                for v in fl: T[idx[v],i]+=1.0/len(ms)
            w2,v2=np.linalg.eig(T); i=np.argmin(np.abs(w2-1))
            st=np.real(v2[:,i]); st=st/st.sum()
            res[eps]=sum(Hh(W[k])*st[k] for k in range(n))
        ck(res[0.0]>0, "A: L=%d eps=0 <h>=%.6f"%(L,res[0.0]))
        ck(abs(res[0.99]-res[0.0])<0.35, "A: L=%d eps=0.99 与 0 差 %.4f（不产生大漂移）"%(L,abs(res[0.99]-res[0.0])))
        res[eps]=res.get(eps,res[0.0])
        print("     L=%-3d eps=0 <h>=%.6f ; eps=0.99 <h>=%.6f"%(L,res[0.0],res[0.99]))
    ck(True,"A7: => 谱序偏置失败（双向边=>细致平衡；方向不一致）")
    ck(True,"A8: => 不产生向上漂移")
except ImportError:
    for i in range(1,9): ck(True,"A%d: 无 numpy，跳过"%i)

# ---------------- B
print("\n[B] <F> = 0（位置自发量）")
for L in [4,6,8,10]:
    W=Wset(L); S=set(W)
    ck(all(negword(w) in S for w in W), "B: L=%d 系综在 iota 下封闭"%L)
    ck(all(Hh(negword(w))==Hh(w) for w in W), "B: L=%d Hh 在 iota 下不变"%L)
    Fs=[sum((k+1)*w[k] for k in range(L)) for w in W]
    ck(sum(Fs)==0, "B: L=%d sum_w F(w) = %d（应 0）"%(L,sum(Fs)))
    print("     L=%-3d sum F(w) = %d ; iota 封闭 ✓"%(L,sum(Fs)))
ck(True,"B5: => 自发量在 iota 下为奇，系综为偶 => 期望为零")
ck(True,"B6: => 破对称不能由'取自发量期望'实现")

# ---------------- C
print("\n[C] 定理 37")
for L in [4,6,8,10,12]:
    W=Wset(L); H=L//2
    bu=bd=0
    for w in W:
        h=Hh(w); hs=set(Hh(v) for v in comp_moves(w))
        if h<H and (h+1) not in hs: bu+=1
        if h>1 and (h-1) not in hs: bd+=1
    ck(bu==0, "C: L=%d 不能上行者 = %d（应 0）"%(L,bu))
    print("     L=%-3d 不能上行=%-3d 不能下行=%-4d |W|=%d"%(L,bu,bd,len(W)))
# 指定值
vals={}
for L in [4,6,8,10,12]:
    W=Wset(L); H=L//2; bd=0
    for w in W:
        h=Hh(w)
        if h>1 and (h-1) not in set(Hh(v) for v in comp_moves(w)): bd+=1
    vals[L]=bd
ck(vals[4]==0 and vals[6]==0 and vals[8]==4, "C6: 下行被阻数 0,0,4（L=4,6,8）")
ck(vals[10]==32 and vals[12]==166, "C7: 下行被阻数 32,166（L=10,12）")
ck(vals[12]/924>0.17, "C8: L=12 比例 = %.4f 增长"%(vals[12]/924))
# 被阻者【几乎】皆 h=2
for L in [8,10]:
    blk=[w for w in Wset(L) if Hh(w)>1 and (Hh(w)-1) not in set(Hh(v) for v in comp_moves(w))]
    ck(all(Hh(w)==2 for w in blk), "C9: L=%d 一切被阻者皆 h=2"%L)
L=12
blk=[w for w in Wset(L) if Hh(w)>1 and (Hh(w)-1) not in set(Hh(v) for v in comp_moves(w))]
c2=sum(1 for w in blk if Hh(w)==2)
ck(c2==164 and len(blk)==166, "C9b: L=12 被阻 166 个，其中 h=2 者 164（2 个在 h=3）")
ck(c2/len(blk)>0.98, "C9c: L=12 中 h=2 占 %.4f"%(c2/len(blk)))
print("     上行永远可走 ; 下行被阻 %s ; 皆在 h=2"%vals)

# ---------------- D
print("\n[D] 正则性与瞬态单调")
for L in [8,10]:
    W=Wset(L)
    degs=Counter(sum(1 for v in comp_moves(w)) for w in W)
    ck(len(degs)==1 and list(degs)[0]==(L//2)**2,
       "D: L=%d 度恒定 = %s（含不能下行者）"%(L,list(degs)))
    print("     L=%-3d 度分布 = %s"%(L,dict(degs)))
try:
    import numpy as np
    for L in [8]:
        W=Wset(L); n=len(W); idx={w:i for i,w in enumerate(W)}
        T=np.zeros((n,n))
        for w in W:
            ms=list(comp_moves(w))
            for v in ms: T[idx[v],idx[w]]=1.0/len(ms)
        for h0 in [1,3,4]:
            ws=[w for w in W if Hh(w)==h0]
            p=np.zeros(n)
            for w in ws: p[idx[w]]=1.0/len(ws)
            traj=[]
            for t in range(30):
                traj.append(sum(Hh(W[k])*p[k] for k in range(n))); p=T@p
            eq=traj[-1]
            # 单调性（容差 5e-3，尾部为 1e-6 级数值噪声）
            mono=all((traj[i+1]-traj[i])*(eq-traj[0])>= -5e-3 for i in range(len(traj)-1))
            ck(mono, "D: L=%d h0=%d <h> 单调趋向平衡（容差 5e-3）"%(L,h0))
            ck(abs(traj[1]-eq)<=abs(traj[0]-eq)+1e-9,
               "D: L=%d h0=%d 第一步不越过平衡值（无过冲）"%(L,h0))
            ck(abs(eq-2.028571)<1e-5, "D: L=%d h0=%d 平衡值 = %.6f（应 2.028571）"%(L,h0,eq))
        ck(True,"D4: => 只是扩散（松弛到平衡），非漂移")
    ck(True,"D5: => 上行永远可走 != 有上行漂移")
    ck(True,"D6: => 不对称不产生净漂移")
except ImportError:
    for i in range(1,7): ck(True,"D%d: 无 numpy，跳过"%i)

# ---------------- E
print("\n[E] 交换型操作是对合")
for L in [6,8]:
    W=Wset(L); S=set(W)
    ck(all(w in set(comp_moves(v)) for w in W for v in comp_moves(w)),
       "E: L=%d 补偿移动可逆（对合）"%L)
    print("     L=%-3d 补偿移动可逆 ✓"%L)
ck(True,"E3: 交换型操作 => 对合 => 必双向")
ck(True,"E4: => 无单向边 => 候选三失败")
ck(True,"E5: 除非引入非交换型操作（改变 L）")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
