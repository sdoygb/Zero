"""z0_wire.py --- L1↔L2 接线：把纯组合的 T 换成 T·diag(w)，w=该模式在**实际行走层**里的实现重数"""
import numpy as np, json, os
from itertools import combinations
from collections import Counter
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"results")
def modes(mp=12):
    out=[]
    for T in range(2,mp+1,2):
        s=set()
        for ones in combinations(range(T),T//2):
            w=[1]*T
            for i in ones: w[i]=-1
            s.add(min(tuple(w[i:]+w[:i]) for i in range(T)))
        out+=sorted(s)
    return out
def cuts(w):
    s=0;o=[]
    for i,x in enumerate(w):
        s+=x
        if s==0 and 0<i+1<len(w): o.append(i+1)
    return o
def mo(s):
    T=len(s); return min(tuple(s[i:]+s[:i]) for i in range(T))
def ch(w,r):
    c=cuts(w)
    if r=="persist": return [w]
    if r=="copy": return [w,w]
    if r=="split_parent":
        o=[w]
        for k in c: o+=[mo(w[:k]),mo(w[k:])]
        return o
    if r=="split_only":
        o=[]
        for k in c: o+=[mo(w[:k]),mo(w[k:])]
        return o
    o=[]
    for mask in range(1<<len(c)):
        ks=[c[i] for i in range(len(c)) if mask>>i&1]
        prev=0
        for k in ks+[len(w)]: o.append(mo(w[prev:k])); prev=k
    return o
M=modes(); ix={m:i for i,m in enumerate(M)}
# w_p = 实际行走层里该模式的实现重数 = 轨道大小（均匀 Γ 上等概率路径 ⟹ 权重=轨道大小）
orb=[]
for w in M:
    T=len(w); orb.append(sum(1 for i in range(T) if mo(w[i:]+w[:i])==w))
orb=np.array(orb,float)
cd=Counter()
for w,o in zip(M,orb): cd[len(cuts(w))]+=int(o)
print(f"模式 {len(M)} 个；轨道大小之和 = {int(orb.sum())}（= 所有平衡词的个数）")
print(f"按内部归零点数 r 的实现重数分布: {dict(sorted(cd.items()))}")
res={}
print(f"\n{'规则':<13} {'ρ(T)':>8} {'ρ(T·diag(w))':>13} {'变化':>7}")
for r in ("persist","copy","split_parent","split_only","split_all"):
    T=np.zeros((len(M),len(M)))
    for p,w in enumerate(M):
        for c in ch(w,r): T[ix[c],p]+=1
    T2=T*orb[None,:]
    r1=float(np.max(np.abs(np.linalg.eigvals(T))))
    r2=float(np.max(np.abs(np.linalg.eigvals(T2))))
    res[r]={"rho_comb":round(r1,4),"rho_wired":round(r2,4)}
    print(f"{r:<13} {r1:>8.3f} {r2:>13.3f} {'×%.2f'%(r2/r1) if r1>0 else '  —':>7}")
json.dump({"orbit_sum":int(orb.sum()),"r_multiplicity":dict(sorted(cd.items())),"rules":res},
          open(os.path.join(OUT,"z0_wire.json"),"w"),ensure_ascii=False,indent=1)
print("\n读法：接线后 ρ 变化 ⟹ L2 的实现重数确实改变了 L1 的生长判决；ρ 不变 ⟹ 接线无效。")
