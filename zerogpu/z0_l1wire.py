"""
z0_l1wire.py --- 接线检验：`copy` 的"2"是双覆盖 M_2 的"2"吗？
关键：模式 = 旋转类 ⟹ **Z_L 已被商掉**。那么作用在模式上的群是 D_L 还是只剩 Z_2？
判据：① 反射 s 在模式上的作用；② T 与 s 是否交换（[T,S]=0 ⟹ 繁殖对符号对称协变 = 原生）。
若 D_L 塌成 Z_2，则 M_2（需要 Z_L 部分）在模式空间上**不可用** ⟹ copy 的 2 不是 M_2。
"""
import numpy as np, json, os
from itertools import combinations
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
def refl(w):                      # G37: 反序 + 变号
    return tuple(-x for x in w[::-1])
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
        ks=[c[i] for i in range(len(c)) if mask>>i&1]; prev=0
        for k in ks+[len(w)]: o.append(mo(w[prev:k])); prev=k
    return o
M=modes(); ix={m:i for i,m in enumerate(M)}; n=len(M)
# ① 反射在模式上是否良定义 + 轨道结构
S=np.zeros((n,n))
ok=True
for p,w in enumerate(M):
    q=mo(refl(w))
    if q not in ix: ok=False; continue
    S[ix[q],p]=1
fixed=int(np.trace(S)); print(f"模式 {n} 个；反射良定义={ok}；不动模式(自镜像)={fixed}；"
                              f"Z_2 轨道数={(n-fixed)//2+fixed}")
# ② 各规则的 T 与 S 是否交换
print(f"\n{'规则':<13} {'[T,S]=0?':>10} {'说明':>34}")
res={"n_modes":n,"reflection_fixed":fixed,"orbits_z2":(n-fixed)//2+fixed}
for r in ("persist","copy","split_parent","split_only","split_all"):
    T=np.zeros((n,n))
    for p,w in enumerate(M):
        for c in ch(w,r): T[ix[c],p]+=1
    comm=float(np.abs(T@S - S@T).max())
    res[r]={"commutator_max":comm}
    print(f"{r:<13} {str(comm<1e-12):>10} {'对符号对称协变 ⟹ 原生' if comm<1e-12 else '不协变':>34}")
json.dump(res, open(os.path.join(OUT,"z0_l1wire.json"),"w"), ensure_ascii=False, indent=1)
print("\n读法：模式=旋转类 ⟹ Z_L 已商掉 ⟹ D_L 只剩 Z_2；若 T 与反射交换，繁殖对符号对称原生。")
