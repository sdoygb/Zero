"""z0_charge.py --- 筛三个候选守恒量：逆序数 / 下降数 / 符号变化数"""
import numpy as np, json, os
import l0_closure as L0, observable_sweep as OS
from zcl import Engine
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"results")
def bits(x,T): return np.array([1 if (int(x)>>i)&1 else -1 for i in range(T)])
def inv(w):   # 逆序数（对环上代表元，标准定义）
    return sum(1 for i in range(len(w)) for j in range(i+1,len(w)) if w[i]>w[j])
def desc(w):  # 下降数（含环绕）
    T=len(w); return sum(1 for i in range(T) if w[i]>w[(i+1)%T])
def flip(w):  # 符号变化数（相邻相异对，环形）
    T=len(w); return sum(1 for i in range(T) if w[i]!=w[(i+1)%T])
eng=Engine(); res={}
for T in (12,14):
    reps,_=L0.enumerate_necklaces(T,eng,chunk_bits=26,verbose=False)
    A=OS.build_sparse(T,reps,eng); D=A.toarray(); N=len(reps)
    W=[bits(int(r),T) for r in reps]
    for name,f in (("逆序数",inv),("下降数",desc),("符号变化数",flip)):
        q=np.array([f(w) for w in W],float)
        Q=np.diag(q)
        comm=float(np.abs(D@Q-Q@D).max())
        nq=len(np.unique(q))
        # 每个本征空间是否落在单一 q 值上
        ev,V=np.linalg.eigh(D); u=np.round(ev,7)
        bad=0; tot=0
        for lam in np.unique(u):
            idx=np.nonzero(u==lam)[0]; tot+=1
            if len(idx)>1 and len(np.unique(q[idx]))>1: bad+=1
        print(f"T={T} {name:<8} 取值数={nq:>3}  [A,Q]max={comm:>8.3f}  "
              f"简并空间跨越多个q值的个数={bad}/{tot}")
        res[f"{T}_{name}"]={"values":nq,"commutator":comm,"degen_split":bad,"degen_total":tot}
json.dump(res,open(os.path.join(OUT,"z0_charge.json"),"w"),ensure_ascii=False,indent=1)
