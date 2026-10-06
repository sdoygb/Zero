"""z0_comm.py --- 交换子代数维数 vs 对称性可解释维数（判定"隐藏守恒量"是否存在）"""
import numpy as np, json, os
import networkx as nx
import l0_closure as L0, observable_sweep as OS
from zcl import Engine
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"results")
eng=Engine(); res={}
print(f"{'T':>3} {'N':>5} {'Σm²(交换子)':>12} {'Σ|orb|²(对称)':>14} {'隐藏荷存在?':>12}")
for T in (10,12,14,16):
    reps,_=L0.enumerate_necklaces(T,eng,chunk_bits=26,verbose=False)
    A=OS.build_sparse(T,reps,eng); N=A.shape[0]
    ev=np.linalg.eigvalsh(A.toarray())
    u,m=np.unique(np.round(ev,6),return_counts=True)
    comm=int((m.astype(np.int64)**2).sum())
    G=nx.from_scipy_sparse_array(A)
    gm=nx.algorithms.isomorphism.GraphMatcher(G,G)
    autos=list(gm.isomorphisms_iter())
    seen=set(); orbits=[]
    for v in range(N):
        if v in seen: continue
        o=frozenset(a[v] for a in autos); orbits.append(o); seen|=set(o)
    sym=int(sum(len(o)**2 for o in orbits))
    extra=comm-sym
    print(f"{T:>3} {N:>5} {comm:>12} {sym:>14} {('是 (差 %d)'%extra) if extra>0 else '否':>12}")
    res[T]={"N":int(N),"commutant":comm,"symmetry":sym,"extra":int(extra),
            "aut_order":len(autos),"n_orbits":len(orbits),
            "max_orbit":int(max(len(o) for o in orbits))}
json.dump(res,open(os.path.join(OUT,"z0_comm.json"),"w"),ensure_ascii=False,indent=1)
print("\n判据: Σm² > Σ|orb|² ⟹ 存在超出对称性的隐藏守恒量（可积）")
