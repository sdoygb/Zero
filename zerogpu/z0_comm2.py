"""z0_comm2.py --- 纠正：对称性不变代数的维数 = Aut 在 V×V 上的**轨道数**（orbitals），不是 Σ|orb|²"""
import numpy as np, json, os
import networkx as nx
import l0_closure as L0, observable_sweep as OS
from zcl import Engine
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"results")
eng=Engine(); res={}
print(f"{'T':>3} {'N':>5} {'Σm²(交换子)':>12} {'orbitals(对称)':>14} {'隐藏荷?':>10}")
for T in (10,12,14,16):
    reps,_=L0.enumerate_necklaces(T,eng,chunk_bits=26,verbose=False)
    A=OS.build_sparse(T,reps,eng); N=A.shape[0]
    ev=np.linalg.eigvalsh(A.toarray()); u,m=np.unique(np.round(ev,6),return_counts=True)
    comm=int((m.astype(np.int64)**2).sum())
    G=nx.from_scipy_sparse_array(A); gm=nx.algorithms.isomorphism.GraphMatcher(G,G)
    autos=[np.array([a[v] for v in range(N)]) for a in gm.isomorphisms_iter()]
    # V×V 上的轨道：用 (i,j) 的陪集代表，直接对全部对做并查集
    parent=list(range(N*N))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for a in autos:
        for i in range(N):
            base=i*N; ai=a[i]*N
            for j in range(N):
                x,y=find(base+j),find(ai+a[j])
                if x!=y: parent[y]=x
    orb=len({find(k) for k in range(N*N)})
    print(f"{T:>3} {N:>5} {comm:>12} {orb:>14} {'是(差%d)'%(comm-orb) if comm>orb else '否':>10}")
    res[T]={"N":int(N),"commutant":comm,"orbitals":orb,"extra":int(comm-orb),"aut_order":len(autos)}
json.dump(res,open(os.path.join(OUT,"z0_comm2.json"),"w"),ensure_ascii=False,indent=1)
print("\n判据: Σm² > orbitals ⟹ 存在超出对称性的隐藏守恒量")
