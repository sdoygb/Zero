"""z0_spec.py --- 条目1：G_T 的简并重数能否被对称群解释（它是"零和振动的能谱"吗）"""
import numpy as np, json, os
import networkx as nx
import l0_closure as L0, observable_sweep as OS
from zcl import Engine
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"results")
eng=Engine(); res={}
for T in (10,12,14,16):
    reps,_=L0.enumerate_necklaces(T,eng,chunk_bits=26,verbose=False)
    A=OS.build_sparse(T,reps,eng); N=A.shape[0]
    ev=np.linalg.eigvalsh(A.toarray())
    u,m=np.unique(np.round(ev,7),return_counts=True)
    G=nx.from_scipy_sparse_array(A)
    try:
        aut=nx.algorithms.isomorphism.GraphMatcher(G,G)

    except Exception: orb=None
    # 自同构群的轨道：用 networkx 的 automorphism_group 若无则用 WL 色当上界
    try:
        gm=nx.algorithms.isomorphism.GraphMatcher(G,G)
        autos=list(gm.isomorphisms_iter())
        order=len(autos); abelian=None
    except Exception as e:
        autos=[]; order=None
    # 顶点轨道（用自同构枚举，若太大则跳过）
    if autos and len(autos)<=200000:
        seen=set(); orbits=[]
        for v in range(N):
            if v in seen: continue
            o={a[v] for a in autos}; orbits.append(o); seen|=o
        norb=len(orbits)
    else:
        norb=None
    print(f"T={T:>2} N={N:>5}  不同本征值={len(u):>4}  最大重数={int(m.max()):>3}  "
          f"|Aut|={order}  轨道数={norb}")
    res[T]={"N":int(N),"distinct":int(len(u)),"max_mult":int(m.max()),
            "mult_hist":{int(k):int(v) for k,v in zip(*np.unique(m,return_counts=True))},
            "aut_order":order,"orbits":norb}
json.dump(res,open(os.path.join(OUT,"z0_spec.json"),"w"),indent=1)
