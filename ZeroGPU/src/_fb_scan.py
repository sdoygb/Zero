#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""反馈扫描：时空→物质反馈强度 → 相互转化振荡。"""
import numpy as np, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from k4_inhomogeneity_scan import InhomogEngine, laplacian, solve_potential
K4=[[1,2,3],[0,2,3],[0,1,3],[0,1,2]]; L4=laplacian(K4); n=4
out={}
for fb in [0.0,0.2,0.4,0.6,0.8,1.0]:
    eng=InhomogEngine(K4,seed=[8,4,1.5,1],T=4,retention="top2",p=1.0,beta=0.2,feedback=fb)
    N=2000; phis=[]; ms=[]
    for _ in range(N):
        eng.run_generation()
        m=eng.mass_dist(eng.record); xs=m-1.0/n
        phis.append(solve_potential(L4,xs)); ms.append(m)
    phis=np.array(phis); ms=np.array(ms)
    tail=float(np.std(phis[-500:],axis=0).max())
    mtail=float(np.std(ms[-500:],axis=0).max())
    verdict="static" if tail<1e-4 else ("dyn" if tail>1e-3 else "micro")
    out[fb]={"phi_tail": round(tail,6), "m_tail": round(mtail,6), "verdict": verdict}
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "fb_scan.json"),"w",encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("saved fb_scan.json")
