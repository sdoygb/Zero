# -*- coding: utf-8 -*-
"""
ZeroGPU · record 回写试算：让完整历史(L1)持续调制 L2 演化，看能否产生新粒子

Z088 六边已有：边⑥(L1→L0 记录累积)、边②(L1→L2 top2 历史作种子)——给稳定
不动点，多粒子谱开放(Z088 §7.2)。

本脚本试算【新回写】：完整 record 作为【演化源项】持续注入 L2（W 每步演化叠加
record 源项，α 控制强度）——超出 Z088 六边，检验能否打破单稳态产生新粒子。

  基线：top2 无源项（Z017 定理7）——预期稳定重复
  回写：record 源项 α∈{0.1,1,10}——观察是否出现新频率(多粒子)

用法：/usr/local/bin/python3 k4_readback.py
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
NORM = 1e-12

def freqs_of(arr, signif=1e-4):
    w = arr / arr.sum()
    mx = w.max()
    sig = [(i, float(w[i])) for i in range(len(w)) if w[i] > signif*mx]
    sig.sort(key=lambda x: -x[1])
    if len(sig) < 2: return []
    return sorted(round(abs(math.log(sig[i][1]/sig[i+1][1])),4) for i in range(len(sig)-1))

def run_config(name, alpha, cycles=3000, L=4, T=4, seed=None):
    V = 4
    A = np.array([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]], float)
    rec = np.zeros(V); W = np.zeros((V,V)); layers = []
    sd = np.array(seed if seed else [8.,4.,1.5,1.], float)
    destructions = 0
    freq_sets = []
    for cyc in range(cycles):
        for _ in range(T):
            Wn = W.copy(); Wn[np.arange(V),np.arange(V)] += sd
            for _ in range(L):
                Wn = Wn @ A
                if alpha > 0:
                    Wn += alpha * np.diag(rec)   # ★ record 源项持续注入
            diag = Wn[np.arange(V),np.arange(V)]
            rec = rec + diag
            W = Wn.copy(); W[np.arange(V),np.arange(V)] = 0
            # 重播种：top2(最近两层闭合)
            layers.append(diag.copy())
            if len(layers) >= 2:
                sd = layers[-2] + diag
            else:
                sd = diag.copy()
            if sd.sum() == 0: sd = np.ones(V)
            # 毁灭
            if (cyc*T + _ + 1) % T == 0:
                W = np.zeros((V,V)); destructions += 1
            # 归一化
            for a in (rec, W, sd):
                t = a.sum()
                if t > NORM: a /= t
        freq_sets.append(tuple(freqs_of(rec)))
    return name, alpha, rec, freq_sets

print("=" * 84)
print("record 回写试算：完整历史调制 L2 演化 → 能否产生新粒子")
print("=" * 84)
configs = [("baseline(top2 无源项)", 0.0), ("回写α=0.1", 0.1),
           ("回写α=1", 1.0), ("回写α=10", 10.0)]
results = {}
for name, a in configs:
    n, alpha, rec, fs = run_config(name, a)
    uniq = len(set(fs))
    final = freqs_of(rec)
    results[name] = {"alpha": alpha, "final_freqs": final,
                     "uniq_freq_sets": uniq,
                     "n_branch": len([x for x in rec if x > 1e-4*rec.max()])}
    print(f"\n[{name}]  α={alpha}")
    print(f"  末 record 频率={final}  跨周期频率集合数={uniq}")
    print(f"  record={[round(float(x),4) for x in rec]}")
print("=" * 84)

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_readback.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=1)
print(f"已存 {os.path.join(RES,'k4_readback.json')}")
