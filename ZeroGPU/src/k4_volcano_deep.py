#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_volcano_deep.py —— 深挖火山形机制：均匀静态 / 适度动态 / 极端锁定
====================================================================================
对 p∈{0,1,2,4,8,16} 分析三区机制：
  源结构 x = m−1/n（粒子弯曲时空源）
  势结构 φ（Lφ=−x 的丰富度）
  反馈驱动 φ−φmin（seed 反馈项, 时空引导粒子）
  seed 集中度（是否单峰锁定）
验证:
  p=0: 源=0 → 势=0 → 无反馈驱动 → 静态(平坦时空)
  p大: 源单峰 → 势单峰 → 反馈锁单点 → 单吸引子静态
  p中: 源多顶点差异 → 势有结构 → 反馈维持多峰 → 动态
"""
from __future__ import annotations
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from k4_inhomogeneity_scan import InhomogEngine, L4, n, K4

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def entropy(w):
    w = w / w.sum()
    w = w[w > 1e-12]
    return float(-np.sum(w * np.log(w)))

print("="*90)
print("深挖火山形机制：源结构 → 势结构 → 反馈驱动 → 动态/锁定")
print("="*90)

ps = [0.0, 1.0, 2.0, 4.0, 16.0]
N = 1500
rows = {}
for p in ps:
    eng = InhomogEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2", p=p)
    for _ in range(N):
        eng.run_generation()
    # 稳态采样
    m = eng.mass_dist(eng.record); xs = m - 1.0/n
    phi = np.linalg.pinv(L4) @ (-(xs - xs.mean())) ; phi -= phi.mean()
    fb = phi - phi.min() + 1e-9                       # 反馈项(seed 驱动)
    seed = eng.seed
    # 指标
    src_norm = float(np.linalg.norm(xs)); src_var = float(np.var(m))
    phi_range = float(phi.max() - phi.min()); phi_norm = float(np.linalg.norm(phi))
    fb_ent = entropy(fb); seed_ent = entropy(seed)
    seed_maxshare = float(seed.max()/seed.sum())
    # 源对称性: 非最大顶点的方差(1+3对称 → 方差~0)
    xs_sorted = np.sort(xs)[::-1]
    nonmax = xs_sorted[1:]
    nonmax_var = float(np.var(nonmax)) if len(nonmax) > 1 else 0.0
    rows[p] = {"src_norm": round(src_norm,5), "src_var": round(src_var,5),
               "phi_range": round(phi_range,5), "phi_norm": round(phi_norm,5),
               "fb_entropy": round(fb_ent,4), "seed_entropy": round(seed_ent,4),
               "seed_maxshare": round(seed_maxshare,4), "nonmax_var": round(nonmax_var,5),
               "src": np.round(xs,4).tolist(), "phi": np.round(phi,4).tolist(),
               "fb": np.round(fb,4).tolist()}
    print(f"\np={p}: 源x={np.round(xs,3)}  |x|={src_norm:.4f}  非最大方差={nonmax_var:.5f}")
    print(f"       势φ={np.round(phi,3)}  范围={phi_range:.4f}  |φ|={phi_norm:.4f}")
    print(f"       反馈项φ−min={np.round(fb,3)}  反馈熵={fb_ent:.3f}")
    print(f"       seed熵={seed_ent:.3f}  seed主占={seed_maxshare:.3f}")

print("\n" + "="*90)
print("火山形机制解读")
print("="*90)
for p in ps:
    r = rows[p]
    if r["src_norm"] < 1e-3:
        verdict = "均匀区: 源≈0 → 势≈0 → 无反馈 → 平坦静态"
    elif r["nonmax_var"] < 1e-4:
        verdict = "对称区: 源1+3对称(单峰+均匀) → 势对称 → 反馈稳定对称 → 锁定静态"
    else:
        verdict = "非对称区: 源多顶点差异 → 势有梯度 → 反馈驱动 → 动态"
    print(f"  p={p:>4}: |源|={r['src_norm']:.4f} 势范围={r['phi_range']:.4f} "
          f"非最大方差={r['nonmax_var']:.5f} 反馈熵={r['fb_entropy']:.3f} → {verdict}")

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_volcano_deep.json"), "w", encoding="utf-8") as f:
    json.dump({"rows": {str(p): rows[p] for p in ps}}, f, ensure_ascii=False, indent=1)
print(f"\n已存 {os.path.join(RES,'k4_volcano_deep.json')}")
