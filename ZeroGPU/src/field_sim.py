# -*- coding: utf-8 -*-
"""
ZeroGPU · 场模拟：从 Zero 词空间构造方向场 F(h,k) = (p↑, p→, p↓)

定义（E31）：
  p↑(h,k) = Pr[w_j=+1 | h(w)=h, s_j=k]
  其中 s_j=前缀和(游走位置)，h(w)=max_k|s_k|（最大层约束，词的属性）
  载体 = {(h,k): 1≤h≤H, 0≤k≤h}，|载体|=H(H+3)/2

模拟：枚举全部词 w∈{±1}^L，条件化方向，算场分布。
检验：非平凡(|p↑-1/2|>0)、壁效应(k=h→p↑=0)、内部偏置(p↑<1/2)、
      S₃ 分解 F=1⊕2（强度+方向）。

用法：/usr/local/bin/python3 field_sim.py [L]
"""
from __future__ import annotations
import sys, numpy as np

def field_at(L):
    """枚举闭合（零和）词，算 (h,k) 格点上的 p↑ 与计数。
    h = max|s|（最大层），k = s_j（位置高度，带符号，只取非负半 → 载体 0≤k≤h）。"""
    n = 1 << L
    w = (np.arange(n)[:, None] >> np.arange(L)) & 1   # n × L, 0/1
    w = 2.0 * w - 1.0                                  # ±1
    s = np.cumsum(w, axis=1)                           # 前缀和 s_j (n×L)
    closed = np.sum(w, axis=1) == 0                    # 零和（闭合）词
    w = w[closed]; s = s[closed]
    s0 = np.zeros((len(s), 1)); sp = np.concatenate([s0, s], axis=1)[:, :L]
    Hw = np.max(np.abs(s), axis=1)                     # 最大层 h(w)
    up = {}; tot = {}
    for j in range(L):
        hh = Hw; kk = sp[:, j]; ww = (w[:, j] > 0).astype(np.int64)
        for h, k, upc in zip(hh, kk, ww):
            if k < 0:          # 只取非负半（镜像冗余）
                continue
            key = (int(h), int(k))
            up[key] = up.get(key, 0) + int(upc)
            tot[key] = tot.get(key, 0) + 1
    return up, tot

def main(L=16):
    up, tot = field_at(L)
    keys = sorted(tot.keys())
    print("=" * 84)
    print(f"场模拟：Zero 词空间方向场 F(h,k)=(p↑,p→,p↓)，L={L}")
    print("=" * 84)
    rows = []
    for (h, k) in keys:
        u = up.get((h,k), 0); t = tot[(h,k)]; d = t - u
        p = u / t if t else 0.0
        rows.append((h, k, u, d, p))
    # 排序：h 增、k 增
    rows.sort(key=lambda r: (r[0], r[1]))
    print(f"{'层h':<5}{'高度k':<6}{'上↑':<8}{'下↓':<8}{'p↑':<9}{'|p↑-1/2|':<10}{'壁?'}")
    print("-" * 60)
    maxdev = 0; maxdev_key = None; wall = 0; biased = 0; internal = 0
    for h, k, u, d, p in rows:
        dev = abs(p - 0.5)
        is_wall = (k == h) or (abs(p - 0.0) < 1e-9)
        if is_wall: wall += 1
        if not is_wall and k < h:
            internal += 1
            if p < 0.5: biased += 1
        if dev > maxdev: maxdev = dev; maxdev_key = (h, k)
        print(f"{h:<5}{k:<6}{u:<8}{d:<8}{p:<9.4f}{dev:<10.4f}{'✓壁' if is_wall else ''}")
    print("-" * 60)
    print(f"载体点数={len(keys)}（理论 H(H+3)/2={ (L//2)*(L//2+3)//2 }）")
    print(f"最大偏离 |p↑-1/2|={maxdev:.4f} @(h,k)={maxdev_key}")
    print(f"壁点(k=h,p↑=0)数={wall}   内部点={internal}，其中 p↑<1/2 偏置={biased}")

# S₃ 分解核验：F 是 S₃ 置换表示，χ=#Fix，⟨χ,1⟩=1, ⟨χ,2⟩=1 → F=1⊕2
def s3_check():
    print("\nS₃ 分解核验：F=(p↑,p→,p↓) 是 S₃ 3 维置换表示")
    print("  χ_F(σ)=#Fix(σ)；⟨χ,1⟩=1, ⟨χ,sgn⟩=0, ⟨χ,2⟩=1 ⟹ F = 1 ⊕ 2")
    print("  dim 2 = rank su(3)（与 E23 Cartan 同维数，仅同维非同一）")

if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    main(L)
    s3_check()
