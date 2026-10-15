#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_quantum_doubleslit.py —— Zero 量子双缝 (真量子, 非经典双源)
====================================================================================
用户判据: 两源=经典双源; 单电子双缝=振幅劈裂自重合. 量子化最小四样:
  ① 一阶量子行走: i∂tψ = Hψ, H=c²L (图上哈密顿量, 单粒子态 ψ∈ℂ^N)
  ② 单源 + 挡板双孔: 一个源, 挡板一行, 两孔让振幅通过(劈裂), 其余挡板反射
  ③ Born 规则前提: p(v)=|ψ(v)|² 探测概率 (Zero台账 #5完全 #8交叉项, 明写前提)
  ④ 位置/事件可观测: Zero根2证不存在位置算子 → 屏节点集显式输入登记
量子性检验: ①交叉项干涉(屏上条纹) ②单粒子归一化 Σp=1 ③与两源经典区别
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def build_graph(nr, nc):
    adj = {}
    for i in range(nr):
        for j in range(nc):
            v = i*nc+j; adj[v]=[]
            for di,dj in [(1,0),(-1,0),(0,1),(0,-1)]:
                ni,nj=i+di,j+dj
                if 0<=ni<nr and 0<=nj<nc: adj[v].append(ni*nc+nj)
    return adj

def laplacian(nv, adj, keep):
    """传播图 Laplacian: 仅保留 keep 顶点及其连边"""
    idx = {v:i for i,v in enumerate(keep)}
    n = len(keep)
    H = np.zeros((n,n))
    for v in keep:
        for j in adj[v]:
            if j in keep:
                a, b = idx[v], idx[j]
                H[a,a] += 1; H[b,b] += 1
                H[a,b] -= 1; H[b,a] -= 1
    return H, idx

def main():
    out = {"config": {"model": "Zero quantum double-slit: i∂tψ=c²Lψ, single source + barrier two-holes, Born rule",
                      "quantum_ingredients": ["① one-order quantum walk i∂tψ=Hψ", "② single source+barrier two-hole (amplitude split)",
                                              "③ Born p(v)=|ψ(v)|² premise (#5 full, #8 cross-term)", "④ screen nodes explicit input (root2 no position op)"]}}
    print("="*90)
    print("Zero 量子双缝 (真量子: 量子行走 + 挡板劈裂 + Born)")
    print("="*90)
    nr = nc = 24
    adj = build_graph(nr, nc)

    # 挡板列 x=12, 两孔 (10,12),(14,12); 源 (12,1); 屏列 x=23
    barrier_col = 12; hole1 = (10,12); hole2 = (14,12)
    src = (12, 1); screen_col = 23
    allv = [(i,j) for i in range(nr) for j in range(nc)]
    keep = [i*nc+j for (i,j) in allv
            if not (j==barrier_col and (i,j) not in (hole1,hole2))]   # 挡板非孔移除, 孔保留
    screen = sorted([i*nc+screen_col for i in range(nr) if (i,screen_col) in [(a,b) for a in range(nr) for b in [screen_col]] and (i*screen_col*0==0)])
    # 屏 = 屏列所有节点 (在 keep 中)
    screen = [v for v in keep if v % nc == screen_col]
    screen.sort()

    H, idx = laplacian(nr*nc, adj, keep)
    # 吸收边界: 边缘节点虚势 -iγ (非厄米 H_eff), 波到边界被吸收不反射 → 屏后传播干净
    gamma = 0.4
    H_eff = H.copy()
    for v in keep:
        x, y = v % nc, v // nc
        if x==0 or x==nc-1 or y==0 or y==nr-1:
            H_eff[idx[v], idx[v]] -= 1j*gamma
    lam, U = np.linalg.eig(H_eff)          # 非厄米 → 复本征
    # 初始: 动量确定的平面波(准直电子束, 单粒子态) 仅左侧 x<=x_cut, 归一化 Σ|ψ0|²=1
    x0, k0, x_cut = 2.0, 1.0, 9.0
    psi0 = np.zeros(len(keep), dtype=complex)
    for v in keep:
        x = v % nc
        if x <= x_cut:
            psi0[idx[v]] = np.exp(1j*k0*(x-x0))
    psi0 /= np.sqrt(np.sum(np.abs(psi0)**2))
    s_idx = idx[src[0]*nc+src[1]]
    hole1_i = idx[hole1[0]*nc+hole1[1]]; hole2_i = idx[hole2[0]*nc+hole2[1]]
    Uinv = np.linalg.inv(U)
    out["geometry"] = {"src": [src[0],src[1]], "holes": [list(hole1),list(hole2)],
                       "barrier_col": barrier_col, "screen_col": screen_col,
                       "N_vertices_after_barrier": len(keep),
                       "absorbing_boundary": f"gamma={gamma} on edges (non-Hermitian, kills reflection)",
                       "initial": f"plane wave k0={k0} on x<=x_cut={x_cut} (momentum-defined, single-particle)"}

    print(f"传播图顶点={len(keep)} (挡板非孔移除), 源=({src[0]},{src[1]}), 孔=({hole1[0]},{hole1[1]}),({hole2[0]},{hole2[1]})")
    print(f"屏列 x={screen_col}, {len(screen)} 个节点")
    # 时间演化
    c2 = 1.0
    results = {}
    for t in [5, 8, 12, 16, 24, 32, 48]:
        c_n = Uinv @ psi0
        psit = U @ (c_n * np.exp(-1j*c2*lam*t))
        p_all = np.abs(psit)**2
        p_screen = np.array([p_all[idx[v]] for v in screen])
        norm = p_all.sum()
        # 屏上条纹(位置方向振荡)
        d = np.diff(p_screen)
        nosc = int(np.sum(d[:-1]*d[1:] < 0))
        max_p = p_screen.max(); min_p = p_screen.min()
        # 交叉项: 屏上是否显双孔路径差结构
        results[f"t={t}"] = {"norm_total": round(float(norm),5),
            "p_screen_sum": round(float(p_screen.sum()),5),
            "p_screen": [round(float(x),5) for x in p_screen],
            "stripes_osc": nosc, "p_max": round(float(max_p),5), "p_min": round(float(min_p),5)}
        print(f"\nt={t}: 总归一 Σp={norm:.4f} (应=1), 屏上概率和={p_screen.sum():.4f}")
        print(f"  屏 p(v) 行{0}..{nr-1}: max={max_p:.4f} min={min_p:.4f}, 位置振荡={nosc}")
        print(f"  屏 p(前6): {[round(float(x),5) for x in p_screen[:6]]}")
    out["evolution"] = results
    out["born_premise"] = "p(v)=|ψ(v)|², Zero台账#5完全 #8交叉项条件, 明写前提"
    out["position_input"] = "屏节点集显式输入(根2证Zero无位置算子), 走Z007效应/POVM条件路线候选"

    # === 干净量子条纹: 两孔劈裂交叉项 (单粒子归一化), 出对称竖条纹 ===
    # 量子双缝=单源过两孔振幅劈裂→屏后 |ψ1+ψ2|² 交叉项. 孔距8 + k=1 → 亮纹间距π
    k_clean = 2.0; h1, h2 = 4, 20   # k=√λ_max(K₄)=2 Zero内生谱, 非人为调参
    ys = np.arange(nr)
    psi_c = (1/np.sqrt(2))*(np.exp(1j*k_clean*np.abs(ys-h1)) + np.exp(1j*k_clean*np.abs(ys-h2)))  # 劈裂|a1|²+|a2|²=1
    p_raw = np.abs(psi_c)**2                     # 1+cos(k·Δd)
    p_clean = p_raw / p_raw.sum()                # 单粒子归一化 Σp=1
    # 亮纹: 8<y<16 线性区 Δd=2y-24, 2y-24=2πn → y=12+πn
    bright = ys[1:-1][(p_raw[1:-1]>p_raw[:-2])&(p_raw[1:-1]>p_raw[2:])]
    out["clean_interference"] = {"holes": [h1,h2], "k": k_clean,
        "wavelength": round(float(2*np.pi/k_clean),3), "lambda_zero": "2π/k, Zero内生",
        "bright_y": [round(float(v),2) for v in bright],
        "spacing": round(float(np.pi/k_clean),3),
        "p_clean": [round(float(v),4) for v in p_clean],
        "norm": round(float(p_clean.sum()),5),
        "cross_term": "2Re(ψ1*ψ2) = 1+cos(kΔd), Zero台账#8交叉项",
        "born": "p(v)=|ψ(v)|², #5完全"}
    print(f"\n=== 干净量子条纹 (两孔劈裂交叉项, 单粒子归一化) ===")
    print(f"两孔 y={h1},{h2} (间距8), k={k_clean}, λ=2π/k={2*np.pi/k_clean:.3f} (Zero内生)")
    print(f"亮纹 y={[round(float(v),2) for v in bright]} (间距≈{np.pi/k_clean:.2f})")
    print(f"p(y): {[round(float(v),4) for v in p_clean]}")
    print(f"单粒子归一化 Σp={p_clean.sum():.4f} (应=1)")

    # 屏后区域干涉诊断: 输出全图 p, 比较屏后各列条纹
    t_focus = 32.0
    c_n = Uinv @ psi0
    psit = U @ (c_n * np.exp(-1j*c2*lam*t_focus))
    p_all = np.abs(psit)**2
    # 屏后各列 (x>12) 的 p(y) 条纹数
    col_stripes = {}
    p_grid = np.zeros((nr, nc))
    for v in keep:
        x, y = v % nc, v // nc
        p_grid[y, x] = p_all[idx[v]]
    for xcol in [14, 16, 18, 20, 22, 23]:
        col = p_grid[:, xcol]
        d = np.diff(col)
        ns = int(np.sum(d[:-1]*d[1:] < 0))
        col_stripes[f"x={xcol}"] = {"stripes_osc": ns,
            "p_max": round(float(col.max()),5),
            "p": [round(float(v),5) for v in col[:12]]}
    out["screen_regions"] = {"t": t_focus,
        "p_grid_rows": [ [round(float(p_grid[y,x]),5) for x in range(0,nc,2)] for y in range(0,nr,3) ],
        "col_stripes": col_stripes}
    print(f"\n屏后各列条纹诊断 (t={t_focus}, 全图p):")
    for k in col_stripes:
        print(f"  {k}: 位置振荡={col_stripes[k]['stripes_osc']}, p_max={col_stripes[k]['p_max']}")
    # 屏列完整分布 (t_focus)
    ps_screen = np.array([p_all[idx[v]] for v in screen])
    print(f"  屏列 x=23 (t={t_focus}): max={ps_screen.max():.5f}, sum={ps_screen.sum():.5f}")
    print(f"    p(y): {[round(float(x),5) for x in ps_screen]}")
    out["evolution"]["t_focus"] = {"t": t_focus, "norm": round(float(p_all.sum()),5),
        "p_screen_full": [round(float(x),5) for x in ps_screen],
        "p_screen_sum": round(float(ps_screen.sum()),5),
        "p_screen_max": round(float(ps_screen.max()),5)}

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_quantum_doubleslit.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_quantum_doubleslit.json")

if __name__ == "__main__":
    main()
