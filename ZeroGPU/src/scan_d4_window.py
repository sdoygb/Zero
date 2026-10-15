# -*- coding: utf-8 -*-
"""
ZeroGPU · 扫描参数：找让 q ∈ (1/2, 3/5) 的配置（D4 窗）

目标：在零和演化的真实参数空间里，找让两条回灌（L1→L0 q_hist、L2→L0 q_live）
落到 D4 窗 (0.5, 0.6) 内的配置。若存在，那个配置就是"三层稳态落在 D4"的候选。

扫描维度：
  图 Γ     ：path / star / 环(C_n) / 环面(T^2) —— 顶点数与度分布影响重数集中度
  seed     ：均匀 / 双枢纽(2 顶点主导) / 三枢纽 / 几何衰减 —— 控制初始集中度
  retention：all(全历史) / top2(近两层) / last(上一层) —— 控制历史累积窗口
  L        ：词长 4 / 6 —— 闭合词的规模
  机制     ：cap+kac（度相关封顶 + 局部毁灭）可选

判定：跑 gens 代后，q_hist 与 q_live 是否 ∈ (0.5, 0.6)。

用法：
    /usr/bin/python3 scan_d4_window.py          # 粗扫 + 精扫，报告进窗配置
"""
from __future__ import annotations
import sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import ZeroCore, torus_adj, path_adj, star_adj
from layer_closure import LayerClosure, QW_MIN, QW_MAX

# ───────── 配置生成 ─────────
def make_seed(V, kind):
    if kind == "all1":
        return [1] * V
    if kind == "hub2":
        s = [1] * V
        for v in (0, 1):
            s[v] = 100
        return s
    if kind == "hub3":
        s = [1] * V
        for v in (0, 1, 2):
            s[v] = 100
        return s
    if kind == "geom":
        return [2 ** (V - 1 - v) for v in range(V)]
    if kind == "one":
        s = [0] * V
        s[0] = 1000
        return s
    raise ValueError(kind)


def graphs():
    g = {}
    for n in (4, 5, 6, 8, 10):
        g[f"path{n}"] = path_adj(n)
    for n in (5, 6, 8, 10):
        g[f"star{n}"] = star_adj(n)
    for n in (6, 8, 12):
        g[f"cycle{n}"] = torus_adj(1, n)
    g["torus2_4"] = torus_adj(2, 4)     # 16 顶点，度 4（均匀）
    g["torus2_3"] = torus_adj(2, 3)     # 9 顶点，度 4（均匀）
    return g


# ───────── 粗扫 ─────────
def coarse_scan(gens=200):
    print("=" * 78)
    print(f"粗扫：跑 {gens} 代，找 q_hist 或 q_live 进窗 (0.5,0.6) 的配置")
    print("=" * 78)
    hits = []
    for gname, adj in graphs().items():
        V = len(adj)
        for skind in ("all1", "hub2", "hub3", "geom", "one"):
            for ret in ("all", "top2", "last"):
                for L in (4, 6):
                    for mech in (False, True):
                        kw = dict(retention=ret, L=L)
                        if mech:
                            kw.update(cap=10 ** 5, kac=True, cap_deg=True)
                        try:
                            lc = LayerClosure(adj, seed=make_seed(V, skind), **kw)
                            lc.run(gens, log_every=0)
                            qh = lc.q_hist(); ql = lc.q_live()
                            in_h = QW_MIN < qh < QW_MAX if qh else False
                            in_l = QW_MIN < ql < QW_MAX if ql else False
                            if in_h or in_l:
                                hits.append((gname, skind, ret, L, mech, qh, ql))
                                print(f"  ★ {gname} {skind} {ret} L{L} mech={mech} "
                                      f"q_hist={qh:.3f} q_live={ql:.3f}")
                        except Exception as e:
                            pass
    print(f"粗扫命中 {len(hits)} 个配置")
    return hits


# ───────── 精扫（放大命中区） ─────────
def fine_scan(base_hits, gens=400):
    print()
    print("=" * 78)
    print(f"精扫：在命中/邻近配置上跑 {gens} 代，稳定进窗者")
    print("=" * 78)
    best = []
    for (gname, skind, ret, L, mech, qh0, ql0) in base_hits:
        adj = graphs()[gname]
        V = len(adj)
        for sk2 in (skind, "hub2", "hub3", "geom"):
            kw = dict(retention=ret, L=L)
            if mech:
                kw.update(cap=10 ** 5, kac=True, cap_deg=True)
            lc = LayerClosure(adj, seed=make_seed(V, sk2), **kw)
            lc.run(gens, log_every=0)
            qh = lc.q_hist(); ql = lc.q_live()
            d = 0
            if qh and QW_MIN < qh < QW_MAX: d |= 1
            if ql and QW_MIN < ql < QW_MAX: d |= 2
            if d:
                best.append((gname, sk2, ret, L, mech, qh, ql, d))
                tag = "hist+live" if d == 3 else ("hist" if d == 1 else "live")
                print(f"  ★ {gname} {sk2} {ret} L{L} mech={mech} 进窗[{tag}] "
                      f"q_hist={qh:.4f} q_live={ql:.4f}")
    return best


if __name__ == "__main__":
    hits = coarse_scan()
    best = fine_scan(hits)
    print()
    print("=" * 78)
    print(f"最终进窗配置：{len(best)} 个")
    for (g, s, r, L, m, qh, ql, d) in best:
        tag = "hist+live" if d == 3 else ("hist" if d == 1 else "live")
        print(f"  {g} seed={s} retention={r} L={L} mech={m} [{tag}] "
              f"q_hist={qh:.4f} q_live={ql:.4f}")
    print("=" * 78)
