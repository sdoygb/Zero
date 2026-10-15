# -*- coding: utf-8 -*-
"""
ZeroGPU · 深化 q=5/9 普适性（CPU float64）

验证几何衰减 seed 的 q→5/9 是否普适，还是特定配置特例：
  ① 图类型依赖：path / cycle(环) / star / torus² —— 5/9 是否与拓扑无关？
  ② L 依赖：闭合词长 L=2,4,6 —— 5/9 是否依赖 L？
  ③ 更大图：path24/32 —— 尺度无关性是否在更大范围成立（需更多代收敛）。

数据保留到 results/。

用法：
    /usr/local/bin/python3 deepen_scan.py
"""
from __future__ import annotations
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z0_engine import path_adj, torus_adj, star_adj
from z0_vectorized import ZeroCoreVec, QW_MIN, QW_MAX

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
FIVE_NINTH = 5 / 9

def geom_seed(V):
    return [float(2 ** (V - 1 - v)) for v in range(V)]


def scan_topology():
    """图类型 × L：几何 seed 的 q 终值。"""
    print("=" * 84)
    print("① 图类型 × L：几何 seed 的 q 是否普适 = 5/9")
    print("=" * 84)
    graphs = {
        "path16": path_adj(16),
        "cycle16": torus_adj(1, 16),
        "star16": star_adj(16),
        "torus2_4x4": torus_adj(2, 4),
    }
    rows = []
    for gname, adj in graphs.items():
        V = len(adj)
        for L in (2, 4, 6):
            vec = ZeroCoreVec(adj, seed=geom_seed(V), retention="last", L=L)
            for _ in range(60000):
                vec.run_generation()
            q = vec.q_hist()
            d59 = abs(q - FIVE_NINTH) if q is not None else float("nan")
            rows.append(dict(graph=gname, L=L, q=round(q, 8), diff_59=round(d59, 8),
                             in_win=QW_MIN < q < QW_MAX))
            tag = "≈5/9" if d59 < 1e-3 else ("窗内" if QW_MIN < q < QW_MAX else "出窗")
            print(f"  {gname:>14} L={L}: q={q:.6f}  差|q-5/9|={d59:.2e}  {tag}")
    return rows


def scan_scale():
    """更大图 path：尺度无关性。"""
    print()
    print("=" * 84)
    print("② 更大 path：path24/32，确认尺度无关（收敛需更多代）")
    print("=" * 84)
    rows = []
    for n in (24, 32):
        V = n
        gens = 80000
        vec = ZeroCoreVec(path_adj(n), seed=geom_seed(V), retention="last", L=4)
        t0 = time.time()
        for _ in range(gens):
            vec.run_generation()
        q = vec.q_hist()
        rows.append(dict(graph=f"path{n}", L=4, q=round(q, 8),
                         diff_59=round(abs(q - FIVE_NINTH), 8), gens=gens,
                         sec=round(time.time() - t0, 1)))
        print(f"  path{n} L=4 {gens}代: q={q:.8f}  差|q-5/9|={abs(q-FIVE_NINTH):.2e}  "
              f"({rows[-1]['sec']}s)")
    return rows


def main():
    os.makedirs(RES, exist_ok=True)
    out = {}
    out["topology"] = scan_topology()
    out["scale"] = scan_scale()

    print()
    print("=" * 84)
    print("总结")
    print("=" * 84)
    top = out["topology"]
    hits = [r for r in top if r["diff_59"] < 1e-3]
    print(f"图类型×L：收敛到 ≈5/9 的 {len(hits)}/{len(top)} 个组合")
    for r in top:
        print(f"  {r['graph']} L={r['L']}: q={r['q']:.6f} ({'≈5/9' if r['diff_59']<1e-3 else '非5/9'})")
    print(f"\n更大 path：5/9 尺度无关性")
    for r in out["scale"]:
        print(f"  {r['graph']}: q={r['q']:.8f} 差={r['diff_59']:.2e}")

    with open(os.path.join(RES, "deepen_59_summary.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {RES}/deepen_59_summary.json")


if __name__ == "__main__":
    main()
