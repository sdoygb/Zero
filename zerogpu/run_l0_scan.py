"""run_l0_scan.py --- 把 L0 闭链图扫到 T=32（文档最远 T=18 的 3.6 倍距离）。

输出 zerogpu/results/l0_scan.json
"""
from __future__ import annotations

import json
import os
import sys
import time
import numpy as np

from zcl import Engine
import l0_closure as L0

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)


def main(Ts):
    eng = Engine()
    info = eng.info()
    print("device:", info, flush=True)
    rows = []
    t_all = time.time()
    for T in Ts:
        t0 = time.time()
        print(f"\n=== T={T} ===  scanning...", flush=True)
        reps, st = L0.enumerate_necklaces(T, eng, chunk_bits=26, verbose=True)
        assert st["necklaces"] == st["necklaces_formula"], f"T={T} 项链数不符"
        print(f"  necklaces = {st['necklaces']:,}  [{time.time()-t0:.1f}s]", flush=True)

        deg, ne = L0.build_graph(T, reps, eng, chunk=1 << 19, verbose=True)
        avg = float(deg.mean())
        hist = np.bincount(deg)
        # 度分布 vs 均值：算二阶矩，看是否收敛到某个分布
        m2 = float((deg.astype(np.float64) ** 2).mean())
        row = {
            "T": T,
            "nodes": int(st["necklaces"]),
            "nodes_formula": int(st["necklaces_formula"]),
            "edges": int(ne),
            "avg_degree": round(avg, 6),
            "avg_degree_minus_half_T": round(avg - T / 2, 6),
            "deg_min": int(deg.min()), "deg_max": int(deg.max()),
            "deg_std": round(float(deg.std()), 4),
            "deg_second_moment": round(m2, 4),
            "deg_hist": hist.tolist(),
            "scan_seconds": round(time.time() - t0, 1),
        }
        rows.append(row)
        print(f"  edges={ne:,}  avg_deg={avg:.4f}  (T/2+0.5={T/2+0.5:.2f})  "
              f"max={deg.max()}  std={deg.std():.3f}  [{row['scan_seconds']}s]", flush=True)
        with open(os.path.join(OUT, "l0_scan.json"), "w") as f:
            json.dump({"device": info, "rows": rows,
                       "total_seconds": round(time.time() - t_all, 1)}, f,
                      indent=2, ensure_ascii=False)
    print("\nDONE", flush=True)
    print(f"{'T':>3} {'nodes':>14} {'edges':>16} {'avg_deg':>9} {'-T/2':>8} {'max':>5} {'std':>7}")
    for r in rows:
        print(f"{r['T']:>3} {r['nodes']:>14,} {r['edges']:>16,} {r['avg_degree']:>9.4f} "
              f"{r['avg_degree_minus_half_T']:>8.4f} {r['deg_max']:>5} {r['deg_std']:>7.3f}")


if __name__ == "__main__":
    Ts = [int(a) for a in sys.argv[1:]] or list(range(8, 34, 2))
    main(Ts)
