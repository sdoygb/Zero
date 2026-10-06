"""
q_L.py --- 类权重 q_L = Σ_C (o_C / N_L)² 的**大 L 推广**（文档只到 L=8）。

R87 H2 的定义（逐位照抄）：
    ω_C = o_C / N_L,      o_C = 旋转类的轨道大小,   N_L = C(L, L/2)
    q_L = Σ_C ω_C²
文档值： q_4 = 5/9, q_6 = 7/25, q_8 = 133/1225。
`R23_dim_desc_results.json` 的 `summary` 里 DIM-COST-Q 的目标窗口是 `3/4 < q < 4/5`。

轨道大小 o_C = 该项链的最小周期 p（因为 w 有最小周期 p ⟺ 它有恰 p 个相异旋转）。
本脚本用向量化周期判定把 L 推到 32。
"""
from __future__ import annotations

import json
import os
import sys
from math import comb, gcd

import numpy as np

from zcl import Engine
import l0_closure as L0


def orbit_sizes(reps: np.ndarray, L: int) -> np.ndarray:
    """向量化：每条项链的最小周期 = 轨道大小 o_C。"""
    mask = np.uint64((1 << L) - 1)
    n = len(reps)
    period = np.full(n, L, dtype=np.int64)
    done = np.zeros(n, dtype=bool)
    divs = sorted(d for d in range(1, L) if L % d == 0)      # 升序，取第一个命中的
    for p in divs:
        rot = ((reps << np.uint64(p)) | (reps >> np.uint64(L - p))) & mask
        m = (rot == reps) & (~done)
        period[m] = p
        done |= m
        if done.all():
            break
    return period


def main(Lmax: int = 32):
    eng = Engine()
    print("device:", eng.info()["name"])
    print(f"\n{'L':>3} {'类数 K(L)':>14} {'N_L=C(L,L/2)':>17} {'q_L':>18} "
          f"{'L/N_L':>14} {'q/(L/N_L)':>10} {'max ω_C':>10} {'1/max ω':>9}")
    rows = []
    for L in range(2, Lmax + 1, 2):
        reps, _ = L0.enumerate_necklaces(L, eng, verbose=False)
        N = comb(L, L // 2)
        sizes = orbit_sizes(reps, L)
        assert int(sizes.sum()) == N, f"L={L}: Σo_C={sizes.sum()} != N={N}"
        w = sizes / N
        q = float((w ** 2).sum())
        rows.append({"L": L, "classes": int(len(reps)), "N_L": N, "q_L": q,
                     "L_over_N_L": L / N, "ratio": q / (L / N),
                     "max_class_weight": float(w.max()),
                     "orbit_size_hist": {int(k): int(v) for k, v in
                                          zip(*np.unique(sizes, return_counts=True))}})
        print(f"{L:>3} {len(reps):>14,} {N:>17,} {q:>18.12g} {L/N:>14.6g} "
              f"{q/(L/N):>10.5f} {w.max():>10.6f} {1/w.max():>9.2f}")

    print("\n--- 与文档/目标的核对 ---")
    doc = {4: 5 / 9, 6: 7 / 25, 8: 133 / 1225}
    for L, v in doc.items():
        got = next(r["q_L"] for r in rows if r["L"] == L)
        print(f"  q_{L}: 本脚本 {got:.12f}  文档 {v:.12f}  "
              f"{'OK' if abs(got - v) < 1e-12 else '*** 不符 ***'}")
    q4plus = [r["q_L"] for r in rows if r["L"] >= 4]
    print(f"  DIM-COST-Q 窗口 (3/4, 4/5) = ({3/4}, {4/5})")
    print(f"  L>=4 的 q_L 全落在 [{min(q4plus):.6g}, {max(q4plus):.6g}]，"
          f"最大值 {max(q4plus):.6f} < 3/4 = 0.75  ->  "
          f"{'窗口不可达' if max(q4plus) < 0.75 else '窗口可达'}")

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "q_L.json"), "w") as f:
        json.dump({"device": eng.info(), "rows": rows,
                   "dim_cost_q_window": [0.75, 0.8]}, f, indent=2)
    print("\nsaved results/q_L.json")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 32)
