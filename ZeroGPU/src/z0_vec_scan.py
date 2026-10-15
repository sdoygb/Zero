# -*- coding: utf-8 -*-
"""
ZeroGPU · 向量化快速扫描（更大图找 D4 窗）

用向量化引擎（22455 代/秒）快速扫更大图的进窗配置。
扫 path_n × seed_kind × gens，报告 q_hist 终值 / 进窗占比 / 漂移。

用法：
    /usr/local/bin/python3 z0_vec_scan.py
"""
from __future__ import annotations
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec, QW_MIN, QW_MAX

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))

def seed_for(V, kind):
    s = [1.0] * V
    if kind == "hub2":
        s[0] = s[1] = 100
    elif kind == "hub3":
        s[0] = s[1] = s[2] = 100
    elif kind == "hub4":
        for v in range(min(4, V)): s[v] = 100
    elif kind == "geom":
        for v in range(V): s[v] = float(2 ** (V - 1 - v))
    elif kind == "one":
        s = [0.0] * V; s[0] = 1000
    return s

def main():
    results = []
    print("=" * 88)
    print("向量化快速扫描：path_n × seed × 更长程 — 找 q ∈ (0.5, 0.6)")
    print("=" * 88)
    for n in (10, 12, 16, 20):
        for sk in ("hub2", "hub3", "hub4", "geom", "one"):
            for gens in (20000, 60000):
                adj = path_adj(n)
                vec = ZeroCoreVec(adj, seed=seed_for(n, sk), retention="last", L=4)
                t0 = time.time()
                SAMP = 2000
                samples = []
                for g in range(gens):
                    vec.run_generation()
                    if (g + 1) % SAMP == 0:
                        qh = vec.q_hist()
                        samples.append((g + 1, qh, QW_MIN < qh < QW_MAX))
                n_in = sum(1 for _, q, w in samples if w)
                qh_first = samples[0][1]; qh_last = samples[-1][1]
                hit = QW_MIN < qh_last < QW_MAX
                row = dict(n=n, seed=sk, gens=gens,
                           q_hist=round(qh_last, 4),
                           drift=round(qh_last - qh_first, 4),
                           frac_in=round(n_in / len(samples), 2),
                           gini=round(vec.gini(), 3),
                           sec=round(time.time() - t0, 1))
                results.append(row)
                mark = " ★" if hit else ""
                print(f"  path{n:>3} {sk:>5} {gens:>6}: q_hist={row['q_hist']:.4f} "
                      f"漂移={row['drift']:+.3f} 进窗={row['frac_in']:.0%} Gini={row['gini']:.2f} {row['sec']}s{mark}")

    # 汇总命中
    hits = [r for r in results if QW_MIN < r["q_hist"] < QW_MAX]
    print("\n" + "=" * 88)
    print(f"进窗配置（q_hist ∈ (0.5,0.6)）：{len(hits)} 个")
    for r in hits:
        print(f"  path{r['n']} {r['seed']} {r['gens']}代  q_hist={r['q_hist']:.4f} 进窗={r['frac_in']:.0%}")
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "d4_vec_scan_summary.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=1)
    print(f"已存 {RES}/d4_vec_scan_summary.json")

if __name__ == "__main__":
    main()
