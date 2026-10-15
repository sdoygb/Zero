# -*- coding: utf-8 -*-
"""
ZeroGPU · 更大规模：确认 q→5/9 在大尺度稳定（CPU float64）

path64 / path128 / path256 几何 seed，按时间预算跑尽可能多代，
确认 q→5/9 尺度无关性在更大图上稳定成立（更大图收敛慢，需更长代）。

数据保留到 results/。

用法：
    /usr/local/bin/python3 scale_up.py [秒预算/配置]
"""
from __future__ import annotations
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
FIVE_NINTH = 5 / 9

def geom_seed(V):
    return [float(2 ** (V - 1 - v)) for v in range(V)]


def run_scale(n, budget_s, L=4):
    V = n
    adj = path_adj(n)
    vec = ZeroCoreVec(adj, seed=geom_seed(V), retention="last", L=L)
    t0 = time.time()
    gen = 0
    samples = []
    SAMP = 5000
    while time.time() - t0 < budget_s:
        vec.run_generation()
        gen += 1
        if gen % SAMP == 0:
            q = vec.q_hist()
            samples.append({"gen": gen, "q_hist": round(q, 8),
                            "diff_59": round(abs(q - FIVE_NINTH), 8)})
    q = vec.q_hist()
    return {
        "graph": f"path{n}", "L": L, "generations": gen,
        "elapsed_s": round(time.time() - t0, 1),
        "rate": round(gen / max(time.time() - t0, 0.01), 1),
        "q_hist_final": round(q, 8),
        "diff_59_final": round(abs(q - FIVE_NINTH), 8),
        "series": samples[-8:],
    }


def main():
    args = sys.argv[1:]
    budget = float(args[0]) if args and args[0].replace(".", "", 1).isdigit() else 90.0
    # 可选指定图：--ns 64,128,256
    ns = [64, 128, 256]
    if "--ns" in args:
        ns = [int(x) for x in args[args.index("--ns") + 1].split(",")]
    print("=" * 88)
    print(f"更大规模：path{'/'.join(map(str,ns))} 几何 seed last L4，每配置预算 {budget}s")
    print("=" * 88)
    out = []
    for n in ns:
        print(f"▶ path{n} ...")
        r = run_scale(n, budget)
        out.append(r)
        print(f"  {r['generations']} 代 ({r['rate']:.0f} 代/秒, {r['elapsed_s']}s)")
        print(f"  q_hist={r['q_hist_final']:.8f}  差|q-5/9|={r['diff_59_final']:.2e}")
        if r["series"]:
            print(f"  末 8 采样: {[(s['gen'], s['q_hist']) for s in r['series']]}")
    print("=" * 88)
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "scale_up_59.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"已存 {RES}/scale_up_59.json")


if __name__ == "__main__":
    main()
