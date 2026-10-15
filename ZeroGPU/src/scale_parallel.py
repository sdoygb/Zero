# -*- coding: utf-8 -*-
"""
ZeroGPU · 多进程并行扫描（吃满 12 线程 CPU + 内存）

单配置的零和演化是顺序的（每代依赖上一代），无法分块并行；
但**多个独立配置天然并行**。本脚本用 ProcessPoolExecutor(10) 同时跑
10 个配置，每个 worker 占一个核 —— 吃满用户要求的 10 线程 CPU。

数据保留到 results/。

用法：
    /usr/local/bin/python3 scale_parallel.py [worker数] [代次]
"""
from __future__ import annotations
import sys, os, time, json
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj, torus_adj, star_adj
from z0_vectorized import ZeroCoreVec, QW_MIN, QW_MAX

RES = os.path.normpath(os.path.join(HERE, "..", "results"))


def geom_seed(V):
    return [float(2 ** (V - 1 - v)) for v in range(V)]

def hub_seed(V, k=2):
    s = [1.0] * V
    for v in range(min(k, V)):
        s[v] = 100.0
    return s

def one_seed(V):
    s = [0.0] * V; s[0] = 1000.0
    return s


# 模块级 worker（可 pickle）
def work_one(cfg):
    n, sk, retention, L, gens = cfg["n"], cfg["seed"], cfg["retention"], cfg["L"], cfg["gens"]
    adj = path_adj(n)
    V = len(adj)
    if sk == "geom":
        seed = geom_seed(V)
    elif sk.startswith("hub"):
        seed = hub_seed(V, int(sk[3:]))
    elif sk == "one":
        seed = one_seed(V)
    else:
        seed = [1.0] * V
    vec = ZeroCoreVec(adj, seed=seed, retention=retention, L=L)
    t0 = time.time()
    for _ in range(gens):
        vec.run_generation()
    q = vec.q_hist()
    return {"config": cfg, "q_hist": q,
            "in_window": QW_MIN < q < QW_MAX if q else False,
            "sec": round(time.time() - t0, 2)}


def build_configs(ns, seeds, rets, Ls, gens):
    cfgs = []
    for n in ns:
        for sk in seeds:
            for ret in rets:
                for L in Ls:
                    cfgs.append(dict(n=n, seed=sk, retention=ret, L=L, gens=gens))
    return cfgs


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 10
    gens = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 20000

    # 配置池：更大图的多种 seed × retention × L
    cfgs = build_configs(
        ns=[64, 128, 256],
        seeds=["geom", "hub2", "hub3", "one", "all1"],
        rets=["last", "all", "top2"],
        Ls=[4],
        gens=gens,
    )
    print("=" * 90)
    print(f"多进程并行：{len(cfgs)} 配置 × {workers} worker（吃满 CPU），每配置 {gens} 代")
    print("=" * 90)
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(work_one, c): c for c in cfgs}
        done = 0
        for f in as_completed(futs):
            r = f.result()
            results.append(r)
            done += 1
            if done % 5 == 0 or done == len(cfgs):
                el = time.time() - t0
                print(f"  [{done}/{len(cfgs)}] {el:.1f}s  吞吐={done*gens/el:.0f} 配置·代/秒")
    wall = time.time() - t0

    hits = [r for r in results if r["in_window"]]
    print("=" * 90)
    print(f"完成 {len(results)} 配置，墙钟 {wall:.1f}s，总吞吐 {len(cfgs)*gens/wall:.0f} 配置·代/秒")
    print(f"进窗配置 {len(hits)}/{len(cfgs)}：")
    for r in hits:
        c = r["config"]
        print(f"  path{c['n']} {c['seed']} {c['retention']} L{c['L']}: q={r['q_hist']:.4f}")
    # 5/9 相关（几何 seed 的所有 retention/L）
    g59 = [r for r in results if r["config"]["seed"] == "geom"]
    print("\n几何 seed 的所有组合（看 5/9 是否跨 retention/L 普适）：")
    for r in g59:
        c = r["config"]
        d = abs(r["q_hist"] - 5 / 9) if r["q_hist"] else float("nan")
        tag = "≈5/9" if d < 1e-2 else ("窗内" if r["in_window"] else "出窗")
        print(f"  path{c['n']} {c['seed']} {c['retention']} L{c['L']}: q={r['q_hist']:.4f} ({tag})")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "scale_parallel.json"), "w", encoding="utf-8") as f:
        json.dump({"workers": workers, "gens": gens, "wall_s": round(wall, 1),
                   "results": results}, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {RES}/scale_parallel.json")


if __name__ == "__main__":
    main()
