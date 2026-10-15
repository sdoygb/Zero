# -*- coding: utf-8 -*-
"""聚焦长代并行：geom seed × retention，path64/128，60 万代（吃满 10 核）。"""
from __future__ import annotations
import sys, os, time, json
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from scale_parallel import work_one, RES


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 10
    gens = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 600000
    cfgs = []
    for n in (64, 128):
        for sk in ("geom", "hub2"):
            for ret in ("last", "all", "top2"):
                cfgs.append(dict(n=n, seed=sk, retention=ret, L=4, gens=gens))
    print(f"聚焦长代并行：{len(cfgs)} 配置 × {workers} worker，每配置 {gens} 代")
    t0 = time.time()
    out = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(work_one, c): c for c in cfgs}
        for f in as_completed(futs):
            r = f.result()
            out.append(r)
            c = r["config"]
            d = abs(r["q_hist"] - 5 / 9) if r["q_hist"] else float("nan")
            tag = "≈5/9" if d < 1e-2 else ("窗内" if r["in_window"] else "出窗")
            print(f"  path{c['n']} {c['seed']} {c['retention']} {gens//1000}k代: "
                  f"q={r['q_hist']:.6f} ({tag})")
    wall = time.time() - t0
    print(f"墙钟 {wall:.1f}s  总吞吐 {len(cfgs)*gens/wall:.0f} 配置·代/秒")
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "scale_parallel_long.json"), "w", encoding="utf-8") as f:
        json.dump({"workers": workers, "gens": gens, "wall_s": round(wall, 1),
                   "results": out}, f, ensure_ascii=False, indent=1)
    print(f"已存 {RES}/scale_parallel_long.json")


if __name__ == "__main__":
    main()
