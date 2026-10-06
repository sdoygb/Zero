"""
evo_long.py --- 判定："时间不够长" 还是 "结构不够"

**判据**：把同一次演化跑 4 个长度（20k/40k/80k/160k 步），对每个观测量算
$$\text{漂移}(t)=\frac{|\langle\text{末 20\% 窗}\rangle-\langle\text{前 20\% 窗}\rangle|}{\langle\text{末 20\% 窗}\rangle}$$
若只是时间不够 ⟹ 漂移应随笔长**衰减**（趋稳，可以说"再跑就会停"）；
若结构不够 ⟹ 漂移**不随 $t$ 衰减**（永远在漂，或按幂律无限漂下去）。
再做**有限尺寸**：同一时长下 $L=16,32,64$ 的对比（区分"时间"与"空间"）。

用法：/usr/bin/python3 evo_long.py
输出：results/evo_long.json
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

import z0_evo as EV

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

OBS = ["record_rate", "regions", "mean_region_size", "occupied_frac",
       "record_gini", "traffic_gini", "ac1", "excursion_mean_len"]


def drift(series, key):
    y = np.array(series[key], float)
    n = len(y)
    if n < 10:
        return None
    q = max(2, n // 5)
    head, tail = y[:q].mean(), y[-q:].mean()
    return {"head": float(head), "tail": float(tail),
            "rel_drift": float(abs(tail - head) / max(abs(tail), 1e-12))}


def main():
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "criterion": "漂移随步长衰减 ⟹ 时间问题；不衰减 ⟹ 结构问题", "time_test": [], "size_test": []}

    A = EV.torus_adj(2, 32)
    print("=" * 100)
    print("① 时间检验：同一系统跑 20k / 40k / 80k / 160k 步（Γ=ℤ² L=32，模式 REC）")
    print(f"{'步数':>8} {'窗口数':>7} " + " ".join(f"{k[:13]:>14}" for k in OBS))
    for steps in (20000, 40000, 80000, 160000):
        ser, info = EV.evolve(A, "REC", steps=steps, win=500, seed=7)
        row = {"steps": steps, "windows": len(ser["t"]), "obs": {}}
        line = []
        for k in OBS:
            d = drift(ser, k)
            row["obs"][k] = d
            line.append(f"{d['rel_drift']:>14.4f}" if d else f"{'-':>14}")
        # 尾部窗口的相对起伏（是否已平稳）
        for k in OBS:
            if row["obs"][k]:
                y = np.array(ser[k], float)[-max(4, len(ser[k]) // 5):]
                row["obs"][k]["tail_rel_std"] = float(y.std() / max(abs(y.mean()), 1e-12))
        out["time_test"].append(row)
        print(f"{steps:>8} {len(ser['t']):>7} " + " ".join(line), flush=True)

    print("\n①b 漂移 vs 步长（同一观测量跨 4 个长度的对比）")
    print(f"{'观测':<20} " + " ".join(f"{s:>10}" for s in ("20k", "40k", "80k", "160k")) + "   判读")
    for k in OBS:
        ds = [r["obs"][k]["rel_drift"] if r["obs"][k] else np.nan for r in out["time_test"]]
        mono = all(ds[i] >= ds[i + 1] * 0.9 for i in range(len(ds) - 1) if np.isfinite(ds[i + 1]))
        tag = "漂移衰减 ⟹ 时间问题" if (np.isfinite(ds[-1]) and ds[-1] < ds[0] * 0.6) else \
              ("漂移不衰减 ⟹ 结构/慢动力学" if np.isfinite(ds[-1]) else "—")
        print(f"{k:<20} " + " ".join(f"{d:>10.4f}" for d in ds) + f"   {tag}")

    print("\n" + "=" * 100)
    print("② 有限尺寸检验：同一时长（40k 步）下 L=16 / 32 / 64（模式 REC）")
    print(f"{'L':>4} {'N':>7} {'区域数':>9} {'区域尺寸':>10} {'记录速率':>10} "
          f"{'记录Gini':>9} {'占据率':>8} {'ac1':>7}")
    for L in (16, 32, 64):
        A2 = EV.torus_adj(2, L)
        ser, info = EV.evolve(A2, "REC", steps=40000, win=500, seed=7)
        row = {"L": L, "N": int(A2.shape[0])}
        vals = {}
        for k in ("regions", "mean_region_size", "record_rate", "record_gini",
                  "occupied_frac", "ac1"):
            y = np.array(ser[k], float)
            vals[k] = float(y[-max(4, len(y) // 5):].mean())
        row["tail_means"] = vals
        out["size_test"].append(row)
        print(f"{L:>4} {row['N']:>7} {vals['regions']:>9.2f} {vals['mean_region_size']:>10.2f} "
              f"{vals['record_rate']:>10.4f} {vals['record_gini']:>9.4f} "
              f"{vals['occupied_frac']:>8.4f} {vals['ac1']:>7.4f}", flush=True)

    # 区域尺寸的 N 依赖：若与 N 无关 ⟹ 有内禀长度（结构性的）；若 ∝N ⟹ 只在粗化
    sizes = [r["tail_means"]["mean_region_size"] for r in out["size_test"]]
    Ns = [r["N"] for r in out["size_test"]]
    sl = float(np.polyfit(np.log(Ns), np.log(sizes), 1)[0])
    out["region_size_scaling"] = {"exponent": round(sl, 4), "sizes": sizes, "Ns": Ns}
    print(f"\n   区域尺寸 vs N：指数 {sl:.3f}  "
          f"（≈0 ⟹ 内禀长度；≈1 ⟹ 只是整块粗化）")

    with open(os.path.join(OUT, "evo_long.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/evo_long.json")


if __name__ == "__main__":
    main()
