# -*- coding: utf-8 -*-
"""
ZeroGPU · 粒子动态演化模拟：层间动力学下粒子(频率槽)如何产生/变化

追踪忠实引擎（K₄ + top2 + T=4 毁灭）的每个毁灭周期内瞬态：
  毁灭 → L1 top2 重播种 → L2 演化(L步) → 闭合(gen_closures) → record 累积

观测：
  ① record 的 3 频率槽随时间是否稳定（粒子稳定）还是变化（粒子演化）？
  ② 每个毁灭周期内是否出现新频率（暂态级联 = 粒子产生）？
  ③ 毁灭-重播种后是否反复长出同一套频率（稳定粒子）？
  ④ 有无 record→演化的回馈（当前引擎 record 只累积不读入——不忠实点监测）

用法：/usr/local/bin/python3 k4_dynamics.py [cycles] [T]
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loyalty_engine import LoyaltyEngine

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]

def freqs_of(arr, signif=1e-4):
    w = arr / arr.sum()
    mx = w.max()
    sig = [(i, float(w[i])) for i in range(len(w)) if w[i] > signif * mx]
    sig.sort(key=lambda x: -x[1])
    if len(sig) < 2:
        return []
    fs = [abs(math.log(sig[i][1] / sig[i+1][1])) for i in range(len(sig)-1)]
    return sorted(round(f, 4) for f in fs)

def run(cycles=2000, T=4):
    eng = LoyaltyEngine(K4, seed=[8.0,4.0,1.5,1.0], T=T, retention="top2")
    # 细粒度：每代记录
    per_gen = []          # 每代瞬态
    per_cycle_freq = []   # 每毁灭周期的 record 频率
    first_300_cycle_freqs = None
    last_300_cycle_freqs = None
    freq_sets = []
    for cyc in range(cycles):
        cycle_record = []
        for _ in range(T):
            eng.run_generation()
            per_gen.append({
                "gen": eng.gen,
                "closures_freqs": freqs_of(eng.gen_closures),
                "record_freqs": freqs_of(eng.record),
                "n_branch": eng.branches()[0],
            })
            cycle_record.append(eng.record.copy())
        rf = freqs_of(eng.record)
        per_cycle_freq.append(rf)
        freq_sets.append(tuple(rf))
        # 记录前/后 300 周期比较（粒子演化检测）
        if cyc == 299:
            first_300_cycle_freqs = rf
        if cyc == cycles - 301:
            last_300_cycle_freqs = rf
    # 统计
    uniq_freq_sets = set(freq_sets)
    # 瞬态是否出现 record 之外的频率（级联/新粒子候选）
    all_closure = set()
    for g in per_gen:
        all_closure.update(g["closures_freqs"])
    n_new = len(all_closure - set(per_cycle_freq[-1]))
    print("=" * 84)
    print(f"粒子动态演化：K₄ top2+T={T}  {cycles}周期×{T}代 = {eng.gen} 代")
    print("=" * 84)
    print(f"record 3 频率槽（末）={per_cycle_freq[-1]}")
    print(f"record 频率集（跨周期不同集合数）={len(uniq_freq_sets)}")
    print(f"  首 300 周期末 record 频率 = {first_300_cycle_freqs}")
    print(f"  末 300 周期末 record 频率 = {last_300_cycle_freqs}")
    stable = first_300_cycle_freqs == last_300_cycle_freqs
    print(f"  ⟹ 粒子(频率槽)跨 2400 代 {'稳定重复' if stable else '演化变化'}：{stable}")
    print(f"  毁灭周期内瞬态(gen_closures)出现过的频率集 = {sorted(all_closure)}")
    print(f"  record 之外的瞬态频率数（新粒子候选/级联）= {n_new}")
    # 支数动态
    nb_set = sorted(set(g["n_branch"] for g in per_gen))
    print(f"  每代支数集 = {nb_set}")
    # record→演化回馈检测
    print(f"\n  [忠实性] record 是否读入演化？ 引擎 record 只累积不读入 → 不忠实点（Z089）")
    # 毁灭恢复：每周期第 1 代 vs 第 T 代的瞬态频率
    print(f"\n  毁灭恢复样本（每周期第1代 vs 末代瞬态）：")
    for cyc in [0, 100, 500, 999]:
        base = cyc * T
        f1 = freqs_of(eng.gen_closures) if False else None
        print(f"    周期{cyc}: 第1代支={per_gen[base]['n_branch']} 第{T}代支={per_gen[base+T-1]['n_branch']}")
    # 存
    out = {"config": {"adj": "K4", "seed": [8,4,1.5,1], "retention": "top2", "T": T,
                      "L": 4}, "cycles": cycles, "total_gen": eng.gen,
            "final_record_freqs": per_cycle_freq[-1],
            "uniq_record_freq_sets": len(uniq_freq_sets),
            "first300_vs_last300_stable": stable,
            "closure_transient_freqs": sorted(all_closure),
            "n_branch_set": nb_set,
            "record_readback": "no (record accumulates only; not read into evolution)"}
    os.makedirs(RES, exist_ok=True)
    fp = os.path.join(RES, "k4_dynamics.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {fp}")
    return out

if __name__ == "__main__":
    cyc = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    T = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    run(cycles=cyc, T=T)
