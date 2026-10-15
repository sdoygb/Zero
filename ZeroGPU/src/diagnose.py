# -*- coding: utf-8 -*-
"""
ZeroGPU · 阻塞诊断（可复现）

三件事：
  ① 复现参考实现的**静默溢出**（在参考实现目录内跑它的引擎，只读不改）
  ② 复现忠实引擎的**三个运行模式**（不动点 / 爆炸 / 自持）
  ③ 量化"稳定性 vs 空间结构"的互斥

用法：
    /usr/bin/python3 diagnose.py
输出：
    results/diagnose.json
"""
from __future__ import annotations
import sys, os, json
from math import log10

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from z0_engine import ZeroCore, torus_adj, star_adj, path_adj

REF = "/Users/oygb/Downloads/Zero-261010/zerogpu"   # 参考实现（只读）
OUT = os.path.join(os.path.dirname(HERE), "results", "diagnose.json")


# ───────── ① 复现参考实现的溢出 ─────────
def ref_overflow():
    print("=" * 74); print("① 参考实现的静默溢出（int64）"); print("=" * 74)
    if not os.path.isdir(REF):
        print("  参考实现目录不在，跳过"); return None
    sys.path.insert(0, REF)
    try:
        import numpy as np
        from z0_core import Z0Core as RefCore, torus_adj as ref_torus
    except Exception as e:
        print("  导入失败：", e); return None
    A = ref_torus(2, 5); V = A.shape[0]
    rng = np.random.default_rng(11)
    eng = RefCore(A, rng.integers(1, 4, size=V).astype(np.int64), L=12, dtype=np.int64)
    rows = []
    frozen = None
    print(f"  int64 上限 = {np.iinfo(np.int64).max}")
    print(f"  {'gen':>4} {'active':>16} {'seed_max':>14}  备注")
    prev = None
    for g in range(10):
        eng.step(); eng.step()
        N = eng.M * eng.seed[:, None]
        act = float(np.sum(N))
        sm = int(eng.seed.max())
        note = ""
        if act < 0: note = "★溢出（负值）"
        if prev is not None and act == prev: note = "★冻结（与前代全同）"
        print(f"  {g:>4} {act:>16.6g} {sm:>14}  {note}")
        rows.append({"gen": g, "active": act, "seed_max": sm, "note": note})
        if prev is not None and act == prev and frozen is None: frozen = g
        prev = act
        if g + 1 < 10: eng.wipe_reseed()
    print(f"\n  → 结论：第 {frozen-1 if frozen else '?'} 代溢出，第 {frozen} 代起冻结")
    return {"rows": rows, "frozen_at": frozen}


# ───────── ② 忠实引擎的三个模式 ─────────
def faithful_modes():
    print(); print("=" * 74); print("② 忠实引擎的三个运行模式"); print("=" * 74)
    out = {}

    # A：传递图 + 均匀种子 ⟹ 精确不动点
    eA = ZeroCore(torus_adj(2, 4), L=4)
    for _ in range(40):
        eA.run_generation()
    aA = [h["active_log10"] for h in eA.hist]
    out["A"] = {"active_log10": aA, "gini": eA.gini(), "bits": eA.bits(),
                "spread": max(aA) - min(aA)}
    print(f"  A 传递图+均匀种子 : 活动 log10 恒定 = {aA[-1]:.4f}（极差 {max(aA)-min(aA):.2e}）"
          f"  Gini = {eA.gini():.4f}  bits = {eA.bits()}")

    # B：非传递图（星图）⟹ 空间结构 + 增长
    eB = ZeroCore(star_adj(9), L=4, seed=[1]*5 + [7]*4)
    tb = []
    for _ in range(10):
        eB.run_generation(); tb.append(eB.hist[-1]["active_log10"])
    out["B"] = {"active_log10": tb, "gini": eB.gini(), "bits": eB.bits()}
    print(f"  B 星图+非均匀种子 : 活动 log10 {tb[0]:.3f} → {tb[-1]:.3f}"
          f"  Gini = {eB.gini():.4f}  bits = {eB.bits()}")

    # C：自持环（闭路径系统的稳定自持）
    seed = [1, 1, 2, 3, 1, 1, 2, 3]
    eC = ZeroCore(torus_adj(1, 8), L=4, seed=seed)
    tc, rc, seeds = [], [], []
    for _ in range(300):
        eC.run_generation()
        tc.append(eC.hist[-1]["active_log10"]); rc.append(eC.hist[-1]["record_log10"])
        seeds.append(tuple(eC.seed))
    repeat = next((i for i, s in enumerate(seeds) if i > 0 and s == seeds[0]), None)
    out["C"] = {"active_log10": tc, "record_log10": rc, "bits": eC.bits(),
                "seed_repeat_at": repeat, "record_monotone":
                all(rc[i] <= rc[i+1] + 1e-12 for i in range(len(rc)-1))}
    print(f"  C 环上自持       : 300 代 活动 log10 极差 = {max(tc)-min(tc):.2e}"
          f"  记录单调 = {out['C']['record_monotone']}  种子首次重现于第 {repeat} 代"
          f"  bits = {eC.bits()}")
    return out


# ───────── ③ 稳定性 vs 空间结构 ─────────
def stability_vs_structure():
    print(); print("=" * 74); print("③ 稳定性 vs 空间结构（封顶强度扫描）"); print("=" * 74)
    print(f"  {'cap':>10} {'activity spread':>18} {'Gini':>10} {'bits':>8}")
    rows = []
    # 在环上扫封顶强度（非均匀种子）
    for cap in [None, 10**3, 10**6, 10**9, 10**12]:
        e = ZeroCore(torus_adj(1, 8), L=4, seed=[1, 1, 2, 3, 1, 1, 2, 3], cap=cap)
        a = []
        for _ in range(120):
            e.run_generation(); a.append(e.hist[-1]["active_log10"])
        spread = max(a) - min(a)
        rows.append({"cap": cap, "spread": spread, "gini": e.gini(), "bits": e.bits()})
        print(f"  {str(cap):>10} {spread:>18.4g} {e.gini():>10.4f} {e.bits():>8}")
    print("\n  → 封顶越低，活动层越稳；但 Gini 不随封顶上升 ⟹ 封顶不是结构的来源")
    return rows


if __name__ == "__main__":
    res = {}
    res["ref_overflow"] = ref_overflow()
    res["faithful"] = faithful_modes()
    res["cap_scan"] = stability_vs_structure()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=2, default=str)
    print(); print("=" * 74); print("已写出", OUT); print("=" * 74)
