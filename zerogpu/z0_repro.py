"""
z0_repro.py --- 把 L1 的**五条子代生成规则**接上来，并**独立复算** 谱半径/幂零指数
（不依赖语料结果；对照 zero_sum_reproduction_transition_theorems.md 的 T2/T3 表）
模式 modes = 周期<=12 的循环平衡词（旋转类，共 123 个 —— 与语料"零特征值数 123"吻合）
"""
from __future__ import annotations
import json, os
import numpy as np
from itertools import combinations

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

def modes(maxp=12):
    out = []
    for T in range(2, maxp + 1, 2):
        seen = set()
        for ones in combinations(range(T), T // 2):
            w = [1] * T
            for i in ones: w[i] = -1
            best = min(tuple(w[i:] + w[:i]) for i in range(T))
            seen.add(best)
        out += sorted(seen)
    return out

def cuts(w):
    """内部归零点：前缀和回到 0 的位置（不含 0 与 T）"""
    s = 0; out = []
    for i, x in enumerate(w):
        s += x
        if s == 0 and 0 < i + 1 < len(w): out.append(i + 1)
    return out

def mode_of(seq):
    T = len(seq)
    return min(tuple(seq[i:] + seq[:i]) for i in range(T))

def children(w, rule):
    c = cuts(w)
    if rule == "persist": return [w]
    if rule == "copy": return [w, w]
    if rule == "split_parent":
        r = [w]
        for k in c: r += [mode_of(w[:k]), mode_of(w[k:])]
        return r
    if rule == "split_only":
        r = []
        for k in c: r += [mode_of(w[:k]), mode_of(w[k:])]
        return r
    if rule == "split_all":
        r = []
        n = len(c)
        for mask in range(1 << n):
            ks = [c[i] for i in range(n) if mask >> i & 1]
            if not ks: continue
            prev = 0
            for k in ks + [len(w)]:
                r.append(mode_of(w[prev:k])); prev = k
        return r
    raise ValueError(rule)

M = modes(); idx = {m: i for i, m in enumerate(M)}
print(f"模式数 = {len(M)}（周期<=12 的循环平衡词旋转类）")
res = {}
for rule in ("persist", "copy", "split_parent", "split_only", "split_all"):
    T = np.zeros((len(M), len(M)), dtype=float)
    for p, w in enumerate(M):
        for ch in children(w, rule):
            T[idx[ch], p] += 1
    ev = np.linalg.eigvals(T)
    rho = float(np.max(np.abs(ev)))
    nz = int(np.sum(np.abs(ev) < 1e-9))
    # 幂零指数
    P = T.copy(); nil = None
    for k in range(1, 40):
        if not P.any(): nil = k; break
        P = P @ T
    res[rule] = {"rho": round(rho, 6), "n_zero_eig": nz, "nilpotent_index": nil}
    print(f"  {rule:<13} ρ={rho:.3f}  零特征值={nz:>3}  幂零指数={nil}")
json.dump(res, open(os.path.join(OUT, "z0_repro.json"), "w"), ensure_ascii=False, indent=1)

# ---- split_parent：novelty 轨迹（新模式出现率 + 成分是否收敛）----
rule = "split_parent"
T = np.zeros((len(M), len(M)))
for p, w in enumerate(M):
    for ch in children(w, rule): T[idx[ch], p] += 1
v = np.zeros(len(M)); v[idx[mode_of([1, 1, -1, -1])]] = 1.0   # T=4 的 ++-- 起
print("\nsplit_parent 的 novelty 轨迹（从 ++-- 起）：")
print(f"{'代':>3} {'总量':>7} {'非零模式数':>9} {'本代新模式':>9} {'最大占比':>9}")
prev = set(np.nonzero(v)[0].tolist()); newmodes = []
for k in range(18):
    v = T @ v
    cur = set(np.nonzero(v)[0].tolist())
    new = len(cur - prev); newmodes.append(new)
    tot = v.sum()
    print(f"{k+1:>3} {tot:>7.0f} {len(cur):>9} {new:>9} {v.max()/max(tot,1e-30):>9.4f}")
    prev = cur
print(f"\n17 代内新模式总数={len(prev)}  每代新增={newmodes}")
print("读法：新增长期为 0 ⟹ novelty 用尽（成分收敛）；长期 >0 ⟹ 持续产生新东西。")
