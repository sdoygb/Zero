# -*- coding: utf-8 -*-
"""Z132 核验：闭合圆正则分划 → 内部多吸引子 2cos 族谱。
A. 2cos 族存在（n=3..12, 比值>1）
B. 每个比值是稳定吸引子（忠实引擎实测频率=logR）
C. 频率谱 = {log|2cos|}
D. K4 三槽对照：log2 在族内, 1.5/2.667 不在（外部输入）
"""
import sys, numpy as np
SRC = "/Users/oygb/Downloads/Zero/ZeroGPU/src"
sys.path.insert(0, SRC)
from loyalty_engine import LoyaltyEngine

def run_freq(adj, seed, gens=1200, T=4):
    eng = LoyaltyEngine(adj, seed=seed, T=T, retention="top2")
    for _ in range(gens): eng.run_generation()
    rec = eng.record / eng.record.sum(); rec = rec[rec > 1e-9]
    if len(rec) < 2: return None
    idx = np.argsort(rec)[-2:]
    return abs(float(np.log(rec[idx[1]] / rec[idx[0]])))

def two_uniform_blocks(R, block=3, weak=1):
    n = 2 * block; adj = [[] for _ in range(n)]
    for i in range(block):
        for j in range(block):
            if i != j: adj[i].append(j)
    for i in range(block, 2 * block):
        for j in range(block, 2 * block):
            if i != j: adj[i].append(j)
    for _ in range(weak): adj[0].append(block); adj[block].append(0)
    return adj, [R] * block + [1.0] * block

print("=" * 84)
print("Z132 核验：闭合圆正则分划 → 内部多吸引子 2cos 族谱")
print("=" * 84)
ok_all = True

# A: 2cos 族存在
fam = {}
for n in range(3, 13):
    for m in range(1, n // 2 + 1):
        v = abs(2 * np.cos(2 * np.pi * m / n))
        if v > 1.0001: fam.setdefault(round(v, 6), (n, m))
distinct = sorted(fam.items())
print(f"[A] 2cos族: n=3..12 独立比值数={len(distinct)}  {'✓' if len(distinct)>=10 else '✗'}")
ok_all &= len(distinct) >= 10

# B+C: 每比值稳定、频率=logR
print("[B/C] 每比值稳定吸引子、频率=logR:")
stable = 0
for v, (n, m) in distinct:
    fr = run_freq(*two_uniform_blocks(v))
    match = fr is not None and abs(fr - np.log(v)) < 0.05
    stable += 1 if match else 0
    print(f"     R={v:.4f}(n={n},m={m}): 实测={fr and round(fr,4)} vs logR={np.log(v):.4f}  "
          f"{'✓' if match else '✗'}")
print(f"     稳定数={stable}/{len(distinct)}  {'✓' if stable==len(distinct) else '✗'}")
ok_all &= stable == len(distinct)

# D: K4 三槽对照
internal = set(round(v, 4) for v, _ in distinct)
internal_f = set(round(float(np.log(v)), 4) for v, _ in distinct)
k4 = [0.4055, 0.6931, 0.9808]
print("\n[D] K4 三槽 vs 2cos族:")
d_ok = True
for f in k4:
    in_int = any(abs(f - g) < 1e-3 for g in internal_f)
    print(f"     {f:.4f} (R={np.exp(f):.3f}) {'✓在族内(内部)' if in_int else '✗不在族(外部)'}")
d_ok = (any(abs(0.6931 - g) < 1e-3 for g in internal_f) and
        not any(abs(0.4055 - g) < 1e-3 for g in internal_f) and
        not any(abs(0.9808 - g) < 1e-3 for g in internal_f))
print(f"     log2在族、1.5/2.667不在:  {'✓' if d_ok else '✗'}")
ok_all &= d_ok

print("=" * 84)
print(f"【判定】{'✓ Z132成立: 闭合圆正则分划=内部多吸引子, 2cos族谱全稳定, K4仅log2内部' if ok_all else '✗'}")
print("=" * 84)
