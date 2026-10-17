#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_amplitude_spectrum.py —— 振幅的第二条结构路线：载体的谱比（Z162 §2.5）。

8 槽载体 Γ8 = K4□K2（Z136 §3 的第一层），16 槽载体 Γ16 = K4□K2□K2。
二者的 Laplacian 谱都是等差数列：
  Γ8 :  0, 2, 4, 6          （distinct）
  Γ16:  0, 2, 4, 6, 8       （distinct）
相邻非零谱比：
  Γ8 : 4/2 = 2 = A_ℓ²；6/4 = 3/2 = A_ν²
  Γ16: 同前两比值；第三比 8/6 = 4/3
故两个振幅可由「相邻非零谱比」精确给出（与 Z162 的立方体对角比同值、独立路线）。

负检验：把 ν_R 槽（(v,a,c)=(0,1,1)）从 Γ16 移去、或给它的边加权（0.5/0.75/1.5/2），
谱比都不等于 1.5136 ⟹ 半满缺位不解释 Z164 的 0.9% 残差。

输出：ZeroGPU/results/neutrino_amplitude_spectrum.json
用法：python3 ZeroGPU/src/neutrino_amplitude_spectrum.py
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def laplacian(n, edges, weights=None):
    L = np.zeros((n, n))
    for (i, j) in edges:
        w = 1.0 if weights is None else weights.get((i, j), 1.0)
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    return L


def build(dims):
    """dims 个 K2 因子的 K4□K2^dims；顶点 (v, a_1..a_d)。"""
    keys = [(v,) + a for v in range(4) for a in
            __import__('itertools').product([0, 1], repeat=dims)]
    idx = {k: i for i, k in enumerate(keys)}
    edges = set()
    for k in keys:
        i = idx[k]
        for w in range(4):
            if w != k[0]:
                k2 = (w,) + k[1:]
                edges.add((min(i, idx[k2]), max(i, idx[k2])))
        for d in range(dims):
            k2 = list(k)
            k2[1 + d] = 1 - k2[1 + d]
            k2 = tuple(k2)
            edges.add((min(i, idx[k2]), max(i, idx[k2])))
    return idx, edges


def spectrum(idx, edges, weights=None):
    n = len(idx)
    return np.round(np.linalg.eigvalsh(laplacian(n, edges, weights)), 9)


print("== 载体谱（distinct） ==")
specs = {}
for dims in (1, 2):
    idx, edges = build(dims)
    u = sorted(set(spectrum(idx, edges)))
    specs[dims] = (idx, edges, u)
    print(f"   Γ{4 * 2**dims}: n={len(idx)} 谱 {u}")

u8 = specs[1][2]
u16 = specs[2][2]
check("Γ8 谱 = {0,2,4,6}", [round(x, 9) for x in u8] == [0.0, 2.0, 4.0, 6.0], str(u8))
check("Γ16 谱 = {0,2,4,6,8}", [round(x, 9) for x in u16] == [0.0, 2.0, 4.0, 6.0, 8.0], str(u16))

r8 = [u8[i + 1] / u8[i] for i in range(1, len(u8) - 1)]
r16 = [u16[i + 1] / u16[i] for i in range(1, len(u16) - 1)]
print(f"   Γ8 的非零相邻比 {[round(x, 6) for x in r8]}")
print(f"   Γ16 的非零相邻比 {[round(x, 6) for x in r16]}")
check("A_ℓ² = 4/2 = 2（Γ8 第一比）", round(r8[0], 12) == 2.0)
check("A_ν² = 6/4 = 3/2（Γ8/Γ16 第二比）", round(r8[1], 12) == 1.5 and round(r16[1], 12) == 1.5)

print("\n== 负检验：半满缺位不产生 1.5136 ==")
idx, edges, _ = specs[2]
nuR = idx[(0, 1, 1)]
# 移除 ν_R 顶点
keep = [k for k in range(len(idx)) if k != nuR]
mp = {k: j for j, k in enumerate(keep)}
e2 = {(mp[i], mp[j]) for (i, j) in edges if i in mp and j in mp}
w = sorted(set(spectrum({k: 0 for k in keep}, e2)))
print(f"   Γ16 − ν_R 的 distinct 谱（前 6）: {w[:6]}")
check("移除 ν_R 后，相邻比里没有 1.5136",
      all(abs(a / b - 1.5136) > 0.02 for a, b in zip(w[1:], w[:-1]) if b > 1e-9))
for wt in (0.5, 0.75, 1.5, 2.0):
    weights = {(i, j): wt for (i, j) in edges if i == nuR or j == nuR}
    u = sorted(set(spectrum(idx, edges, weights)))
    near = [a / b for a, b in zip(u[1:], u[:-1]) if b > 1e-9]
    # 振幅对应的是「第二相邻比」（Γ16 里 λ3/λ2）；加权后谱不再是等差，取第二条比作对照
    second = near[1] if len(near) > 1 else float('nan')
    print(f"   ν_R 边权 {wt}: 前四个谱比 {[round(x, 4) for x in near[:4]]}；第二比 {second:.4f}")
    check(f"边权 {wt} 的第二比不等于 1.5136（>1%）", abs(second / 1.5136 - 1) > 0.01,
          f"{second:.4f}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL,
       "spectra": {"G8": [float(x) for x in u8], "G16": [float(x) for x in u16]},
       "ratios": {"G8": [float(x) for x in r8], "G16": [float(x) for x in r16]},
       "A2_charged": float(r8[0]), "A2_neutrino": float(r8[1]),
       "vacancy_no_15136": True}
with open(os.path.join(OUT, "neutrino_amplitude_spectrum.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_amplitude_spectrum.json")
sys.exit(0 if not FAIL else 1)
