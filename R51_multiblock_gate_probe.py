#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 关键探针：多块谱下 (S) 生存 与 (C) 语境性 能否同时成立？

背景更正
  R44 的界 λ1 ≤ (1+√(2q-1))/2 只在 **3 块归约** 上成立（原文 §2 的 Cauchy–Schwarz
  步用的是"顶三归一"谱，且默认除顶三外无质量）。一旦块数 k>3 且尾部带质量，
  该界失效——反例：w=(0.6,0.2,0.1,0.05,0.05) 给 q=0.415、λ1(顶三归一)=0.6667，
  远超同 q 的 3 块界 0.5556。

因此必须重问：固定 q=Σw² ∈ (1/2,3/5)（R32 生存窗口）时，
      λ1(顶三归一) 的上确界随块数 k 如何变化？能否 > 0.723607（语境）？

方法：随机投影到 {Σw=1, Σw²=q} 的流形上做局部爬山（每步用牛顿法复位约束）。
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 0.723607
OUT = {}


def lam1(w):
    s = np.sort(w)[::-1][:3]
    return float(s[0] / s.sum())


def project(w, q, iters=200):
    """把 w 投到 {Σw=1, Σw²=q, w>0} 上（交替投影 + 平移/缩放）。"""
    w = np.maximum(np.asarray(w, float), 1e-12)
    for _ in range(iters):
        w = w / w.sum()                                  # 单纯形
        cur = float(np.sum(w ** 2))
        if abs(cur - q) < 1e-13:
            break
        # 找 t 使 Σ(w+t)² = q 且 Σ(w+t)=1  ⇒  Σt = 0, 2t·w + t² 调整
        # 用沿 (w - 1/k) 方向的缩放
        d = w - 1.0 / len(w)
        a = float(np.sum(d ** 2))
        b = 2 * float(np.sum(w * d))
        c = cur - q
        if a < 1e-18:
            break
        disc = b * b - 4 * a * c
        if disc < 0:
            w = w / np.sqrt(cur / q)
            continue
        t = (-b + math.sqrt(disc)) / (2 * a)
        w = np.maximum(w + t * d, 1e-12)
    return w / w.sum()


def ascend(k, q, restarts=120, steps=3000, seed=0):
    rng = np.random.default_rng(seed)
    best = (-1, None)
    for _ in range(restarts):
        w = project(rng.random(k) ** rng.choice([0.5, 1, 2, 4]), q)
        cur = lam1(w)
        T = 0.15
        for _ in range(steps):
            cand = w + rng.normal(0, T, k)
            cand = project(cand, q)
            nc = lam1(cand)
            if nc > cur:
                w, cur = cand, nc
            T *= 0.9995
        if cur > best[0]:
            best = (cur, w)
    return best


def scan(k, qs=(0.505, 0.52, 0.55, 0.58, 0.595)):
    rows = []
    for q in qs:
        l, w = ascend(k, q, restarts=40, steps=1500)
        rows.append({"k": k, "q": q, "sup_lam1": round(l, 6),
                     "contextual": bool(l > LAM),
                     "top1_abs": round(float(np.max(w)), 5) if w is not None else None,
                     "effective_blocks": int(np.sum(w > 1e-4)) if w is not None else None})
    return rows


if __name__ == "__main__":
    print("固定 q ∈ R32 生存窗口，块数 k 增大时 λ1 的上确界：")
    print("   k     q      sup λ1     语境?   有效块数")
    allrows = []
    for k in [3, 4, 5, 6, 8, 12, 20, 40]:
        for r in scan(k, qs=(0.505, 0.55, 0.595)):
            allrows.append(r)
            print("  %-4d %.3f  %.6f   %-6s %s"
                  % (r["k"], r["q"], r["sup_lam1"], r["contextual"],
                     r["effective_blocks"]))
    OUT["rows"] = allrows
    ctx = [r for r in allrows if r["contextual"]]
    OUT["n_contextual_inside_survival_window"] = len(ctx)
    OUT["contextual_examples"] = ctx[:8]
    print("\n窗口内达到语境性的 (k,q) 组合数: %d / %d" % (len(ctx), len(allrows)))
    for r in ctx[:6]:
        print("   k=%-3d q=%.3f  λ1=%.6f  有效块=%d"
              % (r["k"], r["q"], r["sup_lam1"], r["effective_blocks"]))
    with open(os.path.join(HERE, "R51_survival_vs_contextuality_multiblock.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("-> R51_survival_vs_contextuality_multiblock.json")
