# -*- coding: utf-8 -*-
"""Z087 核验：稳态二分支 2:1 与内禀频率 log 2（8 条断言）。"""
import sys, os, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "ZeroGPU", "src"))
sys.path.insert(0, SRC)
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

CHECKS = []


def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond), detail))
    print(f"[{'OK' if cond else 'FAIL'}] {name} {detail}")


# A: path8/12/16 稳态 q=5/9、主两支 2:1、log(p/q)=log2
for n in (8, 12, 16):
    adj = path_adj(n); V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]
    vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    for _ in range(60000):
        vec.run_generation()
    w = vec.record / vec.record.sum()
    q = float((w * w).sum())
    nz = sorted(((i, float(w[i])) for i in range(V) if w[i] > 1e-9), key=lambda x: -x[1])
    p, qq = nz[0][1], nz[1][1]
    check(f"path{n}: q=5/9", abs(q - 5 / 9) < 1e-6, f"q={q:.6f}")
    check(f"path{n}: 主两支比=2", abs(p / qq - 2.0) < 1e-6, f"p/q={p/qq:.6f}")
    check(f"path{n}: log(p/q)=log2", abs(math.log(p / qq) - math.log(2)) < 1e-6,
          f"log={math.log(p/qq):.6f}")

# B: path24/32 主两支仍 2:1
for n in (24, 32):
    adj = path_adj(n); V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]
    vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    for _ in range(60000):
        vec.run_generation()
    w = vec.record / vec.record.sum()
    nz = sorted(((i, float(w[i])) for i in range(V) if w[i] > 1e-9), key=lambda x: -x[1])
    p, qq = nz[0][1], nz[1][1]
    check(f"path{n}: 主两支比=2", abs(p / qq - 2.0) < 1e-3, f"p/q={p/qq:.4f}")
    check(f"path{n}: 有泄漏分支(>2支)", len(nz) > 2, f"n_branch={len(nz)}")

n_pass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n断言通过 {n_pass}/{len(CHECKS)}")
assert n_pass == len(CHECKS), "有断言失败"
print("Z087 核验全部通过。")
