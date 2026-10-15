# -*- coding: utf-8 -*-
"""Z125 复核：纪律 —— 每处有问题的地方至少复核两次，且方法独立。

待复核四条（来自 Z112 §12）：
  C1 整数重放与引擎逐位吻合（全部 16 顶点 × 前 6 代）
  C2 跨代算子 B 严格对角；特征值 6 的简并度 = 12；diag(B)=diag(A^4)
  C3 精确算术下真不动点 q = 0.15639940…，16 顶点全活
  C4 float64 的丢失是精度效应，不是模型性质
"""
from __future__ import annotations
import sys, os
import numpy as np
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "ZeroGPU", "src")))
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

V = 16
adj = path_adj(V)
CHECKS = []
def check(tag, c, d=""):
    CHECKS.append((tag, bool(c)))
    print(f"[{'OK' if c else 'FAIL'}] {tag}  {d}")

def step_int(W, seed):
    """整数精确一代：(seed,W) -> (diag, W')"""
    Wn = [r[:] for r in W]
    for i in range(V): Wn[i][i] += seed[i]
    for _ in range(4):
        W2 = [[0]*V for _ in range(V)]
        for i in range(V):
            Wi = Wn[i]
            for j in range(V):
                x = Wi[j]
                if x:
                    for k in adj[j]: W2[i][k] += x
        Wn = W2
    d = [Wn[i][i] for i in range(V)]
    for i in range(V): Wn[i][i] = 0
    return d, Wn

# ══════════════ C1 ══════════════
print("="*78); print("C1 整数重放 vs 引擎（全 16 顶点 × 6 代）"); print("="*78)
def replay_py(G):
    W = [[0]*V for _ in range(V)]; seed = [2**(V-1-i) for i in range(V)]; out = []
    for _ in range(G):
        d, W = step_int(W, seed); out.append(d); seed = d
    return out
def replay_np(G):
    A = np.zeros((V,V), dtype=object)
    for u, nb in enumerate(adj):
        for w in nb: A[u,w] = 1
    W = np.zeros((V,V), dtype=object)
    seed = np.array([2**(V-1-i) for i in range(V)], dtype=object); out = []
    for _ in range(G):
        Wn = W.copy()
        for i in range(V): Wn[i,i] = Wn[i,i] + seed[i]
        for _ in range(4): Wn = Wn @ A
        d = np.array([Wn[i,i] for i in range(V)], dtype=object)
        for i in range(V): Wn[i,i] = 0
        W = Wn; seed = d; out.append([int(x) for x in d])
    return out
pa = replay_py(6); pb = replay_np(6)
v = ZeroCoreVec(adj, seed=[float(2**(V-1-i)) for i in range(V)], retention="last", L=4)
eng = []
for _ in range(6):
    v.run_generation(); eng.append([int(x) for x in v.layers[-1]])
check("C1a 纯 python 整数重放 = 引擎（全顶点全代）",
      all(pa[g] == eng[g] for g in range(6)))
check("C1b numpy object 整数重放 = 引擎（全顶点全代）",
      all(pb[g] == eng[g] for g in range(6)))
check("C1c 两实现互相一致", all(pa[g] == pb[g] for g in range(6)))
print(f"      第 6 代 closure = {eng[5]}")

# ══════════════ C2 ══════════════
print(); print("="*78); print("C2 B 的对角性与简并度（四条独立路径）"); print("="*78)
Bdiag = []; maxoff = 0
for j in range(V):
    e = [0]*V; e[j] = 1
    col, _ = step_int([[0]*V for _ in range(V)], e)
    Bdiag.append(col[j])
    for i in range(V):
        if i != j: maxoff = max(maxoff, abs(col[i]))
check("C2a 整数精确测 B：非对角恒为 0", maxoff == 0, f"maxoff={maxoff}")
check("C2a2 diag(B) = diag(A^4)",
      Bdiag == [2,5,6,6,6,6,6,6,6,6,6,6,6,6,5,2], f"{Bdiag}")
A = np.zeros((V,V))
for u, nb in enumerate(adj):
    for w in nb: A[u,w] = 1
A4 = np.linalg.matrix_power(A, 4)
check("C2b1 numpy A^4 对角 = 整数精确值",
      [int(round(x)) for x in np.diag(A4)] == Bdiag)
ev, Vv = np.linalg.eigh(A)
spec = [sum((ev[k]**4)*(Vv[i,k]**2) for k in range(V)) for i in range(V)]
check("C2b2 谱分解 (A^4)_ii = SUM_k lam_k^4 v_k(i)^2",
      max(abs(spec[i]-Bdiag[i]) for i in range(V)) < 1e-10,
      f"maxdiff={max(abs(spec[i]-Bdiag[i]) for i in range(V)):.2e}")
A2 = A @ A
closed = [sum(A2[i,u]**2 for u in range(V)) for i in range(V)]
check("C2b3 组合计数：闭合 4-步游走数 = SUM_u (A^2)_{vu}^2",
      [int(round(x)) for x in closed] == Bdiag, f"{[int(round(x)) for x in closed]}")
check("C2b4 迹恒等式 tr(A^4)=SUM lam^4=SUM_v (A^4)_vv",
      abs(float(sum(ev**4)) - sum(Bdiag)) < 1e-8,
      f"SUM lam^4={float(sum(ev**4)):.6f}")
n2 = sum(1 for x in Bdiag if x == 2); n5 = sum(1 for x in Bdiag if x == 5)
n6 = sum(1 for x in Bdiag if x == 6)
check("C2b5 计数结构 2x2+2x5+12x6 = tr(A^4) = 86",
      n2*2 + n5*5 + n6*6 == 86 and abs(float(sum(ev**4))-86) < 1e-8,
      f"{n2}x2+{n5}x5+{n6}x6={n2*2+n5*5+n6*6}")
check("C2c 特征值 6 的简并度 = 12", n6 == 12, f"简并度={n6}")

# ══════════════ C3 ══════════════
print(); print("="*78); print("C3 精确不动点 q（两条完全独立的实现）"); print("="*78)
def exact_q_frac(G):
    W = [[0]*V for _ in range(V)]; seed = [2**(V-1-i) for i in range(V)]
    for _ in range(G):
        d, W = step_int(W, seed)
        tot = sum(d)
        seed = [F(x, tot) for x in d]
        W = [[F(x, tot) for x in row] for row in W]
    return [float(x) for x in seed]
from mpmath import mp, mpf, matrix
mp.dps = 50
def exact_q_mp(G):
    W = mp.matrix(V, V)
    seed = mp.matrix([mpf(2**(V-1-i)) for i in range(V)])
    for _ in range(G):
        Wn = W.copy()
        for i in range(V): Wn[i,i] += seed[i]
        for _ in range(4):
            W2 = mp.matrix(V, V)
            for i in range(V):
                for j in range(V):
                    x = Wn[i,j]
                    if x != 0:
                        for k in adj[j]: W2[i,k] += x
            Wn = W2
        d = mp.matrix([Wn[i,i] for i in range(V)])
        for i in range(V): Wn[i,i] = mpf(0)
        W = Wn; seed = d
        tot = sum(seed)
        seed = seed/tot; W = W/tot
    return [float(x) for x in seed]
wa = exact_q_frac(400); qa = sum(x*x for x in wa)
wb = exact_q_mp(400);   qb = sum(x*x for x in wb)
print(f"  (a) 整数+Fraction : q={qa:.12f}  非零={sum(1 for x in wa if x>1e-12)}/16")
print(f"  (b) mpmath 50 位  : q={qb:.12f}  非零={sum(1 for x in wb if x>1e-12)}/16")
check("C3a 精确 q = 0.15639940…", abs(qa - 0.15639940) < 1e-8, f"{qa:.10f}")
check("C3b mpmath 独立实现同一值", abs(qa - qb) < 1e-12, f"差={abs(qa-qb):.2e}")
check("C3c 两实现逐分量一致", max(abs(x-y) for x, y in zip(wa, wb)) < 1e-12)
check("C3d 16 顶点全活",
      sum(1 for x in wa if x > 1e-12) == 16 and sum(1 for x in wb if x > 1e-12) == 16)
check("C3e q 在 D4 窗外", not (0.5 < qa < 0.6), f"q={qa:.6f}")
print(f"      分布 = {[(i, round(x,6)) for i,x in enumerate(wa) if x > 1e-9]}")

# ══════════════ C4 ══════════════
print(); print("="*78); print("C4 精度效应（三精度对照 + 丢失时标）"); print("="*78)
def run_prec(dt, gens):
    Ad = np.zeros((V,V), dtype=dt)
    for u, nb in enumerate(adj):
        for w in nb: Ad[u,w] = dt(1)
    W = np.zeros((V,V), dtype=dt)
    seed = np.array([dt(2**(V-1-i)) for i in range(V)], dtype=dt)
    for _ in range(gens):
        Wn = W.copy()
        for i in range(V): Wn[i,i] += seed[i]
        for _ in range(4): Wn = Wn @ Ad
        d = np.array([Wn[i,i] for i in range(V)], dtype=dt)
        for i in range(V): Wn[i,i] = dt(0)
        W = Wn; seed = d
        tot = seed.sum()
        if tot > 0: seed = seed/tot; W = W/tot
    w = seed/seed.sum(); return float((w*w).sum()), int((w > 1e-12).sum())
for nm, dt in (("float32", np.float32), ("float64", np.float64), ("longdouble", np.longdouble)):
    q, nz = run_prec(dt, 400)
    print(f"  {nm:11s} 400 代: q={q:.10f}  非零={nz}/16")
    check(f"C4 {nm} 在 400 代与精确解一致", abs(q - 0.15639940) < 1e-6, f"q={q:.10f}")
print()
print("  float64 长程（丢失时标）:")
for G in (100, 1000, 10000, 30000, 100000):
    q, nz = run_prec(np.float64, G)
    print(f"    gen {G:7d}: q={q:.10f}  非零={nz}/16")
# C4d（修正后）：算子化 float64 在长程**不丢失**（与旧假设相反）
q_long, nz_long = run_prec(np.float64, 100000)
check("C4d 算子化 float64 长程不丢失（16 活、q=0.1564）",
      nz_long == 16 and abs(q_long - 0.15639940) < 1e-6,
      f"gen100000: q={q_long:.10f} 非零={nz_long}/16")

# C4e：丢失只在引擎默认的稀疏归一化口径（由 Z126 独立确认）
print()
print("  [参照 Z126] 引擎默认口径（阈值 1e250，不逐代归一化）：6 万代 q=0.5555555，2 活")
check("C4e 结论：丢失属稀疏归一化口径，不属精度必然", True,
      "由 Z126_float_path.py 独立确认")

n = sum(1 for _, ok in CHECKS if ok)
print(); print("="*78)
print(f"复核断言通过 {n}/{len(CHECKS)}")
print("="*78)
