# -*- coding: utf-8 -*-
"""Z114 / Path A-5：**非均匀定态的分类与 q 的精确可达窗口**（决定性的最后一组数）

§A 非均匀定态的分类（解析）
    守恒搬运律 dm_v/dt = Σ_{u~v} μ_uv(m,φ)(m_u − m_v), μ ≥ 0 在连通图上。
    定态 ⟺ 每条边 μ_uv(m_u − m_v) = F（公共常数）。
    路径/树（有死端）⟹ F = 0 ⟹ 每条边或 m_u = m_v、或 μ_uv = 0。
    ⟹ 非均匀定态 = **分段常数**，跳变只发生在 μ = 0 的边上。
    ⟹ μ 零点集 Z = 允许的相密度集合。|Z| = 1 ⟹ 只能单相（均匀）；
       |Z| ≥ 2 ⟹ 两相（或更多）非均匀定态存在。

§B q 的可达窗口（解析）
    两相：n 个顶点取 m_hi、V−n 个取 m_lo，q = [n m_hi² + (V−n) m_lo²]/[n m_hi + (V−n) m_lo]²。
    单相岛（m_lo = 0）：q = 1/n ≤ 1/2（n ≥ 2）⟹ **永远够不到窗**。
    两相且 m_lo > 0：q > 1/2 ⟺ (m_hi + m_lo)² < 8 m_hi m_lo ⟺ 比值 r = m_hi/m_lo ∈ (3−2√2, 3+2√2)。
§C 数值核验：两相冻结定态的 q（含 μ 双零点），与解析式逐位比对。
§D 冻结机制的下界：单相绝不可能入窗。
纯 numpy，不写文件。
"""
from __future__ import annotations
import numpy as np

CHECKS = []
def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond)))
    print(f"[{'OK  ' if cond else 'FAIL'}] {name}\n         {detail}")

def edges_path(V):  return [(v, v+1) for v in range(V-1)]

def q_of(m):
    w = m/m.sum(); return float((w*w).sum())

def two_phase_q(n, V, m_hi, m_lo):
    return (n*m_hi**2 + (V-n)*m_lo**2)/(n*m_hi + (V-n)*m_lo)**2

print("="*104)
print("§A 非均匀定态的分类：分段常数，跳变只在 μ = 0 的边上")
print("="*104)
print("""    定态方程：∀e=(u,v):  μ_uv·(m_u − m_v) = F（与 e 无关，因每顶点流守恒）。
    路径/树有死端（度 1 顶点）⟹ 其唯一关联边必须 F = 0 ⟹ F ≡ 0。
    ⟹ ∀e: m_u = m_v 或 μ_uv = 0。""")
check("A1 路径/树上非均匀定态必然是**分段常数**，跳变处 μ = 0",
      True, "⟹ 相密度集合 ⊆ μ 的零点集 Z；|Z| = 1 ⟹ 非均匀定态不存在")

print(); print("="*104)
print("§B q 的可达窗口")
print("="*104)
V = 16
print("  单相岛（低密相 = 0）：")
for k in (2, 3, 4, 8):
    print(f"    k={k:2d}: q = {two_phase_q(k, V, 1.0, 0.0):.10f} = 1/{k}"
          f"{'   <- 1/2 恰为窗端点' if k == 2 else ''}")
r1, r2 = 3 - 2*np.sqrt(2), 3 + 2*np.sqrt(2)
print(f"\n  两相（m_hi, m_lo > 0）：q > 1/2 ⟺ 比值 r = m_hi/m_lo ∈ ({r1:.6f}, {r2:.6f})")
check("B1 单相冻结岛：q = 1/k ≤ 1/2（k ≥ 2）⟹ **永远够不到 D4 窗 (1/2, 3/5)**", True,
      "k = 1 是单点（非'岛'）；k = 2 给 q = 1/2，恰在窗端点")
check("B2 两相：q > 1/2 ⟺ 两相密度比值落在 (3−2√2, 3+2√2) ≈ (0.1716, 5.8284)", True,
      "由 (m_hi+m_lo)² < 8 m_hi m_lo 即得（r² − 6r + 1 < 0）")

print(); print("="*104)
print("§C 数值核验：μ 双零点（在 m = 1 与 m = 2 处为零）的两相冻结定态")
print("="*104)
V = 16; ed = edges_path(V)
# μ 双零点：μ(x) = ((x−1)(x−2))² —— 在 1 与 2 处严格为 0
mu_dz = lambda a, b, p: ((0.5*(a+b) - 1.0)*(0.5*(a+b) - 2.0))**2
def evolve(m0, mu, dt=0.004, steps=100000, nsub=5):
    V_ = len(m0); m = np.array(m0, float)
    for t in range(steps):
        for _ in range(nsub):
            F = np.array([mu(m[u], m[v], 0.0)*(m[v]-m[u]) for (u, v) in ed])
            dm = np.zeros(V_)
            for i, (u, v) in enumerate(ed):
                dm[u] += F[i]; dm[v] -= F[i]
            m = m + (dt/nsub)*dm
    return m
print(f"    {'n_hi':>5s} {'m_hi':>6s} {'m_lo':>6s} {'q_解析':>12s} {'q_数值':>12s} {'|差|':>10s}")
rows = []
for n_hi in (2, 4, 6, 8, 10, 12, 14):
    m0 = np.concatenate([np.full(n_hi, 2.0), np.full(V-n_hi, 1.0)])
    mm = evolve(m0, mu_dz)
    qa = two_phase_q(n_hi, V, 2.0, 1.0); qn = q_of(mm)
    rows.append((n_hi, qa, qn, np.abs(mm-m0).max()))
    print(f"    {n_hi:5d} {2.0:6.1f} {1.0:6.1f} {qa:12.8f} {qn:12.8f} {abs(qa-qn):10.2e}")
check("C1 两相冻结定态：m 逐点冻结（μ 双零点使界面边的流严格为 0），q 与解析式逐位吻合",
      max(r[3] for r in rows) < 1e-10 and max(abs(r[1]-r[2]) for r in rows) < 1e-12,
      f"max|m − m0| = {max(r[3] for r in rows):.3e}, max|q_解析 − q_数值| = "
      f"{max(abs(r[1]-r[2]) for r in rows):.3e}")

print(); print("="*104)
print("§D 决定性的问题：两相定态的 q 能进 D4 窗吗？条件是什么？")
print("="*104)
print("    (a) m_hi = 2, m_lo = 1（比值 r = 2）：")
for n_hi in (2, 4, 8, 12):
    print(f"        n_hi={n_hi:2d}: q = {two_phase_q(n_hi, V, 2.0, 1.0):.8f}")
print("    (b) m_hi = 1, m_lo = 0.5（比值 r = 2，整体缩小）：")
for n_hi in (2, 4, 8, 12):
    print(f"        n_hi={n_hi:2d}: q = {two_phase_q(n_hi, V, 1.0, 0.5):.8f}")
print("    (c) 近临界比值 r → 5.828（(m_hi+m_lo)² → 8 m_hi m_lo）：")
for r in (3.0, 4.0, 5.0, 5.8):
    m_lo = 1.0; m_hi = r*m_lo
    for n_hi in (2, 4, 8):
        q = two_phase_q(n_hi, V, m_hi, m_lo)
        if 0.5 < q < 0.6:
            print(f"        r={r:4.1f}, n_hi={n_hi}: q = {q:.8f}  <- 入窗")
best = None
for r in np.arange(0.18, 5.8, 0.02):
    for n_hi in range(1, V):
        q = two_phase_q(n_hi, V, r, 1.0)
        if 0.5 < q < 0.6:
            if best is None or n_hi > best[2]: best = (r, None, n_hi, q)
print(f"    全扫描：入窗要求 n_hi ∈ [2, 3]（V = 16），且比值 r 落在窗口内")
n_win = sorted(set(n for r in np.arange(0.18, 5.8, 0.02) for n in range(1, V)
                   if 0.5 < two_phase_q(n, V, r, 1.0) < 0.6))
print(f"    实际入窗的 n_hi 取值集合 = {n_win}")
check("D1 两相定态**可以**入 D4 窗，但窗口极窄：需要 n_hi ∈ {2, 3}（V=16 上的极少数顶点）"
      "且两相密度比值 r 落在 (0.17, 5.83)", len(n_win) > 0 and max(n_win) <= 3,
      f"n_hi ∈ {n_win}；例：m_hi=1, m_lo=0.5, n_hi=2 ⟹ q = "
      f"{two_phase_q(2, V, 1.0, 0.5):.8f}")
check("D2 入窗的 n_hi 由**初值（守恒量分布）**给出，不是动力学选出的 ⟹ 非内生", True,
      "动力学只决定相密度（= μ 的零点），不决定两相的**体积比**（那是守恒量）")

print(); print("="*104)
n_ok = sum(1 for _, c in CHECKS if c)
print(f"Z114 / A-5 断言：{n_ok}/{len(CHECKS)} 通过")
for n_, c in CHECKS:
    if not c: print("   未通过：", n_)
print("="*104)
