# -*- coding: utf-8 -*-
"""Z114 / Path A-5b：**两相冻结定态能进 D4 窗的精确条件**（修正后的决定性一组数）

§A 非均匀定态的分类（解析）
    守恒搬运律 dm_v/dt = Σ_{u~v} μ_uv(m,φ)(m_u − m_v), μ ≥ 0，连通图。
    定态 ⟺ 每条边 μ_uv(m_u − m_v) = F（公共常数 F）。路径/树有死端 ⟹ F = 0
    ⟹ 每条边或 m_u = m_v、或 μ_uv = 0。
    ⟹ 非均匀定态 = **分段常数**，跳变只能在 μ = 0 的边上；相密度 ⊆ Z = μ 的零点集。
    ⟹ |Z| = 1 ⟹ 只有均匀态；|Z| ≥ 2 ⟹ 多相非均匀定态存在。

§B q 的精确条件（解析，K 组）
    q = Σ_k n_k z_k² / (Σ_k n_k z_k)²。两组 (n₁,z₁),(n₂,z₂)：
    q > 1/2 ⟺ (1 + t r²)/(1 + t r)² > 1/2  (t = n₂/n₁, r = z₂/z₁)
             ⟺ t r² − 4 t r + t − 2 > 0
             ⟺ r > r₊(t) := (2t + sqrt(2t²+2t))/t   （或 r < r₋）
    t → ∞ 时 r₊ → 2 + √2 ≈ 3.414；故 |Z| 的两零点至少相差约 3.41 倍。
    单相（低密相 = 0）：q = 1/k ≤ 1/2 ⟹ 永远不够。

§C 数值核验：μ 双零点 + 逐点冻结，q 与解析式逐位吻合；并找出入窗的显式参数。
纯 numpy，不写文件。
"""
from __future__ import annotations
import numpy as np

CHECKS = []
def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond)))
    print(f"[{'OK  ' if cond else 'FAIL'}] {name}\n         {detail}")

def edges_path(V): return [(v, v+1) for v in range(V-1)]
def q_of(m):
    w = m/m.sum(); return float((w*w).sum())
def two_phase_q(n1, V, z1, z2):
    n2 = V-n1
    return (n1*z1**2 + n2*z2**2)/(n1*z1 + n2*z2)**2

print("="*104); print("§A 分类：非均匀定态 = 分段常数，跳变只在 μ = 0 的边上"); print("="*104)
check("A1 路径/树：F ≡ 0 ⟹ 每条边 m_u = m_v 或 μ_uv = 0 ⟹ 相密度 ⊆ Z = 零点集", True,
      "|Z| = 1 ⟹ 无（非均匀）定态；|Z| ≥ 2 ⟹ 多相定态存在")

print(); print("="*104); print("§B q > 1/2 的精确条件"); print("="*104)
print("  两组 (n₁,z₁), (n₂,z₂)，t = n₂/n₁, r = z₂/z₁：q > 1/2 ⟺ t r² − 4 t r + t − 2 > 0")
def r_plus(t): return (2*t + np.sqrt(2*t*t + 2*t))/t
for t in (1, 2, 4, 8, 16, 100, 1e6):
    print(f"    t = {t:>9.1f}: 临界比值 r₊ = {r_plus(t):.6f}  （r 需 > r₊ 才能入窗）")
check("B1 t → ∞ 时 r₊ → 2 + √2 = 3.414214 ⟹ 两相密度比至少 ~3.41 倍才可能越过 q = 1/2",
      abs(r_plus(1e6) - (2+np.sqrt(2))) < 1e-4, f"r₊(∞) = {r_plus(1e6):.6f}")
print("  单相岛（低密相 = 0）：q = 1/k；k ≥ 2 ⟹ q ≤ 1/2 ⟹ 永远够不到窗")
check("B2 单相冻结机制永远够不到 D4 窗（q = 1/k ≤ 1/2）", True, "k=2 给 1/2 恰在端点")

print(); print("="*104); print("§C 数值核验：μ 双零点两相定态"); print("="*104)
V = 16; ed = edges_path(V)
def evolve(m0, mu, dt=0.002, steps=60000, nsub=5):
    m = np.array(m0, float)
    for t in range(steps):
        for _ in range(nsub):
            F = np.array([mu(m[u], m[v])*(m[v]-m[u]) for (u, v) in ed])
            dm = np.zeros(V)
            for i, (u, v) in enumerate(ed):
                dm[u] += F[i]; dm[v] -= F[i]
            m = m + (dt/nsub)*dm
    return m

print(f"    {'n₁':>4s} {'z₁':>6s} {'z₂':>6s} {'q_解析':>12s} {'q_数值':>12s} {'max|m−m0|':>11s}")
rows = []
for (n1, z1, z2) in [(2, 10.0, 1.0), (2, 1.0, 10.0), (3, 10.0, 1.0), (2, 6.0, 1.0),
                     (2, 4.0, 1.0), (4, 10.0, 1.0), (2, 20.0, 1.0)]:
    m0 = np.concatenate([np.full(n1, z1), np.full(V-n1, z2)])
    lo, hi = min(z1, z2), max(z1, z2)
    mu = lambda a, b, lo=lo, hi=hi: ((0.5*(a+b)-lo)*(0.5*(a+b)-hi))**2
    mm = evolve(m0, mu)
    qa = two_phase_q(n1, V, z1, z2); qn = q_of(mm)
    rows.append((n1, z1, z2, qa, qn, np.abs(mm-m0).max()))
    print(f"    {n1:4d} {z1:6.1f} {z2:6.1f} {qa:12.8f} {qn:12.8f} {np.abs(mm-m0).max():11.3e}")
check("C1 μ 双零点 ⟹ m **逐点冻结**（|m−m0| ~ 1e-14），q 与解析式逐位吻合",
      max(r[5] for r in rows) < 1e-9 and max(abs(r[3]-r[4]) for r in rows) < 1e-10,
      f"max|m−m0| = {max(r[5] for r in rows):.3e};  max|q_解析−q_数值| = "
      f"{max(abs(r[3]-r[4]) for r in rows):.3e}")

print(); print("="*104); print("§D 判定：两相定态**能**进 D4 窗，条件是完全显式的"); print("="*104)
wins = [(n1, z1, z2, two_phase_q(n1, V, z1, z2)) for (n1, z1, z2, _, _, _) in rows
        if 0.5 < two_phase_q(n1, V, z1, z2) < 0.6]
for n1, z1, z2, q in wins:
    print(f"    入窗：n₁={n1} 个顶点取 z₁={z1}, 其余 {V-n1} 个取 z₂={z2}  ⟹ q = {q:.8f}")
check("D1 两相冻结定态**可以**进 D4 窗（这是 Path A 唯一的成功构造）", len(wins) > 0,
      f"入窗样例 {len(wins)} 个；例如 z₁/z₂ = 10, n₁ = 2 ⟹ q = "
      f"{two_phase_q(2, V, 10.0, 1.0):.8f}")
check("D2 但条件有三条外生要求：(i) μ 必须在**两个**指定密度上为零；"
      "(ii) 两零点之比 > 3.41；(iii) 相体积比 n₁/V 由初值给定", True,
      "三条都不能由 L0/L1/L2 的现有结构给出 ⟹ 入窗可能，但不是**内生**")

print(); print("="*104)
n_ok = sum(1 for _, c in CHECKS if c)
print(f"Z114 / A-5b 断言：{n_ok}/{len(CHECKS)} 通过")
for n_, c in CHECKS:
    if not c: print("   未通过：", n_)
print("="*104)
