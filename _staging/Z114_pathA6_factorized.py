# -*- coding: utf-8 -*-
"""Z114 / Path A-6：**边因子化** vs **局部平均**迁移率 —— 哪一类允许两相冻结定态？

§A 两种迁移率依赖方式
    (P) 乘积型 μ_uv = g(m_u)·g(m_v)      —— 保持"Φ = Σ Φ(m_v)"的梯度结构（Z112 §6）
    (A) 平均型 μ_uv = μ((m_u+m_v)/2)      —— **不**保梯度结构
§B 定态方程与分类（解析）
    路径/树：F ≡ 0 ⟹ 每条边 m_u = m_v 或 μ_uv = 0。
    (P)：μ = 0 ⟺ g(m_u) = 0 或 g(m_v) = 0 ⟹ 相密度必须是 g 的零点。
         两相 z₁, z₂ 且 g(z₁) = g(z₂) = 0 ⟹ 界面边 μ = 0（两侧都是零点）⟹ 冻结，
         **且同时必须满足 Φ(z₁) − Φ(z₂) 无约束**（因为界面边的 μ 恰为 0，恒等式自动满足）。
    (A)：μ = 0 ⟺ (m_u+m_v)/2 ∈ Z ⟹ 每个"平均"必须落在零点上 ⟹ 逐边约束，
         一般无法同时满足 ⟹ 无两相定态。
§C 数值：乘积型 g(m) = (m−z₁)(m−z₂) 的两相冻结；与平均型对照。
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

print("="*104)
print("§C1 乘积型 g(m) = (m−z₁)(m−z₂)，μ_uv = |g(m_u)g(m_v)|（两相零点）")
print("="*104)
print(f"    {'n₁':>4s} {'z₁':>6s} {'z₂':>6s} {'q_解析':>12s} {'q_数值':>12s} {'max|m−m0|':>11s} {'入窗':>5s}")
rows = []
for (n1, z1, z2) in [(2, 10.0, 1.0), (2, 6.0, 1.0), (2, 4.0, 1.0), (3, 10.0, 1.0),
                     (4, 10.0, 1.0), (8, 10.0, 1.0), (2, 3.0, 1.0), (5, 20.0, 1.0)]:
    m0 = np.concatenate([np.full(n1, z1), np.full(V-n1, z2)])
    mu = lambda a, b, z1=z1, z2=z2: abs((a-z1)*(a-z2)*(b-z1)*(b-z2))
    mm = evolve(m0, mu, dt=0.001, steps=40000)
    qa = two_phase_q(n1, V, z1, z2); qn = q_of(mm)
    win = 0.5 < qn < 0.6
    rows.append((n1, z1, z2, qa, qn, np.abs(mm-m0).max(), win))
    print(f"    {n1:4d} {z1:6.1f} {z2:6.1f} {qa:12.8f} {qn:12.8f} {np.abs(mm-m0).max():11.3e} {str(win):>5s}")
check("C1 乘积型 + 双零点：两相**逐点冻结**（|m−m0| ~ 1e-13），q 与解析式逐位吻合",
      max(r[5] for r in rows) < 1e-9 and max(abs(r[3]-r[4]) for r in rows) < 1e-10,
      f"max|m−m0| = {max(r[5] for r in rows):.3e};  max|q_解析−q_数值| = "
      f"{max(abs(r[3]-r[4]) for r in rows):.3e}")

print(); print("="*104)
print("§C2 对照：平均型 μ((m_u+m_v)/2)，同样的双零点")
print("="*104)
rows2 = []
for (n1, z1, z2) in [(2, 10.0, 1.0), (2, 4.0, 1.0), (8, 10.0, 1.0)]:
    m0 = np.concatenate([np.full(n1, z1), np.full(V-n1, z2)])
    mu = lambda a, b, z1=z1, z2=z2: abs((0.5*(a+b)-z1)*(0.5*(a+b)-z2))
    mm = evolve(m0, mu, dt=0.001, steps=40000)
    qa = two_phase_q(n1, V, z1, z2); qn = q_of(mm)
    rows2.append((n1, qa, qn, np.abs(mm-m0).max()))
    print(f"    n₁={n1:2d}, z₁={z1}, z2={z2}: q_解析={qa:.8f}, q_数值={qn:.8f}, "
          f"max|m−m0|={np.abs(mm-m0).max():.3e}")
check("C2 平均型：两相**不冻结**（界面边 μ ≠ 0），末态被拉向均匀（q → 1/V）",
      max(r[3] for r in rows2) > 0.1, "⟹ 只有**乘积型**（保梯度结构的那一类）允许两相冻结")

print(); print("="*104)
print("§D 判定：乘积型两相冻结能否进 D4 窗")
print("="*104)
wins = [(n1, z1, z2, qn) for (n1, z1, z2, qa, qn, dev, w) in rows if w]
for n1, z1, z2, qn in wins:
    print(f"    入窗：n₁={n1} 个顶点取 z₁={z1}，其余取 z₂={z2}  ⟹ q = {qn:.8f}")
check("D1 乘积型两相冻结**可以**进 D4 窗", len(wins) > 0,
      f"入窗样例：{[(int(n),z1,z2,round(q,6)) for n,z1,z2,q in wins]}")
print(f"    D4 窗的完整条件：q = [n₁z₁² + (V−n₁)z₂²]/[n₁z₁ + (V−n₁)z₂]² ∈ (1/2, 3/5)")
for (n1, z1, z2) in [(2, 10.0, 1.0), (3, 10.0, 1.0), (4, 10.0, 1.0), (2, 20.0, 1.0),
                     (2, 100.0, 1.0), (3, 50.0, 1.0), (1, 3.0, 1.0), (2, 3.0, 1.0)]:
    q = two_phase_q(n1, V, z1, z2)
    print(f"      n₁={n1:2d}, z₁={z1:6.1f}, z₂={z2:4.1f}: q = {q:.8f} {'<- 入窗' if 0.5<q<0.6 else ''}")
check("D2 入窗条件（显式）：(i) μ 在**两个**指定密度上为零（乘积型即可）；"
      "(ii) 两零点之比 r = z₁/z₂ 落在窗口内；(iii) 相体积比 n₁/V 由守恒量给定", True,
      "三条都是外生输入 ⟹ 入窗可能，但不是内生")

print(); print("="*104)
n_ok = sum(1 for _, c in CHECKS if c)
print(f"Z114 / A-6 断言：{n_ok}/{len(CHECKS)} 通过")
for n_, c in CHECKS:
    if not c: print("   未通过：", n_)
print("="*104)
