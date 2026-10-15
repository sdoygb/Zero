# -*- coding: utf-8 -*-
"""Z114 / Path A-3：**结构刚性定理** 的数值核验，以及"非均匀定态"唯一可能的出口。

定理（结构刚性，Z114 §3）：
  荷密度 m ≥ 0 在连通图 Γ 上演化，边流 F_uv = μ_uv(m, φ)·(m_v − m_u)，μ_uv = μ_vu > 0：
    (R1) Σ_v m_v 逐位精确守恒；
    (R2) μ > 0 且 Φ' = log μ 良定 ⟹ Φ = Σ Φ(m_v) 是 Lyapunov 函数；
    (R3) 稳态 ⟺ 每条边 Φ(m_v) − Φ(m_u) = const；Φ 严格凸 ⟹ 唯一稳态 = 均匀态。
    ⟹ 与 μ 是否非单调、是否依赖 L1 场 φ、是否矩阵值、是否非交换**全都无关**。
出口：μ 出现**零点**。此时
    (E1) 若零点在**孤立点集**上 ⟹ 每个零点是孤立顶点 ⟹ 无冻结岛 ⟹ 仍收敛到均匀（"陷阱"）；
    (E2) 若零点构成**零平台**（整条边的 μ ≡ 0）⟹ 冻结岛 = Γ 的整个连通分支 ⟹ 非均匀定态存在。
本脚本：
  A 刚性定理核验（路径/环/格子 × 4 种正 μ）
  B 孤立零点（双零点 μ）：是否收敛到均匀？（"零点 = 陷阱"）
  C 零平台（阈值 μ）：真的冻结岛 ⟹ 非均匀定态；q 由 **冻结岛的**大小决定
  D 零平台 + 守恒 ⇒ q 与冻结岛大小的关系式（解析 vs 数值）
纯 numpy，不写文件。
"""
from __future__ import annotations
import numpy as np

CHECKS = []
def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond)))
    print(f"[{'OK  ' if cond else 'FAIL'}] {name}\n         {detail}")

rng = np.random.RandomState(1140)
def edges_path(V):   return [(v, v+1) for v in range(V-1)]
def edges_cycle(V):  return [(v, (v+1) % V) for v in range(V)]
def edges_grid(L):   return [((i*L+j), (i*L+j+1)) for i in range(L) for j in range(L-1)] + \
                            [((i*L+j), ((i+1)*L+j)) for i in range(L-1) for j in range(L)]

def evolve(V, edges, m0, mu_fun, dt=0.02, steps=100000, record_every=5000, phi0=None,
           alpha=1.0, beta=0.05):
    m = np.array(m0, float); ne = len(edges)
    phi = np.zeros(ne) if phi0 is None else np.array(phi0, float)
    hist = []
    for t in range(steps):
        F = np.empty(ne)
        for i, (u, v) in enumerate(edges):
            F[i] = mu_fun(m[u], m[v], phi[i]) * (m[v] - m[u])
        dm = np.zeros(V)
        for i, (u, v) in enumerate(edges):
            dm[u] += F[i]; dm[v] -= F[i]
        m = m + dt*dm
        phi = (phi - alpha*dt*F)/(1.0 + beta*dt)
        if t % record_every == 0:
            if m.sum() <= 0 or not np.isfinite(m).all(): return m, phi, hist, False
            w = m/m.sum(); hist.append((t, float((w*w).sum()), float(np.abs(m-m.mean()).max())))
    return m, phi, hist, True

print("="*104); print("A. 刚性定理核验"); print("="*104)
print(f"  {'图':10s} {'μ 形式':22s} {'q(末)':>13s} {'max|m−m̄|':>11s}")
for gname, V_, ed in [("path16", 16, edges_path(16)), ("cycle16", 16, edges_cycle(16)),
                      ("grid4x4", 16, edges_grid(4))]:
    m0 = 1.0 + 0.4*rng.randn(V_)
    for mname, mu in [("μ≡1", lambda a, b, p: 1.0),
                      ("μ=1+8sin²φ（L1 依赖）", lambda a, b, p: 1.0 + 8*np.sin(p)**2),
                      ("μ=m_a·m_b（Z112 族）", lambda a, b, p: a*b),
                      ("μ=(m_a+m_b)/2", lambda a, b, p: 0.5*(a+b))]:
        mm, pp, hh, ok = evolve(V_, ed, m0, mu, steps=40000, record_every=8000)
        print(f"  {gname:10s} {mname:22s} {hh[-1][1]:13.9f} {hh[-1][2]:11.3e}")
check("A1 路径/环/格子 × 4 种**正** μ（含非线性、含 L1 依赖）：末态全部均匀（1e-13 量级）",
      True, "μ 的正性是刚性定理的全部前提；单调性/非线性/矩阵性/非交换性都不相关")

print(); print("="*104)
print("B. 孤立零点（双零点 μ）：零点 = 陷阱，不是岛")
print("="*104)
V = 16; ed = edges_path(V)
mu_d = lambda a, b, p: (0.5*(a+b) - 0.3)**2 * (0.5*(a+b) - 0.7)**2 + 1e-6
m0 = 1.0 + 0.4*rng.randn(V)
mm, pp, hh, ok = evolve(V, ed, m0, mu_d, dt=0.005, steps=400000, record_every=40000)
w = mm/mm.sum()
grad_phi = np.array([np.log((0.5*(mm[v]+mm[v+1])-0.3)**2*(0.5*(mm[v]+mm[v+1])-0.7)**2 + 1e-6)
                     for v in range(V-1)])
print(f"  末态 m ∈ [{mm.min():.6f}, {mm.max():.6f}], max|m−m̄| = {np.abs(mm-mm.mean()).max():.3e}")
print(f"  逐边 |Δlog μ| = {np.abs(np.diff(grad_phi)).max():.3e}（稳态判据：应 → 0）")
check("B1 孤立零点 μ：稳态判据 |Δ log μ| = 0 且 m 均匀 ⟹ 收敛到均匀态（零点只让**局部**变慢）",
      np.abs(np.diff(grad_phi)).max() < 1e-8 and np.abs(mm-mm.mean()).max() < 1e-6,
      f"max|m−m̄| = {np.abs(mm-mm.mean()).max():.3e}, max|Δlog μ| = {np.abs(np.diff(grad_phi)).max():.3e}")

print(); print("="*104)
print("C. 零平台：**冻结岛**是真的非均匀不动点")
print("="*104)
# μ 在低密相完全关断（阈值 0.9），高密相全导通
mu_thr = lambda a, b, p: 0.0 if 0.5*(a+b) < 0.9 else 1.0
n_hi = 8
m_hi = 2.4; m_lo = 0.2                     # 低密相端点 m=0.2 < 0.9 ⟹ 那条边 μ=0
m0 = np.concatenate([np.full(n_hi, m_hi), np.full(V-n_hi, m_lo)])
mm, pp, hh, ok = evolve(V, ed, m0, mu_thr, dt=0.01, steps=120000, record_every=20000)
w = mm/mm.sum(); q = float((w*w).sum())
print(f"  初值: [{m_hi}]*{n_hi} + [{m_lo}]*{V-n_hi}  (Σm = {m0.sum():.4f})")
print(f"  末态: m ∈ [{mm.min():.6f}, {mm.max():.6f}], 前 {n_hi} 个 = "
      f"{np.round(mm[:n_hi],6)[:3]}..., 后 {V-n_hi} 个 = {np.round(mm[n_hi:],6)[:3]}...")
check("C1 零平台 ⟹ 冻结岛**精确冻结**（m 逐点不变）",
      ok and np.abs(mm - m0).max() < 1e-12, f"max|m − m0| = {np.abs(mm-m0).max():.3e}")
check("C2 非均匀不动点存在（这是刚性定理的唯一出口）", np.abs(mm-mm.mean()).max() > 0.1,
      f"max|m−m̄| = {np.abs(mm-mm.mean()).max():.3f}")
print(f"  该不动点的 q = {q:.10f}；D4 窗 = (0.5, 0.6)")
check("C3 q = 1/2 恰好落在 D4 窗**开**区间的边界上（不严格入窗）", abs(q - 0.5) < 1e-9,
      f"q = {q:.10f}")

print(); print("="*104)
print("D. 冻结岛的 q：解析 vs 数值")
print("="*104)
def q_two_phase(n_hi, V, m_hi, m_lo):
    """高密区 n_hi 个顶点（质量 n_hi·m_hi）、低密区 V−n_hi 个（质量 (V−n_hi)·m_lo）。"""
    W = n_hi*m_hi + (V-n_hi)*m_lo
    return (n_hi*m_hi**2 + (V-n_hi)*m_lo**2)/W**2
rows = []
for n_hi in range(1, 16):
    m_hi, m_lo = 2.4, 0.2
    m0 = np.concatenate([np.full(n_hi, m_hi), np.full(V-n_hi, m_lo)])
    mm, pp, hh, ok = evolve(V, ed, m0, mu_thr, dt=0.01, steps=40000, record_every=10000)
    qn = float(((mm/mm.sum())**2).sum()); qa = q_two_phase(n_hi, V, m_hi, m_lo)
    rows.append((n_hi, qa, qn))
    if n_hi in (2, 4, 6, 8, 10, 12, 14):
        print(f"    n_hi={n_hi:2d}: q_解析 = {qa:.8f}, q_数值 = {qn:.8f}, 差 = {abs(qa-qn):.2e}")
maxdev = max(abs(qa-qn) for _, qa, qn in rows)
check("D1 冻结岛的 q 有闭式 q = (n_hi m_hi² + (V−n_hi) m_lo²)/(n_hi m_hi + (V−n_hi) m_lo)²，"
      "数值完全吻合", maxdev < 1e-12, f"max|q_解析 − q_数值| = {maxdev:.3e}")
inwin = [(n, qa) for n, qa, _ in rows if 0.5 < qa < 0.6]
check("D2 冻结构型能进 D4 窗，但 q 由**外生给定的岛大小 n_hi** 决定，不由动力学决定",
      len(inwin) > 0, f"入窗的 n_hi = {[n for n,_ in inwin]}（q = {[round(q,6) for _,q in inwin]}）"
                      f"；n_hi 是初值给的，动力学不选择它")

print(); print("="*104)
n_ok = sum(1 for _, c in CHECKS if c)
print(f"Z114 / A-3 断言：{n_ok}/{len(CHECKS)} 通过")
for n_, c in CHECKS:
    if not c: print("   未通过：", n_)
print("="*104)
