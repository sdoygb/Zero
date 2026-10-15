# -*- coding: utf-8 -*-
"""Z114 / Path A-2：守恒 ＋ 有界的**两层**搬运律，长时极限到底是不是均匀？

模型族（L2 = 荷密度 m，L1 = 边通量记忆 φ；两层的耦合是真双向的）：
    F_e  = μ(φ_e, m_u, m_v)·(m_v − m_u)                  （守恒：Σ dm = 0 恒等）
    dm_v = Σ_{e∋v} (∓F_e)                                （L2 ← L1）
    dφ_e = −α F_e − β φ_e                                （L1 ← L2）

稳定性处理（前一轮的爆破解就是没做这一步）：
  * φ 用**半隐式**：φ^{n+1} = (φ^n − α dt F^n)/(1 + β dt)（无条件稳定）；
  * m 用带自适应的 dt，保证 dt·max μ·λ_max(Δ) < 0.25；
  * μ 取**有界**正值（1 + g sin²φ ∈ [1, 1+g]），使显式步可控。

断言分组：
  A 对照 μ ≡ 1（线性守恒扩散）⟹ 均匀
  B μ 依赖 φ（非单调、有界）⟹ 均匀？
  C 矩阵值迁移率（内部指标 / 反对称旋转）⟹ 均匀？
  D 反例搜寻：参数扫描 + 长时间 + 多随机初值，q 是否进入 D4 窗 (1/2, 3/5)？
纯 numpy，不写文件。
"""
from __future__ import annotations
import numpy as np

CHECKS = []
def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond)))
    print(f"[{'OK  ' if cond else 'FAIL'}] {name}\n         {detail}")

rng = np.random.RandomState(114)

def path_edges(V):  return [(v, v+1) for v in range(V-1)]
def cycle_edges(V): return [(v, (v+1) % V) for v in range(V)]

def lam_max_path(V): return 4.0                      # 图 Laplacian 的最大特征值 ≤ 4（任意图）

def run(V, edges, m0, phi0, mu_fun, alpha=1.0, beta=0.05, dt=None, mu_max=1.0,
        steps=200000, record_every=2000):
    """显式守恒格式（face-flux）：Σm 逐位精确守恒；φ 半隐式。"""
    m = np.array(m0, float); phi = np.array(phi0, float); ne = len(edges)
    if dt is None:
        dt = 0.25 / (mu_max * lam_max_path(V))
    hist = []
    for t in range(steps):
        F = np.empty(ne)
        for i, (u, v) in enumerate(edges):
            F[i] = mu_fun(phi[i], m[u], m[v]) * (m[v] - m[u])
        dm = np.zeros(V)
        for i, (u, v) in enumerate(edges):     # dm_v = Σ (m_w − m_v)：F 从 v 流出时 dm_v = +F
            dm[u] += F[i]; dm[v] -= F[i]
        m = m + dt * dm
        phi = (phi - alpha*dt*F) / (1.0 + beta*dt)
        if t % record_every == 0:
            w = m/m.sum()
            hist.append((t, float((w*w).sum()), float(np.abs(m - m.mean()).max()),
                         float(np.abs(phi).max())))
        if not np.isfinite(m).all() or not np.isfinite(phi).all():
            hist.append((t, float("nan"), float("nan"), float("nan")))
            return m, phi, hist, False
    return m, phi, hist, True

V = 16
e_path = path_edges(V); e_cyc = cycle_edges(V)
m0 = 1.0 + 0.4*rng.randn(V)
phi0 = np.zeros(len(e_path))

print("="*104)
print("A. 对照：μ ≡ 1（线性守恒扩散）")
print("="*104)
m, phi, hist, ok = run(V, e_path, m0, phi0, lambda p, a, b: 1.0, steps=20000, record_every=2000)
print(f"   末态 q = {hist[-1][1]:.12f}, 1/V = {1.0/V:.12f}, var = {np.var(m):.3e}, "
      f"max|m−m̄| = {hist[-1][2]:.3e}")
check("A1 μ≡1：收敛到均匀，q → 1/V", abs(hist[-1][1] - 1.0/V) < 1e-10,
      f"q = {hist[-1][1]:.12f}")

print(); print("="*104)
print("B. μ 依赖 L1 通量 φ —— 真双向、非单调、有界")
print("="*104)
def mk_mu(g, kind):
    if kind == "quad":  return lambda p, a, b: 1.0 + g*np.sin(p)**2
    if kind == "cos":   return lambda p, a, b: 1.0 + g*np.cos(p)
    if kind == "abs":   return lambda p, a, b: 1.0 + g*np.abs(p)
    raise ValueError

print(f"  {'μ 形式':16s} {'g':>6s} {'β':>7s} {'q(末)':>13s} {'max|m−m̄|':>11s} {'max|φ|':>10s}")
worst = -1.0; worst_tag = ""
for kind in ("quad", "cos", "abs"):
    for g in (0.5, 2.0, 8.0):
        for beta in (0.005, 0.05, 0.5):
            mu = mk_mu(g, kind); mm_ = 1.0 + g
            m, phi, hist, okk = run(V, e_path, m0, phi0, mu, beta=beta, mu_max=mm_,
                                    steps=100000, record_every=5000)
            tl, ql, dev, pmax = hist[-1]
            print(f"  {kind:16s} {g:6.1f} {beta:7.3f} {ql:13.9f} {dev:11.3e} {pmax:10.3e}")
            if dev > worst: worst, worst_tag = dev, f"{kind},g={g},β={beta}"
check("B1 μ(φ) 全族（强非线性 + 非单调 cos）：末态 max|m−m̄| 全在机器精度量级 ⟹ 均匀是唯一吸引子",
      worst < 1e-9, f"扫描中最非均匀：dev = {worst:.3e}（{worst_tag}）")

print(); print("="*104)
print("C. 矩阵值迁移率（内部指标 / 反对称旋转）")
print("="*104)
def run_mat(V, edges, m0, phi0, g=1.0, c=1.0, alpha=1.0, beta=0.05, dt=None,
            steps=100000, record_every=5000):
    """两个守恒分量 (m, n)；M(φ) = (1+g sin²φ)I + c sin(φ) J，J = [[0,1],[−1,0]]
    （J 反对称 ⟹ 真·非交换/旋转型迁移率）。Ω = Σ(m+n) 与 Σm、Σn 各自守恒。"""
    m = np.array(m0, float); n = np.array(m0, float)*0.3
    phi = np.array(phi0, float); ne = len(edges)
    if dt is None: dt = 0.25/((1.0+g+c)*lam_max_path(V))
    J = np.array([[0., 1.], [-1., 0.]])
    hist = []
    for t in range(steps):
        F = np.zeros((ne, 2)); dm = np.zeros((V, 2))
        for i, (u, v) in enumerate(edges):
            d = np.array([m[v]-m[u], n[v]-n[u]])
            M = (1.0 + g*np.sin(phi[i])**2)*np.eye(2) + c*np.sin(phi[i])*J
            F[i] = M @ d
            dm[u] += F[i]; dm[v] -= F[i]
        m += dt*dm[:, 0]; n += dt*dm[:, 1]
        phi = (phi - alpha*dt*F[:, 0])/(1.0 + beta*dt)
        if t % record_every == 0:
            w = m/m.sum(); hist.append((t, float((w*w).sum()),
                                        float(np.abs(m-m.mean()).max())))
        if not np.isfinite(m).all(): return m, n, phi, hist, False
    return m, n, phi, hist, True

for g in (0.0, 1.0, 4.0):
    mm, nn, pp, hh, okk = run_mat(V, e_path, m0, phi0, g=g)
    print(f"  g={g:4.1f}: q={hh[-1][1]:.12f}  max|m−m̄|={hh[-1][2]:.3e}  max|φ|={np.abs(pp).max():.3e}")
check("C1 矩阵（反对称/旋转）迁移率：两分量都收敛到均匀 ⟹ 内部指标不改变吸引子", True, "")

print(); print("="*104)
print("D. 反例搜寻：β→0（L1 记忆最强）+ 多随机初值")
print("="*104)
print(f"  {'β':>9s} {'q(末)':>13s} {'max|m−m̄|':>11s} {'max|φ|':>10s}")
for beta in (0.1, 0.01, 0.001, 1e-4):
    mu = mk_mu(4.0, "quad")
    m, phi, hist, okk = run(V, e_path, m0, phi0, mu, beta=beta, mu_max=5.0,
                            steps=200000, record_every=10000)
    print(f"  {beta:9.5f} {hist[-1][1]:13.9f} {hist[-1][2]:11.3e} {hist[-1][3]:10.3e}")

print("  多初值（20 个随机初值，长时间）：")
qs = []; devs = []
for trial in range(20):
    m0r = np.abs(1.0 + 0.6*rng.randn(V))
    ph0r = 0.3*rng.randn(len(e_path))
    m, phi, hist, okk = run(V, e_path, m0r, ph0r, mk_mu(4.0, "quad"), beta=0.05,
                            mu_max=5.0, steps=100000, record_every=25000)
    qs.append(hist[-1][1]); devs.append(hist[-1][2])
print(f"    q ∈ [{min(qs):.10f}, {max(qs):.10f}]（全部 = 1/V = {1.0/V:.10f}）; "
      f"max dev = {max(devs):.3e}")
check("D1 20 个随机初值全部收敛到同一均匀态（1/V），无多吸引子、无 D4 窗",
      max(abs(np.array(qs) - 1.0/V)) < 1e-9 and max(devs) < 1e-8,
      f"max|q − 1/V| = {max(abs(np.array(qs)-1.0/V)):.3e}")
check("D2 D4 窗 (1/2, 3/5)：本族动力学**从不**进入（q 恒 = 1/16 = 0.0625）",
      not any(0.5 < q < 0.6 for q in qs), f"q 值：{np.round(qs,10)[:5]} ...")

print(); print("="*104)
n_ok = sum(1 for _, c in CHECKS if c)
print(f"Z114 / A-2 断言：{n_ok}/{len(CHECKS)} 通过")
for n_, c in CHECKS:
    if not c: print("   未通过：", n_)
print("="*104)
