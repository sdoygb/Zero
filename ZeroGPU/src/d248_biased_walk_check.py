#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
d248_biased_walk_check.py —— 能否同时拿到「D248 形式 + 守恒 + 稳定」
====================================================================
背景（规矩 1 先查库）：
  · **Z002 引理 4【条件】**：无偏好 ⟹ 均匀导纳（路径图 P4/P5 给不出 ⟹ 非导出）
  · **Z093 §四**：力 = 势的梯度（Zero 有 E=fᵀLf）＋ 时间演化律（Zero 无，Z005）⟹"给了方向，没给力"
  · **D65**：新库忠实 L2 的输运由密度驱动（R²=1），不由读出势 φ 驱动（R²<1）

本脚本测三条路，判据是「D248 形式（净流 ∝ φ_i−φ_j，K 对称唯一）」＋「总荷精确守恒」＋「稳定」：
  (i)  偏置随机步（Boltzmann / 线性梯度）—— 保持"走者"图像的候选
  (ii) 确定性交换更新 dq_i/dt = −Σ_j J_{i→j}，J_{i→j} = K_ij(φ_i−φ_j)  —— D248 原文的更新
  (iii) (ii) 的自洽闭环稳定性（吸引 vs 排斥）

运行：cd ZeroGPU/src && /usr/bin/python3 d248_biased_walk_check.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import torus_adj, path_adj      # noqa: E402

RES = os.path.join(HERE, "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


def lap(adj):
    V = len(adj)
    L = np.zeros((V, V))
    for i, nb in enumerate(adj):
        L[i, i] = len(nb)
        for j in nb:
            L[i, j] -= 1.0
    return L


def solve_potential(L, x):
    """Z033 约定：Lφ = −x ⟹ φ = −L⁺x（零和）"""
    phi = -np.linalg.pinv(L) @ x
    return phi - phi.mean()


def biased_A(adj, phi, beta, mode):
    V = len(adj)
    A = np.zeros((V, V))
    for i, nb in enumerate(adj):
        w = np.array([np.exp(-beta * phi[j]) if mode == "boltzmann"
                      else 1.0 + beta * (phi[i] - phi[j]) for j in nb])
        if (w <= 0).any():
            return None
        w = w / w.sum()
        for k, j in enumerate(nb):
            A[i, j] = w[k]
    return A


def edges(adj):
    out = []
    for i in range(len(adj)):
        for j in adj[i]:
            if j > i:
                out.append((i, j))
    return out


def flow_stats(adj, phi, Jvals):
    """Jvals: 每条边 (i,j) 的净流 J_{i→j}。回归 J ~ (φ_i−φ_j)，返回 (R², 斜率, K 相对离散度)"""
    e = edges(adj)
    dphi = np.array([phi[i] - phi[j] for i, j in e])
    Jn = np.array(Jvals)
    ss_tot = float(((Jn - Jn.mean()) ** 2).sum())
    if ss_tot == 0:
        return 1.0, 0.0, 0.0
    slope = float((dphi @ Jn) / (dphi @ dphi)) if (dphi @ dphi) > 0 else 0.0
    r2 = 1.0 - float(((Jn - slope * dphi) ** 2).sum()) / ss_tot
    Kimp = Jn / dphi
    spread = float(np.std(Kimp) / max(1e-30, abs(np.mean(Kimp))))
    return r2, slope, spread


# ============================================================
print("=" * 76)
print("第 0 步：台账定位")
print("=" * 76)
print("  Z002 引理 4【条件】：无偏好 ⟹ 均匀导纳（P4/P5 给不出 ⟹ 非导出）")
print("  Z093 §四：力 = 势的梯度（有）＋ 时间演化律（无，Z005）⟹『给了方向，没给力』")
print("  D65：忠实 L2 输运由密度驱动（R²=1），不由读出势 φ 驱动（R²<1）")

# ============================================================
print()
print("=" * 76)
print("第 1 步：(i) 偏置随机步 —— 净流是否满足 D248 形式？")
print("=" * 76)
print("  均匀密度 q≡1 ⟹ J_ij = p_{i→j} − p_{j→i}")
for gname, adj in [("cycle8", torus_adj(1, 8)), ("torus4x4", torus_adj(2, 4))]:
    V = len(adj)
    rng = np.random.default_rng(31)
    x = rng.normal(size=V); x -= x.mean()
    phi = solve_potential(lap(adj), x)
    print(f"\n  ── {gname} ──   {'mode':<10}{'β':>6} | {'R²':>10} | {'斜率':>9} | {'K 离散度':>9} | {'守恒偏差':>9}")
    for mode in ("boltzmann", "linear"):
        for beta in (0.05, 0.2, 0.5):
            A = biased_A(adj, phi, beta, mode)
            if A is None:
                continue
            Jv = [A[i, j] - A[j, i] for i, j in edges(adj)]
            r2, slope, spread = flow_stats(adj, phi, Jv)
            cons = float(np.abs(A.sum(axis=1) - 1).max())
            print(f"                 {mode:<10}{beta:>6} | {r2:>10.6f} | {slope:>9.6f} | {spread:>9.4f} | {cons:>9.1e}")
            if gname == "cycle8" and mode == "linear" and abs(beta - 0.2) < 1e-9:
                ck("F1a 偏置随机步的净流**不**满足 D248 形式（R² 明显 < 1）",
                   r2 < 0.99, f"R² = {r2:.6f}，K 相对离散度 = {spread:.3f}")
                note("R² 不随 β→0 恢复", "cycle8 上 β=0.05 与 0.5 的 R² 几乎相同（0.660 vs 0.677）",
                     "说明破坏不是高阶项，而是**行归一化的 site-dependent 因子**"
                     "（p_{i→j} ∝ 1/(deg_i(1+β(φ_i−φ̄_i)))，多出邻域均值项）")

# ============================================================
print()
print("=" * 76)
print("第 2 步：守恒（两种偏置都是行随机矩阵）")
print("=" * 76)
adj = torus_adj(1, 8); V = len(adj)
rng = np.random.default_rng(41)
x = rng.normal(size=V); x -= x.mean()
phi = solve_potential(lap(adj), x)
ok_cons = True
for mode, beta in [("boltzmann", 0.2), ("linear", 0.2)]:
    A = biased_A(adj, phi, beta, mode)
    q = rng.random(V) * 3 + 0.1
    qa = q @ A
    print(f"  {mode:<10}: max|Σ_j A_ij − 1| = {np.abs(A.sum(axis=1)-1).max():.2e}；"
          f"总荷 {q.sum():.8f} → {qa.sum():.8f}（差 {abs(qa.sum()-q.sum()):.2e}）")
    ok_cons &= abs(qa.sum() - q.sum()) < 1e-13
ck("F2a 偏置输运总荷精确守恒", ok_cons, "行随机矩阵 ⟹ 精确（与偏置形式无关）")

# ============================================================
print()
print("=" * 76)
print("第 3 步：(ii) 确定性交换更新 —— D248 原文的更新，是否精确满足？")
print("=" * 76)
print("  dq_i/dt = −Σ_j J_{i→j}，J_{i→j} = K_ij(φ_i−φ_j)，K = A（对称、≥0、零对角）")
for gname, adj in [("cycle8", torus_adj(1, 8)), ("torus4x4", torus_adj(2, 4))]:
    V = len(adj)
    rng = np.random.default_rng(53)
    x = rng.normal(size=V); x -= x.mean()
    phi = solve_potential(lap(adj), x)
    Jv = []
    for i, j in edges(adj):
        Jv.append(phi[i] - phi[j])                 # K_ij = 1（边存在）
    r2, slope, spread = flow_stats(adj, phi, Jv)
    # 反对称 + 守恒 + 与 −Lφ 一致
    Jmat = np.zeros((V, V))
    for i, j in edges(adj):
        Jmat[i, j] = phi[i] - phi[j]
        Jmat[j, i] = phi[j] - phi[i]
    asym = float(np.abs(Jmat + Jmat.T).max())
    dq = Jmat.sum(axis=0)                          # Σ_i J_{i→j}
    Lphi = lap(adj) @ phi
    err = float(np.abs(dq - (-Lphi)).max())        # dq_j/dt = Σ_i J_{i→j} = −(Lφ)_j
    print(f"\n  {gname}: R² = {r2:.12f}；K 相对离散度 = {spread:.2e}；"
          f"反对称 max|J+Jᵀ| = {asym:.2e}；|Σ_i J_{{i→j}} + (Lφ)_j| = {err:.2e}")
    if gname == "cycle8":
        ck("F3a 交换型更新**精确**满足 D248 形式（R² = 1，K 唯一＝A）",
           r2 > 1 - 1e-12 and spread < 1e-12, f"R² = {r2:.12f}，K 离散度 = {spread:.2e}")
        ck("F3b 交换型更新精确反对称且精确守恒（代数恒等）",
           asym < 1e-14 and err < 1e-14, f"max|J+Jᵀ| = {asym:.2e}；连续性方程差 = {err:.2e}")

# ============================================================
print()
print("=" * 76)
print("第 4 步：(iii) 交换闭环的稳定性 —— 吸引 vs 排斥（自洽势）")
print("=" * 76)
print("  自洽：φ 由当前荷解得。**Z033 约定 Lφ = −x ⟹ φ = −L⁺x（势井＝吸引）**；φ = +L⁺x（势丘＝排斥）")
print("  更新：q ← q + dt·(−Lφ)。因 q 线性，偏差按常数因子逐代变化 ⟹ 测 log std 的斜率")
print("  预期：吸引（Lφ=−x）⟹ 斜率 = log(1+dt) > 0（增长）；排斥（Lφ=+x）⟹ log(1−dt) < 0（衰减）")

for gname, adj in [("cycle8", torus_adj(1, 8)), ("torus4x4", torus_adj(2, 4))]:
    V = len(adj); L = lap(adj)
    dt = 0.05
    print(f"\n  ── {gname} ──   dt = {dt}")
    slopes = {}
    for tag, sign in [("吸引 φ=−L⁺x（Z033 约定）", +1.0), ("排斥 φ=+L⁺x", -1.0)]:
        rng = np.random.default_rng(67)
        q = rng.random(V) + 1.0
        devs = []
        for t in range(60):
            xv = q - q.mean()
            phi = sign * solve_potential(L, xv)      # sign=+1 ⟹ φ = −L⁺x（吸引）
            q = q - dt * (L @ phi)
            devs.append(float(np.std(q)))
        t_arr = np.arange(len(devs))
        slope = float(np.polyfit(t_arr, np.log(devs), 1)[0])
        slopes[tag] = slope
        pred = np.log(1 + dt) if sign > 0 else np.log(1 - dt)
        print(f"    {tag:<26}: std {devs[0]:.6f} → {devs[-1]:.6f}；"
              f"log-std 斜率 = {slope:+.6f}（预期 {pred:+.6f}）")
    a = slopes["吸引 φ=−L⁺x（Z033 约定）"]
    r = slopes["排斥 φ=+L⁺x"]
    if gname == "cycle8":
        ck("F4a 吸引（Z033 约定 Lφ=−x）⟹ 自洽交换**不稳定**（偏差指数增长）",
           a > 0 and abs(a - np.log(1 + dt)) < 5e-3,
           f"log-std 斜率 = {a:+.6f}，预期 log(1+dt) = {np.log(1+dt):+.6f}")
        ck("F4b 排斥（φ=+L⁺x）⟹ 偏差指数**衰减到均匀**（归零）",
           r < 0 and abs(r - np.log(1 - dt)) < 5e-3,
           f"log-std 斜率 = {r:+.6f}，预期 log(1−dt) = {np.log(1-dt):+.6f}")

# ============================================================
print()
print("=" * 76)
print("第 5 步：判定")
print("=" * 76)
print("  (i)  偏置随机步：守恒 ✓ 精确；D248 形式 ✗（R² = 0.66–0.94，K 不唯一）")
print("  (ii) 交换型更新：D248 形式 ✓ 精确（K = A）；守恒 ✓ 精确 —— 但**放弃走者图像**")
print("  (iii) 自洽交换闭环：吸引 ⟹ 偏差增长（聚类/塌缩）；排斥 ⟹ 衰减归零")
print("        ⟹ **没有稳定非平凡稳态**（与 D57「势反馈终态归零」、Z025「源须外部输入」一致）")
print()
print("  ⟹ 三条路不能同时拿到：守恒＋φ 驱动＋稳定。要稳定，必须有**外部固定的源/势**")
print("     （或压强支撑——对接 D250 的最小类 p = ρ 刚性物态、GR 的 Jeans 不稳定性）")

n_fail = sum(1 for x in CK if not x["pass"])
print()
print("=" * 76)
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 76)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "d248_biased_walk_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "verdict": {
                   "biased_walk": "守恒精确，但净流不满足 D248 形式（行归一化带 site-dependent 因子）",
                   "exchange_update": "精确满足 D248 形式（K=A）且精确守恒；但非随机步",
                   "self_consistent_loop": "吸引⟹增长（不稳定）；排斥⟹归零；无稳定非平凡稳态",
                   "conclusion": "守恒＋φ驱动＋稳定 三者不可兼得；要稳定需外部固定源/势（或压强，D250 p=ρ）"}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/d248_biased_walk_check.json")
sys.exit(1 if n_fail else 0)