#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
d248_orientation_check.py —— 冻结图案的"取向"是被动力学选出的，还是从初值继承的？
=================================================================================
D68 发现：容量机制下图案冻结在 Fiedler 子空间，但子空间内的取向（相位）依赖初值。
本脚本判定该取向的**来源**：

逻辑要点：**λ_min 若简并，则简并模的增长速率完全相同** ⟹ 它们在子空间内的**比值是守恒量**
          ⟹ 若取向由初值继承，则"冻结相位 = 初态在 Fiedler 子空间上的投影相位"（严格相等）

三测：
  ① 相位继承：θ_frozen 与 θ_initial（初态的 Fiedler 投影相位）之差（mod 180°）
  ② 微扰灵敏度：初值加 ε·δ，看相位偏移是否 ∝ ε（连续依赖 ⟹ 非"随机选择"）
  ③ 非简并对照：取一个 **非简并** 的模式（如 cycle8 的 λ_min 不简并 / 或构造），看是否不同

运行：cd ZeroGPU/src && /usr/bin/python3 d248_orientation_check.py
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


def run(adj, kappa, Kc, dt=0.02, max_steps=200000, n_cycles=20, amp0=0.05, x0=None):
    V = len(adj); L = lap(adj)
    x = x0.copy(); x = x - x.mean(); x = x / np.linalg.norm(x) * amp0
    pats = []; last2 = []; steps = 0
    while len(pats) < n_cycles and steps < max_steps:
        x = x + dt * (x - kappa * (L @ x)); x -= x.mean(); steps += 1
        if np.max(np.abs(x)) > Kc:
            pat = x / np.linalg.norm(x)
            pats.append(pat.copy()); last2.append(pat.copy())
            if len(last2) > 2:
                last2.pop(0)
            s = np.sum(last2, axis=0); s -= s.mean()
            nn = np.linalg.norm(s)
            x = (s / nn) * amp0 if nn > 0 else np.zeros(V)
    return np.array(pats)


def phase_of(v, Q, fid):
    a, b = Q[:, fid].T @ v
    return float(np.degrees(np.arctan2(b, a))) % 180.0     # 180° 制（±p 同一取向）


# ============================================================
print("=" * 78)
print("第 0 步：逻辑要点")
print("=" * 78)
print("  dx/dt = x − κLx；若 λ_min 简并（两个模同速率），则子空间内**比值守恒**")
print("  ⟹ 若取向由初值继承，则冻结相位应**严格等于**初态在该子空间上的投影相位")

# ============================================================
print()
print("=" * 78)
print("第 1 步：相位继承检验（cycle16，λ_min 简并度 2）")
print("=" * 78)
adj = torus_adj(1, 16); V = len(adj); L = lap(adj)
w, Q = np.linalg.eigh(L)
nzm = w > 1e-9
lmin, lmax = w[nzm].min(), w[nzm].max()
fid = np.where(nzm & (np.abs(w - lmin) < 1e-9))[0]
kappa = 0.5 * (1.0 / lmax + 1.0 / lmin)
print(f"  λ_min = {lmin:.6f}（简并度 {len(fid)}），κ = {kappa:.4f}")
print(f"\n  {'#':>2} | {'初态投影相位':>12} | {'冻结相位':>10} | {'差(mod 180)':>12}")
diffs = []
for k in range(8):
    rng = np.random.default_rng(2000 + k)
    x0 = rng.normal(size=V); x0 -= x0.mean()
    th0 = phase_of(x0 / np.linalg.norm(x0), Q, fid)
    pats = run(adj, kappa, 1.0, x0=x0)
    thf = phase_of(pats[-1], Q, fid)
    d = min(abs(thf - th0), 180 - abs(thf - th0))
    diffs.append(d)
    print(f"  {k:>2} | {th0:>12.4f} | {thf:>10.4f} | {d:>12.6f}")
ck("K1a 冻结相位**严格继承**初态在 Fiedler 子空间上的投影相位",
   max(diffs) < 1e-6, f"8 个初值的最大相位差 = {max(diffs):.2e}°")

# ============================================================
print()
print("=" * 78)
print("第 2 步：微扰灵敏度（相位偏移是否 ∝ ε）")
print("=" * 78)
rng = np.random.default_rng(777)
x_base = rng.normal(size=V); x_base -= x_base.mean(); x_base /= np.linalg.norm(x_base)
th_base = phase_of(x_base, Q, fid)
print(f"  基准相位 = {th_base:.6f}°")
print(f"\n  {'ε':>10} | {'冻结相位':>12} | {'|Δθ|(°)':>12} | {'|Δθ|/ε':>12}")
rows = []
# 注意：**同一个 δ 必须贯穿所有 ε**（否则每次扰动方向不同，比值不可比）—— 这是我的首版 bug
delta = rng.normal(size=V); delta -= delta.mean(); delta /= np.linalg.norm(delta)
for eps in (1e-2, 1e-4, 1e-6, 1e-9):
    x0 = x_base + eps * delta
    x0 = x0 - x0.mean()
    pats = run(adj, kappa, 1.0, x0=x0)
    thf = phase_of(pats[-1], Q, fid)
    d = min(abs(thf - th_base), 180 - abs(thf - th_base))
    rows.append((eps, d))
    print(f"  {eps:>10.0e} | {thf:>12.6f} | {d:>12.3e} | {d/eps:>12.4f}")
# 线性依赖检验：|Δθ| ∝ ε（同一 δ）
if len(rows) >= 3:
    dd = np.array([d for _, d in rows])
    ee = np.array([e for e, _ in rows])
    ratio = dd / ee
    ck("K2a 相位偏移 ∝ ε（同一 δ；连续依赖 ⟹ 非'随机选择'，而是初值数据的连续函数）",
       float(ratio.std() / ratio.mean()) < 0.05,
       f"|Δθ|/ε = {np.round(ratio,4)}（相对散布 {ratio.std()/ratio.mean():.4f}）")
note("含义", "相位是初值的**连续函数**（不是被动力学随机选出，也不是离散量子化的）",
     "⟹ 简并子空间内的取向 = **未被动力学选择的方向** = 需要额外输入的数据")

# ============================================================
print()
print("=" * 78)
print("第 3 步：非简并对照（path8，λ_min 简并度 1）")
print("=" * 78)
adj2 = path_adj(8); V2 = len(adj2); L2 = lap(adj2)
w2, Q2 = np.linalg.eigh(L2)
nz2 = w2 > 1e-9
lmin2, lmax2 = w2[nz2].min(), w2[nz2].max()
fid2 = np.where(nz2 & (np.abs(w2 - lmin2) < 1e-9))[0]
print(f"  path8: λ_min = {lmin2:.6f}（简并度 {len(fid2)}）")
kappa2 = 0.5 * (1.0 / lmax2 + 1.0 / lmin2)
phs = []
for k in range(6):
    rng = np.random.default_rng(3000 + k)
    x0 = rng.normal(size=V2); x0 -= x0.mean()
    pats = run(adj2, kappa2, 1.0, x0=x0)
    # 简并度 1 ⟹ "取向"只剩 ±，用与唯一 Fiedler 向量的重叠绝对值测
    ov = abs(float(Q2[:, fid2[0]] @ pats[-1]))
    phs.append(ov)
print(f"  6 个初值的 |⟨冻结图案, Fiedler⟩| = {np.round(phs,6)}")
ck("K3a 非简并情形：冻结图案**唯一**（± 之外无自由度）",
   min(phs) > 0.9999, f"最小重叠 = {min(phs):.6f}")

# ============================================================
print()
print("=" * 78)
print("第 4 步：判定")
print("=" * 78)
print("  ⟹ 简并子空间内的取向 **不是被动力学选出的**：")
print("     简并模同速率增长 ⟹ 子空间内比值守恒 ⟹ 冻结相位 = 初态投影相位（严格相等）")
print("  ⟹ 因此准确的说法是：")
print("     · **波长/模式**：被动力学**选择**（其他模衰减，Fiedler 存活）——真正的选择")
print("     · **取向（简并子空间内）**：**不被选择**，是**守恒量**，继承自初值")
print("  ⟹ 结构性判词：**简并 = 动力学未选择的方向 = 需要额外输入的数据**")
print("     （与 D66「守恒＋势驱动＋稳定三者不可兼得」、Z025「源须外部输入」同一族结论）")

n_fail = sum(1 for x in CK if not x["pass"])
print()
print("=" * 78)
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 78)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "d248_orientation_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "phase_inherit_maxdiff": float(max(diffs)),
               "perturbation_ratios": [[float(e), float(d)] for e, d in rows],
               "verdict": {
                   "orientation": "不被动力学选择（简并模同速率 ⟹ 比值守恒 ⟹ 严格继承初值）",
                   "wavelength": "被动力学选择（非 Fiedler 模衰减）",
                   "structural": "简并 = 动力学未选择的方向 = 需要额外输入的数据"}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/d248_orientation_check.json")
sys.exit(1 if n_fail else 0)