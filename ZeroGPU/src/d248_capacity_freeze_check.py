#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
d248_capacity_freeze_check.py —— 把 Z022 §277 的容量饱和接进交换闭环：图案会被冻结吗？
=====================================================================================
链条：
  · D66：守恒交换＋自洽势（吸引）⟹ 指数增长（无稳态）
  · D67：加压强 p=ρ（离散代理 κ(q_i−q_j)）⟹ 仍无静态解；但有**选择律**：渐近形状恒为 Fiedler（λ_min）
  · **Z022 §277**（台账已有的 L2 非线性）：**容量 K 是总量上限；触顶 ⟹ 清空活动层、写终端层、重播种**

本脚本问：把 Z022 的容量机制接上去后，是否得到**冻结的稳定结构**？它是什么？是否由选择律（Fiedler）定？

模型（守恒交换 + 压强 + 容量）：
    dx/dt = x − κ L x                         （D67 的线性部分；x = q − q̄，φ = −L⁺x）
    容量      ：max|x_i| > K_c ⟹ ① 写终端（累积归一化图案）② 清空活动层 ③ 重播种
    重播种    ：按 Z017 定理 7「保留最高两层」⟹ seed = 最近两个周期图案之和（归一化）× 小幅度

运行：cd ZeroGPU/src && /usr/bin/python3 d248_capacity_freeze_check.py
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


def fidler_ids(L, tol=1e-9):
    w, Q = np.linalg.eigh(L)
    nzm = w > tol
    lmin = w[nzm].min()
    return Q, w, np.where(nzm & (np.abs(w - lmin) < tol))[0], lmin


def run(adj, kappa, Kc, dt=0.02, max_steps=200000, n_cycles=24, amp0=0.05, seed_vec=None):
    """跑交换+压强+容量；返回每个周期末的图案（归一化）序列"""
    V = len(adj); L = lap(adj)
    rng = np.random.default_rng(2026)
    if seed_vec is None:
        x = rng.normal(size=V); x -= x.mean(); x /= np.linalg.norm(x); x *= amp0
    else:
        x = seed_vec.copy()
    terminal = np.zeros(V)
    pats = []
    last2 = []
    steps = 0
    while len(pats) < n_cycles and steps < max_steps:
        x = x + dt * (x - kappa * (L @ x))
        x -= x.mean()
        steps += 1
        if np.max(np.abs(x)) > Kc:
            pat = x / np.linalg.norm(x)
            pats.append(pat.copy())
            terminal = terminal + pat                      # ② 写终端层
            last2.append(pat.copy())
            if len(last2) > 2:
                last2.pop(0)
            seedn = np.sum(last2, axis=0)                  # ③ 重播种：保留最高两层
            seedn = seedn - seedn.mean()
            nn = np.linalg.norm(seedn)
            x = (seedn / nn) * amp0 if nn > 0 else np.zeros(V)   # ① 清空并重播种
    return np.array(pats), terminal / max(1e-30, np.linalg.norm(terminal)), steps


# ============================================================
print("=" * 78)
print("第 0 步：模型（D67 的交换 + Z022 §277 的容量）")
print("=" * 78)
print("  dx/dt = x − κ L x ；容量：max|x_i| > K_c ⟹ 写终端＋清空＋重播种（保留最高两层）")

# ============================================================
print()
print("=" * 78)
print("第 1 步：图案序列是否收敛（冻结）？")
print("=" * 78)
for gname, adj in [("cycle16", torus_adj(1, 16)), ("torus4x4", torus_adj(2, 4))]:
    L = lap(adj)
    Q, w, fid, lmin = fidler_ids(L)
    lmax = w[w > 1e-9].max()
    kappa = 0.5 * (1.0 / lmax + 1.0 / lmin)          # 窗口内
    Kc = 1.0
    pats, term, steps = run(adj, kappa, Kc)
    print(f"\n  ── {gname} ──  κ = {kappa:.4f}（窗口内），K_c = {Kc}，步数 = {steps}，周期数 = {len(pats)}")
    if len(pats) < 3:
        print("    周期太少")
        continue
    # 相邻周期图案差
    diffs = [float(np.linalg.norm(pats[i + 1] - pats[i])) for i in range(len(pats) - 1)]
    print(f"    相邻周期图案差（前 5）= {np.round(diffs[:5],4)}")
    print(f"    相邻周期图案差（末 5）= {np.round(diffs[-5:],6)}")
    ovs = [float(np.linalg.norm(Q[:, fid].T @ p)) for p in pats]
    print(f"    Fiedler 子空间投影（末 5）= {np.round(ovs[-5:],6)}")
    conv = diffs[-1] < 1e-6
    print(f"    ⟹ 图案序列{'收敛（冻结）' if conv else '未收敛'}")
    if gname == "cycle16":
        ck("H1a 容量机制使图案序列**收敛**（末两周期图案差 → 0）",
           conv, f"末次差 = {diffs[-1]:.3e}")
        ck("H1b 冻结图案 = Fiedler 模（选择律，D67）",
           ovs[-1] > 0.9999, f"末周期 Fiedler 投影 = {ovs[-1]:.6f}")

# ============================================================
print()
print("=" * 78)
print("第 2 步：选择律是否是**吸引子**（换初值给同一冻结图案）？")
print("=" * 78)
adj = torus_adj(1, 16); L = lap(adj)
Q, w, fid, lmin = fidler_ids(L)
lmax = w[w > 1e-9].max(); kappa = 0.5 * (1.0 / lmax + 1.0 / lmin)
frozen = []
for k in range(6):
    rng = np.random.default_rng(1000 + k)
    x0 = rng.normal(size=len(adj)); x0 -= x0.mean(); x0 /= np.linalg.norm(x0); x0 *= 0.05
    pats, term, steps = run(adj, kappa, 1.0, seed_vec=x0)
    frozen.append(pats[-1])
    ov = float(np.linalg.norm(Q[:, fid].T @ pats[-1]))
    print(f"  初值 #{k}: 周期数={len(pats):>3}  末图案 Fiedler 投影 = {ov:.6f}")
print(f"  末图案之间的最大差 = {max(float(np.linalg.norm(frozen[i]-frozen[0])) for i in range(len(frozen))):.3e}")
ck("H2a 冻结图案**不依赖初值**（选择律是吸引子）",
   all(float(np.linalg.norm(Q[:, fid].T @ f)) > 0.9999 for f in frozen),
   "6 个随机初值全部冻结到 Fiedler 子空间")

# ============================================================
print()
print("=" * 78)
print("第 3 步：对照 —— 不同 κ 区（选择律是否随 κ 变）")
print("=" * 78)
print(f"  {'κ/κ_窗口中点':>12} | 冻结图案(末周期) Fiedler 投影 | 周期数")
for frac, tag in [(0.2, "远小于窗口"), (0.5, "窗口中"), (0.9, "接近上界"), (1.5, "超上界(全衰减)")]:
    kappa = frac * (1.0 / lmin) if frac != 0.2 else 0.5 * (1.0 / lmax)
    pats, term, steps = run(adj, kappa, 1.0, max_steps=100000)
    if len(pats) == 0:
        print(f"  {tag:<12} κ={kappa:<8.4f} | 无周期（衰减到零，归零）")
        continue
    ov = float(np.linalg.norm(Q[:, fid].T @ pats[-1]))
    print(f"  {tag:<12} κ={kappa:<8.4f} | {ov:.6f} | {len(pats)}")

# ============================================================
print()
print("=" * 78)
print("第 4 步：判定")
print("=" * 78)
print("  ⟹ 容量机制（Z022 §277）把指数增长**截断成周期**；每个周期末的图案**归一化后收敛到同一模式**")
print("  ⟹ 冻结图案 = **Fiedler 模**（λ_min，最长波长）—— 即 D67 的选择律在非线性下**仍是吸引子**")
print("  ⟹ 与 GR 对照：这不是 TOV 型平衡剖面，而是**『谱隙选出的冻结图案』**")

n_fail = sum(1 for x in CK if not x["pass"])
print()
print("=" * 78)
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 78)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "d248_capacity_freeze_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "verdict": {
                   "capacity": "Z022 §277（触顶⟹清空+写终端+重播种）把线性指数行为截断成周期",
                   "frozen_pattern": "周期末图案收敛到 Fiedler 模（λ_min，最长波长）",
                   "attractor": "选择律是吸引子（6 个随机初值同结果）",
                   "contrast_GR": "不是 TOV 平衡剖面，而是谱隙选出的冻结图案"}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/d248_capacity_freeze_check.json")
sys.exit(1 if n_fail else 0)