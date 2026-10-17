#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
d248_pressure_jeans_check.py —— 把 p=ρ（旧库 D250）接进交换闭环：是否给出稳定非平凡稳态？
=========================================================================================
背景：
  · D66 已证：守恒交换＋自洽势（吸引，Z033 约定 Lφ=−x）⟹ 偏差指数增长（无稳态）
  · 旧库 **D250**：最小标量类给出 **p = ρ**（刚性物态）；GR 里静态自引力解需**压强支撑**（TOV）
  · 本脚本把"压强"实现为**密度梯度驱动的交换流** $J^p_{i→j} = κ(q_i − q_j)$（从高密度流向低密度），
    与引力流 $J^g_{i→j} = K(φ_i−φ_j)$ 叠加（两者都是**守恒交换**，沿用 D66 路 (ii) 的形式）：

        J_{i→j} = K(φ_i − φ_j) + κ(q_i − q_j),   dq_i/dt = −Σ_j J_{i→j}
        ⟹ 线性化（x = q − q̄，φ = −L⁺x）：dx/dt = x − κ L x
        ⟹ 逐模生长率  γ_k = 1 − κ λ_k   （λ_k 是 L 的本征值）

预测（线性的、可证伪的）：
    小 κ  ⟹ 全部非零模增长（最不稳定）
    窗口  1/λ_max < κ < 1/λ_min ⟹ **选择性图案**：只有 λ_k < 1/κ 的模增长，主导模＝最小非零 λ（Fiedler）
    大 κ  ⟹ 全部模衰减 ⟹ 归零（均匀，平凡）

运行：cd ZeroGPU/src && /usr/bin/python3 d248_pressure_jeans_check.py
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
    phi = -np.linalg.pinv(L) @ x
    return phi - phi.mean()


# ============================================================
print("=" * 78)
print("第 0 步：模型的物理与线性分析")
print("=" * 78)
print("  引力流  J^g_{i→j} = K(φ_i − φ_j)      （K = A；φ = −L⁺x，势井＝吸引）")
print("  压强流  J^p_{i→j} = κ(q_i − q_j)      （从高密度流向低密度；刚性物态 p=ρ 的最简离散代理）")
print("  合计    dq_i/dt = −Σ_j J_{i→j}       ⟹ 线性化：dx/dt = x − κ L x")
print("  ⟹ 逐模生长率  γ_k = 1 − κ λ_k       （x 在 L 的本征基下对角）")

# ============================================================
print()
print("=" * 78)
print("第 1 步：谱与三区边界")
print("=" * 78)
spec = {}
for gname, adj in [("cycle16", torus_adj(1, 16)), ("torus4x4", torus_adj(2, 4)),
                   ("path8", path_adj(8))]:
    L = lap(adj)
    w = np.linalg.eigvalsh(L)
    nz = w[w > 1e-9]
    lmin, lmax = float(nz.min()), float(nz.max())
    spec[gname] = {"lam_min": lmin, "lam_max": lmax,
                   "k_lo": 1.0 / lmax, "k_hi": 1.0 / lmin}
    print(f"  {gname:<10}: λ_min = {lmin:.6f}  λ_max = {lmax:.6f}  ⟹ "
          f"窗口 κ ∈ ({1.0/lmax:.6f}, {1.0/lmin:.6f})，宽度 {1.0/lmin - 1.0/lmax:.6f}")
note("窗口的物理含义", "κ 大于 1/λ_min ⟹ 全部模衰减 ⟹ 归零；小于 1/λ_max ⟹ 全部模增长 ⟹ 塌缩",
     "窗口内 ⟹ 只长波模增长 ⟹ **选择性图案**（Turing 型）")

# ============================================================
print()
print("=" * 78)
print("第 2 步：逐模生长率 —— 解析 vs 实跑（完整矩阵更新，再投影到本征基）")
print("=" * 78)
adj = torus_adj(1, 16)
L = lap(adj)
w, V = np.linalg.eigh(L)
nz_mask = w > 1e-9
nz_idx = np.where(nz_mask)[0]
lam_nz = w[nz_mask]
dt = 0.02
print(f"  cycle16：dt = {dt}；非零模数 = {len(lam_nz)}")

rng = np.random.default_rng(101)
x0 = rng.normal(size=len(adj)); x0 -= x0.mean()
x0 /= np.linalg.norm(x0)

ok_all = True
for kappa in (0.0, 0.05, 0.3, 1.0, 3.0):
    # 完整矩阵更新（不分解）跑 200 步
    x = x0.copy()
    for _ in range(200):
        x = x + dt * (x - kappa * (L @ x))          # dx/dt = x − κ L x
    # 投影到非零模
    c = V[:, nz_idx].T @ x
    amp = np.abs(c)
    # 主导模
    dom = int(np.argmax(amp))
    lam_dom = float(lam_nz[dom])
    # 解析：γ_k = 1 − κλ_k，逐步放大率 = 1 + dt·γ_k
    gam_pred = 1.0 - kappa * lam_nz
    # 实测：末 50 步的逐模放大率（用解析式的对数残差代替逐步追踪，直接比"是否增长"）
    grow_pred = gam_pred > 0
    print(f"\n  κ = {kappa:<5}: 主导模 λ = {lam_dom:.6f}（解析最不稳定模 λ_min = {lam_nz.min():.6f}）；"
          f"增长模数（解析）= {int(grow_pred.sum())}/{len(lam_nz)}")
    # 检出：末态模幅是否与解析一致（用解析放大 200 步的结果对比）
    amp_theory = np.abs(V[:, nz_idx].T @ x0) * np.abs(1 + dt * gam_pred) ** 200
    rel = float(np.abs(amp / max(1e-300, amp.max()) - amp_theory / max(1e-300, amp_theory.max())).max())
    print(f"            末态模幅（归一化后）与解析预言的最大偏差 = {rel:.3e}")
    ok_all &= rel < 1e-6
    if abs(kappa - 0.3) < 1e-9:
        ck("G2a 逐模生长率 γ_k = 1 − κλ_k 被实跑精确确证（完整矩阵更新）",
           rel < 1e-6, f"归一化模幅偏差 = {rel:.3e}")
        ck("G2b 窗口内主导模 = 最小非零 λ（Fiedler）",
           abs(lam_dom - lam_nz.min()) < 1e-9,
           f"主导 λ = {lam_dom:.6f}，λ_min = {lam_nz.min():.6f}")

# ============================================================
print()
print("=" * 78)
print("第 3 步：三区行为（振幅趋势 ＋ 渐近形状）")
print("=" * 78)
print("  注意：**渐近形状在一切 κ 下都是最慢模（λ_min 的模）** —— 因为它是'最不稳定'（增长最快）")
print("        或'最慢衰减'（衰减最慢）的那一支。κ 只决定**振幅是长还是消**、以及**参与模数**。")
print("  判据：① 振幅对数斜率（>0 增长 / <0 衰减）② 末态在 **Fiedler 本征子空间**（λ=λ_min，可能简并）上的投影范数")
for gname, adj in [("cycle16", torus_adj(1, 16)), ("torus4x4", torus_adj(2, 4)), ("path8", path_adj(8))]:
    L = lap(adj); w, V = np.linalg.eigh(L)
    nzm = w > 1e-9
    lmin, lmax = float(w[nzm].min()), float(w[nzm].max())
    k_lo, k_hi = 1.0 / lmax, 1.0 / lmin
    fid = np.where(nzm & (np.abs(w - lmin) < 1e-9))[0]     # Fiedler 子空间（含简并）
    print(f"\n  ── {gname} ──  窗口 κ ∈ ({k_lo:.4f}, {k_hi:.4f})；"
          f"λ_min = {lmin:.6f}（简并度 {len(fid)}）")
    for tag, kappa in [("小 κ（全增长）", 0.5 * k_lo),
                       ("窗口内", 0.5 * (k_lo + k_hi)),
                       ("大 κ（全衰减）", 2.0 * k_hi)]:
        rng = np.random.default_rng(103)
        x = rng.normal(size=len(adj)); x -= x.mean(); x /= np.linalg.norm(x)
        dtx = 0.02
        # ① 振幅趋势（前 60 步，不归一化）
        norms = []
        xx = x.copy()
        for _ in range(60):
            xx = xx + dtx * (xx - kappa * (L @ xx))
            norms.append(float(np.linalg.norm(xx)))
        slope = float(np.polyfit(np.arange(len(norms)), np.log(np.maximum(norms, 1e-300)), 1)[0])
        # ② 渐近形状（每步归一化 ＋ 强制零和，跑长）
        #    注：dx/dt = x − κLx 精确保持零和（L 杀常数），但浮点舍入会在常数方向漏分量，
        #    而常数方向的因子是 1+dt（≈增长）⟹ 必须每步强制零和，否则大 κ 下会被舍入噪声淹没
        xa = x.copy()
        for _ in range(3000):
            xa = xa + dtx * (xa - kappa * (L @ xa))
            xa -= xa.mean()
            n_ = np.linalg.norm(xa)
            if n_ < 1e-300:
                break
            xa /= n_
        ov = float(np.linalg.norm(V[:, fid].T @ (xa / max(1e-300, np.linalg.norm(xa)))))
        trend = "增长" if slope > 0 else "衰减"
        print(f"    {tag:<16} κ={kappa:<7.4f}: log‖x‖ 斜率 = {slope:+.6f}（{trend}）；"
              f"Fiedler 子空间投影 = {ov:.6f}")
        if gname == "cycle16" and tag == "窗口内":
            ck("G3a 窗口内 ⟹ 振幅增长 且 形状收敛到 Fiedler 子空间（图案放大、非归零）",
               slope > 0 and ov > 0.9999, f"斜率 {slope:+.6f}；投影 {ov:.6f}")
        if gname == "cycle16" and tag == "大 κ（全衰减）":
            ck("G3b 大 κ ⟹ 振幅衰减（归零到均匀；形状仍沿最慢模 Fiedler）",
               slope < 0 and ov > 0.9999, f"斜率 {slope:+.6f}；投影 {ov:.6f}")

# ============================================================
print()
print("=" * 78)
print("第 4 步：判定")
print("=" * 78)
print("  ⟹ 加上压强项 p=ρ（离散代理：κ(q_i−q_j)）后：")
print("     小 κ   ：全部非零模增长（塌缩型）")
print("     窗口内 ：**选择性图案** —— 收敛到最小非零 λ 的模（Fiedler），波长由**谱隙**选择")
print("     大 κ   ：全部模衰减 ⟹ 归零（均匀）")
print("  ⟹ **仍没有「稳定非平凡静态解」**；有的是**图案选择**（Turing 型），选择律由图的谱隙给：")
print("       最不稳定模 ⟺ λ_k 最小 ⟺ 最长波长（cut-off 是**图本身的有限尺寸/谱隙**，不是 Jeans 长度）")
print()
print("  与 GR 对照：GR 的静态自引力解 = TOV（压强梯度 vs 引力，边界固定）；")
print("              本模型里对应物 = **谱隙截断的图案**，而非 TOV 型平衡剖面")

n_fail = sum(1 for x in CK if not x["pass"])
print()
print("=" * 78)
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 78)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "d248_pressure_jeans_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT, "spectra": spec,
               "linear_law": "γ_k = 1 − κ λ_k（dx/dt = x − κ L x）",
               "verdict": {
                   "no_static_solution": "仍无稳定非平凡静态解",
                   "selected_pattern": "窗口内收敛到 Fiedler 模（最小非零 λ）；波长由谱隙选择",
                   "contrast_GR": "GR 的静态解＝TOV 剖面；此处对应物＝谱隙截断的图案"}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/d248_pressure_jeans_check.json")
sys.exit(1 if n_fail else 0)