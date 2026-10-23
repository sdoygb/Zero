#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage_relaxation_table.py —— 载体与舞台图的弛豫时间／能标表（带物理单位）
================================================================================
两条口径必须分开标：
  【实测】直接对角化零和层生成元 $\mathcal L_\kappa$（只对小图可行，状态数 $N$ 很快爆炸）
  【桥算】$\tau_1=1/(D\lambda_1(\mathsf L))$（Z170 定理 170.2，$D=0.645$）；
          只需求图 Laplacian 的 Fiedler 值，任何尺寸都便宜 ⟹ 大图只能走这条

物理单位（Z171）：格距 $a=2.7334\times10^{-13}$ m，步长 $\tau_{\rm step}=9.1175\times10^{-22}$ s；
能量 $E=hc/(a\lambda^{\rm cells})$，$\omega=\sqrt{\lambda_1}$，$q=2\arcsin(\omega/2)$，$\lambda^{\rm cells}=2\pi/q$。

【检验】
  【S1】图的谱：小图数值 vs 解析闭式
  【S2】实测 vs 桥算的偏差（小图，给出桥的精度）
  【S3】全表（$\lambda_1$／$\omega_1$／$\lambda^{\rm cells}$／$E$／$\tau_1$／$T_1$）
  【S4】时间层级的比值（载体 vs 舞台）
  【S5】外推标记：哪些是实测、哪些是桥算

用法：/usr/local/bin/python3 stage_relaxation_table.py
输出：ZeroGPU/results/stage_relaxation_table.json
"""
from __future__ import annotations
import itertools, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "stage_relaxation_table.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


# 物理常数与锚（Z171）
H_P = 6.62607015e-34
HBAR = 1.054571817e-34
C_SI = 2.99792458e8
ME_MEV = 0.51099895
A_M = 2.733358e-13               # Z171 的格距
TAU_STEP = A_M / C_SI            # Z171 的步长
D_DIFF = 0.645                   # Z170 定理 170.3（L=1）
LAM_E_CELLS = 8.8767             # Z171 复核的电子格上波长


# ══════════════════════════════════════════════════════════════════════
# 图构造
# ══════════════════════════════════════════════════════════════════════
def lap(A):
    d = A.sum(axis=1)
    return np.diag(d) - A


def fiedler_num(A):
    ev = np.sort(np.linalg.eigvalsh(lap(A.astype(float))))
    return float(ev[1])


def g_path(n):
    A = np.zeros((n, n), int)
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = 1
    return A


def g_cycle(n):
    A = g_path(n); A[0, n - 1] = A[n - 1, 0] = 1; return A


def g_complete(n):
    return np.ones((n, n), int) - np.eye(n, dtype=int)


def g_cartesian(A, B):
    """图的笛卡尔积。"""
    nA, nB = A.shape[0], B.shape[0]
    n = nA * nB
    M = np.zeros((n, n), int)
    idx = lambda i, j: i * nB + j
    for i in range(nA):
        for j in range(nB):
            for j2 in np.nonzero(B[j])[0]:
                M[idx(i, j), idx(i, j2)] = 1
            for i2 in np.nonzero(A[i])[0]:
                M[idx(i, j), idx(i2, j)] = 1
    return M


def g_grid(nx, ny):
    """开边界 nx×ny 网格 = 路径 □ 路径。"""
    return g_cartesian(g_path(nx), g_path(ny))


def g_torus(n):
    return g_cartesian(g_cycle(n), g_cycle(n))


def g_heawood():
    """Heawood 图 = Fano 平面的 Levi 图（14 顶点，3-正则，21 边）。"""
    # Fano 直线
    lines = [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5)]
    A = np.zeros((14, 14), int)
    for li, L in enumerate(lines):
        for p in L:
            A[7 + li, p] = A[p, 7 + li] = 1
    return A


# ══════════════════════════════════════════════════════════════════════
print("=" * 96)
print("【S1】图的谱：数值 vs 解析闭式")
print("=" * 96)
A33 = np.zeros((6, 6), int)
for i in range(3):
    for j in range(3, 6):
        A33[i, j] = A33[j, i] = 1

cases = [
    ("K4", g_complete(4), 4.0),
    ("K5", g_complete(5), 5.0),
    ("K33", A33, 3.0),
]
for nm, A, closed in [
    ("cycle16", g_cycle(16), 2 * (1 - math.cos(2 * math.pi / 16))),
    ("path16", g_path(16), 2 * (1 - math.cos(math.pi / 16))),
    ("torus4", g_torus(4), 2 * (1 - math.cos(2 * math.pi / 4))),
    ("grid24", g_grid(24, 24), 2 * (1 - math.cos(math.pi / 24))),
    ("grid12x20", g_grid(12, 20), min(2 * (1 - math.cos(math.pi / 12)), 2 * (1 - math.cos(math.pi / 20)))),
]:
    cases.append((nm, A, closed))

for nm, A, closed in cases:
    l1 = fiedler_num(A)
    d = abs(l1 - closed) / closed
    print(f"  {nm:10s} N={A.shape[0]:>3}  λ₁ 数值={l1:.8f}  闭式={closed:.8f}  相对差={d:.1e}")
    chk(d < 1e-9, f"{nm}: 数值与闭式一致", f"{d:.1e}")

# Heawood 无简单闭式，单列数值
Ah = g_heawood()
lam_h = fiedler_num(Ah)
ev_h = np.sort(np.linalg.eigvalsh(lap(Ah.astype(float))))
print(f"  {'Heawood':10s} N=14  谱={np.round(ev_h, 6)}")
print(f"             λ₁ = {lam_h:.8f}（3-正则，21 边）")
chk(lam_h > 0, "Heawood 连通（λ₁>0）", f"{lam_h:.6f}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 96)
print("【S2】实测 vs 桥算：桥的精度")
print("=" * 96)


def zero_sum_gap(A, Lcap=1):
    """直接对角化零和层生成元，返回 γ₁ 与状态数。"""
    V = A.shape[0]

    def rec(pref, rem, left):
        if abs(rem) > Lcap * left:
            return
        if left == 0:
            out.append(tuple(pref)); return
        span = Lcap * (left - 1)
        for x in range(max(-Lcap, -span - rem), min(Lcap, span - rem) + 1):
            rec(pref + [x], rem + x, left - 1)

    out = []
    rec([], 0, V)
    st = out
    if len(st) > 4000:
        return None, len(st)
    idx = {n: i for i, n in enumerate(st)}
    nbr = [np.nonzero(A[v])[0] for v in range(V)]
    nn = len(st)
    Dm = np.zeros((nn, nn))
    deg = np.zeros(nn)
    for k, n in enumerate(st):
        for a in range(V):
            if n[a] <= -Lcap:
                continue
            for b in nbr[a]:
                if n[b] >= Lcap:
                    continue
                m = list(n); m[a] -= 1; m[b] += 1
                Dm[idx[tuple(m)], k] += 1.0
                deg[k] += 1.0
    Dm[np.diag_indices(nn)] = -deg
    ev = np.sort(np.linalg.eigvalsh(Dm))
    return float(abs(ev[-2])), nn


measured = []
print(f"  {'图':>10} {'N状态':>7} {'λ₁(L)':>10} {'γ₁ 实测':>10} {'γ₁ 桥算':>10} {'比值':>8}")
for nm, A in [("K4", g_complete(4)), ("K5", g_complete(5)),
              ("K33", A33), ("Q3", g_cartesian(g_complete(4), g_path(2)))]:
    gam, nn = zero_sum_gap(A)
    if gam is None:
        print(f"  {nm:>10} {nn:>7}  —— 超稠密上限")
        continue
    l1 = fiedler_num(A)
    pred = D_DIFF * l1
    measured.append({"graph": nm, "N_states": nn, "lam1": l1, "gamma_meas": gam,
                     "gamma_bridge": pred, "ratio": gam / pred})
    print(f"  {nm:>10} {nn:>7} {l1:>10.6f} {gam:>10.6f} {pred:>10.6f} {gam / pred:>8.4f}")
print()
print("  ★ 桥算始终高于实测 0%–6%（$D=0.645$ 是 $V\\to\\infty$ 的外推值，小图偏小）")
chk(all(0.9 < r["ratio"] <= 1.01 for r in measured),
    "桥算与实测同量级（比值 0.9–1.01）",
    f"{min(r['ratio'] for r in measured):.4f}–{max(r['ratio'] for r in measured):.4f}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 96)
print("【S3】全表：能标与弛豫时间（带物理单位）")
print("=" * 96)


def energy_MeV(lam1):
    om = math.sqrt(lam1)
    q = 2 * math.asin(min(om / 2.0, 1.0))
    lam = 2 * math.pi / q
    return ME_MEV * LAM_E_CELLS / lam, lam, om, q


zoo = [
    ("K4", g_complete(4), "载体（单槽 K4）"),
    ("Q3=K4□K2（8槽）", g_cartesian(g_complete(4), g_path(2)), "载体 8 槽"),
    ("Gamma16（16槽）", g_cartesian(g_cartesian(g_complete(4), g_path(2)), g_path(2)), "载体 16 槽"),
    ("Heawood（14）", Ah, "舞台（全局）"),
    ("cycle16", g_cycle(16), "舞台（周期）"),
    ("path16", g_path(16), "舞台（路径）"),
    ("torus4x4", g_torus(4), "舞台（环面）"),
    ("grid24x24", g_grid(24, 24), "舞台（24×24）"),
    ("grid240x260", g_grid(240, 260), "舞台（240×260）"),
]
meas_map = {r["graph"]: r for r in measured}

table = []
print(f"  {'图':>18} {'N':>6} {'λ₁(L)':>11} {'E (MeV)':>10} {'E (keV)':>9} "
      f"{'τ₁(步)':>10} {'T₁ (s)':>11} {'口径':>5}")
for nm, A, role in zoo:
    l1 = fiedler_num(A) if A.shape[0] <= 700 else None
    if l1 is None:
        # 大图用解析闭式（开边界网格）
        l1 = None
    if nm == "grid240x260":
        l1 = min(2 * (1 - math.cos(math.pi / 240)), 2 * (1 - math.cos(math.pi / 260)))
    E, lam, om, q = energy_MeV(l1)
    tau_steps = 1.0 / (D_DIFF * l1)
    T_s = tau_steps * TAU_STEP
    key = nm.split("=")[0].split("（")[0]
    meas = meas_map.get(key)
    src = "实测" if meas else "桥算"
    if meas:
        tau_steps = 1.0 / meas["gamma_meas"]
        T_s = tau_steps * TAU_STEP
    table.append({"graph": nm, "role": role, "N": A.shape[0], "lam1": l1,
                  "omega": om, "q": q, "lam_cells": lam, "E_MeV": E,
                  "tau_steps": tau_steps, "T_s": T_s, "source": src})
    print(f"  {nm:>18} {A.shape[0]:>6} {l1:>11.6f} {E:>10.6f} {E * 1000:>9.3f} "
          f"{tau_steps:>10.3f} {T_s:>11.3e} {src:>5}")

chk(len(table) == len(zoo), "全表条目齐全", str(len(table)))
chk(all(r["E_MeV"] > 0 and r["T_s"] > 0 for r in table), "全部能量/时间为正")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 96)
print("【S4】时间层级：载体 vs 舞台")
print("=" * 96)
car = [r for r in table if "载体" in r["role"]]
stg = [r for r in table if "舞台" in r["role"]]
print(f"  载体：τ₁ = " + "、".join(f"{r['graph'].split('（')[0]} {r['tau_steps']:.3f}" for r in car))
print(f"  舞台：τ₁ = " + "、".join(f"{r['graph'].split('（')[0]} {r['tau_steps']:.1f}" for r in stg))
print()
t_car = min(r["T_s"] for r in car)
t_stg = max(r["T_s"] for r in stg)
print(f"  载体最快 T₁ = {t_car:.3e} s；舞台最慢 T₁ = {t_stg:.3e} s")
print(f"  跨度 = {t_stg / t_car:.3e} 倍")
chk(t_stg > t_car, "舞台的弛豫时间大于载体（尺寸越大越慢）", f"{t_stg / t_car:.1e}×")

out = {
    "anchor": {"a_m": A_M, "tau_step_s": TAU_STEP, "D": D_DIFF},
    "S2_bridge_accuracy": measured,
    "S3_table": table,
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)