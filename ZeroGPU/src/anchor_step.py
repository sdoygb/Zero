#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
anchor_step.py —— 锚定唯一锚：把"一步"换成物理时间，算出绝对时间/能量
================================================================================
链条（全部取自库内既有结果，不新增输入）：
  ① 唯一量纲锚 = $v$（Z001／Z013 定理 16／Z153）。$m_e$ 由 $v$ **算出**：
     $A^2=v/784.5=313.856$ MeV，$M_e=1+\sqrt2\cos(2\pi/3+2/9)$，$m_e=A^2M_e^2$。
  ② 电子在格上的波长 $\lambda_e^{\rm cells}=2\pi/q$，$q=2\arcsin(\omega/2)$，$\omega=\log2$，
     实测 $8.874$ 格（Z118 §2；Z152 定理 152.1 称该式精确）。
  ③ 物理上 $m_ec^2 = hc/\lambda_{\rm phys}$，而 $\lambda_{\rm phys}=a\,\lambda^{\rm cells}$
     ⟹ 格距 $a=\lambda_C/\lambda_e^{\rm cells}$，步长 $\tau_{\rm step}=a/c$。

【纪律备忘】print 里不要嵌 ASCII 双引号（本脚本作者在此摔过多次）：中文引用一律用「」。

【检验】
  【A1】$m_e$ 由 $v$ 算出（复核相律，确认锚链没断）
  【A2】锚定 $a$ 与 $\tau_{\rm step}$
  【A3】两条路线互检：波长路线 vs $\hbar\omega$ 路线；差应为色散效应（Z118 的 ratio 0.979）
  【A4】载体与舞台图的绝对能量表
  【A5】绝对时间表（$\tau_1$（步）× $\tau_{\rm step}$）
  【A6】两个"能量"的区分：模能 $\hbar\omega$ vs 线宽 $\hbar\gamma$，各自标度

用法：/usr/local/bin/python3 anchor_step.py
输出：ZeroGPU/results/anchor_step.json
"""
from __future__ import annotations
import json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "anchor_step.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


# ══════════════════════════════════════════════════════════════════════
# 物理常数（CODATA）
# ══════════════════════════════════════════════════════════════════════
h_P = 6.62607015e-34
hbar = 1.054571817e-34
c_SI = 2.99792458e8
eV = 1.602176634e-19
me_kg = 9.1093837015e-31

print("=" * 92)
print("【A1】$m_e$ 由 $v$ 算出（复核锚链）")
print("=" * 92)
v_GeV = 246.2205
A2_MeV = v_GeV * 1000 / 784.5
M_e = 1 + math.sqrt(2) * math.cos(2 * math.pi / 3 + 2 / 9)
me_MeV = A2_MeV * M_e ** 2
print(f"  A² = v/784.5 = {A2_MeV:.4f} MeV")
print(f"  M_e = 1+√2cos(2π/3+2/9) = {M_e:.10f}")
print(f"  m_e = A²·M_e² = {me_MeV:.6f} MeV   （CODATA {me_kg * c_SI ** 2 / eV / 1e6:.6f} MeV）")
chk(abs(me_MeV - 0.51099895) / 0.51099895 < 3e-3, "m_e 由 v 算出（相对差 < 0.3%）",
    f"{me_MeV:.6f} vs 0.510999 (差 {100 * (me_MeV - 0.51099895) / 0.51099895:+.4f}%)")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A2】锚定：格距 a 与步长 τ_step")
print("=" * 92)
omega_e = math.log(2.0)
q_e = 2 * math.asin(omega_e / 2.0)
lam_e_cells = 2 * math.pi / q_e
print(f"  电子：ω = log2 = {omega_e:.6f}，q = 2arcsin(ω/2) = {q_e:.6f}")
print(f"        λ_e = 2π/q = {lam_e_cells:.4f} 格   （Z118 实测 8.874）")
chk(abs(lam_e_cells - 8.874) / 8.874 < 1e-3, "复现 Z118 的 8.874 格",
    f"{lam_e_cells:.4f}")

lambdabar_e = hbar / (me_kg * c_SI)          # 约化康普顿波长
lambda_C = 2 * math.pi * lambdabar_e          # 全长
a_m = lambda_C / lam_e_cells
tau_s = a_m / c_SI
print()
print(f"  约化康普顿 λ̄_e = {lambdabar_e:.6e} m")
print(f"  全波长   λ_C   = {lambda_C:.6e} m")
print(f"  ⟹ 格距 a = λ_C/λ_e = {a_m:.6e} m")
print(f"  ⟹ 步长 τ_step = a/c = {tau_s:.6e} s")
out_anchor = {"a_m": a_m, "tau_step_s": tau_s, "lam_e_cells": lam_e_cells}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A3】两条路线互检：波长路线 vs ħω 路线")
print("=" * 92)
tau_hbar = hbar * omega_e / (me_kg * c_SI ** 2)
print(f"  波长路线（精确）：τ = a/c = {tau_s:.6e} s")
print(f"  ħω 路线（连续极限）：τ = ħω/(m_ec²) = {tau_hbar:.6e} s")
ratio = tau_s / tau_hbar
print(f"  比值 = {ratio:.5f}   （Z118 记 λ_e/(2π/f) = 0.979 = 1/{1 / 0.979:.5f}）")
chk(abs(ratio - 1 / 0.979) < 0.01, "两路线之差 = 已登记的色散效应（≈2.2%）",
    f"{ratio:.5f} vs {1 / 0.979:.5f}")
print()
print("  ★ 采用波长路线（Z152 定理 152.1 称 λ=2π/q 精确）；ħω 路线是 q→0 极限。")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A4】绝对能量表（E = hc/(a·λ_cells)，λ_cells = 2π/q，ω=√λ₁，c=1）")
print("=" * 92)
D = 0.645
me_MeV_codata = 0.51099895


def energy_MeV(lam1):
    """由 Fiedler 值给绝对能量（MeV）。"""
    om = math.sqrt(lam1)
    q = 2 * math.asin(min(om / 2.0, 1.0))
    lam = 2 * math.pi / q
    return me_MeV_codata * lam_e_cells / lam, lam, q, om


def lam1_path(V):
    return 2 * (1 - math.cos(math.pi / V))


def lap_eigs_path(V):
    return np.sort(np.linalg.eigvalsh(2 * np.eye(V) - np.eye(V, k=1) - np.eye(V, k=-1)))


rows = []
print(f"  {'图':>14} {'λ₁(L)':>11} {'ω₁':>9} {'q₁':>9} {'λ_cells':>9} {'E (MeV)':>10}")
for nm, lam1 in [("Q3 (8槽)", 2.0), ("Γ16 (16槽)", 2.0),
                 ("path8", lam1_path(8)), ("path16", lam1_path(16)),
                 ("path32", lam1_path(32)), ("path128", lam1_path(128))]:
    E, lam, q, om = energy_MeV(lam1)
    rows.append({"graph": nm, "lam1": lam1, "omega": om, "q": q,
                 "lam_cells": lam, "E_MeV": E})
    print(f"  {nm:>14} {lam1:>11.6f} {om:>9.6f} {q:>9.6f} {lam:>9.4f} {E:>10.6f}")
chk(all(r["E_MeV"] > 0 for r in rows), "全部能量为正")
chk(rows[0]["E_MeV"] == rows[1]["E_MeV"], "8 槽与 16 槽能量相同（λ₁ 冻结在 2）", )

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A5】绝对时间表：τ₁（步）× τ_step")
print("=" * 92)
print(f"  {'图':>14} {'λ₁':>11} {'τ₁(步)':>10} {'物理寿命 (s)':>16}")
rows_t = []
for r in rows:
    tau_steps = 1.0 / (D * r["lam1"])
    T_s = tau_steps * tau_s
    rows_t.append({**r, "tau_steps": tau_steps, "T_s": T_s})
    print(f"  {r['graph']:>14} {r['lam1']:>11.6f} {tau_steps:>10.4f} {T_s:>16.6e}")
chk(all(r["T_s"] > 0 for r in rows_t), "全部寿命为正")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A6】两个能量必须分开：模能 ħω vs 线宽 ħγ")
print("=" * 92)
print(f"  {'图':>14} {'ħω (MeV)':>11} {'ħγ (MeV)':>11} {'比值 γ/ω':>10}")
for r in rows:
    E_mode = r["E_MeV"]                                     # = ħω（由波长路线标定）
    gam = D * r["lam1"]
    # 同一条标定：ħγ/(ħω) = γ/ω
    E_width = E_mode * gam / r["omega"]
    print(f"  {r['graph']:>14} {E_mode:>11.6f} {E_width:>11.6f} {gam / r['omega']:>10.4f}")
print()
print("  标度：ħω ∝ √λ₁ ∝ V⁻¹；  ħγ ∝ λ₁ ∝ V⁻²")
print("  ⟹ 线宽随体系**快于**模能下降（比值为 γ/ω = D·ω，随 V 减半每次翻倍）")

out = {
    "anchor": out_anchor,
    "me_MeV_from_v": me_MeV,
    "tau_two_routes": {"wavelength_s": tau_s, "hbar_omega_s": tau_hbar,
                       "ratio": ratio},
    "energy_table": rows,
    "time_table": rows_t,
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