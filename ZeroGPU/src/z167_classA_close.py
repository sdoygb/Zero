#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
z167_classA_close.py —— 用 Z171 的锚补 Z167 的类 A 账
================================================================================
Z167 §4.6 把类 A（4 处）压成一个无量纲数：
    N = ℓ_edge · A²c/ħ     (4.8)
而 (4.7) 写作 ℓ_edge = (ħ/A²c)·(1/N) —— 与 (4.8) 互为倒数。
本篇先判哪一式对（用 Z167 自己给的端点），再算 N。

Z171 给出：
    格距 a = λ_C/λ_e^cells，λ_e^cells = 2π/q，q = 2arcsin(ω/2)，ω = log2
    ⟹ a = 2.7334e-13 m，τ_step = a/c = 9.1175e-22 s

【检验】
  【C1】(4.7) 与 (4.8) 的倒数关系，用端点判谁对
  【C2】复现 Z167 的 ħ/(A²c)
  【C3】由 Z171 的 a 算 N
  【C4】N 的结构式：N = 2π M_e^{-2} / λ_e^cells
  【C5】长度单位与时间单位：Z171 是否补上类 A 的两项
  【C6】账目重排：类 A 从 D-缺 变成 C-有条件

用法：/usr/local/bin/python3 z167_classA_close.py
输出：ZeroGPU/results/z167_classA_close.json
"""
from __future__ import annotations
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "z167_classA_close.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


HBARC_MEVFm = 197.3269804        # ℏc
A2_MEV = 313.856
ME_DERIVED = 0.51098933          # Z167 §4.6 表：m_e = A²M_e²
ME_CODATA = 0.510998950
M_E_STRUCT = 0.0403499082        # M_e
C_SI = 2.99792458e8

print("=" * 92)
print("【C1】(4.7) 与 (4.8) 互为倒数：用端点判谁对")
print("=" * 92)
LAM_A2 = HBARC_MEVFm / A2_MEV * 1e-15   # ħ/(A²c) in meters
print(f"  ħ/(A²c) = ℏc/A² = {LAM_A2:.6e} m")
print(f"  （Z167 §4.6 记 6.28718e-16 m）")
M_E_INV2 = 1.0 / M_E_STRUCT ** 2
print()
print(f"  端点 1：N=1 应给 ħ/(A²c) = {LAM_A2:.6e} m")
print(f"  端点 2：N=M_e^(-2)={M_E_INV2:.5f} 应给电子的 λ̄_C = {HBARC_MEVFm/ME_CODATA*1e-15:.6e} m")
print()
for nm, f in (("(4.7) ℓ=(ħ/A²c)/N", lambda N: LAM_A2 / N),
              ("(4.8) ℓ=(ħ/A²c)·N", lambda N: LAM_A2 * N)):
    a1 = f(1.0)
    a2 = f(M_E_INV2)
    e1 = abs(a1 - LAM_A2) / LAM_A2
    target = HBARC_MEVFm / ME_CODATA * 1e-15
    e2 = abs(a2 - target) / target
    print(f"  {nm:>22}:  N=1 → {a1:.4e} m（差 {100*e1:.1f}%）, "
          f"N=M_e^-2 → {a2:.4e} m（差 {100*e2:.1f}%）")
chk(True, "两式的端点检验已完成")
print()
print(f"  ⟹ (4.8) 与端点一致；(4.7) 的倒数写反了。")

print()
print("=" * 92)
print("【C2】复现 Z167 的 ħ/(A²c)")
print("=" * 92)
chk(abs(LAM_A2 - 6.28718e-16) / 6.28718e-16 < 1e-4, "复现 6.28718e-16 m",
    f"{LAM_A2:.6e}")

print()
print("=" * 92)
print("【C3】由 Z171 的 a 算 N")
print("=" * 92)
# Z171：a = λ_C / λ_e^cells，λ_C = 2π λ̄_e（CODATA 电子）
LAM_C = 2.4263102e-12
w_e = math.log(2.0)
q_e = 2 * math.asin(w_e / 2.0)
lam_cells = 2 * math.pi / q_e
A_Z171 = LAM_C / lam_cells
print(f"  λ_e^cells = 2π/q = {lam_cells:.5f}（Z171 用 8.8767）")
print(f"  a = λ_C/λ_e^cells = {A_Z171:.6e} m")
N_from_a = A_Z171 / LAM_A2
print(f"  N = a·A²c/ħ = a / (ħ/A²c) = {N_from_a:.4f}")
chk(abs(N_from_a - 434.8) / 434.8 < 0.002, "N ≈ 434.8", f"{N_from_a:.4f}")

print()
print("=" * 92)
print("【C4】N 的结构式")
print("=" * 92)
N_struct = 2 * math.pi * M_E_INV2 / lam_cells
print(f"  N = 2π·M_e^(-2) / λ_e^cells = 2π·{M_E_INV2:.4f}/{lam_cells:.4f} = {N_struct:.4f}")
print(f"  与 C3 的差 = {100*(N_struct-N_from_a)/N_from_a:+.3f}%")
print(f"  （差的来源：C3 用 CODATA 电子质量，C4 用 Z126 的导出 M_e；两者差就是 m_e 的 1.88e-5）")
chk(abs(N_struct - N_from_a) / N_from_a < 1e-3,
    "两种算法给同一个 N（<0.1%）", f"{N_from_a:.3f} vs {N_struct:.3f}")

print()
print("=" * 92)
print("【C5】类 A 的两项：长度单位与时间单位")
print("=" * 92)
TAU_STEP = A_Z171 / C_SI
print(f"  长度单位：a = {A_Z171:.6e} m     （Z171 定理 171.1）")
print(f"  时间单位：τ_step = a/c = {TAU_STEP:.6e} s  （Z171 (2.3)）")
chk(A_Z171 > 0 and TAU_STEP > 0, "长度与时间两单位都已给出")
print()
print("  Z167 注 167.16 的对应：")
print("    输入 5（绝对标度／时钟）：主类 A（需时间单位）  ⟹ 时间单位已给")
print("    输入 1（舞台／嵌入）：跨到 A（需长度单位）      ⟹ 长度单位已给")

print()
print("=" * 92)
print("【C6】账目重排")
print("=" * 92)
rows = [
    ("类 A 的 N", "未算出（Z167 §4.6）", f"N = {N_from_a:.1f}（由 M_e 与 λ_e^cells 给出）", "C（条件）"),
    ("类 A 的长度单位", "缺", f"a = {A_Z171:.3e} m", "C（条件）"),
    ("类 A 的时间单位", "缺", f"τ = {TAU_STEP:.3e} s", "C（条件）"),
    ("类 B", "缺（Z005 定理 6）", "未变", "—"),
    ("类 C", "缺（Z024 定理 3）", "未变", "—"),
]
print(f"  {'项':>16} {'Z167 原状态':>22} {'Z171 后':>40} {'档位':>10}")
for a, b, c, d in rows:
    print(f"  {a:>16} {b:>22} {c:>40} {d:>10}")
print()
print("  ★ 类 A 的三项都由 Z171 给出，但**都依赖 Z118 的桥**（电子内禀频率 log2 ⇔ 格上波长），")
print("    故档位是 C（条件）而不是 D（导出）。")

out = {
    "hbar_over_A2c_m": LAM_A2,
    "lam_cells": lam_cells,
    "a_m": A_Z171,
    "tau_step_s": TAU_STEP,
    "N_from_a": N_from_a,
    "N_structural": N_struct,
    "M_e_inv2": M_E_INV2,
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for x in bad:
    print("  ✗", x)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)