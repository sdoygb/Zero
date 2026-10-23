#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage_ratio_decomposition.py —— 10^39 的分解：舞台图能解释多少
================================================================================
问题：a_cosmo/a_micro = 1.07e39 是"舞台图未定"造成的吗？攻舞台图能把它变成预言吗？

分解（无新输入，全部取自既有篇）：
  a = L_obs / N      （a 是格距，L_obs 是该支的观测长度，N 是该支的格数）
  ⟹ a_cosmo/a_micro = (L_cosmo/L_micro) × (N_micro/N_cosmo)

两处 N：
  N_cosmo = ⌈|C|^{1/3}⌉ ，|C|=14（Heawood 顶点数）        ← Z023 L121「|C|=14 的立方排布」
  N_micro = 2π/arcsin(ω/2)，ω=log2                        ← Z118 §2／Z171 (2.1)

【检验】
  【D1】N_cosmo 由 14 的立方排布给 3
  【D2】N_micro 复现 8.8767
  【D3】两个 N 的比 ≈ 2.96（量级、非拟合）
  【D4】观测长度比 ≈ 3.6e38
  【D5】乘积复现 1.07e39
  【D6】"攻舞台图"的收益上界 = log10(2.96) ≈ 0.47 dex，而总跨度 38.6 dex

用法：/usr/local/bin/python3 stage_ratio_decomposition.py
输出：ZeroGPU/results/stage_ratio_decomposition.json
"""
from __future__ import annotations
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "stage_ratio_decomposition.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


C_SI = 2.99792458e8

print("=" * 92)
print("【D1】宇宙支的格数：|C|=14 的立方排布")
print("=" * 92)
V_HEA = 14
n_cosmo = math.ceil(V_HEA ** (1.0 / 3.0))
print(f"  |C| = {V_HEA}（Heawood 顶点数，Z021 §8 申报的最小性）")
print(f"  14^(1/3) = {V_HEA ** (1/3):.4f}  ⟹  ⌈·⌉ = {n_cosmo}")
print(f"  即 3×3×3 = {3**3} 的立方里放 14 个站 ⟹ 直径 3 格（Z023 L121）")
chk(n_cosmo == 3, "N_cosmo 由 14 的立方排布给 3", f"{n_cosmo}")

print()
print("=" * 92)
print("【D2】微观支的格数：电子的格上波长")
print("=" * 92)
w_e = math.log(2.0)
q_e = 2 * math.asin(w_e / 2.0)
n_micro = 2 * math.pi / q_e
print(f"  ω = log2 = {w_e:.6f}；q = 2arcsin(ω/2) = {q_e:.6f}")
print(f"  N_micro = 2π/q = {n_micro:.4f} 格（Z118 实测 8.874）")
chk(abs(n_micro - 8.8767) < 0.01, "复现 Z171 的 8.8767", f"{n_micro:.4f}")

print()
print("=" * 92)
print("【D3】格数比")
print("=" * 92)
ratio_N = n_micro / n_cosmo
print(f"  N_micro/N_cosmo = {n_micro:.4f}/{n_cosmo} = {ratio_N:.4f}")
chk(2.5 < ratio_N < 3.5, "格数比在 3 倍量级", f"{ratio_N:.4f}")

print()
print("=" * 92)
print("【D4】观测长度比")
print("=" * 92)
LY = 9.4607304725808e15
L_cosmo = 9.3e10 * LY                       # 可观测宇宙直径
ME = 9.1093837015e-31
HBAR = 1.054571817e-34
L_micro = 2 * math.pi * HBAR / (ME * C_SI)  # 电子全康普顿波长
ratio_L = L_cosmo / L_micro
print(f"  L_cosmo = 9.3e10 光年 = {L_cosmo:.4e} m")
print(f"  L_micro = 电子全康普顿波长 = {L_micro:.4e} m")
print(f"  比 = {ratio_L:.4e}  (= 10^{math.log10(ratio_L):.2f})")
chk(3.0e38 < ratio_L < 4.0e38, "观测长度比 ≈ 3.6e38", f"{ratio_L:.3e}")

print()
print("=" * 92)
print("【D5】乘积复现 10^39")
print("=" * 92)
total = ratio_L * ratio_N
print(f"  观测长度比 × 格数比 = {total:.4e}")
print(f"  直接由 a 算         = {(L_cosmo/n_cosmo)/(L_micro/n_micro):.4e}")
chk(abs(total / ((L_cosmo / n_cosmo) / (L_micro / n_micro)) - 1) < 1e-9,
    "两条路算出同一个比", f"{total:.4e}")

print()
print("=" * 92)
print("【D6】攻舞台图的收益上界")
print("=" * 92)
payoff = math.log10(ratio_N)
span = math.log10(total)
print(f"  总跨度          = 10^{span:.2f}")
print(f"  格数（舞台）贡献 = 10^{payoff:.2f} = {ratio_N:.2f} 倍")
print(f"  观测长度贡献     = 10^{math.log10(ratio_L):.2f}")
print(f"  ⟹ 攻舞台图的收益上界 = {payoff:.2f} dex，占总跨度 {100*payoff/span:.1f}%")
chk(payoff / span < 0.02, "舞台图贡献 < 总跨度的 2%", f"{100*payoff/span:.2f}%")

print()
print("  分解式：")
print(f"    a_cosmo/a_micro = (10^{math.log10(ratio_L):.2f}) × ({ratio_N:.2f})")
print(f"                      = 观测 × 舞台")
print()
print("  ★ 结论：舞台图（格数）这一侧**已由 |C|=14 的立方排布定住**；")
print("     剩下的 10^%.2f 是两个观测量之比，与舞台图无关。" % math.log10(ratio_L))

out = {
    "N_cosmo": n_cosmo, "N_micro": n_micro, "ratio_N": ratio_N,
    "L_cosmo_m": L_cosmo, "L_micro_m": L_micro, "ratio_L": ratio_L,
    "total_ratio": total,
    "payoff_dex": payoff, "span_dex": span,
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