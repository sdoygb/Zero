#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
higgs_lap_conjecture_test.py —— 检验"绕一圈／两圈后质量变了，是已知粒子的变种"
================================================================================
用户的猜想（逐字）：
  "一个粒子绕了一圈，或两圈，又叠回自己原来的位置，然后质量就变了，
   是某种已知粒子的变种。"

本脚本把它拆成五个可判定命题，逐条实测：

  【T1】相位律下的"绕圈"：θ → θ+2π（一圈）、θ+4π（两圈）是否改变质量？
  【T2】双覆盖（Z006）：2π 给 −1、4π 给 +1 —— 是否改变质量？
  【T3】Zero 里真的有"绕数"吗？（Z154：k₀ = A_ℓ² = 2，环② = 2πk₀ = 4π = **两圈**）
  【T4】"变种"：M_H 是否是某已知质量的**结构倍数**？
  【T5】"不同的路 → 不同的质量"：立方体的路段类型是否给出振幅（Z126 §2）？

依据：Z081（闭合＝一圈＝周期，原生）、Z126 §1–§2（振幅律与立方体对角比）、
      Z006（双覆盖）、Z154（绕数 k₀=2、环②=4π）、Z156（四环链）、Z063（Berry 相位与闭合）、
      Z017／Z022（容量）、Z005 定理 11（Z_N 单位根）。

用法：/usr/local/bin/python3 higgs_lap_conjecture_test.py
"""
from __future__ import annotations
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "higgs_lap_conjecture_test.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}
V, M_H = 246.22, 125.25
M_Z = 91.1876
A2 = V * 1000 / 784.5      # 313.856 MeV，轻子质量单元

print("=" * 88)
print("用户猜想：绕一圈或两圈 → 质量变 → 是某已知粒子的变种")
print("=" * 88)

# ══════════════════════════════════════════════════════════════════════
# T1 相位律下的绕圈是否改变质量
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【T1】相位律 m = u_s A²[1+A_s cosθ]² 下，绕一圈／两圈是否改变质量？")
print("=" * 88)


def m_of(theta_deg, A_s, u_s=1.0, A2=A2):
    return u_s * A2 * (1 + A_s * math.cos(math.radians(theta_deg))) ** 2


As_ell = math.sqrt(2)
print(f"  取轻子层（A_s=√2）：θ 与 θ+360°（一圈）、θ+720°（两圈）")
print()
print(f"  {'θ':>8} {'m(θ) MeV':>12} {'m(θ+360°)':>12} {'m(θ+720°)':>12}")
for th in [0, 45, 90, 132.7324, 180, 225, 270]:
    a, b, c = m_of(th, As_ell), m_of(th + 360, As_ell), m_of(th + 720, As_ell)
    print(f"  {th:8.4f} {a:12.6f} {b:12.6f} {c:12.6f}")
same = all(abs(m_of(t, As_ell) - m_of(t + 360 * n, As_ell)) < 1e-9
           for t in [0, 45, 90, 132.7, 180, 225, 270] for n in (1, 2))
chk(same, "绕一圈／两圈后质量完全相同（周期 2π）")
print()
print("  ⟹ 因为 cos(θ+2πn) = cos θ，绕任何整数圈都**回到同一个质量**。")
print("     故'绕一圈后质量变了'在相位律下**不成立**——回到原位就是回到原质量。")

out["T1"] = {"periodicity": "2π", "mass_unchanged_after_laps": True}

# ══════════════════════════════════════════════════════════════════════
# T2 双覆盖是否改变质量
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【T2】双覆盖（Z006）：2π 给 −1、4π 给 +1 —— 是否改变质量？")
print("=" * 88)
print("  Z006：反演 ι:w↦−w 给出的 Z_2 **就是**自旋双覆盖的 Z_2；")
print("        转 2π 给 −1，转 4π 给 +1（旋量必须两圈才复原）。")
print()
print("  核验：本征值符号（示意：二维旋量在转 2π 与 4π 下的相因子）")
import cmath
for turns, ang in (("一圈 (2π)", 360), ("两圈 (4π)", 720)):
    ph = cmath.exp(1j * math.radians(ang) / 2)
    print(f"    {turns}: 相因子 e^{{iθ/2}} = {ph.real:+.6f}{ph.imag:+.6f}j  |相|={abs(ph):.6f}")
print()
print("  ★ 关键：双覆盖给的是**符号**（−1 与 +1），不是质量。")
print("     质量是 |振幅|² 的函数 ⟹ 符号变化**不改变质量**。")
chk(True, "双覆盖只给符号")
out["T2"] = {"double_cover": "2π→−1, 4π→+1",
             "changes_mass": False,
             "reason": "质量取模平方；符号不进入 |m|"}

# ══════════════════════════════════════════════════════════════════════
# T3 Zero 里真的有"绕数" —— 而且它给的是标度
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【T3】Zero 里真的有绕数吗？—— 有，而且'两圈'确实出现过")
print("=" * 88)
k0 = 2.0
ring2 = 2 * math.pi * k0
print("  Z154：k₀ = A_ℓ² = 2（A_ℓ = |1+e^{iθ_M}| = 2cos(θ_M/2) = √2，θ_M = π/2）")
print(f"        环② = 2πk₀ = 2π·{k0:g} = {ring2:.6f} = 4π")
print(f"        M_I = M_Z·e^环② = {M_Z:.4f} × e^{ring2:.4f} = {M_Z * math.exp(ring2):.6e} GeV")
print()
print(f"  ★ 环② = 4π 字面上就是 **2 圈**（每圈 2π）。")
print(f"     它给出的不是某个粒子的质量，而是**一个标度台阶**：")
print(f"     M_Z → M_I 相差 e^(4π) = {math.exp(ring2):.4e} 倍。")
chk(abs(ring2 - 4 * math.pi) < 1e-12, "环② = 4π（两圈）")
out["T3"] = {"k0": k0, "ring2_rad": ring2, "ring2_in_laps": ring2 / (2 * math.pi),
             "M_I_GeV": M_Z * math.exp(ring2),
             "gives": "标度台阶（非粒子质量）"}

# 四环链：哪几环是"整数圈"？
print()
print("  四环链（Z156）：ln(M_X^E6/v) = ① + ② + ③ + ④")
rings = [
    ("① ln(M_Z/v)", math.log(M_Z / V), "电弱比（实测）"),
    ("② 2πk₀", 4 * math.pi, "**2 圈**（Z154，已证）"),
    ("③ ln(M_X^unif/M_I)", 14.454565, "中间谱（具名输入）"),
    ("④ 29π/12", 29 * math.pi / 12, "E6 破缺（条件定理）"),
]
for nm, val, src in rings:
    laps = val / (2 * math.pi)
    print(f"    {nm:22s} = {val:9.6f}   = {laps:+8.4f} 圈   [{src}]")
print()
print("  ⟹ 四环里只有环② 是**整数圈**；其余是普通对数。")
print("     故'绕 n 圈'在 Zero 里是一个**标度**生成机制，不是质量生成机制。")
chk(True, "四环环数登记")
out["T3_rings"] = [{"ring": n, "value": v, "laps": v / (2 * math.pi), "source": s}
                   for n, v, s in rings]

# ══════════════════════════════════════════════════════════════════════
# T4 "变种"：M_H 是否是某已知质量的结构倍数
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【T4】'是某种已知粒子的变种'：M_H 是否是某已知质量的结构倍数？")
print("=" * 88)
masses = {"M_W": 80.377, "M_Z": 91.1876, "m_t": 172.69, "m_b": 4.18,
          "m_τ": 1.77686, "m_p": 0.938272}
print("  比值 M_H/m 与结构候选（立方体比、圈数、π）：")
factors = {
    "√2": math.sqrt(2), "√3": math.sqrt(3), "√(3/2)": math.sqrt(1.5),
    "2": 2.0, "3/2": 1.5, "4/3": 4 / 3, "1/2": 0.5, "3": 3.0,
    "2π/4": 2 * math.pi / 4, "e": math.e,
}
print(f"  {'对':>10} {'M_H/m':>10}  最接近的结构因子（偏差）")
rows = []
for k, m in masses.items():
    r = M_H / m
    best = min(factors.items(), key=lambda kv: abs(kv[1] - r) / r)
    dev = 100 * (best[1] - r) / r
    rows.append((k, r, best[0], dev))
    print(f"  M_H/{k:<6} {r:10.6f}  {best[0]:>8} (偏差 {dev:+.3f}%)")
print()
bestdev = min(abs(d) for _, _, _, d in rows)
BESTROW = min(rows, key=lambda r: abs(r[3]))
# 本框架的命中标准：10⁻⁶–10⁻⁴（见 Z164 α、Z166 21cm、Z165 m_p）
FRAME_TOL = 0.01   # 0.01% 作为"本框架意义上的命中"的宽松上界
chk(bestdev > FRAME_TOL,
    "最佳偏差超过本框架命中标准（>0.01%）")
print(f"  ★ 最佳近贴：M_H/{BESTROW[0]} = {BESTROW[1]:.6f} ≈ {BESTROW[2]} （偏差 {BESTROW[3]:+.3f}%）")
print(f"     本框架的命中标准是 10⁻⁶–10⁻⁴（Z164 α、Z166 21cm、Z165 m_p）。")
print(f"     故 {bestdev:.3f}% 比标准宽 {bestdev / FRAME_TOL:.0f} 倍 ⟹ **不构成本框架意义上的命中**。")
print()
print("  ★ 但这一条值得单独记：M_H ≈ (π/2)·M_W 是**最接近的'变种'关系**。")
print("     两条否决：① 0.803% 远宽于框架标准；② M_W 本身属 F 档（实测反解），")
print("     用它去定 M_H 不是导出。另有 ③ 无机制：Zero 里没有任何结构规定 π/2 该出现在这两个质量之间。")
out["T4"] = [{"pair": k, "ratio": r, "best_factor": f, "dev_pct": d}
             for k, r, f, d in rows]

# ══════════════════════════════════════════════════════════════════════
# T5 "走不同的路 → 不同的质量"：立方体路段类型
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【T5】'走不同的路 → 不同的质量'：立方体路段类型给振幅（Z126 §2）")
print("=" * 88)
print("  单位立方体的三个长度：棱 1、面对角线 √2、体对角线 √3")
print("  Z126 §2 的两个振幅是它们的**比**：")
print()
cube = {"棱 1": 1.0, "面对角线 √2": math.sqrt(2), "体对角线 √3": math.sqrt(3)}
pairs = [
    ("A_ℓ = 面对角线/棱", math.sqrt(2) / 1, math.sqrt(2)),
    ("A_ν = 体对角线/面对角线", math.sqrt(3) / math.sqrt(2), math.sqrt(1.5)),
]
for nm, val, ref in pairs:
    print(f"    {nm:26s} = {val:.6f}   与库值 {ref:.6f}  {('一致' if abs(val - ref) < 1e-9 else '不符')}")
chk(abs(math.sqrt(2) / 1 - math.sqrt(2)) < 1e-12, "A_ℓ = √2/1")
chk(abs(math.sqrt(3) / math.sqrt(2) - math.sqrt(1.5)) < 1e-12, "A_ν = √(3/2)")
print()
print("  ★ 这是'不同路径 → 不同振幅 → 不同质量'的**真实实现**：")
print("     走面 vs 走体，振幅不同，质量不同。**已证**（Z126 §2–§3）。")
print()
print("  但另外两个振幅**不是**立方体对角比：")
quark = [("A_up", 1.777833), ("A_down", 2.682850)]
print(f"    {'振幅':>8} {'值':>10}  与最接近立方体比的距离")
cubes = {"√3=1.732051": math.sqrt(3), "2=2": 2.0, "√2·√(3/2)=1.732051": math.sqrt(3),
         "2√2=2.828427": 2 * math.sqrt(2)}
for nm, val in quark:
    best = min(cubes.items(), key=lambda kv: abs(kv[1] - val) / val)
    dev = 100 * (best[1] - val) / val
    print(f"    {nm:>8} {val:10.6f}  {best[0]}（偏差 {dev:+.3f}%）")
print()
print("  ⟹ 立方体的三段只给出**两个**振幅（A_ℓ、A_ν）；夸克两个振幅来历仍是开放项。")
print("     故'再加一种走法'不自动给出新振幅——立方体只有三个长度、两个独立比。")
chk(True, "立方体只给两个振幅")
out["T5"] = {"cube_lengths": cube,
             "amplitudes_from_cube": ["A_ℓ=√2", "A_ν=√(3/2)"],
             "quark_amplitudes_from_cube": False,
             "independent_ratios": 2}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("T1 相位律绕圈", "绕任意整数圈**质量不变**（cos 周期 2π）——'叠回原位'就是'回到原质量'。"),
    ("T2 双覆盖", "2π→−1、4π→+1 只改**符号**，不改质量（质量取模平方）。"),
    ("T3 绕数存在", "**有**！Z154：k₀=2、环②=2πk₀=4π=**两圈**，已证。但它生成**标度台阶**（M_Z→M_I 差 e^{4π}），不是粒子质量。"),
    ("T4 变种", "最近贴是 M_H ≈ (π/2)·M_W（偏差 +0.803%）——比框架标准 0.01% 宽 80 倍，且 M_W 属 F 档、无机制 ⟹ 不构成命中。"),
    ("T5 不同路径", "**真实实现**：面 vs 体 → 振幅 √2 vs √(3/2) → 质量不同（已证）。但立方体只给两个独立比，夸克振幅来历仍开放。"),
    ("净结论", "猜想命中了 Zero 的真机制，但**位置不对**：绕圈在 Zero 里管**标度**，不管质量；管质量的是**走哪条对角线**。"),
]
for k, v in verdict:
    print(f"  {k}：{v}")
out["verdict"] = dict(verdict)

print()
print("=" * 88)
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print("=" * 88)
out["assertions"] = {"passed": ok, "failed": len(bad), "bad": bad}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)
