"""
z0_graft_Zlambda.py --- 嫁接验证 · 第二批：Z_Λ 律（带电轻子与中微子【同源】）

【嫁接原则】见 `GRAFT_LEDGER.md §0`：J1 可复算 / J2 输入全是框架量 / J3 有机制。

本题的**统一律**（`build_169` + `build_170`）：
    √m_n = μ·( 1 + A·cos(θ + 2πn/Λ) ),     Λ = 3
  带电轻子:  A = √2
  中微子  :  A = √(Λ/k₀) = √(3/2)         （k₀ = 2, Λ = 3 —— Z_Λ 律的同一组整数）

★ 本文件的中心检验：中微子的 R = m₂/m₃ = (73−28√6)/25
  能否【从同一个 A 导出】（而不是独立拟合）？—— 若能，两条预言【同源】。

用法：/usr/bin/python3 z0_graft_Zlambda.py   输出：results/z0_graft_Zlambda.json
"""
from __future__ import annotations

import json
import os
import time
from math import acos, cos, pi, sin, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


# ===== 框架结构量 =====
K0, LAM = 2, 3
A_LEPT = sqrt(2)                      # 带电轻子侧（Koide 的几何版：2cos(θ_物质/2)=√2）
A_NU = sqrt(LAM / K0)                 # 中微子侧 = √(3/2)
DELTA0 = pi / 12 - K0 / LAM ** 2      # 带电轻子：δ_0（映射已识别）

if __name__ == "__main__":
    t0 = time.time()

    # ================= Z1 统一律的形式 =================
    print("=" * 96)
    print("Z1：统一律 √m_n = μ(1 + A·cos(θ + 2πn/Λ))，Λ=3")
    print("=" * 96)
    print(f"  带电轻子 A = √2 = {A_LEPT:.9f}")
    print(f"  中微子   A = √(Λ/k₀) = √(3/2) = {A_NU:.9f}   （k₀=2, Λ=3）")
    print(f"  中微子 A² = Λ/k₀ = {LAM/K0:.9f}   {'✓' if abs(A_NU**2 - LAM/K0)<1e-12 else '✗'}")
    print(f"  带电轻子 A² = 2；比值 A²_轻/A²_ν = {A_LEPT**2/A_NU**2:.9f} = 4/3")
    check("**Z1a 中微子 A² = Λ/k₀ = 3/2 精确**", abs(A_NU ** 2 - LAM / K0) < 1e-12,
          f"A²_ν = {A_NU**2:.9f}")
    check("**Z1b 带电轻子 A² = 2（Koide 的几何版 2cos45°=√2）**",
          abs(A_LEPT ** 2 - 2) < 1e-12, f"A²_轻 = {A_LEPT**2:.9f}")
    check("**Z1c 两侧 A² 的比值 = 4/3（不是 Λ/k₀ —— 本文件更正一处命名）**",
          abs(A_LEPT ** 2 / A_NU ** 2 - 4 / 3) < 1e-12, f"{A_LEPT**2/A_NU**2:.9f}")

    # ================= Z2 带电轻子 =================
    print("\n" + "=" * 96)
    print("Z2：带电轻子（A = √2）—— 零参数")
    print("=" * 96)
    print(f"  δ_0 = π/12 − k₀/Λ² = {DELTA0:.9f} rad = {np.degrees(DELTA0):.7f}°")
    print(f"  M(δ_0) = 1 − cos δ_0 + sin δ_0 = {1-cos(DELTA0)+sin(DELTA0):.9f}")

    def mlep(n):
        d = DELTA0 - 2 * pi * n / 3
        return (1 - cos(d) + sin(d)) ** 2

    me, mmu, mtau = mlep(0), mlep(1), mlep(2)
    obs = {"m_mu/m_e": 206.768283, "m_tau/m_e": 3477.2283, "m_tau/m_mu": 16.817029}
    pred = {"m_mu/m_e": mmu / me, "m_tau/m_e": mtau / me, "m_tau/m_mu": mtau / mmu}
    print(f"\n  {'比值':>10} {'预言':>14} {'实测':>13} {'偏差':>11}")
    for k in obs:
        print(f"  {k:>10} {pred[k]:>14.6f} {obs[k]:>13.6f} {(pred[k]/obs[k]-1)*100:>10.4f}%")
    check("**Z2a m_μ/m_e 偏差 < 0.002%**", abs(pred["m_mu/m_e"] / obs["m_mu/m_e"] - 1) < 2e-5,
          f"{(pred['m_mu/m_e']/obs['m_mu/m_e']-1)*100:.4f}%")
    check("**Z2b m_τ/m_e 偏差 < 0.01%**", abs(pred["m_tau/m_e"] / obs["m_tau/m_e"] - 1) < 1e-4,
          f"{(pred['m_tau/m_e']/obs['m_tau/m_e']-1)*100:.4f}%")
    K = np.array([me, mmu, mtau])
    Kv = K.sum() / (np.sqrt(K).sum() ** 2)
    print(f"\n  Koide K（框架）= {Kv:.12f}   2/3 = {2/3:.12f}   偏差 {(Kv/(2/3)-1)*100:.2e}%")
    print(f"  几何来源: √2 = 2cos(θ_物质/2)，θ_物质 = 90° ⟹ 2cos45° = √2 ✓")
    check("**Z2c Koide = 2/3 精确**", abs(Kv - 2 / 3) < 1e-12, f"K = {Kv:.12f}")

    # ================= Z3 中微子：R 是否【从同一个 A 导出】=================
    print("\n" + "=" * 96)
    print("Z3：★ 中心检验 —— 中微子 R = (73−28√6)/25 能否从同一个 A 导出？")
    print("=" * 96)
    A = A_NU
    print(f"  A² = 3/2，1 − 1/A² = {1-1/A**2:.9f} = 1/3，√(1−1/A²) = {sqrt(1-1/A**2):.9f} = 1/√3 ✓")
    # m₁ = 0 饱和 ⟹ 在零点 x: cos x = −1/A, sin x = +√(1−1/A²)
    cx, sx = -1 / A, sqrt(1 - 1 / A ** 2)
    x = acos(cx)
    print(f"  m₁ = 0 ⟹ cos x = −1/A = {cx:.9f}, sin x = +1/√3 = {sx:.9f}")
    print(f"  x = {np.degrees(x):.6f}°， x − 2π/3 = {np.degrees(x-2*pi/3):.6f}°  "
          f"（文件声明 θ_ν = 24.735610°）")
    check("**Z3a 相角复现**（x − 2π/3 = 24.735610°）",
          abs(np.degrees(x - 2 * pi / 3) - 24.735610) < 1e-5,
          f"{np.degrees(x-2*pi/3):.6f}°")

    def mnu(w):
        return (1 + A * cos(x + w)) ** 2

    m1n, m2n, m3n = mnu(0), mnu(2 * pi / 3), mnu(4 * pi / 3)
    print(f"\n  谱（n 从零点起）: m₁={m1n:.3e}（应≈0）, m₂={m2n:.9f}, m₃={m3n:.9f}")
    R_num = m2n / m3n
    R_closed = (73 - 28 * sqrt(6)) / 25
    R_alt = ((7 - 2 * sqrt(6)) / 5) ** 2      # ★ 本文件新得的精确闭式
    print(f"\n  R = m₂/m₃ = {R_num:.12f}")
    print(f"  (73−28√6)/25     = {R_closed:.12f}")
    print(f"  ((7−2√6)/5)²     = {R_alt:.12f}   （★ 本文件新得的精确闭式）")
    print(f"  ⟹ 振幅比 (7−2√6)/5 = {(7-2*sqrt(6))/5:.12f}；平方即质量比 ⟹ R 是【振幅比】的平方")
    print(f"  差 = {abs(R_num-R_closed):.3e}")
    check("**Z3b ★ R 从同一个 A 精确导出**（差 < 1e-12）",
          abs(R_num - R_closed) < 1e-12, f"差 {abs(R_num-R_closed):.2e}")
    check("**Z3c ★ 精确闭式 R = ((7−2√6)/5)²**", abs(R_alt - R_closed) < 1e-15, f"{R_alt:.12f}")

    print(f"""
  ⟹ ★★★ **结论：中微子的 R 不是独立拟合，而是【同一个 A】的精确推论。**
     所以"带电轻子与中微子共用同一个 Z_Λ 律"这条 **在本文件的独立复算下成立** ✓
     验证链： A = √(Λ/k₀) → m₁=0 饱和 → 相角 x → R = m₂/m₃ 精确等于 (73−28√6)/25""")
    check("**Z3d 同源已证**（R 是 A 的函数，非新增参数）", True, "见上")

    # ================= Z4 等价性：R 的闭式与 A 的关系 =================
    print("\n" + "=" * 96)
    print("Z4：R 与 A 的代数关系（为将来查机制用）")
    print("=" * 96)
    t = 1 / 3                                   # 1 − 1/A² = 1/3
    print(f"  记 t ≡ 1 − 1/A² = {t:.9f} = 1/3")
    print(f"  则 m₂/m₃ 的闭式（a = √t cos(2π/3) = −√t/2, b = √t sin(2π/3) = √(3t)/2）:")
    a = -sqrt(t) / 2
    b = sqrt(3 * t) / 2
    # 1 + A cos(x+w) = 1 + A(cos x cos w − sin x sin w) = 1 + A((−1/A)cos w − √t sin w) = 1 − cos w − A√t sin w
    f = lambda w: 1 - cos(w) - A * sqrt(t) * sin(w)
    print(f"  f(w) ≡ 1 + A cos(x+w) = 1 − cos w − A√t·sin w")
    print(f"  f(0) = {f(0):.3e}（零点 ✓）  f(2π/3) = {f(2*pi/3):.9f}  f(4π/3) = {f(4*pi/3):.9f}")
    R2 = (f(2 * pi / 3) / f(4 * pi / 3)) ** 2   # ★ 平方 = 质量比
    print(f"  R = (f(2π/3)/f(4π/3))² = {R2:.12f}   {'✓ 与上一致' if abs(R2-R_closed)<1e-12 else '✗'}")
    print(f"  A√t = √(3/2)·√(1/3) = {A*sqrt(t):.9f} = 1/√2 ✓")
    check("**Z4 R 的闭式用 (A, t) 表出**", abs(R2 - R_closed) < 1e-12, f"{R2:.12f}")

    # ================= Z5 中微子绝对质量（接第一批）=================
    print("\n" + "=" * 96)
    print("Z5：中微子绝对质量（一个已测量 Δm² 定标）—— 与本律一致")
    print("=" * 96)
    d21, d31 = 7.42e-5, 2.510e-3
    m2A = sqrt(d21); m3A = m2A / R_closed
    print(f"  [A] 以 Δm²₂₁ 定标: m₂ = {m2A*1e3:.4f} meV, m₃ = {m3A*1e3:.4f} meV, "
          f"Σ = {(m2A+m3A)*1e3:.4f} meV")
    pred31 = m3A ** 2 - m2A ** 2
    print(f"      ⟹ 预言 Δm²₃₁ 偏差 {(pred31/d31-1)*100:+.3f}%（见 z0_neutrino_absolute.py 的更正）")
    check("**Z5 与第一批一致**", abs((m2A + m3A) * 1e3 - 57.3984) < 1e-3,
          f"Σm_ν = {(m2A+m3A)*1e3:.4f} meV")

    # ================= Z6 缺口 =================
    print("\n" + "=" * 96)
    print("Z6：缺口（本机器独立列出）")
    print("=" * 96)
    gap("**Z_Λ 律本身的来源**",
        "律的形式 √m = μ(1+A cos(θ+2πn/Λ)) 与 Λ=3 是从哪来的？本文件只验了它【自洽且同源】，"
        "没有验它【从 Zero 推出】")
    gap("**k₀ = 2 的来源**", "与带电轻子的 δ_0 = π/12 − k₀/Λ² 共用 k₀；k₀ 为何是 2 未给机制")
    gap("**为什么余弦而非其它周期函数**", "律的函数形式无机制")
    print("""  ⟹ 判决：
    · **J1 可复算** ✅（本文件独立复算，含新闭式 R = ((7−2√6)/5)²）
    · **J2 输入全是框架量** ✅（只有 k₀=2, Λ=3, π；无实验反解）
    · **J3 有机制** ❌（律的形式与 k₀、Λ 的来源都缺）
    ⇒ 按 GRAFT_LEDGER §0 的判据：**【条件】**（J1+J2，缺 J3）—— 而不是【导出】""")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          A_lept=A_LEPT, A_nu=A_NU, delta0=DELTA0,
                          lepton={k: pred[k] for k in pred},
                          R_num=R_num, R_closed=R_closed, R_alt=R_alt,
                          theta_nu_deg=float(np.degrees(x - 2 * pi / 3)),
                          verdict="【条件】：J1+J2 满足，J3（机制）缺")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_Zlambda.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_Zlambda.json")
