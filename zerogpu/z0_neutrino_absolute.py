"""
z0_neutrino_absolute.py --- 中微子绝对质量：独立复现

复现 `cosmos-construct/NEUTRINO_ABSOLUTE_MASSES.md`（脚本 `build/build_171`）的整套计算，
并补上原脚本缺的两项：**m_β**（KATRIN）与 **m_ββ**（0νββ），以及**比值张力的准确表述**。

框架给的两条（零参数）：
  m₁ = 0                                  ← 「最轻态饱和」
  R = m₂/m₃ = (73 − 28√6)/25 = 0.176571488 ← Z_Λ 律，A = √(3/2)

定标（**不需要输入**，由任何一条已测量的 Δm² 定死）：
  [A] 以 Δm²₂₁ 定标 ⟹ m₂ = √(Δm²₂₁), m₃ = m₂/R ⟹ **Δm²₃₁ 成为预言**
  [B] 以 Δm²₃₁ 定标 ⟹ m₃ = √(Δm²₃₁), m₂ = R·m₃ ⟹ **Δm²₂₁ 成为预言**

★ 本文件更正的一处（原脚本 `build_171` 第 119 行）：
     p31 = m3A**2          ✗ 漏了 − m₂²
     p31 = m3A**2 − m2A**2 ✓
  后果：原报告 NuFIT 5.2 的偏差为 −5.18% / +5.47%（看起来"两边都差 ~5%"），
  正确值是 **−8.139% / +5.466%**（**不对称**）。

用法：/usr/bin/python3 z0_neutrino_absolute.py   输出：results/z0_neutrino_absolute.json
内存纪律：纯标量计算，峰值可忽略。
"""
from __future__ import annotations

import json
import math
import os
import time
from math import sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []

# 框架（零参数）
R = (73 - 28 * sqrt(6)) / 25          # m₂/m₃
# 全球拟合（正序 NO）
SETS = [("NuFIT 5.2", 7.42e-5, 2.510e-3),
        ("NuFIT 6.0", 7.49e-5, 2.513e-3),
        ("较早的值",  7.53e-5, 2.453e-3)]
# 宇宙学上限（95%）
BOUNDS = {"Planck 2018+BAO": 0.12, "DESI DR2 2024+CMB": 0.072}


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def report(label, d21, d31):
    """两种定标"""
    m2A = sqrt(d21); m3A = m2A / R
    m3B = sqrt(d31); m2B = R * m3B
    p31 = m3A ** 2 - m2A ** 2          # ★ 正确式
    p21 = m2B ** 2                     # m₁ = 0 ⟹ Δm²₂₁ = m₂²
    return dict(label=label, d21=d21, d31=d31,
                A=dict(m2=m2A, m3=m3A, S=m2A + m3A, pred_d31=p31,
                       dev31=p31 / d31 - 1),
                B=dict(m2=m2B, m3=m3B, S=m2B + m3B, pred_d21=p21,
                       dev21=p21 / d21 - 1))


if __name__ == "__main__":
    t0 = time.time()

    print("=" * 96)
    print("N1：框架的两条零参数输入")
    print("=" * 96)
    print(f"  m₁ = 0")
    print(f"  R = m₂/m₃ = (73 − 28√6)/25 = {R:.12f}")
    check("**N1 零参数比 R**（与 cosmos-construct 声明一致）",
          abs(R - 0.176571488) < 1e-8, f"{R:.9f}")

    print("\n" + "=" * 96)
    print("N2：两种定标给出的绝对质量（带单位）")
    print("=" * 96)
    rows = []
    for lab, d21, d31 in SETS:
        r = report(lab, d21, d31)
        rows.append(r)
        print(f"\n  --- {lab}:  Δm²₂₁={d21:.3e}  Δm²₃₁={d31:.3e} eV² ---")
        A, B = r["A"], r["B"]
        print(f"    [A] 以 Δm²₂₁ 定标: m₂={A['m2']*1e3:7.4f} meV  m₃={A['m3']*1e3:7.4f} meV  "
              f"Σ={A['S']*1e3:7.4f} meV")
        print(f"        ⟹ 预言 Δm²₃₁ = {A['pred_d31']:.6e} eV²   偏差 = {A['dev31']*100:+.3f}%")
        print(f"    [B] 以 Δm²₃₁ 定标: m₂={B['m2']*1e3:7.4f} meV  m₃={B['m3']*1e3:7.4f} meV  "
              f"Σ={B['S']*1e3:7.4f} meV")
        print(f"        ⟹ 预言 Δm²₂₁ = {B['pred_d21']:.6e} eV²   偏差 = {B['dev21']*100:+.3f}%")

    r0 = rows[0]                       # NuFIT 5.2
    check("**N2a 代表值（NuFIT 5.2，[A]）**",
          abs(r0["A"]["m2"] * 1e3 - 8.6139) < 1e-3 and abs(r0["A"]["m3"] * 1e3 - 48.7844) < 1e-3
          and abs(r0["A"]["S"] * 1e3 - 57.3984) < 1e-3,
          f"m₂={r0['A']['m2']*1e3:.4f} m₃={r0['A']['m3']*1e3:.4f} Σ={r0['A']['S']*1e3:.4f} meV")
    check("**N2b [B] 代表值（NuFIT 5.2）**",
          abs(r0["B"]["m3"] * 1e3 - 50.0999) < 1e-3 and abs(r0["B"]["S"] * 1e3 - 58.9461) < 1e-3,
          f"m₃={r0['B']['m3']*1e3:.4f} Σ={r0['B']['S']*1e3:.4f} meV")

    print("\n" + "=" * 96)
    print("N3：★ 更正 —— 原脚本的 Δm²₃₁ 偏差是错的")
    print("=" * 96)
    m3A = r0["A"]["m3"]; m2A = r0["A"]["m2"]; d31 = r0["d31"]
    wrong = m3A ** 2 / d31 - 1
    right = (m3A ** 2 - m2A ** 2) / d31 - 1
    print(f"  原脚本（p31 = m3²）        : {wrong*100:+.3f}%   ← 报告值 −5.18%")
    print(f"  正确（p31 = m3² − m2²）    : {right*100:+.3f}%   ← 本文件")
    print(f"  差 = {(right-wrong)*100:+.3f} 个百分点")
    check("**N3 更正确认**（原式漏 −m₂²，正确偏差 −8.139%）",
          abs(right * 100 + 8.139) < 0.01, f"{right*100:+.3f}%")
    check("**N3' 两种定标的偏差【不对称】**（原脚本称『两边都差 ~5%』）",
          abs(abs(r0["A"]["dev31"]) - abs(r0["B"]["dev21"])) > 0.02,
          f"|{r0['A']['dev31']*100:+.3f}%| vs |{r0['B']['dev21']*100:+.3f}%|")

    print("\n" + "=" * 96)
    print("N4：比值张力（定标无关的准确表述）")
    print("=" * 96)
    print("  实测 Δm²₃₁/Δm²₂₁（= 1/R² 因为 m₁ = 0）应等于 1/R²")
    print(f"  1/R² = {1/R**2:.6f}")
    print(f"  {'数据':>12} {'实测 Δm²₃₁/Δm²₂₁':>18} {'偏差':>10}")
    for lab, d21, d31 in SETS:
        obs = d31 / d21
        print(f"  {lab:>12} {obs:>18.6f} {(obs/(1/R**2)-1)*100:>9.3f}%")
    print("""
  ⟹ 这是**定标无关**的张力：形状 R 与两条 Δm² 的比值不一致。
     由于两种定标给出【不同】的预言，理论在这一点上不是自洽的 —— 而是有一个 ~8% 的张力。""")
    obs = SETS[0][2] / SETS[0][1]
    check("**N4 比值张力**（= 1/R² 与实测之比，定标无关）",
          abs(obs / (1 / R ** 2) - 1) > 0.05, f"{SETS[0][0]}: {(obs/(1/R**2)-1)*100:+.3f}%")

    print("\n" + "=" * 96)
    print("N5：新算的两项（原脚本缺）—— m_β 与 m_ββ")
    print("=" * 96)
    # 混合角（NuFIT 5.2 NO 代表值）
    s122, s132 = 0.303, 0.02225
    c132 = 1 - s132
    s13 = sqrt(s132); c13 = sqrt(c132)
    U_e1 = c13 * sqrt(1 - s122)          # |U_e1|
    U_e2 = c13 * sqrt(s122)              # |U_e2|
    U_e3 = s13                           # |U_e3|（c13 已含）
    # m_β² = Σ |U_ei|² m_i²
    for tag, r in (("[A]", r0["A"]), ("[B]", r0["B"])):
        m1, m2, m3 = 0.0, r["m2"], r["m3"]
        mbeta2 = (U_e1 ** 2) * m1 ** 2 + (U_e2 ** 2) * m2 ** 2 + (U_e3 ** 2) * m3 ** 2
        mbb_max = (U_e1 ** 2) * m1 + (U_e2 ** 2) * m2 + (U_e3 ** 2) * m3   # 同相上限
        print(f"  {tag} m_β  = {sqrt(mbeta2)*1e3:7.3f} meV      "
              f"m_ββ(同相上限) = {mbb_max*1e3:7.3f} meV")
    print("""
  对照：
    KATRIN (2024): m_β < 0.45 eV        (90% CL)
    0νββ（下一代，如 nEXO/LEGEND-1000）目标灵敏度 ~ 10–20 meV 量级 ⟹ 可判。""")
    r = r0["A"]
    mbeta2 = (U_e1 ** 2) * 0 + (U_e2 ** 2) * r["m2"] ** 2 + (U_e3 ** 2) * r["m3"] ** 2
    mbb = (U_e1 ** 2) * 0 + (U_e2 ** 2) * r["m2"] + (U_e3 ** 2) * r["m3"]
    check("**N5 m_β 与 m_ββ 已给出**（原脚本缺）",
          sqrt(mbeta2) < 0.45 and 1e-3 < mbb < 20e-3,   # KATRIN 上限 0.45 **eV**
          f"m_β={sqrt(mbeta2)*1e3:.3f} meV, m_ββ≤{mbb*1e3:.3f} meV")

    print("\n" + "=" * 96)
    print("N6：宇宙学上限检验")
    print("=" * 96)
    S = r0["A"]["S"]
    for k, v in BOUNDS.items():
        print(f"  Σm_ν = {S*1e3:.3f} meV = {S:.4f} eV   vs  {k} < {v} eV  "
              f"{'✅ 通过' if S < v else '❌ 超限'}")
    check("**N6 Σm_ν 通过当前所有宇宙学上限**",
          all(S < v for v in BOUNDS.values()), f"Σm_ν = {S:.4f} eV")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    print("=" * 96)
    print("""结论：
  · 框架给**零参数形状**（m₁=0, R=(73−28√6)/25）；
  · **一个已测量的 Δm² 定标** ⟹ **绝对质量**（不需『输入』）；
  · 代表值（NuFIT 5.2, [A]）：m₂ = 8.614 meV, m₃ = 48.784 meV, Σm_ν = 57.40 meV；
  · Σm_ν = 0.0574 eV **通过 Planck(<0.12) 与 DESI DR2(<0.072)** ⟹ **可否证预言**；
  · ⚠ **张力**：另一个 Δm² 偏离 **−8.14%**（[A]）或 **+5.47%**（[B]）——
    本文件更正了原脚本在此处的算式错误（原报 −5.18%，因漏 −m₂²）。

诚实边界：
  · R = (73−28√6)/25 的**几何推导**在 cosmos-construct 的 build_170；本文件只复现其下游；
  · 混合角（s122=0.303, s132=0.02225）是**外部输入**，用于 m_β/m_ββ；
  · 0νββ 的 m_ββ 依赖 Majorana 相位（本文件取同相上限）。""")

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL,
                          R=R, rows=rows,
                          correction=dict(wrong_pct=wrong * 100, right_pct=right * 100),
                          m_beta_meV=float(sqrt(mbeta2) * 1e3),
                          m_bb_max_meV=float(mbb * 1e3),
                          bounds=BOUNDS)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_neutrino_absolute.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_neutrino_absolute.json")
