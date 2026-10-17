#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
anchor_ledger_check.py —— 锚点账（补 #1）
=========================================
处理 D33/D35/D36 的「量纲归约与定标」。判定其数值成果的地基：

  A1 色散对**任意** c 成立 ⟹ "相关 0.9989 测出 c=0.8944" 是重言式（测的是你放进去的 c）
  A2 G_disc 随舞台变 ⟹ 单个数字不可能是台账常数（cubic 1.0000 vs diamond 0.2887）
  A3 D35 §5.2 的 (α,β,γ) 与 SI 的 c、G 互不自洽（文档却称"误差 0.01%"）
  A4 用 (α,β,γ) 反算 c、G 是**恒等式**（不携带信息）
  A5 正确的账：1 个量纲锚 + N 个无量纲数（Buckingham π）

运行：cd ZeroGPU/src && /usr/bin/python3 anchor_ledger_check.py
"""
import json
import os
import sys

import numpy as np

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


# ============================================================
print("=" * 72)
print("A1  色散测量 c：重言式检验（Z104/Z105 的 ω²=c²λ_k）")
print("=" * 72)
N = 256
k_mode = 16
lam = 2 - 2 * np.cos(2 * np.pi * k_mode / N)          # 环图上该本征模式的 λ
j = np.arange(N)
eig = np.cos(2 * np.pi * k_mode * j / N)
print(f"  环图 N={N}，本征模式 k={k_mode}：λ = {lam:.6f}（ω_理论 = c·√λ）")
rows = []
for c in (0.5, 0.8944, 1.2):
    dt, steps = 0.005, 40000                       # T = 200
    x = eig.copy()
    v = np.zeros(N)
    lap = np.roll(x, 1) + np.roll(x, -1) - 2 * x
    v += c * c * lap * dt                          # 半步（∂²x=-c²Lx ⟺ +c²·lap）
    ts, a = [], np.empty(steps)
    for t in range(steps):
        lap = np.roll(x, 1) + np.roll(x, -1) - 2 * x
        v += c * c * lap * dt
        x += v * dt
        a[t] = float(x @ eig) / float(eig @ eig)
        ts.append(t * dt)
    ts = np.array(ts)
    cr = np.where((a[:-1] < 0) & (a[1:] >= 0))[0]   # 正斜率过零点
    tz = ts[cr] - a[cr] * (ts[cr + 1] - ts[cr]) / (a[cr + 1] - a[cr])
    per = float(np.mean(np.diff(tz)))               # 正斜率过零点相隔一个整周期
    w_fit = 2 * np.pi / per
    rows.append((c, w_fit, c * np.sqrt(lam), w_fit / (c * np.sqrt(lam))))
    print(f"  c={c:<6} ω_拟合={w_fit:.6f}  ω_理论=c√λ={c*np.sqrt(lam):.6f}  比={w_fit/(c*np.sqrt(lam)):.6f}")
ratios = [r[3] for r in rows]
ck("A1a 『色散相关 0.9989』构成对 c 的测量（c=0.8944 顶点/代）",
   any(abs(r - 1) > 1e-3 for r in ratios),
   f"三个不同的 c 都给出 ω/(c√λ) = {[round(r,6) for r in ratios]} ⟹ 拟合总是返回你放进去的 c；"
   "该「测量」不携带 c 的信息（c 是方程系数，Z105 §5.3 自认由 CFL 自适应 c²=3.5/λ_max 反推）")
note("c 的真实地位", "方程系数（可任取）", "物理 c 只能由外部锚定，不能由色散「测出」")

# ============================================================
print()
print("=" * 72)
print("A2  G_disc 随舞台变（vs D33 的『G = 2.6842 离散单位』）")
print("=" * 72)


def sigma_of(neighbors):
    m = np.array(neighbors, float)
    c = float(np.trace(m.T @ m) / 3)
    lb = float(np.linalg.norm(m, axis=1).mean())
    return c / 2.0, lb


CUB = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
DIA = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
for lab, nb in (("cubic(z=6)", CUB), ("diamond(z=4)", DIA)):
    s, lb = sigma_of(nb)
    print(f"  {lab:<14} σ={s:.4f} ℓ_b={lb:.4f} ⟹ A = 1/(σℓ_b) = {1/(s*lb):.4f}")
s_c, lb_c = sigma_of(CUB)
s_d, lb_d = sigma_of(DIA)
A_c, A_d = 1 / (s_c * lb_c), 1 / (s_d * lb_d)
ck("A2a 『G = 2.6842（离散单位）』是 Zero 内生的常数（与舞台无关）",
   abs(A_d / A_c - 1) < 0.01,
   f"同一方程在两个舞台上给 A = {A_c:.4f}（cubic）与 {A_d:.4f}（diamond），比 {A_d/A_c:.3f}；"
   "G = 1/Q 是舞台量（见 stage_ledger_check.py B 组），单个数字无脚本不可复现")

# ============================================================
print()
print("=" * 72)
print("A3  D35 §5.2 的 (α,β,γ) 与 SI 的 c、G")
print("=" * 72)
c_SI, G_SI = 2.99792458e8, 6.67430e-11
alpha, beta, gamma = 2.01e-33, 5.99e-42, 9.08e-6          # D35 §5.2 / D36 §1.10 的打印值
c_from, G_from = alpha / beta, alpha * (alpha / beta) ** 2 / gamma
print(f"  α/β = {c_from:.6e} m/s（SI c = {c_SI:.6e}，偏 {c_from/c_SI-1:+.2%}）")
print(f"  αc²/γ = {G_from:.6e}（SI G = {G_SI:.6e}，偏 {G_from/G_SI-1:+.2%}）")
print(f"  由 α 与 SI c、G 反推应有：β = {alpha/c_SI:.4e} s，γ = {alpha*c_SI**2/G_SI:.4e} kg")
ck("A3a D35『验证：c 误差 0.0024%、G 误差 0.0071%』",
   abs(c_from / c_SI - 1) < 1e-4 and abs(G_from / G_SI - 1) < 1e-4,
   f"与它自己打印的 (α,β,γ) 矛盾：c 偏 {c_from/c_SI-1:+.2%}、G 偏 {G_from/G_SI-1:+.2%}；"
   f"β 偏 {beta/(alpha/c_SI):.4f} 倍（≈ 文档自己的 c_disc=0.8944），γ 偏 {gamma/(alpha*c_SI**2/G_SI):.4f} 倍")

# ============================================================
print()
print("=" * 72)
print("A4  『用 (α,β,γ) 反算 c、G』是不是独立检验")
print("=" * 72)
alpha2 = 2.01e-33
beta2 = alpha2 / c_SI            # 按 SI 一致性定义
gamma2 = alpha2 * c_SI ** 2 / G_SI
print(f"  按 SI 一致性定义 β、γ 后：α/β = {alpha2/beta2:.10e}（= SI c，残差 {abs(alpha2/beta2/c_SI-1):.2e}）")
print(f"                            αc²/γ = {alpha2*c_SI**2/gamma2:.10e}（= SI G，残差 {abs(alpha2*c_SI**2/gamma2/G_SI-1):.2e}）")
ck("A4a 『验证 c、G』是对定标的独立检验",
   abs(alpha2 / beta2 / c_SI - 1) > 1e-9,
   "β、γ 若由 SI 的 c、G 定义，则反算必为恒等式（残差 ~1e-16）⟹ 该『验证』只检验了算术，不检验理论")

# ============================================================
print()
print("=" * 72)
print("A5  正确的账：1 个量纲锚 + N 个无量纲数（Buckingham π）")
print("=" * 72)
LEDGER = {
    "锚（有量纲，必须实测）": {"v = 246.22 GeV": "电弱真空期望值（D35 的选择）"},
    "无量纲数（模型内，可算/须付账）": {
        "33.6197 = ln(M_X/v)": "D34 §12 自认『从旧理论（cosmos-construct）嫁接』",
        "c_disc = 0.8944 = √0.8": "Z105 §5.3 自认的 CFL 自适应值（c²=3.5/λ_max）",
        "G_disc = 2.6842": "1/r 拟合读数；随舞台变（A2），无脚本不可复现",
    },
}
for k, v in LEDGER.items():
    print(f"  {k}：")
    for kk, vv in v.items():
        print(f"      {kk:<26} —— {vv}")
print(f"  由锚链：M_X = v·e^33.6197 = {246.22*np.exp(33.6197):.4e} GeV；"
      f"r_int = ħc/M_X = {1.973269804e-16/(246.22*np.exp(33.6197)):.4e} m ✓（与文档 2.01e-33 相符）")
ck("A5a 『框架恰好只需要一个锚』",
   False,
   "实为 1 个量纲锚（v）＋ 3 个无量纲数（33.6197、c_disc、G_disc）；其中 2 个不是台账常数"
   "（c_disc 由 CFL 反推、G_disc 随舞台变），故『只差 0.01%』不是理论精度而是单位换算的算术")

# ============================================================
print()
print("=" * 72)
n_fail = sum(1 for x in CK if not x["pass"])
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 72)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "anchor_ledger_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "c_tautology": [{"c": r[0], "omega_fit": r[1], "omega_theory": r[2], "ratio": r[3]} for r in rows],
               "A_platform": {"cubic": A_c, "diamond": A_d},
               "D35_recompute": {"c_from": c_from, "G_from": G_from,
                                 "beta_should": alpha / c_SI, "gamma_should": alpha * c_SI ** 2 / G_SI}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/anchor_ledger_check.json")
sys.exit(1 if n_fail else 0)