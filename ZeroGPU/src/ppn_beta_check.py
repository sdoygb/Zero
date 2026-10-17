#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ppn_beta_check.py —— 补 #4：PPN β 的范畴错误

【口径】本脚本里的 FAIL ＝ **该篇要报告的判定**（「预期形态不成立」），不是脚本错误；
        退出码 1 表示存在此类判定。读退出码前先读本行。
=============================================
对象：D33 §2–§4 与 D34 §10「非线性已经内生在体系里——来自 freq_mass 函数（sort+log）」，
      以及 D36 §3.8.2「近日点进动：空间弯曲(γ) + 非线性修正(β)」。

要钉死的一件事：**非线性发生在哪一端**
  GR 的 PPN β 由**场方程**的非线性定义（引力场自身的能量也产生引力 ⟹ 自相互作用）；
  D33/D34 的非线性在**源映射**（mass → x）上 —— 位置不同，不能充当 β。

机器核验：
  P1 台账场方程 Lφ=−x 对源是**线性**的（叠加原理精确成立）
  P2 台账里 φ 不进入源 ⟹ **自相互作用项恒为 0**（结构性）
  P3 给场方程加自相互作用 Lφ = −x(1+κφ) 后，叠加原理才被破坏 —— 演示"β 要在哪一端加"
  P4 该机制能否给出水星近日点进动（需要 β 与 γ 都在场方程层）

运行：cd ZeroGPU/src && /usr/bin/python3 ppn_beta_check.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import torus_adj      # noqa: E402

RES = os.path.join(HERE, "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


# ============================================================
# 构台：3D 立方格上的图 Laplacian（台账 torus_adj）
# ============================================================
D, Ls = 3, 6
adj = torus_adj(D, Ls)
V = Ls ** D
L = np.zeros((V, V))
for u, nb in enumerate(adj):
    L[u, u] = len(nb)
    for w in nb:
        L[u, w] -= 1.0
print(f"  舞台：torus_adj({D},{Ls})，|V|={V}（图 Laplacian，L ≥ 0，ker L = 常量）")


def solve_poisson(x, kappa=0.0):
    """解 Lφ = −x（kappa=0）或 Lφ = −x(1+κφ)（自相互作用，用不动点迭代）"""
    x = x - x.mean()
    if kappa == 0.0:
        return -np.linalg.lstsq(L, x, rcond=None)[0]
    phi = -np.linalg.lstsq(L, x, rcond=None)[0]
    for _ in range(200):                     # 不动点迭代解非线性方程
        rhs = x * (1 + kappa * phi)
        rhs = rhs - rhs.mean()
        phi_new = -np.linalg.lstsq(L, rhs, rcond=None)[0]
        if np.max(np.abs(phi_new - phi)) < 1e-14:
            phi = phi_new
            break
        phi = 0.5 * phi + 0.5 * phi_new       # 阻尼
    return phi


# ============================================================
print("=" * 72)
print("P1  台账场方程的线性性（叠加原理）")
print("=" * 72)
rng = np.random.default_rng(20261016)
x1 = rng.normal(size=V); x1 -= x1.mean()
x2 = rng.normal(size=V); x2 -= x2.mean()
p1, p2, p12 = solve_poisson(x1), solve_poisson(x2), solve_poisson(x1 + x2)
err = float(np.max(np.abs(p12 - (p1 + p2))))
scale = float(np.max(np.abs(p12)))
print(f"  φ[x1+x2] − (φ[x1]+φ[x2]) 的最大偏差 = {err:.3e}（相对 {err/scale:.2e}）")
ck("P1a 台账场方程 Lφ=−x 是**非线性**的（源映射的非线性会传进场方程）",
   err > 1e-8,
   f"叠加原理精确成立（相对偏差 {err/scale:.1e}）⟹ 场方程端**没有**非线性；"
   "D33/D34 的『非线性已内生』不成立")

# ============================================================
print()
print("=" * 72)
print("P2  自相互作用项：φ 是否进入源")
print("=" * 72)
print("  台账：Lφ = −x，其中 x = 源（L2 物质/非中心态分布），φ **不出现**在 x 里")
print("  ⟹ 自相互作用项 ≡ 0 ⟹ 引力场自身的能量不产生引力")
ck("P2a 台账里存在『引力自能也产生引力』这一项（PPN β 的来源）",
   False,
   "源 x 只由物质分布给出，不含 φ；GR 的 β 正来自 T_μν 里的引力场贡献 ⟹ 台账无此机制（结构性 no-go）")

# ============================================================
print()
print("=" * 72)
print("P3  演示：β 要加在**场方程端**才改变叠加性")
print("=" * 72)
print(f"  {'κ（自相互作用强度）':>20} | {'叠加偏差':>12} | {'相对':>10}")
for kappa in (0.0, 1e-3, 1e-2, 5e-2):
    q1, q2 = solve_poisson(x1, kappa), solve_poisson(x2, kappa)
    q12 = solve_poisson(x1 + x2, kappa)
    e = float(np.max(np.abs(q12 - (q1 + q2))))
    print(f"  {kappa:>20.4g} | {e:>12.3e} | {e/float(np.max(np.abs(q12))):>10.2e}")
q_small = solve_poisson(x1, 0.0)
q_large = solve_poisson(x1, 0.05)
print(f"  κ=0.05 时同一源 x1 的势变化幅度 = {float(np.max(np.abs(q_large-q_small)))/float(np.max(np.abs(q_small))):.2%}")
note("自相互作用位置", "必须加在场方程端", "台账的场方程端是线性的；源端的非线性不改变叠加性")

# ============================================================
print()
print("=" * 72)
print("P4  水星近日点进动是否可给出")
print("=" * 72)
print("  进动率 ∝ (2γ − β + 2)/3 形式的组合，其中 γ、β **都**在场方程层定义：")
print("    γ（空间曲率）—— 见补 #3：台账无 c_eff(Φ)/h_xx 与源的耦合；")
print("    β（非线性）—— 本节 P1/P2：台账场方程线性、无自能项。")
ck("P4a 该机制可给出水星近日点进动（43″/世纪）",
   False,
   "γ 与 β 两个参数在台账里都没有定义（γ 缺 h_xx 耦合、β 缺自相互作用）⟹ 不是『算得不准』，是"
   "『两个参数都不存在』；D36 §3.8.2 的『空间弯曲 + 非线性 → 进动』在台账里没有对应物")

# ============================================================
print()
print("=" * 72)
print("P5  源端的非线性是什么（把范畴错误写清楚）")
print("=" * 72)
print("  D33/D34 的非线性：mass → x 的映射（freq_mass = sort + log）")
print("    · 它是**读出端**的非线性：同一物质分布，读出方式不同则 x 不同")
print("    · 它不改变场方程的线性性（P1 已测：对给定的 x，φ 线性）")
print("  GR 的 PPN β：**场方程端**的非线性（引力能产生引力）")
print("    · 它使『两个源的总势 ≠ 各自势之和』（引力自能）")
print("  ⟹ 两者位置不同：前者是 x(mass) 的曲率，后者是 φ(x) 的曲率。前者无法替代后者。")
ck("P5a 源映射的非线性（sort+log）可以充当 PPN β",
   False,
   "范畴错误：β 定义在 φ(x) 端（自相互作用），台账该端线性；x(mass) 端的非线性是读出约定，"
   "不产生引力自能")

# ============================================================
print()
print("=" * 72)
n_fail = sum(1 for x in CK if not x["pass"])
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 72)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "ppn_beta_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "superposition_err": err, "superposition_rel": err / scale,
               "kappa_scan": [{"kappa": k} for k in (0.0, 1e-3, 1e-2, 5e-2)],
               "config": {"D": D, "L": Ls, "V": V}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/ppn_beta_check.json")
sys.exit(1 if n_fail else 0)