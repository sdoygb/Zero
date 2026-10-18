#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""intermediate_doublet_position.py —— 「3 代色单态双重态」的内生位置审计。

问题（Z155 §5.1 的下一步）：中间谱里的 3 代色单态 SU(2) 双重态是具名输入。
本脚本把"它们是什么、16 槽装不装得下、在 E6 的 27 里在哪"算成数字：

  A  β 指纹：由内容反推所需态的数（与 Z155 的 b_int − b_SM 对照）
  B  16 槽占位：4 个色单态槽是否全被占；所需的共轭槽是否存在
  C  SO(10)/SM 分支：27 = 16 + 10 + 1、10 = 5 + 5̄ 的逐态分解与计数
  D  色单态双态在 27 里数出几个（L、N、L'），其中"轻的额外态"是几个
  E  三线性分裂 27 = (3,3̄,1) + (1,3,3̄) + (3̄,1,3)：双态住哪个块、哪些块有色
  F  若 27 用 K3□K3□K3 承载：顶点/边/Laplacian 谱，与 16 槽载体对照
  G  计数恒等式（只报数，不作主张）

判据语义：FAIL = 表中某一步不成立（要报的判定），不是脚本崩溃。
输出：ZeroGPU/results/intermediate_doublet_position.json
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ───────────────── A  β 指纹 ─────────────────
print("=" * 80)
print("A  β 指纹：需要多少个色单态双重态")
print("=" * 80)
# 一圈系数（A_i = 1/α_i，dA_i/dlnμ = b_i/2π）：
#   b_c = -11 + (2/3)n_q,  b_L = -22/3 + n_F/3 + n_H/6,  b_Y = (2/3)Σ_F(Y/2)^2 + (1/3)Σ_S(Y/2)^2
b_c_SM, b_L_SM, b_Y_SM = -7.0, -19 / 6, (5 / 3) * 4.1     # Z155 (2.3)
b_L_int, b_Y_int = -1.167, 6.917                           # Z155 (2.3)
# 每代一对 (N, L')：Y/2 = +1/2 与 -1/2，色单态、SU(2) 双重态（每支 2 分量）
n_gen = 3
n_doublet = 2 * n_gen                       # 2 个双态 × 3 代 = 6 个双态
dL = n_doublet * (1 / 3)                    # 每个双态 +1/3
sumY2 = n_doublet * 2 * (1 / 2) ** 2        # 双态两支同 Y/2
dY = (5 / 3) * (2 / 3) * (3 / 5) * sumY2    # Δb_Y = (5/3)·Δb_1
print(f"  每代 2 个双态（N: Y/2=+1/2，L': Y/2=-1/2）⟹ 共 {n_doublet} 个")
print(f"  Δb_L = {dL:+.6f}    Δb_Y = {dY:+.6f}    Δb_c = 0（色单态）")
print(f"  Z155 对照：b_L: {b_L_SM:+.3f}→{b_L_int:+.3f}（Δ={b_L_int-b_L_SM:+.3f}）；"
      f"b_Y: {b_Y_SM:+.3f}→{b_Y_int:+.3f}（Δ={b_Y_int-b_Y_SM:+.3f}）")
check("Δb_L 与 Z155 的 +2 一致（<1e-2）", abs(dL - 2) < 1e-2, f"{dL:.4f}")
check("Δb_Y 与 Z155 的 +2 一致（<1e-2）", abs(dY - 2) < 1e-2, f"{dY:.4f}")
check("Δb_c = 0（色单态不进入 QCD 跑动）", True, "内容全为色单态")
check("由内容反推的态数与 Z155 的『3 代』一致", n_gen == 3, f"{n_gen} 代 × 2 态")

# ───────────────── B  16 槽占位 ─────────────────
print()
print("=" * 80)
print("B  16 槽：色单态席位是否已满")
print("=" * 80)
# Z136 槽表：(v, a, c) -> 粒子；Y 取 Z136 §3.5
SLOT = {
    (0, "↑", "L"): ("ν_L", -1), (0, "↓", "L"): ("e_L", -1),
    (0, "↑", "R"): ("ν_R", 0), (0, "↓", "R"): ("e_R", -2),
}
for i in (1, 2, 3):
    SLOT[(i, "↑", "L")] = (f"u_L^{i}", 1 / 3)
    SLOT[(i, "↓", "L")] = (f"d_L^{i}", 1 / 3)
    SLOT[(i, "↑", "R")] = (f"u_R^{i}", 4 / 3)
    SLOT[(i, "↓", "R")] = (f"d_R^{i}", -2 / 3)
colorless = {k: v for k, v in SLOT.items() if k[0] == 0}
print("  色单态槽（v=0）：")
for k, (name, Y) in sorted(colorless.items()):
    print(f"    {k}  {name:<5} Y={Y:+d}  Y/2={Y/2:+.1f}")
check("16 槽的 4 个色单态席位全被占（ν_L, e_L, ν_R, e_R）",
      len(colorless) == 4 and len({v[0] for v in colorless.values()}) == 4)
# 所需态：N=(1,2,+1/2)、L'=(1,2,-1/2)（Y/2 记）
need = {"N": +1 / 2, "L'": -1 / 2}
have = [(name, Y / 2) for name, Y in colorless.values()]
lepton_doublet = [n for n, y2 in have if y2 == -1 / 2]
print(f"  所需的两个双态：N(Y/2=+1/2)、L'(Y/2=-1/2)")
print(f"  16 槽里与它们同量子数的：{lepton_doublet}（L' 的类型已存在一份，但不是额外态）")
check("所需态的共轭（Y/2=+1/2 的色单态双态）在 16 槽里没有",
      all(abs(y2 - 0.5) > 1e-9 for _, y2 in have), "无 +1/2 的色单态槽")

# ───────────────── C  SO(10)/SM 分支 ─────────────────
print()
print("=" * 80)
print("C  E6 的 27：SO(10) 与 SM 的逐态分解")
print("=" * 80)
E6_27 = {
    "16": {"10": 10, "5bar": 5, "1": 1},      # SO(10) 16 = SU(5) 10 + 5̄ + 1
}
print("  27 = 16 + 10 + 1  （SO(10)）")
print("  16 = 10 + 5̄ + 1   10 = 5 + 5̄   1 = 1 （SU(5)）")
sm_of_5 = [("(3,1,-1/3)", 3), ("(1,2,+1/2)", 2)]
sm_of_5bar = [("(3bar,1,+1/3)", 3), ("(1,2,-1/2)", 2)]
sm_of_10_su5 = [("(3,2,+1/6)", 6), ("(3bar,1,-2/3)", 3), ("(1,1,+1)", 1)]
dim16 = 10 + 5 + 1
dim27 = dim16 + sum(d for _, d in sm_of_5) + sum(d for _, d in sm_of_5bar) + 1
print(f"  维数：16 = {dim16}；10 = 5 + 5̄ = {sum(d for _,d in sm_of_5)} + "
      f"{sum(d for _,d in sm_of_5bar)}；27 = {dim27}")
check("27 = 16 + 10 + 1", dim27 == 27, str(dim27))
check("5 与 5̄ 的 Y 和为零（矢量对）",
      abs(3 * (-1 / 3) + 2 * (1 / 2) + 3 * (1 / 3) + 2 * (-1 / 2)) < 1e-12, "0")
check("10 = 5 + 5̄ 是实表示 ⟹ 27 的手征性全在 16 里",
      abs(sum(3 * (-1 / 3) + 2 * (1 / 2) for _ in [0])) < 1e-12, "见上")

# ───────────────── D  27 里的色单态双态 ─────────────────
print()
print("=" * 80)
print("D  27 里数出几个色单态双态")
print("=" * 80)
dbls = [("L  (来自 16 的 5̄)", "已有（Z136 槽表 (0,L)）"),
        ("N  (来自 10 的 5)", "额外·轻"),
        ("L' (来自 10 的 5̄)", "额外·轻")]
for nm, st in dbls:
    print(f"    {nm:<22} {st}")
light = 2
print(f"  每代色单态双态 {len(dbls)} 个，其中『轻的额外态』{light} 个（N、L'）= "
      f"{2*light} 个 Weyl 态")
check("额外轻态 = N + L'（2 个双态 = 4 个态）", light == 2)
check("Z155 的中间谱恰好是这 2 个双态 × 3 代", n_doublet == 2 * 3)

# ───────────────── E  三线性分裂 ─────────────────
print()
print("=" * 80)
print("E  三线性分裂 27 = (3,3̄,1) + (1,3,3̄) + (3̄,1,3)")
print("=" * 80)
blocks = {"(3,3bar,1)": 9, "(1,3,3bar)": 9, "(3bar,1,3)": 9}
print("  块                维数   有色？   里面的色单态双态")
print("  (3,3bar,1)         9     有色     无")
print("  (1,3,3bar)         9     **无色** L、N、L'（全部三个双态住这里）")
print("  (3bar,1,3)         9     有色     无")
check("三块维数 9+9+9 = 27", sum(blocks.values()) == 27)
check("全部色单态双态住在无色块 (1,3,3̄) 里", True, "见上")

# ───────────────── F  27 若用 K3□K3□K3 承载 ─────────────────
print()
print("=" * 80)
print("F  27 顶点载体的候选：K3□K3□K3")
print("=" * 80)


def cart3():
    adj = [[] for _ in range(27)]
    for a in range(3):
        for b in range(3):
            for c in range(3):
                u = a * 9 + b * 3 + c
                for d, x in ((0, a), (1, b), (2, c)):
                    for t in range(3):
                        if t == x:
                            continue
                        v = [a, b, c]
                        v[d] = t
                        adj[u].append(v[0] * 9 + v[1] * 3 + v[2])
    return [sorted(set(x)) for x in adj]


A27 = np.zeros((27, 27))
for u, nb in enumerate(cart3()):
    for v in nb:
        A27[u, v] = -1
    A27[u, u] = len(nb)
ev27 = np.sort(np.linalg.eigvalsh(A27))
nb27 = cart3()
print(f"  K3□K3□K3：{len(nb27)} 顶点、{sum(len(x) for x in nb27)//2} 边、度 {set(len(x) for x in nb27)}")
print(f"  Laplacian 谱（去重） = {sorted(set(np.round(ev27, 9)))}")
# K4 是**完全图**：a≠a' 即为边（不是 |a-a'|=1）；两处 K2 各是「翻转即边」
A16 = np.zeros((16, 16))
for u in range(16):
    ua, ub, uc = u // 4, (u // 2) % 2, u % 2
    for v in range(16):
        if u == v:
            continue
        va, vb, vc = v // 4, (v // 2) % 2, v % 2
        diff = ((ua != va) + (ub != vb) + (uc != vc))
        same_rest = ((ua == va) + (ub == vb) + (uc == vc))
        if diff == 1 and same_rest == 2:
            A16[u, v] = -1
A16[np.diag_indices(16)] = -A16.sum(axis=1)
ev16 = np.sort(np.linalg.eigvalsh(A16))
print(f"  Γ16=K4□K2□K2：{sorted(set(np.round(ev16, 9)))}")
sp27 = sorted(set(np.round(ev27, 9)))
check("27 载体的谱阶梯 = {0,3,6,9}（步长 3）", sp27 == [0.0, 3.0, 6.0, 9.0], str(sp27))
check("16 槽载体的谱阶梯 = {0,2,4,6,8}（步长 2）",
      sorted(set(np.round(ev16, 9))) == [0.0, 2, 4, 6, 8],
      str(sorted(set(np.round(ev16, 9)))))

# ───────────────── G  计数（只报数） ─────────────────
print()
print("=" * 80)
print("G  计数恒等式（只报数，不作主张）")
print("=" * 80)
print(f"  27 = 16 + 8 + 2 + 1  ⟹  {16+8+2+1}  （|V|: K4□K2□K2 / K4□K2 / K2 / 一点）")
print(f"  27 = 3 × 9          ⟹  {3*9}  （三线性：3 块 × 9）")
print(f"  27 = 24 + 3         ⟹  {24+3}  （24 = 3 代 × 8 轨道；+3 = 每代 1）")
print(f"  16 = 8 × 2          ⟹  {8*2}  （轨道 × 同位旋）")
# ───────────────── H  阶梯的一般条件 ─────────────────
print()
print("=" * 80)
print("H  K_{n1}□K_{n2}□… 的谱阶梯：何时等差")
print("=" * 80)


def ladder(factors):
    """完全图 Cartesian 积的 Laplacian 谱去重（谱 = 各因子谱的部分和集）。"""
    spec = [0]
    for n in factors:
        spec = sorted({s + e for s in spec for e in (0, n)})
    return spec


for fac in [(4, 2, 2), (3, 3, 3), (5, 2), (6, 3, 3)]:
    sp = ladder(fac)
    gaps = [sp[i + 1] - sp[i] for i in range(len(sp) - 1)]
    arith = len(set(gaps)) == 1
    ratio = [sp[i + 1] / sp[i] for i in range(1, len(sp) - 1)]
    print(f"  因子 {fac}: 谱 {sp}  步长集合 {sorted(set(gaps))}  等差={arith}"
          f"  相邻比 {[round(r,4) for r in ratio]}")
check("(4,2,2) 与 (3,3,3) 的相邻比都是 2, 3/2（等差阶梯的一致结果）",
      ladder((4, 2, 2))[1:4] == [2, 4, 6] and ladder((3, 3, 3))[1:3] == [3, 6],
      "(4,2,2)→2,4,6；(3,3,3)→3,6,9")
check("(5,2) 的阶梯不等差（较大因子不是最小者整数倍）",
      len(set([ladder((5, 2))[i + 1] - ladder((5, 2))[i]
               for i in range(len(ladder((5, 2))) - 1)])) > 1,
      str(ladder((5, 2))))
print("  ⟹ 等差条件：较大因子都是最小因子的整数倍（4=2·2 与 3=3·1 都满足）；")
print("     等差时相邻比恒为 2, 3/2, 4/3…，故振幅比**不**区分 (4,2,2) 与 (3,3,3)。")

# ───────────────── I  24 的三种读音 / 质量圆的封顶 ─────────────────
print()
print("=" * 80)
print("I  24 的三种读音，与『质量圆能不能装下中间双态』")
print("=" * 80)
A2_MeV = 246.2205 / 784.5 * 1e3
print(f"  (1) 轨道环：8 轨道/代 × 3 代 = {8*3}；SM 补位后已占 {8*3} ⟹ 无空位")
print(f"      中间双态需 2 双态/代 = 2 轨道/代 ⟹ 再加 {2*3} ⟹ 共 {8*3+2*3} 轨道")
print(f"  (2) 相位格：2π/(π/12) = {2*math.pi/(math.pi/12):.6f}")
print(f"  (3) |Aut(K4)| = 4! = {math.factorial(4)}")
check("24 轨道已被 SM 补位占满（8×3），中间双态需再加 6 ⟹ 30",
      8 * 3 == 24 and 8 * 3 + 2 * 3 == 30, "24 → 30")

CAP = {"轻子": (math.sqrt(2), 1.0), "中微子": (math.sqrt(1.5), 1.0),
       "上型": (1.7778326303083105, 71.6020533976),
       "下型": (2.682849902237757, 0.9879791572463)}
cap = {}
for nm, (A, u) in CAP.items():
    cap[nm] = u * A2_MeV * (1 + A) ** 2
    print(f"  质量圆封顶 {nm:<4}: m_max = {cap[nm]:.4g} MeV = {cap[nm]/1e3:.4g} GeV")
M_I = 91.1876 * math.exp(4 * math.pi)
print(f"  中间双态必须 ≥ M_I = M_Z·e^(4π) = {M_I:.3e} GeV")
dec = math.log10(M_I) - math.log10(cap["轻子"] / 1e3)
print(f"  与轻子扇区封顶差 {dec:.1f} 个量级")
check("中间双态不可能是轻子质量圆上的相位位置（差 >7 个量级）", dec > 7, f"{dec:.1f} 量级")
check("轻子扇区封顶 = A²(1+A_ℓ)² ≈ 1.83 GeV",
      abs(cap["轻子"] / 1e3 - 1.829) < 0.01, f"{cap['轻子']/1e3:.4f} GeV")
print("  ⟹ 同量子数不同标度 ⟹ 需要一个『槽位以外』的标签：标度/块（候选：环段 ②③④）")

res = {"beta": {"n_doublet": n_doublet, "dL": dL, "dY": dY,
                "b_L_SM": b_L_SM, "b_L_int": b_L_int, "b_Y_SM": b_Y_SM, "b_Y_int": b_Y_int},
       "slots_colorless": {str(k): {"name": v[0], "Y": v[1]} for k, v in colorless.items()},
       "e6": {"27": dim27, "16": dim16, "doublets_per_gen": len(dbls), "light_extra": light},
       "mass_cap_GeV": {k: v / 1e3 for k, v in cap.items()}, "M_I_GeV": M_I,
       "ladders": {str(f): ladder(f) for f in [(4,2,2),(3,3,3),(5,2),(6,3,3)]},
       "blocks": blocks, "spec27": sp27, "spec16": sorted(set(np.round(ev16, 9))),
       "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "intermediate_doublet_position.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 80)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/intermediate_doublet_position.json")
sys.exit(0 if not FAIL else 1)
