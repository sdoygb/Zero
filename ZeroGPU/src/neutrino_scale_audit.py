#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_scale_audit.py —— 中微子绝对标度「(2/9)^15」的审计（推导稿 §18 / 审计台账）。

四件事，逐件独立核验：
  A  无自由参数的比值式：m_νk/m_e = (2/9)^15·tan²(kπ/8) / M_e²，
     M_e = 1+√2·cos(2π/3+2/9)。整条预测不含任何拟合标度（A² 消掉）。
  B  引擎实测：2/9 在忠实引擎里是什么角色（稳态块的概率乘积），
     以及引擎真正的「每步／每代」因子是多少（每步 2 的支化、每代闭合计数 6）。
  C  多重性审计：把库里出现过的原生数（每个带出处）当底、整数当幂，
     看 (2/9)^15 在命中排名里排第几；报告全部 0.5% 内命中。
  D  零模型：把目标质量换成对数均匀随机值，重跑同一搜索，
     量出「0.3% / 0.1% 级命中」的假阳性率。
  E  判别力：竞争候选的 m3 之差与 JUNO 精度比较。

判据语义：本脚本的 FAIL 是「预期的核验事项不成立」= 该审计要报的判定，
不是脚本崩溃；正常结束一律 exit 0（除真异常）。

输出：ZeroGPU/results/neutrino_scale_audit.json
用法：python3 ZeroGPU/src/neutrino_scale_audit.py
"""
from __future__ import annotations

import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ───────────────────────────── 常量（库内出处随行） ─────────────────────────────
V_EW = 246.2205              # GeV。由 A²·784.5 反推：0.313856 GeV × 784.5
A2 = V_EW / 784.5           # 轻子质量单元（Z121 质量律 m=A²M²；c_ℓ=1/784.5 见 Z037/入口稿）
ME = 0.51099895e-3          # GeV（PDG）
M3_OBS = 50.10e-12          # GeV（Δm²31 定标，正序 m1=0）
DM21, DM31 = 7.42e-5, 2.510e-3   # eV²（NuFIT 5.3 附近）
DM21b, DM31b = 7.49e-5, 2.513e-3  # eV²（NuFIT 6.0，版本稳健性）
TH12, TH13 = math.radians(33.4), math.radians(8.6)

TWO_NINE = 2 / 9
N_STEP = 15
TAN2_PI8 = math.tan(math.pi / 8) ** 2          # 3−2√2
TARGET = M3_OBS / A2                            # = 1.5963e-10


def tan2(k):
    return math.tan(k * math.pi / 8) ** 2


def mu_e(two_nine=TWO_NINE):
    """M_e = 1+√2cos(2π/3+2/9)（Z121／Z123 的电子振幅）。"""
    return 1 + math.sqrt(2) * math.cos(2 * math.pi / 3 + two_nine)


# ────────────────────────── A  无自由参数的比值式 ──────────────────────────
print("=" * 78)
print("A  无自由参数的比值式（预测 m_ν3 不引入任何拟合标度）")
print("=" * 78)
Me = mu_e()
ratio = TWO_NINE ** N_STEP / Me ** 2
m3_pred = ME * ratio
print(f"  M_e = 1+√2cos(2π/3+2/9) = {Me:.9f}   M_e² = {Me**2:.9e}")
print(f"  m_ν3/m_e = (2/9)^15/M_e² = {ratio:.6e}")
print(f"  m_ν3 = {m3_pred*1e12:.4f} meV   （观测 50.10：{100*(m3_pred/M3_OBS-1):+.3f}%）")
check("比值式 m_ν3/m_e=(2/9)^15/M_e² 与 50.10 meV 差 <0.5%",
      abs(m3_pred / M3_OBS - 1) < 0.005, f"{m3_pred*1e12:.3f} meV")
check("PDG m_e 与 A²·M_e² 一致到 0.02%（c_ℓ=1/784.5 的取整）",
      abs(A2 * Me ** 2 / ME - 1) < 2e-4, f"A²M_e²/m_e−1 = {100*(A2*Me**2/ME-1):+.4f}%")

masses = [ME * ratio * tan2(k) * 1e12 for k in range(3)]
R_model = masses[1] / masses[2]
X_model = R_model ** 2
R_obs = math.sqrt(DM21 / DM31)
print(f"  三代 = {[round(x,4) for x in masses]} meV   Σ = {sum(masses):.3f} meV")
print(f"  R = {R_model:.6f}（tan²(π/8)=3−2√2={3-2*math.sqrt(2):.6f}，观测 {R_obs:.6f}："
      f"{100*(R_model/R_obs-1):+.3f}%）")
print(f"  X = R² = {X_model:.6f}（观测 {DM21/DM31:.6f}：{100*(X_model/(DM21/DM31)-1):+.3f}%）")
check("R = tan²(π/8) 与观测差 <0.5%", abs(R_model / R_obs - 1) < 0.005,
      f"{100*(R_model/R_obs-1):+.3f}%")
Ue = [math.cos(TH13) * math.cos(TH12), math.cos(TH13) * math.sin(TH12), math.sin(TH13)]
mbb = abs(sum(Ue[i] ** 2 * masses[i] for i in range(3)))
print(f"  m_ββ（正序，m1=0）= {mbb:.3f} meV")
check("m1 = 0（tan 0）", masses[0] == 0.0, "0")

# 版本稳健性（NuFIT 6.0）
R6 = math.sqrt(DM21b / DM31b)
print(f"  [稳健性] NuFIT6.0：R_obs={R6:.6f}，模型偏差 {100*(R_model/R6-1):+.3f}%；"
      f"X 偏差 {100*(X_model/(DM21b/DM31b)-1):+.3f}%")

# ────────────────────────── B  引擎实测：2/9 的角色 ──────────────────────────
print()
print("=" * 78)
print("B  忠实引擎实测：2/9 是什么角色，引擎的每步／每代因子是多少")
print("=" * 78)
import numpy as np
sys.path.insert(0, HERE)
from z0_engine import path_adj                      # noqa: E402
from z0_vectorized import ZeroCoreVec               # noqa: E402

V = 16
adj = path_adj(V)
A = np.zeros((V, V))
for u, nb in enumerate(adj):
    for w in nb:
        A[u][w] = 1.0
A4 = np.linalg.matrix_power(A, 4)
closed = [int(round(A4[u, u])) for u in range(V)]
print(f"  每代闭合计数 C_4(u) = (A^4)_uu = {closed}")
print(f"  内部顶点 C_4 = {set(closed[2:14])}；端点 = {closed[0]},{closed[-1]}")
print(f"  每步支化 = 路径图内点度数 2 ⟹ 每步因子 2（log2 = {math.log(2):.4f} rad = "
      f"{math.degrees(math.log(2)):.1f}°）")

seed = [2.0 ** (15 - v) for v in range(V)]
vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
for _ in range(60000):
    vec.run_generation()
w = np.array(vec.record, dtype=float)
w = w / w.sum()
support = [(i, float(w[i])) for i in range(V) if w[i] > 1e-6]
print(f"  6 万代稳态 ω（支撑） = {[(i, round(x, 8)) for i, x in support]}")
q = float((w ** 2).sum())
prod = float(np.prod([x for _, x in support]))
ratio_br = support[0][1] / support[1][1]
print(f"  q = Σω² = {q:.9f}（5/9 = {5/9:.9f}）；块内两支比 = {ratio_br:.6f}；"
      f"块概率乘积 = {prod:.9f}（2/9 = {TWO_NINE:.9f}）")
check("稳态支撑恰 2 支（2:1，7 与 8）", len(support) == 2 and abs(ratio_br - 2) < 1e-6,
      f"{len(support)} 支")
check("稳态 q = 5/9（Z088）", abs(q - 5 / 9) < 1e-6, f"{q:.9f}")
check("两支概率乘积 = 2/9（精确）", abs(prod - TWO_NINE) < 1e-6, f"{prod:.9f}")
check("引擎每代因子是放大（6），不是 2/9", max(closed) == 6 and all(c > 1 for c in closed[2:14]),
      f"内部 {max(closed)}")
print("  ⟹ 2/9 是「稳态块两支概率的乘积」（每代一次），不是每步乘子；")
print("     引擎的每步因子是 2（支化）、每代因子是 6（闭合词数 C(4,2)）。")

# ────────── B2  同一引擎换个图：2:1 块在闭圆与载体上是否还在 ──────────
print()
print("  B2 跨图对照（同一几何 seed、同一 last/L=4、同一 6 万代）")


def cycle_adj(n):
    return [[(v - 1) % n, (v + 1) % n] for v in range(n)]


def carrier():
    """K₄□K₂□K₂：16 顶点 (a,b,c)，a∈Z₄、b,c∈Z₂；邻接 = 换一个坐标。"""
    def idx(a, b, c):
        return a * 4 + b * 2 + c
    adj = [[] for _ in range(16)]
    for a in range(4):
        for b in range(2):
            for c in range(2):
                u = idx(a, b, c)
                adj[u] += [idx(a2, b, c) for a2 in range(4) if a2 != a]
                adj[u] += [idx(a, 1 - b, c), idx(a, b, 1 - c)]
    return [sorted(set(x)) for x in adj]


def run_graph(adj):
    V = len(adj)
    seed = [2.0 ** (15 - v) for v in range(V)]
    vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    for _ in range(60000):
        vec.run_generation()
    w = np.array(vec.record, dtype=float)
    w = w / w.sum()
    return [float(x) for x in w]


GRAPH = {"path16（开路径）": path_adj(16), "C16（闭圆）": cycle_adj(16),
         "K₄□K₂□K₂（16 槽载体）": carrier()}
B2 = {}
for nm, adj in GRAPH.items():
    A_ = np.zeros((len(adj), len(adj)))
    for u, nb in enumerate(adj):
        for v_ in nb:
            A_[u][v_] = 1.0
    A4_ = np.linalg.matrix_power(A_, 4)
    c4 = sorted({int(round(A4_[u, u])) for u in range(len(adj))})
    w_ = run_graph(adj)
    sup = [(i, x) for i, x in enumerate(w_) if x > 1e-6]
    qq = sum(x * x for x in w_)
    B2[nm] = {"deg": sorted({len(x) for x in adj}), "C4_values": c4,
              "n_support": len(sup), "q": qq,
              "block_product": (sup[0][1] * sup[1][1]) if len(sup) == 2 else None}
    print(f"    {nm:<22} 度={sorted({len(x) for x in adj})}  C4={c4}  支撑={len(sup)}  "
          f"q={qq:.6f}")
    if len(sup) == 2:
        print(f"      → 两支 {sup[0][1]:.4f}:{sup[1][1]:.4f}，块乘积 = {sup[0][1]*sup[1][1]:.9f}")
    else:
        print(f"      → 权重 = 几何 seed 本身（{sup[0][1]:.4f},{sup[1][1]:.4f},"
              f"{sup[2][1]:.4f},…），无 2:1 两支块")
check("2:1 两支块只出现在开路径（闭圆与载体上退化为几何继承）",
      B2["path16（开路径）"]["n_support"] == 2 and
      B2["C16（闭圆）"]["n_support"] > 2 and
      B2["K₄□K₂□K₂（16 槽载体）"]["n_support"] > 2,
      "见上三行")

# B3  载体上确有哈密顿圈（独立构造）
def hamilton_cycle(adj):
    n = len(adj)
    path_ = [0]
    used = [False] * n
    used[0] = True

    def dfs(u):
        if len(path_) == n:
            return 0 in adj[u]
        for v_ in adj[u]:
            if not used[v_]:
                used[v_] = True
                path_.append(v_)
                if dfs(v_):
                    return True
                path_.pop()
                used[v_] = False
        return False

    return path_ if dfs(0) else None


cyc = hamilton_cycle(carrier())
all_adj = all(cyc[i + 1] in carrier()[cyc[i]] for i in range(len(cyc) - 1)) if cyc else False
closed_ok = cyc is not None and cyc[-1] in carrier()[cyc[0]]
print(f"    B3 载体上的哈密顿圈：{'构造成功' if cyc else '未找到'}；"
      f"相邻全成立={all_adj}；首尾闭合={closed_ok}；步数={len(cyc) if cyc else 0}")
check("16 槽载体上存在哈密顿圈（16 步全相邻、含闭合）", bool(cyc) and all_adj and closed_ok,
      f"{len(cyc) if cyc else 0} 步")

# B4  载体的代数谱（供"15 个因子的乘积"这条路线用）
Ac = np.zeros((16, 16))
for u, nb in enumerate(carrier()):
    for v_ in nb:
        Ac[u][v_] = 1.0
ev = np.sort(np.linalg.eigvalsh(Ac))
uniq = []
for x in ev:
    if not uniq or abs(x - uniq[-1][0]) > 1e-9:
        uniq.append([x, 1])
    else:
        uniq[-1][1] += 1
prod15 = float(np.prod(ev[:-1]))       # 去掉最大本征值 λmax=5 后的 15 个
print(f"    B4 载体邻接谱 = " + "、".join(f"{round(x,6)}^{{{m}}}" for x, m in uniq))
print(f"       det A = {np.linalg.det(Ac):.6f}；15 个非平凡本征值之积 = {prod15:.6f}"
      f"（= −3^5 = −243）；λmax = {ev[-1]:.1f}（= 度数）")
check("载体 15 个非平凡本征值之积 = −3^5", abs(prod15 + 243) < 1e-6, f"{prod15:.6f}")
B4 = {"spectrum": [[float(x), int(m)] for x, m in uniq], "det": float(np.linalg.det(Ac)),
      "prod_15": prod15, "lambda_max": float(ev[-1])}

# ────────────────────── C  多重性审计（原生数当底、整数当幂） ──────────────────────
print()
print("=" * 78)
print("C  多重性审计：原生底 × 整数幂，谁命中 50.10 meV")
print("=" * 78)
# 每个底都带库内出处（出处写在注释里，生成台账时照录）
BASES = {
    "2/9": (2 / 9, "(2/3)(1/3) 非中心两支积 Z123"),
    "2/3": (2 / 3, "稳态主支权重 Z087/Z088"),
    "1/3": (1 / 3, "稳态次支权重 Z087/Z088"),
    "5/9": (5 / 9, "q=4/9+1/9 Z088"),
    "3/8": (3 / 8, "sin²θ_w|统一 精确 Z151；零点 3/8 圈 Z121"),
    "7/12": (7 / 12, "K_ν（嫁接台账）"),
    "3/4": (3 / 4, "λ_max/|V| = 3/4（K₄ Perron）"),
    "1/2": (1 / 2, "半满槽 1/2（Z136）；二部 λ₂"),
    "1/4": (1 / 4, "(3/4)(1/3) 组合"),
    "1/8": (1 / 8, "八分圆 1/8 圈（Z121）"),
    "1/16": (1 / 16, "16 槽（Z136）"),
    "1/12": (1 / 12, "12 格（π/12 量子 Z123）"),
    "1/15": (1 / 15, "15 = 16−1"),
    "1/24": (1 / 24, "|S₄| = |Aut K₄| Z021；24 相位格"),
    "1/28": (1 / 28, "binom(8,2)=28"),
    "1/56": (1 / 56, "2·binom(8,2)=56（Lean）"),
    "1/784.5": (1 / 784.5, "c_ℓ = 1/784.5（轻子单元）"),
    "√2/3": (math.sqrt(2) / 3, "A_ℓ/3 = 面对角比/代"),
    "√3/2": (math.sqrt(3) / 2, "体对角/棱·1/2 = cos30°"),
    "1/√2": (1 / math.sqrt(2), "反对称投影 1/√2"),
    "1/√3": (1 / math.sqrt(3), "立方体 1/√3"),
    "tan(π/8)": (math.tan(math.pi / 8), "八分圆半角正切 √2−1"),
    "tan²(π/8)": (math.tan(math.pi / 8) ** 2, "形状比 3−2√2（本篇）"),
    "4/9": (4 / 9, "q 的主支 (2/3)²"),
    "8/27": (8 / 27, "(2/3)³"),
    "2/27": (2 / 27, "(2/3)(1/3)²/…"),
    "4/81": (4 / 81, "(2/9)²"),
    "1/9": (1 / 9, "(1/3)²"),
    "1/6": (1 / 6, "|E(K₄)|=6 的倒数"),
    "1/5": (1 / 5, "K₄□K₂□K₂ 度数 5 的倒数"),
    "3/16": (3 / 16, "3 代/16 槽"),
    "5/16": (5 / 16, "度数/槽数"),
    "2/7": (2 / 7, "二支/七边形（§6 已否）"),
    "π/24": (math.pi / 24, "π/12 的半格"),
}
rows = []
for name, (b, src) in BASES.items():
    if not (0 < b < 1):
        continue
    n = math.log(TARGET) / math.log(b)
    ni = round(n)
    if ni < 1:
        continue
    dev = b ** ni / TARGET - 1
    rows.append({"base": name, "sym": b, "n": ni,
                 "n_implied": n, "dev_pct": 100 * dev,
                 "m3_meV": A2 * b ** ni * 1e12, "source": src})
rows.sort(key=lambda r: abs(r["dev_pct"]))
print(f"  {'底':<11}{'幂':>4}{'m3(meV)':>11}{'偏差':>10}   出处")
for r in rows[:10]:
    print(f"  {r['base']:<11}{r['n']:>4}{r['m3_meV']:>11.3f}{r['dev_pct']:>9.3f}%   {r['source']}")
for thr in (0.5, 0.3, 0.2, 0.1):
    k = sum(1 for r in rows if abs(r["dev_pct"]) < thr)
    print(f"  命中 |偏差|<{thr}% 的 (底,幂) 对数：{k} / {len(rows)}")
best = rows[0]
print(f"  ⟹ 最优命中：({best['base']})^{best['n']}，偏差 {best['dev_pct']:+.3f}%"
      f"（(2/9)^15 为 {[r for r in rows if r['base']=='2/9'][0]['dev_pct']:+.3f}%）")
check("最优单底幂命中 ≤0.3%（模型自身水平）", abs(best["dev_pct"]) < 0.3,
      f"{best['base']}^{best['n']} {best['dev_pct']:+.3f}%")

# 幂必须是库内计数（严格词表）与「计数−1」（本篇自己的约定）两种口径
COUNTS = {2, 3, 4, 6, 8, 12, 15, 16, 24, 28, 56}
strict = [r for r in rows if r["n"] in COUNTS]
minus1 = [r for r in rows if r["n"] in {c - 1 for c in COUNTS}]
print(f"  口径①（幂∈库内计数 {sorted(COUNTS)}）：命中 {[ (r['base'],r['n'],round(r['dev_pct'],3)) for r in strict if abs(r['dev_pct'])<0.5]}")
print(f"  口径②（幂∈计数−1，本篇 15=16−1 用的就是这条）："
      f"{[ (r['base'],r['n'],round(r['dev_pct'],3)) for r in minus1 if abs(r['dev_pct'])<0.5]}")

# ────────────── C2  形状的多重性（正 m 分之一圈的正切） ──────────────
print()
print("  形状侧：R = tan²(以 1/m 圈为步长的角)，看 m 是否被数据唯一选中")
R_obs_v = math.sqrt(DM21 / DM31)
cells = []
for m in range(3, 65):
    ang = 2 * math.pi / m
    R = math.tan(ang) ** 2
    cells.append((abs(R / R_obs_v - 1), m, R))
cells.sort()
for e, m, R in cells[:5]:
    print(f"    m={m:2d} 一 步 = 2π/{m}   R={R:.6f}   偏差 {100*(R/R_obs_v-1):+.3f}%")
h_shape = sum(1 for e, _, _ in cells if e < 0.005)
print(f"    命中 |偏差|<0.5% 的 m 数：{h_shape} / {len(cells)}（m=3..64）")
check("形状侧命中唯一（1/m 圈的 m 只有 16 落在 0.5% 内）", h_shape == 1,
      f"{h_shape} 个：{[(m, round(100*(r/R_obs_v-1),3)) for e,m,r in cells if e<0.005]}")

# ────────────────────────── D  零模型：假阳性率 ──────────────────────────
print()
print("=" * 78)
print("D  零模型：同一搜索在随机目标上的命中率（假阳性率）")
print("=" * 78)
logs = [math.log(b) for b, _ in BASES.values() if 0 < b < 1]


def any_hit(ln_target, eps):
    """给定目标对数，是否存在某底 + 整数幂使 |b^n/T−1|<eps。"""
    for L in logs:
        if L == 0:
            continue
        n = ln_target / L
        if abs(n - round(n)) < eps / abs(L) and round(n) >= 1:
            return True
    return False


lnT = math.log(TARGET)
ln_lo, ln_hi = math.log(M3_OBS / 3 / A2), math.log(M3_OBS * 3 / A2)   # 质量窗口 ×/÷3
N = 20000
random.seed(20261018)
P_NULL = {}

for eps, label in [(0.003, "0.30%（模型水平）"), (0.001, "0.10%（最优命中水平）")]:
    cnt = 0
    for _ in range(N):
        cnt += any_hit(random.uniform(ln_lo, ln_hi), eps)
    p = cnt / N
    pa = sum(2 * eps / abs(L) for L in logs)      # 并集上界（各底相加）
    P_NULL[label] = p
    print(f"  阈值 {label}: 随机目标命中率 = {100*p:.1f}%（MC {N} 次）；相加上界 Σ2ε/|ln b| = {100*pa:.1f}%")
p279 = 2 * 0.003 / abs(math.log(2 / 9))
p38 = 2 * 0.001 / abs(math.log(3 / 8))
print(f"  [口径①] 底与幂都事先定好（不许挑）：单组合 0.3% 级命中概率 = {100*p279:.2f}%（底 2/9 自身）")
print(f"  [口径②] 允许在 {len(logs)} 个底 × 整数幂里挑：0.3% 级命中概率 = {100*P_NULL['0.30%（模型水平）']:.1f}%")
check("0.3% 级命中在随机目标上不是必然（假阳性率 <50%）", P_NULL["0.30%（模型水平）"] < 0.5,
      f"{100*P_NULL['0.30%（模型水平）']:.1f}%")

# ────────────────────────── E  判别力 ──────────────────────────
print()
print("=" * 78)
print("E  判别力：模型 vs 最优竞争候选 vs JUNO")
print("=" * 78)
m_model = A2 * TWO_NINE ** N_STEP * 1e12
m_comp = best["m3_meV"]
dmev = abs(m_model - m_comp)
sigma_juno = 0.0026 * M3_OBS * 1e12      # JUNO 对 m3 的 0.26% 精度（0.13 meV）
print(f"  (2/9)^15  → m3 = {m_model:.3f} meV")
print(f"  {best['base']}^{best['n']} → m3 = {m_comp:.3f} meV")
print(f"  两者差 = {dmev:.3f} meV = {100*dmev/(M3_OBS*1e12):.3f}% ；JUNO σ(m3) ≈ {sigma_juno:.3f} meV")
print(f"  ⟹ 分离显著性 = {dmev/sigma_juno:.2f}σ（<3σ ⟹ 近期不可判别）")
print(f"  对照：稿 §5 的 A/B 两式差 0.43 meV = {0.43/sigma_juno:.1f}σ（可判别）")
check("模型与最优竞争候选在 JUNO 精度下不可判别", dmev / sigma_juno < 3,
      f"{dmev/sigma_juno:.2f}σ")

# ────────────────────────── 汇总 ──────────────────────────
print()
print("=" * 78)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
res = {
    "A_ratio": {"M_e": Me, "m_nu3_over_me": ratio, "m3_meV": m3_pred * 1e12,
                "dev_pct": 100 * (m3_pred / M3_OBS - 1), "masses_meV": masses,
                "Sigma_meV": sum(masses), "mbb_meV": mbb, "R_model": R_model,
                "R_obs": R_obs, "X_model": X_model, "X_obs": DM21 / DM31},
    "B_engine": {"closed_counts": closed, "steady_state": support, "q": q,
                 "block_product": prod, "branch_ratio": ratio_br},
    "C_multiplicity": {"rows": rows, "n_bases": len(rows),
                       "hits_0p5": sum(1 for r in rows if abs(r["dev_pct"]) < 0.5),
                       "hits_0p3": sum(1 for r in rows if abs(r["dev_pct"]) < 0.3),
                       "hits_0p1": sum(1 for r in rows if abs(r["dev_pct"]) < 0.1),
                       "best": best},
    "B2_graphs": B2,
    "B4_carrier_spectrum": B4,
    "D_null": {"p_0p3": P_NULL["0.30%（模型水平）"], "p_0p1": P_NULL["0.10%（最优命中水平）"],
                  "p_single_2_9": p279, "p_single_3_8": p38, "n_bases": len(logs)},
    "E_discrimination": {"m_model": m_model, "m_comp": m_comp, "delta_meV": dmev,
                         "sigma_juno_meV": sigma_juno, "sigma_sep": dmev / sigma_juno},
    "pass": PASS, "fail": FAIL,
}
with open(os.path.join(OUT, "neutrino_scale_audit.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("  → results/neutrino_scale_audit.json")
sys.exit(0 if not FAIL else 1)
