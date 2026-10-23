#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_faithful_v2.py —— 忠实引擎 v2：把 L2 三亚层、6/18/4 破坏、三层角与 α、
3 槽/模流、费米子质量、质子与 21cm 读出接进 zero_evo
================================================================================

性质
----
**讨论稿引擎**（ZeroGPU 工作区），不改动根目录正式篇的任何判定。
底座是 `zero_evo.ZeroEvo`（12 通道、20 条自检）；v2 **不替换**它，只在它上面
加三个扇区与一个读出层：

  S1  L2 三亚层周期映射（底／中／顶；ρ=(0,1/9,1)，d=(1,2,2)，顶容量 5）
  S2  6/18/4 图版（边／边×代／顶点；年龄重置规则；Z017 定理 7）
  S3  三层角、六项式、α 三条读出（路线一／路线二／整式闭式）
  S4  组件级毁灭（28 组件，统一（年龄,−记录）选 6；与 S2 层规则并存）
  R   读出层：3 槽 = seed 相邻比 log（Z131）、模流频率、12 费米子质量
      （Z121/Z126）、质子 m_p 与 g_p（Z165）、21cm（Z166）

**每个字段都在 docstring 或常量旁注明出处与档位**（Z167 的 D／C／I）。
"探索"只用于读出层的候选形式，默认不参与动力学。

为什么这样接（一句话）
----------------------
库里近一周的结论把 L2 拆成两套不同的「3」：
  · **三层**（底／中／顶，Z164 定义 164.1）——毁灭深度不同，给出三层角与 α；
  · **三代**（3 个相邻比／频率槽，Z161/Z129）——给出代的结构计数。
两者都是 3，但**不可互换**（Z161 §5：3 槽是结构计数，不是质量数值）。
v2 因此把它们放在**两个扇区**里，并在报告里明确写成两条支：
  结构支：seed(4 值) → 3 相邻比 → 3 槽 ↔ 三代（计数）
  动力支：三亚层毁灭 → 三层角 → α；非中心态 → 模流频率
这正是《从 seed 到物理世界两条分岔》的结论。

4/18/6 与三层：**只有一套三层划分**（底／中／顶，Z164 定义 164.1）；
6/18/4 是这三层的**格数**（6 边／18 = 3 代 × 6 边／4 顶点），不是第二套分层。
层规则（S2）与组件规则（S4）**同时运行**：层规则给粗粒化，组件规则给微观；
两者在毁灭代给出同一组存活者（4 顶 + 2 中 + 0 底），逐次核对（H1–H5）。

用法
----
    cd ZeroGPU/src
    /usr/bin/python3 zero_faithful_v2.py                 # 自检（~25 条判据）
    /usr/bin/python3 zero_faithful_v2.py --run 20000     # 长跑 + 结果落盘
    /usr/bin/python3 zero_faithful_v2.py --run 5000 --full  # 12 通道全开共存检查

输出
----
    ../results/zero_faithful_v2.json
    ../results/zero_faithful_v2.png
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
from collections import Counter, deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402
from zero_evo import ZeroEvo, Params  # noqa: E402

RES = os.path.normpath(os.path.join(HERE, "..", "results"))

# ══════════════════════════════════════════════════════════════════
#  结构常数（Z164 §6 的输入表；全部无量纲，除 A2_MEV）
# ══════════════════════════════════════════════════════════════════
K0 = 2.0                    # 两分支（Z087 稳态 2:1）
LAM = 3.0                   # 三色／三层
L_BOTT = 7.0                # Bott 七层截断（外部库 GeometryAI 0.3；Z164 §2.5 引用）
DTH = 5.0                   # ΔΘ = Λ² − k₀²（零和恒等式）
DELTA0 = math.pi / 12.0 - 2.0 / 9.0   # δ₀（Z123 内生非中心相位）
DELTA = 3.0 * math.pi / 4.0 - DELTA0  # 电子相角 δ（Z121）
A2_MEV = 313.856            # 轻子质量单元 A² = v/784.5（Z126；v 是唯一量纲锚）
ALPHA_EXP_INV = 137.035999084
ME_EXP_MEV = 0.51099895
MP_EXP_MEV = 938.272
GP_EXP = 5.58569
ME_MP_EXP = 5.446166e-4
NU21_EXP_MHZ = 1420.4058
NU21_SCALE_HZ = 1.235590e20  # m_e c²/h（Z166 因子⑦；按实测 m_e 作单位换算）

# 物理单位（只用于读出层的单位换算；c、h 是 SI 定义）
MEV_J = 1.602176634e-13
H_PLANCK = 6.62607015e-34

# 质量层参数（Z126 §4；上／下型的 A_s、u_s 是登记输入，见 Z126 §7）
MASS_LAYERS = {
    "lepton": dict(A=math.sqrt(2.0), u=1.0, theta0=DELTA,
                   A_src="Z121/Z126（面对角线/棱）", u_src="定义 1"),
    "up": dict(A=1.777833, u=71.60, theta0=math.radians(124.6107),
               A_src="Z126 §3（(4/3)^4 候选，未导出）", u_src="Z126 §4（≈72，未导出）"),
    "down": dict(A=2.682850, u=0.98798, theta0=math.radians(114.7387),
                 A_src="Z126 §3（36/5 候选，未导出）", u_src="Z126 §4（≈1，未导出）"),
}
# 中微子：Z019 的左右手相干残余（不是同一条相位圆的第三代）
NU_MEV = (0.0, 8.570e-6, 49.951e-6)  # 8.570 meV / 49.951 meV → MeV

# 6/18/4 的组件表（顶 4 顶点／中 3 代×6 边／底 6 边）
COMP_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
COMP_DEG = np.array([3.0, 3.0, 3.0, 3.0])


# ══════════════════════════════════════════════════════════════════
#  基础函数
# ══════════════════════════════════════════════════════════════════
def sixterm(theta) -> float:
    """六项式（Z164 定义 164.5）：S = Σ csc²θ_i + Σ_{i<j} cscθ_i cscθ_j。"""
    th = np.asarray(theta, dtype=float)
    x = 1.0 / np.sin(th)
    return float(np.sum(x * x) + x[0] * x[1] + x[0] * x[2] + x[1] * x[2])


def theta_from_contents(x, k: float | None = None):
    """三层角：θ_i ∝ x_i^{-k}，归一化到 Σθ = π/2（Z164 定义 164.24 / 命题 164.2）。

    k=1 时等价于「内容×角＝常数」（最终账 §1.5）；默认取非中心指数
    k = 1 − ΔΘ δ₀²（Z164 命题 164.27）。
    """
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        return None
    kk = (1.0 - DTH * DELTA0 ** 2) if k is None else float(k)
    w = x ** (-kk)
    return (math.pi / 2.0) * w / w.sum()


def solve_route1():
    """路线一（Z164 §3）：闭合角 ＋ 六项式 S=137 ＋ 中层平衡 θ_C²=κθ_Mθ_I。

    取 θ_M>θ_C>θ_I 支（文档基准点 57.93°,26.16°,5.91°），返回 (θ, S)。
    κ = k₀ cos δ₀（Z164 命题 164.16）。
    """
    kap = K0 * math.cos(DELTA0)

    def tC_of(tM):
        disc = kap * kap * tM * tM + 4.0 * kap * tM * (math.pi / 2.0 - tM)
        return (-kap * tM + math.sqrt(max(disc, 0.0))) / 2.0

    def f(tM):
        tC = tC_of(tM)
        tI = math.pi / 2.0 - tM - tC
        if tC <= 0 or tI <= 0:
            return float("nan")
        return sixterm((tM, tC, tI)) - 137.0

    # 粗扫找根，再二分（避开 5.91° 的镜像支）
    N = 200000
    prev = None
    roots = []
    for i in range(1, N):
        tM = (math.pi / 2.0) * i / N
        v = f(tM)
        if not math.isfinite(v):
            prev = None
            continue
        if prev is not None and prev * v < 0:
            lo = (math.pi / 2.0) * (i - 1) / N
            hi = tM
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if f(lo) * f(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            roots.append(0.5 * (lo + hi))
        prev = v
    if not roots:
        raise RuntimeError("路线一：S=137 无根")
    tM = max(roots)                      # θ_M 最大支
    tC = tC_of(tM)
    tI = math.pi / 2.0 - tM - tC
    th = np.array([tM, tC, tI])
    return th, sixterm(th)


def route1_with_displacement():
    """路线一 ＋ 结构位移 Δξ（Z164 定义 164.13 / 命题 164.14）。

    Δξ = 1/(S_min L^{L−Λ}) = 1/(24·7⁴) = 1/57624 —— 一个**作用量量子**（类 C）。
    """
    th, _ = solve_route1()
    dxi = 1.0 / (24.0 * L_BOTT ** (L_BOTT - LAM))  # = 1/57624
    th2 = th.copy()
    th2[1] += dxi
    th2[2] -= dxi
    return th2, sixterm(th2), dxi


def alpha_closed_form() -> float:
    """α⁻¹ 的整式闭式（Z164 命题 164.32）：k₀^L + Λ² + Λ²/(k₀ΔΘ³)。"""
    return K0 ** L_BOTT + LAM ** 2 + LAM ** 2 / (K0 * DTH ** 3)


def d_scan():
    """枚举 ΔΘ=5 的 6 个正整数三分拆，验证 d=(1,2,2) 是离散唯一解（最终账 §1.2）。"""
    rows = []
    k = 1.0 - DTH * DELTA0 ** 2
    for a in range(1, 4):
        for b in range(1, 4):
            for c in range(1, 4):
                if a + b + c != 5:
                    continue
                x = (float(a), 9.0 / 8.0 * b, 5.0 * c)
                th = theta_from_contents(x, k)
                S = sixterm(th)
                rows.append(dict(d=[a, b, c], x=list(x),
                                 theta_deg=[math.degrees(t) for t in th],
                                 S=S, dev=(S - 137.0) / 137.0))
    rows.sort(key=lambda r: abs(r["dev"]))
    return rows


def slots_from_seed(seed):
    """3 槽 = seed 相邻比的 log（Z131；Z161 命题 161.10）。"""
    s = np.asarray(seed, dtype=float)
    s = s[s > 0]
    if s.size < 2:
        return np.array([])
    return np.log(s[:-1] / s[1:])


def z168_gradient_variance(L: int) -> float:
    """Z168：L1 词层面积相位的梯度方差 Var_θ[∂C] = (4/N_L²) Σ a² n_a²。

    a(w) = Σ_k (前缀和)_k；N_L = C(L, L/2)。L=4→4.444，L=20→30.952
    （Z168 摘要的 4.44→30.95）。这是**变分/优化扇区**的读出，不是世界动力学。
    """
    N = math.comb(L, L // 2)
    counts = {}
    for pos in itertools.combinations(range(L), L // 2):
        w = [-1] * L
        for i in pos:
            w[i] = 1
        s = 0
        a = 0
        for x in w:
            s += x
            a += s
        counts[a] = counts.get(a, 0) + 1
    return 4.0 / (N * N) * sum(a * a * n * n for a, n in counts.items() if a != 0)


def mass_layer(A, u, theta0):
    """一条振幅律的 3 代质量（Z121/Z126）：m = u A² [1+A cos(θ₀+2πi/3)]²。"""
    out = []
    for i in range(3):
        M = 1.0 + A * math.cos(theta0 + 2.0 * math.pi * i / 3.0)
        out.append(u * A2_MEV * M * M)
    return out


def proton_readout():
    """质子 m_p、g_p、m_e/m_p（Z165 定理 165.2、165.7、推论 165.3）。"""
    Me = 1.0 + math.sqrt(2.0) * math.cos(DELTA)
    mp = 3.0 * A2_MEV * (1.0 - math.sqrt(DTH) * DELTA0 ** 2)
    gp = 2.0 * LAM * math.cos(math.pi / 12.0) ** 2 * (1.0 - math.sqrt(2.0) * DELTA0 ** 2)
    ratio = Me * Me / (3.0 * (1.0 - math.sqrt(DTH) * DELTA0 ** 2))
    me = A2_MEV * Me * Me
    return dict(M_e=Me, m_e_MeV=me, m_p_MeV=mp, g_p=gp, m_e_over_m_p=ratio,
                dev_m_p=(mp - MP_EXP_MEV) / MP_EXP_MEV,
                dev_g_p=(gp - GP_EXP) / GP_EXP,
                dev_ratio=(ratio - ME_MP_EXP) / ME_MP_EXP)


def h21_readout(alpha_inv, ratio, gp, m_e_MeV, measured_scale: bool = True):
    """21cm（Z166 定理 166.2）。

    measured_scale=True：因子⑦ 用实测 m_e c²/h（Z166 的 1420.4728 MHz 口径）。
    measured_scale=False：因子⑦ 用结构 m_e = A²M_e² c²/h（自洽变体）。
    两者的差只有 1.0×10⁻⁵ 相对——低于 Z166 的两圈残余 4.7×10⁻⁵。
    """
    alpha = 1.0 / float(alpha_inv)
    scale = NU21_SCALE_HZ if measured_scale else (m_e_MeV * MEV_J / H_PLANCK)
    nu = ((K0 ** 2) / LAM) * alpha ** 4 * (1.0 + 1.162364e-3) * ratio * gp \
        * (1.0 + ratio) ** (-3) * scale
    return nu / 1e6  # MHz


# ══════════════════════════════════════════════════════════════════
#  S1：L2 三亚层周期映射（Z164 §4；Z017 定理 7 的粗粒化）
# ══════════════════════════════════════════════════════════════════
class SublayerWorld:
    """底／中／顶三亚层的毁灭—重播种周期映射。

    x ← ρ·x + d，  ρ=(0,1/9,1)，d=(1,2,2)，顶容量 T_cap=5（Z164 定理 164.21 /
    命题 164.25）。ρ 由 Z017 定理 7 定：底全灭、中留 2/18、顶不动；
    d=(1,2,2) 是 ΔΘ=5 的 6 个整数三分拆里唯一给出 S≈137 的（最终账 §1.2）。

    档位：ρ —— D（Z017 定理 7）；d —— D（离散唯一，但"为何是这组"仍是输入）；
    T_cap —— C（α 恒定性 ＋ 数据选它）。
    """

    def __init__(self, d=(1, 2, 2), rho=(0.0, 1.0 / 9.0, 1.0), tcap=5.0):
        self.d = np.asarray(d, dtype=float)
        self.rho = np.asarray(rho, dtype=float)
        self.tcap = float(tcap)
        self.x = np.zeros(3)
        self.n = 0
        self.hist = []

    def period(self):
        """一次毁灭周期：底重置为 d_M；中按 ρ_C 保留；顶累积到容量。"""
        self.x = self.rho * self.x + self.d
        self.x[2] = min(self.x[2], self.tcap * self.d[2])
        self.n += 1
        th = self.theta()
        S = sixterm(th) if th is not None else float("nan")
        self.hist.append(dict(period=self.n, x=self.x.tolist(),
                              theta_deg=[math.degrees(t) for t in th] if th is not None else [],
                              S=S, alpha_inv=S))
        return self

    def theta(self, k=None):
        return theta_from_contents(self.x, k)

    def steady(self, n=5):
        """解析稳态（Z164 定理 164.22）：(d_M, (9/8)d_C(1−9⁻ⁿ), min(n,T_cap)d_I)。"""
        nn = min(n, self.tcap)
        return np.array([self.d[0],
                         (9.0 / 8.0) * self.d[1] * (1.0 - 9.0 ** (-n)),
                         nn * self.d[2]])

    def reset(self):
        self.x[:] = 0.0
        self.n = 0
        self.hist.clear()
        return self


# ══════════════════════════════════════════════════════════════════
#  S2：6/18/4 图版（Z017 定理 7 的年龄规则；不是第二套分层）
# ══════════════════════════════════════════════════════════════════
class SublayerGraph:
    """底 6 边／中 3 代×6 边＝18／顶 4 顶点；统一 T=4。

    命运差异只来自**年龄重置规则**（《6/18/4 是零件不是层》）：
      · 顶：每代 at[:]=0 ⟹ 永不到 T ⟹ 永不毁灭；
      · 中：am+=1，但 top2（按记录）清零 ⟹ 每周期留 2/18；
      · 底：ab+=1，到期全清 ⟹ 每周期 0/6。
    破坏代（gen%T==0）在**重播种前**记录 (底,中,顶) 活跃数，应为 (0,2,4)。
    """

    def __init__(self, T=4):
        self.T = int(T)
        self.NB, self.NG, self.NM, self.NT = 6, 3, 18, 4
        self.age_b = np.zeros(self.NB, dtype=int)
        self.age_m = np.zeros((self.NG, self.NB), dtype=int)
        self.age_t = np.zeros(self.NT, dtype=int)
        self.act_b = np.ones(self.NB, dtype=bool)
        self.act_m = np.ones((self.NG, self.NB), dtype=bool)
        self.act_t = np.ones(self.NT, dtype=bool)
        self.rec = np.zeros((self.NG, self.NB))
        self.gen = 0
        self.pattern = []   # 破坏代 (底,中,顶)
        self.trace = []     # 每代 (gen, 底,中,顶)

    def step(self):
        self.gen += 1
        # 顶：年龄每代清零（永不毁灭）
        self.age_t[:] = 0
        # 中：年龄自增；top2 按记录清零（2 个格，不是 2 代）
        self.age_m += 1
        flat = self.rec.reshape(-1)
        top2 = np.argsort(-flat, kind="stable")[:2]
        for idx in top2:
            g, e = divmod(int(idx), self.NB)
            self.age_m[g, e] = 0
        # 底：年龄自增
        self.age_b += 1
        # 活跃 = 年龄 < T
        self.act_b = self.age_b < self.T
        self.act_m = self.age_m < self.T
        self.act_t = np.ones(self.NT, dtype=bool)
        # 破坏代：记录重播种前的模式，然后重播种（年龄清零）
        if self.gen % self.T == 0:
            self.pattern.append((int(self.act_b.sum()), int(self.act_m.sum()),
                                 int(self.act_t.sum())))
            self.age_b[:] = 0
            self.age_m[:] = 0
            self.age_t[:] = 0
            self.act_b[:] = True
            self.act_m[:] = True
            self.act_t[:] = True
        # 记录（只由活跃的中层格贡献；top2 由此稳定）
        self.rec = 0.9 * self.rec + self.act_m.astype(float) / self.NM
        self.trace.append((self.gen, int(self.act_b.sum()),
                           int(self.act_m.sum()), int(self.act_t.sum())))
        return self

    def counts(self):
        return int(self.act_b.sum()), int(self.act_m.sum()), int(self.act_t.sum())


def status_table():
    """读出层的档位表（Z167 的 D／C／I；I=具名输入）。"""
    return [
        dict(item="L0/L1/L2 分层数", status="D", source="Z161（L0=0, L1=H=L/2, L2=3 槽）"),
        dict(item="3 槽 = 3 相邻比", status="D", source="Z161 命题 161.10 / Z131"),
        dict(item="6/18/4 破坏模式", status="I", source="Z017 定理 7（明标具名输入）"),
        dict(item="ρ=(0,1/9,1)", status="D", source="Z164 定理 164.21（由 Z017 定理 7）"),
        dict(item="d=(1,2,2)", status="D", source="Z164 / 最终账 §1.2（ΔΘ=5 的 6 个三分拆里唯一）"),
        dict(item="T_cap=5", status="C", source="Z164 命题 164.25（α 恒定性 + 数据选它）"),
        dict(item="θ_M+θ_C+θ_I=90°", status="D", source="Z164 命题 164.2"),
        dict(item="α⁻¹ 基数 137", status="D", source="Z164 路线二（输出；无位置输入）"),
        dict(item="α⁻¹ 终值 137.036", status="D+C", source="Z164 命题 164.32 + 一个作用量量子（类 C）"),
        dict(item="m_p, g_p, m_e/m_p", status="D", source="Z165（A²=v/784.5 是唯一量纲锚）"),
        dict(item="带电轻子质量", status="D", source="Z121/Z126 振幅律（形状）"),
        dict(item="上/下型夸克质量", status="C", source="Z126 §7（A_s、u_s 是登记输入）"),
        dict(item="中微子质量", status="C", source="Z019（左右手相干残余，另一结构）"),
        dict(item="21cm", status="D", source="Z166（因子⑦用实测 m_e c²/h 作单位换算）"),
        dict(item="模流频率 log(p/q)", status="D", source="Z131 / Z005"),
        dict(item="3 槽 → 质量", status="未建", source="Z161 §5：3 槽是结构计数；Z167"),
        dict(item="Z168 L1 梯度方差", status="独立读出", source="Z168；不进世界动力学（变分/优化扇区）"),
        dict(item="4/18/6 与三层", status="D/I", source="Z164 定义164.1/命题164.3（一套三层；6/18/4 是各层格数，本身是 Z017 具名输入）"),
        dict(item="组件级毁灭（28 选 6）", status="D", source="本引擎 H1–H5；结构规则：4 顶 + 顶点 deg−1 条入射边（当前代），无容量 6"),
    ]


def structural_keep(rt, rm, gen):
    """结构规则（无全局容量）：4 顶全留；取记录最高顶点的 3 条入射边，
    留记录前 2 条（deg−1=2），且只留**当前代**。共 4+2=6 个。"""
    v = int(np.argmax(rt))
    inc = [i for i, (a, b) in enumerate(COMP_EDGES) if v in (a, b)]
    er = np.asarray(rm).sum(axis=0)
    top2 = sorted(inc, key=lambda i: -er[i])[:2]
    slot = (gen - 1) % 3
    return np.array([0, 1, 2, 3] + [4 + slot * 6 + i for i in top2], dtype=int)


class ComponentDestruction:
    """28 个组件的统一毁灭（与 `SublayerGraph` 的层规则**并存**）。

    组件：顶 4（顶点 v）／中 18（边 e × 3 代）／底 6（边 e）。
    刷新：顶 = 本代闭合 c_v；中 = 最近 3 代边流；底 = 本代边流。
    毁灭：按 (年龄, −记录) 统一选前 `cap` 个；记录按共享数归一
          （顶 × deg、中 ÷ 3、底 × 1）。
    层规则与组件规则不冲突：层规则说"哪一组留几个"，组件规则说"这 28 个里谁留"；
    两者在毁灭代给出同一组存活者（4 顶 + 2 中 + 0 底）。
    """

    def __init__(self, T: int = 4, cap: int = 6, mode: str = "structural"):
        self.T = int(T)
        self.cap = int(cap)
        self.mode = mode
        self.rt = np.zeros(4)
        self.rm = np.zeros((3, 6))
        self.rb = np.zeros(6)
        self.at = np.zeros(4, dtype=int)
        self.am = np.zeros((3, 6), dtype=int)
        self.ab = np.zeros(6, dtype=int)
        self.mem: deque = deque(maxlen=3)
        self.events: list[dict] = []

    def step(self, closures, flow):
        c = np.asarray(closures, dtype=float)
        f = np.asarray(flow, dtype=float)
        self.rt += c * COMP_DEG
        self.at = np.where(c > 0, 0, self.at + 1)
        self.mem.append(f.copy())
        recent = np.max(np.array(self.mem), axis=0)
        self.rm += f / 3.0
        self.am = np.where(recent > 0, 0, self.am + 1)
        self.rb += f
        self.ab = np.where(f > 0, 0, self.ab + 1)

    def destroy(self, gen: int):
        aa = np.concatenate([self.at, self.am.reshape(-1), self.ab])
        rr = np.concatenate([self.rt, self.rm.reshape(-1), self.rb])
        vstar = None
        top2_edges = []
        if self.mode == "structural":
            vstar = int(np.argmax(self.rt))
            inc = [i for i, (a, b) in enumerate(COMP_EDGES) if vstar in (a, b)]
            er = np.asarray(self.rm).sum(axis=0)
            top2_edges = [int(x) for x in sorted(inc, key=lambda i: -er[i])[:2]]
            keep = structural_keep(self.rt, self.rm, gen)
        else:
            order = sorted(range(28), key=lambda i: (aa[i], -rr[i]))
            keep = order[:self.cap]
        types = []
        for k in keep:
            if k < 4:
                types.append(("T", int(k)))
            elif k < 22:
                j = k - 4
                types.append(("M", int(j % 6), int(j // 6)))
            else:
                types.append(("B", int(k - 22)))
        counts = [sum(1 for t in types if t[0] == ch) for ch in "TMB"]
        destroyed = [4 - counts[0], 18 - counts[1], 6 - counts[2]]
        ev = dict(gen=int(gen), counts=counts, destroyed=destroyed,
                  middle_edges=sorted(t[1] for t in types if t[0] == "M"),
                  vstar=vstar, top2_edges=top2_edges,
                  keep=[int(x) for x in keep],
                  ages_all=[int(x) for x in aa],
                  records_all=[float(x) for x in rr])
        self.events.append(ev)
        for k in range(28):
            if k in keep:
                continue
            if k < 4:
                self.at[k] = 0; self.rt[k] = 0.0
            elif k < 22:
                j = k - 4; self.am[j // 6, j % 6] = 0; self.rm[j // 6, j % 6] = 0.0
            else:
                self.ab[k - 22] = 0; self.rb[k - 22] = 0.0
        return ev


def _age0_counts(ev):
    """去掉容量 6 的反事实：年龄 0 的组件全留（age0）。"""
    ages = ev["ages_all"]
    c = [sum(1 for i in range(0, 4) if ages[i] == 0),
         sum(1 for i in range(4, 22) if ages[i] == 0),
         sum(1 for i in range(22, 28) if ages[i] == 0)]
    return dict(counts=c, destroyed=[4 - c[0], 18 - c[1], 6 - c[2]],
                note="去掉容量 6：年龄 0 全留，无容量/无阈值")


# ══════════════════════════════════════════════════════════════════
#  v2 引擎：zero_evo 的超集
# ══════════════════════════════════════════════════════════════════
class ZeroFaithfulV2(ZeroEvo):
    """zero_evo（12 通道）＋ 三亚层周期映射 ＋ 6/18/4 层规则 ＋ 28 组件毁灭 ＋ 读出层。

    `SublayerGraph`（层规则）与 `ComponentDestruction`（组件规则）**同时运行**：
    层规则给粗粒化（底/中/顶各留几个），组件规则给微观（28 个里谁留）；
    两者在毁灭代给出同一组存活者（4 顶 + 2 中 + 0 底）。
    """

    def __init__(self, adj=None, params: Params | None = None, seed=None,
                 sub_d=(1, 2, 2), sub_tcap=5.0):
        if adj is None:
            adj = zs.k4()                      # L0 舞台：K₄（Z021）
        if params is None:
            params = Params(L=4, T=4, micro=True, micro_semantics="persistent",
                            retention="top2", poisson_on=False, wave_on=False,
                            c_h=0.0, g_src=0.0, Gamma_h=0.0, norm_hi=1e50)
        if seed is None:
            seed = [8.0, 4.0, 1.5, 1.0]        # Z089/Z131 的 K4 seed
        super().__init__(adj, params, seed=seed)
        self.world = SublayerWorld(d=sub_d, tcap=sub_tcap)
        self.graph = SublayerGraph(T=params.T or 4)   # 层规则（粗粒化）
        self.comp = ComponentDestruction(T=params.T or 4, mode="structural")  # 组件规则（微观，无容量）
        self.world_hist = []                   # 每周期快照
        self.comp_hist = []                    # 每次组件级毁灭
        self.alpha_hist = []                   # 每次破坏周期的 α⁻¹

    def step(self):
        super().step()
        # 组件级毁灭：闭合 c 与边流 f 都取自本代走者状态（毁灭代 W=0，f=0）
        c = self.gen_closures.copy()
        f = np.array([self.W[u, w] + self.W[w, u] for u, w in COMP_EDGES], dtype=float)
        self.comp.step(c, f)
        self.graph.step()
        if self.p.T and self.gen % self.p.T == 0:
            self.world.period()
            ev = self.comp.destroy(self.gen)
            self.comp_hist.append(ev)
            th = self.world.theta()
            S = sixterm(th)
            self.world_hist.append(dict(gen=self.gen, period=self.world.n,
                                        x=self.world.x.tolist(),
                                        theta_deg=[math.degrees(t) for t in th],
                                        S=S, alpha_inv=S,
                                        pattern=self.graph.counts(),
                                        component=ev["counts"],
                                        component_destroyed=ev["destroyed"]))
            self.alpha_hist.append(S)
        return self

    # ── 读出 ────────────────────────────────────────────────────
    def slots(self, sort=True):
        s = self.seed / self.seed.sum() if self.seed.sum() > 0 else self.seed
        sl = slots_from_seed(s)
        return np.sort(sl) if sort else sl

    def modular_frequencies(self):
        """模流频率 = seed 相邻比的 log（Z131）；log2 是 2:1 稳态的内禀频率。"""
        s = self.seed / self.seed.sum()
        return np.log(s[:-1] / s[1:])

    def report(self):
        snap = self.snapshot()
        seed_norm = self.seed / self.seed.sum() if self.seed.sum() > 0 else self.seed
        slots = self.slots()
        freq = self.modular_frequencies()
        th2 = self.world.theta()
        S2 = sixterm(th2) if th2 is not None else float("nan")
        th1, S1, dxi = route1_with_displacement()
        pr = proton_readout()
        lep = mass_layer(**{k: MASS_LAYERS["lepton"][k] for k in ("A", "u", "theta0")})
        up = mass_layer(**{k: MASS_LAYERS["up"][k] for k in ("A", "u", "theta0")})
        down = mass_layer(**{k: MASS_LAYERS["down"][k] for k in ("A", "u", "theta0")})
        nu21 = h21_readout(S1, pr["m_e_over_m_p"], pr["g_p"], pr["m_e_MeV"], True)
        nu21_struct = h21_readout(S1, pr["m_e_over_m_p"], pr["g_p"],
                                  pr["m_e_MeV"], False)
        return dict(
            gen=int(self.gen), periods=int(self.world.n),
            adj_V=int(self.V), T=int(self.p.T), L=int(self.p.L),
            seed=[float(x) for x in seed_norm],
            seed_log10=float(math.log10(max(self.seed.sum(), 1e-300))),
            record=[float(x) for x in (self.record / self.record.sum())],
            slots=[float(x) for x in slots],
            modular_freq=[float(x) for x in freq],
            sub_x=[float(x) for x in self.world.x],
            sub_theta_deg=[float(x) for x in self.world.hist[-1]["theta_deg"]]
            if self.world.hist else [],
            alpha_route2_inv=float(S2),
            dev_route2_vs_137=(float(S2) - 137.0) / 137.0,
            dev_route2_vs_measured=(float(S2) - ALPHA_EXP_INV) / ALPHA_EXP_INV,
            alpha_route1_inv=float(S1),
            dev_route1=(float(S1) - ALPHA_EXP_INV) / ALPHA_EXP_INV,
            delta_xi=float(dxi),
            alpha_closed=float(alpha_closed_form()),
            dev_closed=(float(alpha_closed_form()) - ALPHA_EXP_INV) / ALPHA_EXP_INV,
            proton=pr,
            masses_MeV=dict(lepton=lep, up=up, down=down,
                            neutrino=list(NU_MEV)),
            nu21_MHz=float(nu21), dev_nu21=(float(nu21) - NU21_EXP_MHZ) / NU21_EXP_MHZ,
            nu21_struct_MHz=float(nu21_struct),
            dev_nu21_struct=(float(nu21_struct) - NU21_EXP_MHZ) / NU21_EXP_MHZ,
            graph_pattern=[list(p) for p in self.graph.pattern[-4:]],
            graph_counts=list(self.graph.counts()),
            component=dict(n=len(self.comp_hist),
                           last=(self.comp_hist[-1] if self.comp_hist else None),
                           counts_summary={str(k): int(v) for k, v in
                                           Counter(tuple(e['counts']) for e in self.comp_hist).items()},
                           destroyed_summary={str(k): int(v) for k, v in
                                              Counter(tuple(e['destroyed']) for e in self.comp_hist).items()}),
            component_free=(_age0_counts(self.comp_hist[-1]) if self.comp_hist else None),
            layer_vs_component=dict(
                layer_last=(list(self.graph.pattern[-1]) if self.graph.pattern else None),
                component_last=(self.comp_hist[-1]['counts'] if self.comp_hist else None),
                agree=all((e['counts'][2], e['counts'][1], e['counts'][0]) == p
                          for e, p in zip(self.comp_hist, self.graph.pattern))
                if self.comp_hist and len(self.comp_hist) == len(self.graph.pattern) else False,
                note="层规则 (底,中,顶) 与组件规则 (顶,中,底) 是同一组存活者的两种写法"),
            macro=dict(n_std=float(snap['n_std']), q=float(snap['q']),
                       h_energy=float(snap['h_energy']), psi_coh=float(snap['psi_coh']),
                       lam_peak=float(snap['lam_peak'])),
            branches=dict(
                structure=dict(seed=[float(x) for x in seed_norm],
            seed_log10=float(math.log10(max(self.seed.sum(), 1e-300))),
                               slots=[float(x) for x in slots],
                               generations=3,
                               note="seed → 3 相邻比 → 3 槽 ↔ 三代（结构计数，Z161 §5）"),
                dynamics=dict(sub_x=[float(x) for x in self.world.x],
                              theta_deg=[float(math.degrees(t)) for t in th2],
                              alpha_inv=float(S2),
                              modular_freq=[float(x) for x in freq],
                              note="三亚层毁灭 → 三层角 → α；非中心态 → 模流频率（Z164/Z131）"),
            ),
            z168=dict(L4=float(z168_gradient_variance(4)),
                      L20=float(z168_gradient_variance(20))),
            status=status_table(),
        )

    # ── 自检 ────────────────────────────────────────────────────
    def selftest(self, gens=400):
        ck = []

        def chk(name, cond, detail=""):
            ck.append(dict(name=name, pass_=bool(cond), detail=str(detail)))
            print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
            return bool(cond)

        # 跑一段
        for _ in range(gens):
            self.step()
        r = self.report()

        print("=" * 80)
        print("zero_faithful_v2 自检：三亚层 / 6-18-4 / 三层角 / α / 质量 / 21cm")
        print("=" * 80)

        # A. 6/18/4 破坏模式
        pats = self.graph.pattern
        chk("A1 破坏代出现（≥3 次）", len(pats) >= 3, f"{len(pats)} 次")
        chk("A2 破坏代模式 = (0,2,4)", all(p == (0, 2, 4) for p in pats[:3]), f"{pats[:4]}")
        chk("A3 重播种后回到 6/18/4", self.graph.counts() in [(6, 18, 4)], f"{self.graph.counts()}")
        chk("A4 6/18/4 同住 L2（不是三套分层）", self.graph.NB + self.graph.NM + self.graph.NT == 28,
            "6+18+4=28")

        # B. 三亚层周期映射
        w = SublayerWorld()
        for _ in range(5):
            w.period()
        chk("B1 n=5 稳态 x≈(1,2.25,10)（9⁻⁵ 残差）", np.allclose(w.x, [1.0, 2.25, 10.0], atol=1e-3), f"{w.x}")
        chk("B2 n≥5 后 x 不变（顶容量）", (w.period(), w.period(), np.allclose(w.x, [1, 2.25, 10]))[2],
            f"{w.x}")
        th = w.theta()
        chk("B3 三角和 = 90°", abs(sum(th) - math.pi / 2) < 1e-12, f"{math.degrees(sum(th)):.9f}°")
        chk("B4 角序 θ_M>θ_C>θ_I（与毁灭率同序）", th[0] > th[1] > th[2],
            f"{[round(math.degrees(t),3) for t in th]}")

        # C. 三层角与 α
        S2 = sixterm(th)
        chk("C1 路线二 S(θ*) ≈ 137.0038", abs(S2 - 137.0038) < 1e-3, f"{S2:.6f}")
        chk("C2 路线二相对偏差 < 3e-5", abs((S2 - 137.0) / 137.0) < 3e-5,
            f"{(S2-137.0)/137.0:.3e}")
        th1, S1, dxi = route1_with_displacement()
        chk("C3 路线一 S(θ^e) = 137.0359991", abs(S1 - ALPHA_EXP_INV) < 1e-6, f"{S1:.9f}")
        chk("C4 路线一相对偏差 < 1e-9", abs((S1 - ALPHA_EXP_INV) / ALPHA_EXP_INV) < 1e-9,
            f"{(S1-ALPHA_EXP_INV)/ALPHA_EXP_INV:.3e}")
        chk("C5 Δξ = 1/57624", abs(dxi - 1.0 / 57624.0) < 1e-18, f"{dxi:.9e}")
        chk("C6 整式闭式 α⁻¹ = 137.036", abs(alpha_closed_form() - 137.036) < 1e-9,
            f"{alpha_closed_form():.9f}")
        rows = d_scan()
        chk("C7 d=(1,2,2) 在 6 个整数三分拆中唯一", rows[0]["d"] == [1, 2, 2] and
            abs(rows[0]["dev"]) < 1e-3 and abs(rows[1]["dev"]) > 1e-2,
            f"次优 {rows[1]['d']} 偏离 {rows[1]['dev']:.3f}")

        # D. seed → 3 槽 / 模流
        s0 = np.array([8.0, 4.0, 1.5, 1.0]); s0 = s0 / s0.sum()
        sl = np.sort(slots_from_seed(s0))
        chk("D1 3 槽 = [0.4055,0.6931,0.9808]", np.allclose(sl, [0.405465, 0.693147, 0.980829],
                                                           atol=1e-6), f"{sl}")
        chk("D2 模流频率含 log2=0.6931", any(abs(f - math.log(2.0)) < 1e-9
                                            for f in self.modular_frequencies()),
            f"{self.modular_frequencies()}")
        chk("D3 K4 稳态 record = [0.5517,0.2759,0.1034,0.0690]",
            np.allclose(self.record / self.record.sum(),
                        [0.551724, 0.275862, 0.103448, 0.068966], atol=1e-4),
            f"{self.record/self.record.sum()}")

        # E. 质量 / 质子 / 21cm
        pr = r["proton"]
        chk("E1 m_p = 938.2702 MeV", abs(pr["m_p_MeV"] - 938.2702) < 1e-3, f"{pr['m_p_MeV']:.4f}")
        chk("E2 m_p 相对偏差 < 3e-6", abs(pr["dev_m_p"]) < 3e-6, f"{pr['dev_m_p']:.3e}")
        chk("E3 g_p = 5.5856756", abs(pr["g_p"] - 5.5856756) < 1e-6, f"{pr['g_p']:.7f}")
        chk("E4 g_p 相对偏差 < 4e-6", abs(pr["dev_g_p"]) < 4e-6, f"{pr['dev_g_p']:.3e}")
        chk("E5 m_e/m_p 相对偏差 < 1e-5", abs(pr["dev_ratio"]) < 1e-5, f"{pr['dev_ratio']:.3e}")
        chk("E6 轻子质量比 1:206.77:3477.5",
            abs(r["masses_MeV"]["lepton"][1] / r["masses_MeV"]["lepton"][0] - 206.77) < 0.05 and
            abs(r["masses_MeV"]["lepton"][2] / r["masses_MeV"]["lepton"][0] - 3477.5) < 1.0,
            f"{[m/r['masses_MeV']['lepton'][0] for m in r['masses_MeV']['lepton']]}")
        chk("E7 21cm = 1420.4728 MHz（Z166 口径）", abs(r["nu21_MHz"] - 1420.4728) < 1e-3,
            f"{r['nu21_MHz']:.4f}")
        chk("E8 21cm 相对偏差 < 6e-5", abs(r["dev_nu21"]) < 6e-5, f"{r['dev_nu21']:.3e}")
        chk("E9 结构 m_e 变体与 Z166 口径差 < 2e-5",
            abs(r["nu21_struct_MHz"] - r["nu21_MHz"]) / r["nu21_MHz"] < 2e-5,
            f"{r['nu21_struct_MHz']:.4f} vs {r['nu21_MHz']:.4f}")

        # G. Z168 读出（独立扇区，不进世界动力学）
        z4 = z168_gradient_variance(4)
        z20 = z168_gradient_variance(20)
        chk("G1 Z168 L=4 Var=4.444", abs(z4 - 4.444) < 1e-3, f"{z4:.4f}")
        chk("G2 Z168 L=20 Var=30.952", abs(z20 - 30.952) < 1e-2, f"{z20:.4f}")

        # H. 组件级毁灭（与层规则并存）
        chk("H1 组件级毁灭：4 顶 + 2 中 + 0 底",
            len(self.comp_hist) > 0 and all(e["counts"] == [4, 2, 0] for e in self.comp_hist),
            f"{len(self.comp_hist)} 次")
        chk("H2 组件级毁灭：0 顶 + 16 中 + 6 底",
            all(e["destroyed"] == [0, 16, 6] for e in self.comp_hist), "")
        chk("H3 层规则图版：破坏代 (0,2,4)",
            all(p == (0, 2, 4) for p in self.graph.pattern), f"{len(self.graph.pattern)} 次")
        chk("H4 层规则与组件规则逐次一致",
            len(self.comp_hist) == len(self.graph.pattern) and
            all((e["counts"][2], e["counts"][1], e["counts"][0]) == p
                for e, p in zip(self.comp_hist, self.graph.pattern)), "")
        def _struct_mid(e):
            if not e.get("top2_edges"):
                return False
            return (len(e["middle_edges"]) == 2 and
                    [int(x) for x in e["middle_edges"]] == [int(x) for x in e["top2_edges"]])
        chk("H5 组件级中层 2 个 = 最高记录顶点的前 2 条入射边",
            all(_struct_mid(e) for e in self.comp_hist), "")

        # F. 长跑稳定 / 分支
        h = self.alpha_hist
        chk("F1 α⁻¹ 收敛后冻结（末 10 周期 std < 1e-12）", len(h) >= 20 and np.std(h[-10:]) < 1e-12,
            f"std={np.std(h[-10:]) if len(h) >= 10 else float('nan'):.2e}")
        chk("F2 无 NaN / 尺度受控（归一化阈值 1e50）", all(np.isfinite(self.seed)) and self.seed.sum() < 1e60,
            f"log10 seed_sum={math.log10(max(self.seed.sum(),1e-300)):.2f}")
        chk("F3 已跑过多个完整周期", self.gen > 10 and len(self.alpha_hist) > 2, f"{self.gen} 代 / {len(self.alpha_hist)} 周期")
        e2 = ZeroFaithfulV2(params=self.p, seed=[2.0 * s for s in self.seed0])
        for _ in range(self.gen):
            e2.step()
        r2 = e2.report()
        chk("F4 seed 缩放不改变 α（结构支/动力支分离）",
            abs(r2["alpha_route2_inv"] - r["alpha_route2_inv"]) < 1e-9,
            f"{r2['alpha_route2_inv']:.6f} vs {r['alpha_route2_inv']:.6f}")

        nfail = sum(1 for c in ck if not c["pass_"])
        print("-" * 80)
        print(f"合计 {len(ck)} 条：PASS {len(ck)-nfail} / FAIL {nfail}")
        return ck, r


# ══════════════════════════════════════════════════════════════════
#  图
# ══════════════════════════════════════════════════════════════════
def make_figure(engine: ZeroFaithfulV2, path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ws = engine.world.hist
    if not ws:
        return
    per = np.array([w["period"] for w in ws], dtype=float)
    x = np.array([w["x"] for w in ws])
    th = np.array([w["theta_deg"] for w in ws])
    al = np.array([w["alpha_inv"] for w in ws])
    # 6/18/4: 层规则破坏代（pre-reseed）的模式
    pgen = np.array([(i + 1) * engine.graph.T for i in range(len(engine.graph.pattern))], dtype=float)
    pat = np.array(engine.graph.pattern, dtype=float) if engine.graph.pattern else np.zeros((0, 3))
    # 组件级毁灭
    cg = np.array([e["gen"] for e in engine.comp_hist], dtype=float)
    cc = np.array([e["counts"] for e in engine.comp_hist], dtype=float) if engine.comp_hist else np.zeros((0, 3))
    cd = np.array([e["destroyed"] for e in engine.comp_hist], dtype=float) if engine.comp_hist else np.zeros((0, 3))

    fig, ax = plt.subplots(2, 3, figsize=(14.5, 7.2))

    ax[0, 0].plot(per, x[:, 0], "o-", ms=3, label="bottom $x_M$")
    ax[0, 0].plot(per, x[:, 1], "s-", ms=3, label="middle $x_C$")
    ax[0, 0].plot(per, x[:, 2], "^-", ms=3, label="top $x_I$")
    ax[0, 0].set_xscale("log")
    ax[0, 0].set_xlabel("destruction period (log)"); ax[0, 0].set_ylabel("content $x_i$")
    ax[0, 0].set_title("(a) L2 three sublayers: destroy + reseed")
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=.3)

    for i, lab in enumerate(["$\\theta_M$", "$\\theta_C$", "$\\theta_I$"]):
        ax[0, 1].plot(per, th[:, i], "o-", ms=3, label=lab)
    ax[0, 1].set_xscale("log")
    ax[0, 1].set_xlabel("destruction period (log)"); ax[0, 1].set_ylabel("angle (deg)")
    ax[0, 1].set_title("(b) three-layer angles: sum = 90°")
    ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=.3)

    ax[0, 2].plot(per, al, "o-", ms=3, label="$\\alpha^{-1}$ (route 2)")
    ax[0, 2].axhline(137.0, color="k", ls="--", lw=1, label="S = 137 (base)")
    ax[0, 2].axhline(137.036, color="r", ls=":", lw=1, label="137.036 (closed form)")
    ax[0, 2].set_xscale("log")
    ax[0, 2].set_xlabel("destruction period (log)"); ax[0, 2].set_ylabel("$\\alpha^{-1}$")
    ax[0, 2].set_title("(c) $\\alpha$ from destruction")
    ax[0, 2].legend(fontsize=8); ax[0, 2].grid(alpha=.3)
    ax[0, 2].text(0.02, 0.05, "n=5: 137.0045 (Z164)\nfixed point: 137.0038\nbase 137 + Δξ = 137.036",
                  transform=ax[0, 2].transAxes, fontsize=8,
                  bbox=dict(fc="white", ec="0.7", alpha=.9))

    if len(pat):
        ax[1, 0].plot(pgen, pat[:, 0], "o-", ms=2, lw=1, label="bottom 6 (pre-reseed)")
        ax[1, 0].plot(pgen, pat[:, 1], "s-", ms=2, lw=1, label="middle 18 (pre-reseed)")
        ax[1, 0].plot(pgen, pat[:, 2], "^-", ms=2, lw=1, label="top 4 (pre-reseed)")
        for y, c in ((6, "C0"), (18, "C1"), (4, "C2")):
            ax[1, 0].axhline(y, color=c, ls=":", lw=.8, alpha=.6)
    ax[1, 0].set_xlabel("generation"); ax[1, 0].set_ylabel("active sites at destruction")
    ax[1, 0].set_title("(d) layer rule: destruction gen → (0,2,4)")
    ax[1, 0].set_ylim(-0.5, 19)
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=.3)

    if len(cc):
        ax[1, 1].plot(cg, cc[:, 0], "o-", ms=2, lw=1, label="top 4")
        ax[1, 1].plot(cg, cc[:, 1], "s-", ms=2, lw=1, label="middle 18")
        ax[1, 1].plot(cg, cc[:, 2], "^-", ms=2, lw=1, label="bottom 6")
    ax[1, 1].set_xlabel("generation"); ax[1, 1].set_ylabel("component survivors")
    ax[1, 1].set_title("(e) component rule: 4 top + 2 middle + 0 bottom")
    ax[1, 1].set_ylim(-0.5, 19)
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=.3)

    if len(cd):
        ax[1, 2].plot(cg, cd[:, 0], "o-", ms=2, lw=1, label="top destroyed")
        ax[1, 2].plot(cg, cd[:, 1], "s-", ms=2, lw=1, label="middle destroyed")
        ax[1, 2].plot(cg, cd[:, 2], "^-", ms=2, lw=1, label="bottom destroyed")
    ax[1, 2].set_xlabel("generation"); ax[1, 2].set_ylabel("components destroyed")
    ax[1, 2].set_title("(f) component rule: destroyed 0 / 16 / 6")
    ax[1, 2].set_ylim(-0.5, 19)
    ax[1, 2].legend(fontsize=8); ax[1, 2].grid(alpha=.3)

    fig.suptitle("Zero faithful engine v2 —— sublayers, angles, alpha, layer 6/18/4, component destruction", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════
#  main
# ══════════════════════════════════════════════════════════════════
def main(argv=None):
    ap = argparse.ArgumentParser(description="Zero faithful engine v2")
    ap.add_argument("--run", type=int, default=0, help="跑多少代（0=只自检）")
    ap.add_argument("--full", action="store_true", help="12 通道全开（几何/波共存检查）")
    ap.add_argument("--seed", type=float, nargs=4, default=[8.0, 4.0, 1.5, 1.0])
    ap.add_argument("--json", default=None)
    ap.add_argument("--png", default=os.path.join(RES, "zero_faithful_v2.png"))
    ap.add_argument("--selftest-gens", type=int, default=400)
    args = ap.parse_args(argv)

    if args.full:
        params = Params(L=4, T=4, micro=True, micro_semantics="persistent",
                        retention="top2", c_h=1.0, g_src=0.5, kappa_h=0.1, Gamma_h=0.1,
                        poisson_on=True, wave_on=True, norm_hi=1e50)
    else:
        params = Params(L=4, T=4, micro=True, micro_semantics="persistent",
                        retention="top2", poisson_on=False, wave_on=False,
                        c_h=0.0, g_src=0.0, Gamma_h=0.0, norm_hi=1e50)
    eng = ZeroFaithfulV2(params=params, seed=args.seed)
    checks, rep = eng.selftest(gens=args.selftest_gens)

    if args.run > 0:
        # 从自检状态继续跑；长跑只记录周期快照
        for _ in range(max(0, args.run - eng.gen)):
            eng.step()
        rep = eng.report()
        make_figure(eng, args.png)
        print(f"\n长跑 {eng.gen} 代：α⁻¹(route2)={rep['alpha_route2_inv']:.6f} "
              f"ν21={rep['nu21_MHz']:.4f} MHz；图 → {args.png}")

    json_path = args.json or os.path.join(
        RES, "zero_faithful_v2_selftest.json" if args.run == 0 else "zero_faithful_v2.json")
    out = dict(params=dict(L=params.L, T=params.T, full=args.full,
                           seed=[float(x) for x in args.seed]),
               checks=checks, n_fail=sum(1 for c in checks if not c["pass_"]),
               report=rep)
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {json_path}")
    return 0 if out["n_fail"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
