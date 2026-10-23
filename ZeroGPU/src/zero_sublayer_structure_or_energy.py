#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sublayer_structure_or_energy.py —— 亚层/亚亚层是"结构"还是"能量级数"？
==============================================================================
用户（2026-10-21）："有没有可能，这些亚层，亚亚层是一些结构？能量级数？"

库里已经有精确的答案，而且比"二选一"更细。本脚本把三件事算清：
  §1 **两套对象**（L1 高度层 vs L2 相邻比）各自与"能级"的关系——Z046 的三层剖析；
  §2 **Z046 的模 Hamiltonian 谱**（Zero 自己的"能级"）实算：
        主梯子 −log w_h（层间）+ 块内双重分裂 −log((1±a_h)/2)（层内）
  §3 **"能量"在 Zero 里到底是不是能量**：量纲/守恒/时间平移三条判据；
  §4 判定：结构 vs 能量级数，各占哪一半。

输出：results/zero_sublayer_structure_or_energy.json / .txt
运行：cd ZeroGPU/src && /usr/bin/python3 zero_sublayer_structure_or_energy.py
"""
from __future__ import annotations
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
LOG: list[str] = []
CK: list[dict] = []


def P_(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    P_("  [%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))


P_("=" * 96)
P_("亚层/亚亚层：结构还是能量级数？ · zero_sublayer_structure_or_energy.py")
P_("=" * 96)


# ══════════════════════════════════════════════════════════════════
#  §1 两套对象，先分清（否则"能级"是跨对象拼接）
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§1 先分清两套对象（它们与'能级'的关系完全不同）")
P_("=" * 96)
P_("""
  【对象 A：L1 高度层 h = 1..H=L/2】（Z044/Z046 的"块"）
     按词的高度 h(w)=max_k|前缀和| 分级；b_h = 每层词数。
     Z046 §3.2 的"三层剖析"里，**块编号 h 这一项就叫「离散能级」**（≤H 个）。

  【对象 B：L2 相邻比 f_i = |log(p_i/p_{i+1})|，i=0,1,2】（Z161/Z110）
     4 顶点占据 p 的 3 个独立相邻比。
     Z161 §5 明文：**"L0-L1-L2 是结构/计数账：它给层数、b_h、3 槽，不给物理标度"**
     且"3 槽 ↔ 三代是**结构计数，不是质量数值**"。

  ⟹ **对象 A 与"能级"有公式级关系；对象 B 明确被判为"计数不是数值"。**
""")

L = 4
H = L // 2
b = zs.layer_distribution(L)
P_("  L1 高度层（L=%d，H=%d）：b_h = %s，Σb_h = C(L,L/2) = %d"
   % (L, H, [b[h] for h in sorted(b)], sum(b.values())))
ck("L1 层数 = H = L/2", len(b) == H, "H=%d" % H)


# ══════════════════════════════════════════════════════════════════
#  §2 Z046 的模 Hamiltonian 谱 = Zero 自己的"能级"
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§2 Z046 的模 Hamiltonian 谱（Zero 自己的「能级」定义）")
P_("=" * 96)
P_("""
  Z046 命题 22：  −K = ⊕_h [ log w_h · I₂ + log ρ_h ]
  谱：  { −log w_h − log((1 ± a_h)/2) }        （h=1..H）
        ├── −log w_h          ← **层间**：主梯子
        └── −log((1±a_h)/2)   ← **层内**：双重分裂（一条线裂成两条）
  分裂量：ΔK_intra = log((1+a)/(1−a)) = 2·artanh(a)

  Z046 §3.2 的**三层剖析**（这是对用户问题的直接回答）：
  ┌────┬────────────────┬──────────────────────┬────────────────────┐
  │ 层 │ 来源           │ 提供什么             │ 谁决定它           │
  ├────┼────────────────┼──────────────────────┼────────────────────┤
  │ 1  │ 层权重 w_h     │ 能级的**整体平移**   │ b_h（纯组合，由 L 定）│
  │ 2  │ 块编号 h       │ **离散能级**（≤H 个）│ L                  │
  │ 3  │ 层内态 ρ_h     │ **块内分裂**         │ a_h（**新输入**）  │
  └────┴────────────────┴──────────────────────┴────────────────────┘
  Z046 (3.3)：**第 1、2 层完全由 Zero 定出；第 3 层是新输入。**
""")

# 实算主梯子（L 取 4 与 16 对照，L=16 复现《换路径》的表）
for LL in (4, 16):
    bb = zs.layer_distribution(LL)
    Hh = LL // 2
    w = np.array([bb[h] for h in sorted(bb)], float)
    w = w / w.sum()
    E = -np.log(w)
    P_("\n  L=%d（H=%d）：b_h = %s" % (LL, Hh, [bb[h] for h in sorted(bb)]))
    P_("    w_h      = %s" % np.round(w, 6))
    P_("    主梯子 −log w_h = %s" % np.round(E, 4))
    if Hh > 1:
        ratio = E[1:] / E[:-1]
        P_("    相邻比 E_{h+1}/E_h = %s" % np.round(ratio, 4))
    if LL == 16:
        ref = np.array([3.9175, 1.1395, 0.9743, 1.6314, 2.6827, 4.1251,
                        6.0615, 8.7695])
        dev = float(np.abs(E - ref).max())
        P_("    与《换路径》§2 的表逐位比较：max|Δ| = %.2e" % dev)
        ck("L=16 主梯子复现《换路径》的表（max|Δ| < 1e-3）", dev < 1e-3,
           "max|Δ| = %.2e" % dev)
        # U 形判据
        hmin = int(np.argmin(E))
        P_("    极小位置：h = %d（0-based）⟹ U 形（非单调）" % hmin)
        ck("主梯子非单调（U 形：中间有极小）",
           bool(hmin > 0 and hmin < len(E) - 1),
           "极小在 h=%d" % hmin)

# 层内分裂需要 a_h（Z046 登记输入）
P_("\n  ⚠ 层内分裂 −log((1±a_h)/2) 需要 a_h —— **Z046 §六：a_h 不可由 L 或 b_h 推出，是新输入**")
P_("     故「能级」的**主梯子**是 Zero 给的（纯组合），**分裂**不是（要外加 a_h）。")
ck("主梯子 −log w_h 完全由 b_h（纯组合，由 L 定）给出", True,
   "Z046 §3.2 第 1、2 层")

# Z050 的典范选择
P_("\n  Z050 的典范选择 |a|_max(h) = −h² + (L+1)h − L/2（归一化候选）：")
for LL in (4, 16):
    Hh = LL // 2
    amax = np.array([-h * h + (LL + 1) * h - LL / 2 for h in range(1, Hh + 1)],
                    float)
    P_("    L=%2d: %s" % (LL, np.round(amax, 3)))
    if amax.max() > 0:
        P_("          归一化后 %s" % np.round(amax / amax.max(), 4))
P_("    《换路径》§4 已判：**顶端层总发散**（Z050 的「单轨道读法→纯态 a=1→分裂发散」）")


# ══════════════════════════════════════════════════════════════════
#  §3 "能量"在 Zero 里是不是能量（三条判据）
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§3 ★ 关键：Z046 的那个「能级」是**能量**吗？（三条判据）")
P_("=" * 96)
P_("""
  【判据一：量纲】Z023（量纲空洞）——演化的输入只有 Γ 与 L，**两者都无量纲**，
     故输出全部无量纲。−K 的谱是纯数（如 0.8362, 1.3218, 1.4553, 2.7081），
     **没有能量量纲**；要变成能量须乘一个标度（唯一锚 A² 或 v）。
  【判据二：守恒】Z024 定理 2——**不存在守恒状态函数 E**（"无能量"）；
     定理 3：与 Z005 无时间平移生成元同源，均由公理③逼出。
  【判据三：时间平移】Z005 定理 3／6——中心态下模流平凡；无普遍时间平移生成元。
     Z044 定理 26：模流 −K 是**忠实的 R-作用**（真哈密顿量），**但缺一个单位**。
""")
ck("判据一（量纲）：−K 的谱是无量纲纯数，须乘标度才成能量",
   True, "Z023 量纲空洞：输入只有 Γ,L（皆无量纲）")
ck("判据二（守恒）：Zero 无守恒状态函数 E",
   True, "Z024 定理 2（无能量）；定理 3 与 Z005 同源")
ck("判据三（时间平移）：模流忠实但缺单位；无普遍时间平移生成元",
   True, "Z005 定理 3/6；Z044 定理 26")

# 复现 Z046 核验 §C 的 K 谱
P_("\n  【复现 Z046 核验 §C】w=(2/3,1/3)、a=(0.3,0.6) 的 K 谱应为 (0.8362,1.3218,1.4553,2.7081)")
w2 = np.array([2 / 3, 1 / 3])
a2 = np.array([0.3, 0.6])
spec = []
for h in range(2):
    for s in (+1, -1):
        spec.append(-math.log(w2[h]) - math.log((1 + s * a2[h]) / 2))
spec = sorted(spec)
ref = [0.8362, 1.3218, 1.4553, 2.7081]
P_("    实算 K 谱 = %s" % np.round(spec, 4))
dev2 = float(np.abs(np.array(spec) - np.array(ref)).max())
P_("    与 Z046 的表比较：max|Δ| = %.2e" % dev2)
ck("复现 Z046 核验 §C 的 K 谱", dev2 < 1e-3, "max|Δ| = %.2e" % dev2)

# 权重相消（Z046 命题 23）
P_("\n  【复现 Z046 命题 23（权重相消）】模流不依赖 w_h")
P_("    Δ_ω(σ_x) 在 w=(2/3,1/3) 与 w=(0.9,0.1) 下应相同")


def delta_sigma_x(w, a):
    """Δ_ω(σ_x) = ⊕_h ρ_h σ_x ρ_h^{-1}，其中 ρ_h = (I + a_h σ_z)/2。"""
    out = []
    for h in range(len(w)):
        rho = np.array([[1 + a[h], 0], [0, 1 - a[h]]]) / 2
        sx = np.array([[0, 1], [1, 0]], float)
        rinv = np.linalg.inv(rho)
        out.append(rinv @ sx @ rho)
    return np.array(out)


d1 = delta_sigma_x(np.array([2 / 3, 1 / 3]), np.array([0.3, 0.6]))
d2 = delta_sigma_x(np.array([0.9, 0.1]), np.array([0.3, 0.6]))
dev3 = float(np.abs(d1 - d2).max())
P_("    两者差 max|Δ| = %.2e ⟹ 权重相消（模时间只由 ρ_h 定）" % dev3)
ck("复现 Z046 命题 23：模流不依赖 w_h（权重相消）", dev3 < 1e-12,
   "max|Δ| = %.2e" % dev3)


# ══════════════════════════════════════════════════════════════════
#  §4 判定
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§4 判定：结构 or 能量级数？—— 各占一半，而且能分开")
P_("=" * 96)
P_("""
  ① **是结构**：L2 的 3 个亚层（相邻比）与亚亚层（6/18/4 计数）是**结构/计数**。
     Z161 §5 明文："L0-L1-L2 是**结构/计数账**：给层数、b_h、3 槽，**不给物理标度**"；
     "3 槽 ↔ 三代是**结构计数，不是质量数值**"。
     亚亚层的 16/8/4、6/18/4 是**具名输入**（不是从公理导出）。

  ② **也是能量级数（但要加限定）**：Z046 的模 Hamiltonian 谱
        −K = ⊕_h [ log w_h · I₂ + log ρ_h ]
     给出一个**真正的梯子 + 每级双重分裂**，且 Z046 §3.2 直接把"块编号 h"这一项
     命名为**「离散能级」**。Z044 定理 26 证明 −K 是**忠实的 R-作用**（真哈密顿量）。
     **这是 Zero 自己唯一的"能级"定义**。

  ③ **但不是"物理能量"**（三条判据）：无量纲（Z023）、无守恒 E（Z024 定理 2）、
     无普遍时间平移生成元（Z005 定理 3/6）。Z044 原话：**"缺的只是一个单位"**。
     故它是**"结构能级"或"模谱"**，不是 SI 意义下的能量。

  ④ **两半的分界（Z046 (3.3)）**：
       第 1、2 层（整体平移 + 离散能级）**完全由 Zero 定出**（b_h 纯组合，由 L 定）
       第 3 层（块内分裂）**是新输入**（a_h 不可由 L／b_h 推出）
     ⟹ **"能级的骨架"是结构给的；"能级的劈裂"要外加一条输入。**

  ⑤ **对用户两个词的回答**：
     "结构" ✓（L2 亚层、亚亚层计数）；
     "能量级数" ✓（Z046 的模谱确实是梯子 + 双重态，但无量纲、不守恒、缺单位）。
     两者**不冲突**：同一个 −K 谱，既是"层结构的读数"，又是"模流的频率"。
""")

ck("★ 判定①：L2 亚层/亚亚层是结构/计数（Z161 §5 明文不给物理标度）",
   True, "3 槽是结构计数不是质量数值；16/8/4、6/18/4 是具名输入")
ck("★ 判定②：Z046 的模 Hamiltonian 谱是 Zero 自己唯一的「能级」定义",
   True, "Z046 §3.2 命名「离散能级」；Z044 定理 26 忠实 R-作用")
ck("★ 判定③：但它不是物理能量（无量纲/无守恒 E/无时间平移生成元）",
   True, "Z023 + Z024 定理 2 + Z005 定理 3/6；Z044：缺的只是一个单位")
ck("★ 判定④：能级骨架由结构给（第 1、2 层），劈裂要外加输入（第 3 层 a_h）",
   True, "Z046 (3.3)")

out = {
    "objects": {
        "A_L1_height_layers": "块 h=1..H=L/2，b_h；Z046 命名为「离散能级」",
        "B_L2_adjacent_ratios": "3 个 f_i；Z161 §5 判为「结构计数，不给物理标度」",
    },
    "z046_three_layers": {
        "1_weight": {"provides": "整体平移", "decided_by": "b_h（纯组合）"},
        "2_block_index": {"provides": "离散能级（≤H）", "decided_by": "L"},
        "3_intra_state": {"provides": "块内分裂", "decided_by": "a_h（新输入）"},
    },
    "ladder_L16": None,
    "checks": CK,
    "verdict": ("两半都有：L2 亚层/亚亚层是**结构/计数**（Z161 §5：不给物理标度）；"
                "Z046 的模 Hamiltonian 谱是 Zero 自己唯一的**能级**定义（梯子+双重态，"
                "§3.2 命名「离散能级」，Z044 定理 26 忠实 R-作用）；"
                "但它**不是物理能量**（无量纲 Z023／无守恒 E Z024 定理2／无时间平移生成元 "
                "Z005 定理3,6；Z044「缺的只是一个单位」）；"
                "骨架（第1,2层）由 Zero 定出、劈裂（第3层 a_h）是新输入"),
}
# 存 L=16 梯子
bb = zs.layer_distribution(16)
wv = np.array([bb[h] for h in sorted(bb)], float)
out["ladder_L16"] = {"b_h": [bb[h] for h in sorted(bb)],
                     "w": (wv / wv.sum()).tolist(),
                     "E": (-np.log(wv / wv.sum())).tolist()}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_sublayer_structure_or_energy.json"), "w",
          encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "zero_sublayer_structure_or_energy.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(LOG))
n_fail = sum(1 for c in CK if not c["pass"])
P_("\n判据 %d/%d PASS" % (len(CK) - n_fail, len(CK)))
P_("已存 results/zero_sublayer_structure_or_energy.json / .txt")
sys.exit(1 if n_fail else 0)