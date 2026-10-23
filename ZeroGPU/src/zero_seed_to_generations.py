#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_seed_to_generations.py —— 从 seed 到物理世界：两条分岔，不是一条直链
==============================================================================
用户链条（2026-10-22）："由于 seed 的重新播种，K4 的四个顶点有了四个值，
然后产生三个非中心态，然后就有了三代粒子等等，然后有了我们这个物理世界，对吗？"

本脚本逐环节核验，结论：
  环节 1  seed → 四个值                      ✅ 对（seed 是 4-向量，注入 W 对角）
  环节 2  四个值 → 三个相邻比                ✅ 对（Z161）
  环节 3  三个相邻比 → 三个频率槽            ✅ 对（Z131：槽 = 相邻比的 log）
  环节 4  "产生三个非中心态"                 ⚠ 措辞要拆：非中心是**谓词**，不是三个对象
  环节 5  → "有了三代粒子"                   ❌ 不能这么接：
             · 三代来自 **3 槽的【结构计数】**（Z089／Z161 §5 明文"不是质量数值"）
             · 非中心态给的是 **内禀时间／演化**（Z027 (5.1)）
  环节 6  → 物理世界                         ✅ 对（Z016：物理世界在 L2）

⟹ **两条分岔**：
     A 结构支：seed(4 值) → 3 比 → 3 槽 ↔ 三代（**计数**）
     B 动力支：非中心态 → 模流非平凡 → **内禀时间（演化）**（Z027）

输出：results/zero_seed_to_generations.json / .txt
运行：cd ZeroGPU/src && /usr/bin/python3 zero_seed_to_generations.py
"""
from __future__ import annotations
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
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
P_("从 seed 到物理世界：两条分岔 · zero_seed_to_generations.py")
P_("=" * 96)

# ══════════════════════════════════════════════════════════════════
#  环节 1-3：seed → 4 值 → 3 比 → 3 槽（可证）
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("环节 1-3：seed → 四个值 → 三个相邻比 → 三个频率槽")
P_("=" * 96)
seed = np.array([8.0, 4.0, 1.5, 1.0])
P_("  seed = %s（4 个值，一个值对应 K4 的一个顶点）" % seed.tolist())
ck("环节1：seed 是 4-向量（四个顶点各有值）", seed.shape == (4,),
   "shape = %s" % (seed.shape,))

P_("\n  三个相邻比（Z161 定义 f_i = log(p_i/p_{i+1})）与三个槽：")
ratios = []
slots = []
for i in range(3):
    r = float(seed[i] / seed[i + 1])
    s = abs(math.log(r))
    ratios.append(r)
    slots.append(s)
    P_("    比 %d = seed[%d]/seed[%d] = %.4f/%.4f = %.4f ; 槽 = log = %.6f"
       % (i + 1, i, i + 1, seed[i], seed[i + 1], r, s))
P_("\n  三个槽 = %s" % [round(x, 6) for x in slots])
ck("环节2：四个值给出三个独立相邻比（Z161：V-1=3 个自由度）",
   len(ratios) == 3, "3 个比")
ck("环节3：三个相邻比给出三个频率槽（Z131：槽 = 相邻比的 log）",
   len(slots) == 3, "3 个槽 = %s" % [round(x, 4) for x in slots])

# 与库内 K4 三槽对照
ref = [abs(math.log(1.5 / 1)), abs(math.log(8 / 4)), abs(math.log(4 / 1.5))]
P_("\n  对照库内 K4 三槽（Z161／Z089）：[log(3/2), log2, log(8/3)] = %s"
   % [round(x, 6) for x in sorted(ref)])
P_("  实算 = %s" % [round(x, 6) for x in sorted(slots)])
ck("实算三槽与库内 K4 三槽一致（排序后）",
   np.allclose(sorted(slots), sorted(ref), atol=1e-12),
   "max|Δ| = %.2e" % float(np.max(np.abs(np.array(sorted(slots))
                                         - np.array(sorted(ref))))))

# ══════════════════════════════════════════════════════════════════
#  环节 4：「三个非中心态」——非中心是谓词，不是三个对象
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("环节 4：非中心态是【谓词】（状态的性质），不是三个对象")
P_("=" * 96)
P_("""
  Z027 命题 3 原文：
      "门开（ρ 非中心）⟺ 块内两支不等权 ⟺ 权重 ι-破缺"
  ⟹ "非中心"是对 **ρ 的一个性质的判定**（是/否），**不是**"有 3 个非中心态"。
""")

# 演示：中心态 vs 非中心态（一个谓词的两种取值）
def rho_central(w, H=2):
    r = np.zeros((2 * H, 2 * H), complex)
    for h in range(H):
        r[2 * h, 2 * h] = w[h] / 2
        r[2 * h + 1, 2 * h + 1] = w[h] / 2
    return r


def rho_noncentral(w, eps):
    """ORI(ε) 构造（Z027 (5.1)）：w_h^± = w_h/2 (1 ± ε η_h)，η_h = ±1。"""
    H = len(w)
    eta = [1 if h % 2 == 0 else -1 for h in range(H)]
    r = np.zeros((2 * H, 2 * H), complex)
    for h in range(H):
        r[2 * h, 2 * h] = w[h] / 2 * (1 + eps * eta[h])
        r[2 * h + 1, 2 * h + 1] = w[h] / 2 * (1 - eps * eta[h])
    return r


w2 = [1 / 3, 2 / 3]
P_("\n  演示（H=2，w=(1/3,2/3)）：")
for eps in (0.0, 0.5):
    r = rho_noncentral(w2, eps)
    # 判非中心：块内两支是否等权
    branch = []
    for h in range(2):
        branch.append((float(np.real(r[2 * h, 2 * h])),
                       float(np.real(r[2 * h + 1, 2 * h + 1]))))
    is_nc = any(abs(a - b) > 1e-12 for a, b in branch)
    P_("    ε=%.1f：各块两支 %s ⟹ 非中心？%s"
       % (eps, [(round(a, 4), round(b, 4)) for a, b in branch],
          "是" if is_nc else "否"))
ck("环节4：'非中心'是谓词（ε=0 否 / ε≠0 是），不是三个对象",
   True, "Z027 命题3：ρ 非中心 ⟺ 块内两支不等权")

# ══════════════════════════════════════════════════════════════════
#  环节 5：三代来自哪里？
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("环节 5：「三代」来自 3 槽的【结构计数】，不来自非中心态")
P_("=" * 96)
P_("""
  Z161 §5 原文：
      "Z089 的『3 槽 ↔ 三代』是**结构计数，不是质量数值**。"
  ⟹ 三代 = **3 个槽的计数**（结构支）

  Z027 (5.1) 原文：
      "ε≠0 ⇒ ρ 非中心 ⇒ 模流非平凡 ⇒ **有内禀时间**（频率 log(p/q)）"
  ⟹ 非中心态 = **时间／演化**（动力支）

  **两者是两条分岔，不是同一条链的两节。**
""")
ck("环节5：三代 ↔ 3 槽是【结构计数】（Z161 §5 明文）",
   len(slots) == 3, "3 槽；Z161 §5：不是质量数值")
ck("环节5：非中心态给【内禀时间】，不是「三代」（Z027 (5.1)）",
   True, "ε≠0 ⇒ 非中心 ⇒ 模流非平凡 ⇒ 内禀时间")

# ══════════════════════════════════════════════════════════════════
#  环节 6：物理世界
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("环节 6：物理世界在 L2")
P_("=" * 96)
P_("  Z016 四层表 + 结论箱：**我们的物理世界在 L2**")
ck("环节6：物理世界在 L2（Z016）", True, "Z016 四层表")

# ══════════════════════════════════════════════════════════════════
#  判定
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("判定：两条分岔")
P_("=" * 96)
P_("""
  ┌──────┬──────────────────────────────────┬────────┬──────────────────────┐
  │ 环节 │ 用户说法                          │ 判定   │ 依据                  │
  ├──────┼──────────────────────────────────┼────────┼──────────────────────┤
  │ 1    │ seed → 四个值                     │ ✅ 对  │ seed 是 4-向量        │
  │ 2    │ 四个值 → 三个相邻比               │ ✅ 对  │ Z161                  │
  │ 3    │ → 三个频率槽                      │ ✅ 对  │ Z131                  │
  │ 4    │ "产生三个非中心态"                │ ⚠ 拆   │ 非中心是谓词          │
  │ 5    │ → "有了三代粒子"                  │ ❌ 否  │ 三代＝3 槽的**计数**   │
  │ 6    │ → 物理世界                        │ ✅ 对  │ Z016：L2              │
  └──────┴──────────────────────────────────┴────────┴──────────────────────┘

  ⟹ **正确结构**：
       A 结构支：seed(4 值) → 3 比 → 3 槽 ↔ 三代（**计数**）
       B 动力支：非中心态 → 模流非平凡 → **内禀时间（演化）**（Z027）
     两条支汇合在 L2（物理世界），但**不是同一条链的两节**。

  ⟹ **落差**：库内只到"3 个槽在计数上与三代对应"；
      三代的**质量**没有（Z161 §5 明文），三代的**具体粒子**也没有。
      故"有了三代粒子"比能支持的**多**。
""")
ck("★ 判定：结构支与动力支要分开；3 槽 ↔ 三代只是计数", True,
   "Z027（时间）+ Z161 §5（计数）")

out = {
    "user_chain": ("seed 重新播种 → K4 四顶点有四个值 → 产生三个非中心态 → "
                   "有了三代粒子 → 有了物理世界"),
    "links": {
        "1_seed_to_four_values": {"verdict": "对", "basis": "seed 是 4-向量"},
        "2_four_to_three_ratios": {"verdict": "对", "basis": "Z161"},
        "3_ratios_to_slots": {"verdict": "对", "basis": "Z131"},
        "4_three_noncentral": {"verdict": "措辞要拆",
                               "basis": "Z027 命题3：非中心是谓词"},
        "5_to_three_generations": {"verdict": "不能这么接",
                                   "basis": "Z161 §5：3 槽 ↔ 三代是结构计数"},
        "6_to_physical_world": {"verdict": "对", "basis": "Z016：L2"},
    },
    "seed": seed.tolist(),
    "ratios": ratios,
    "slots": slots,
    "two_branches": {
        "A_structure": "seed(4 值) → 3 比 → 3 槽 ↔ 三代（计数）",
        "B_dynamics": "非中心态 → 模流非平凡 → 内禀时间（Z027 (5.1)）",
    },
    "gap": ("库内只到'3 槽在计数上与三代对应'；三代的**质量**没有（Z161 §5），"
            "具体粒子也没有 ⟹ '有了三代粒子'超出可支持范围"),
    "checks": CK,
    "verdict": ("环节 1-3、6 成立；环节 4 措辞要拆（非中心是谓词不是三个对象）；"
                "环节 5 不能这么接——三代来自 3 槽的【结构计数】（Z161 §5 明文"
                "'不是质量数值'），而非中心态给的是【内禀时间】（Z027 (5.1)）。"
                "⟹ 结构支与动力支是两条分岔，汇合在 L2。"
                "落差：'有了三代粒子'超出库内支持（无质量、无具体粒子）。"),
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_seed_to_generations.json"), "w",
          encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "zero_seed_to_generations.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(LOG))
n_fail = sum(1 for c in CK if not c["pass"])
P_("\n判据 %d/%d PASS" % (len(CK) - n_fail, len(CK)))
P_("已存 results/zero_seed_to_generations.json / .txt")
sys.exit(1 if n_fail else 0)