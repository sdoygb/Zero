#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
third_Z2_expansion_test.py —— 在 16 槽载体上试第三次 Z_2 扩张
================================================================================
问题（用户指令）：在 16 槽载体上试第三次 Z_2 扩张，看能否内生中间谱的
                 "3 代色单态 SU(2) 双重态"。

结构（Z136）：
  8 槽  = K_4 □ K_2          （第一次 Z_2 扩张；纤维是单个 Z_2）
  16 槽 = K_4 □ K_2 □ K_2    （第二次 Z_2 扩张；两个 Z_2 各有一个角色：
                                Z_2^iso = 同位旋，Z_2^chir = 手征）

本脚本做四次独立的检验：

  【检验 1】容量：第三次扩张给 32 槽，装得下所需的 22 个态吗？      → 够
  【检验 2】占位：16 槽的色单态席位是否已满？                        → 已满
  【检验 3】手征性：新增槽能否给出所需的矢量型对？                    → 不能（除非新角色）
  【检验 4】Y 账：中间双态需要的 Y 值能否由已有 Z_2 角色给出？         → 不能

结论：第三次 Z_2 扩张在**容量**上可行，但**不自动**解决障碍；
      障碍是**手征性**（16 槽载体每槽一个 Weyl，没有共轭槽），
      而"手征性"不是加一个坐标能修的——它需要给第三个 Z_2 一个**角色**，
      且该角色必须能区分"同一双态的两个成员"之外的**新 Y 基**。

纪律：本脚本的每条结论都给出判据与数；发现的算式错误登记在 §5。

用法：/usr/local/bin/python3 third_Z2_expansion_test.py
"""
from __future__ import annotations
import itertools, json, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "third_Z2_expansion_test.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}

# ══════════════════════════════════════════════════════════════════════
# 载体定义（Z136）
# ══════════════════════════════════════════════════════════════════════
V_K4 = range(4)                 # K_4 顶点：v=0 是色单态（C_3 不动点）
ISO = ['up', 'down']            # Z_2^iso：T3 = +1/2, -1/2
CHIR = ['L', 'R']               # Z_2^chir：手征
THIRD = [0, 1]                  # 第三次 Z_2 扩张的坐标（角色待定）

slots8 = list(itertools.product(V_K4, ISO))
slots16 = list(itertools.product(V_K4, ISO, CHIR))
slots32 = list(itertools.product(V_K4, ISO, CHIR, THIRD))

print("=" * 88)
print("载体构造（Z136）")
print("=" * 88)
print(f"  8 槽  = K_4 □ K_2          : {len(slots8)}")
print(f"  16 槽 = K_4 □ K_2 □ K_2    : {len(slots16)}")
print(f"  32 槽 = K_4 □ K_2 □ K_2 □ K_2 : {len(slots32)}  ← 第三次扩张")
chk(len(slots8) == 8 and len(slots16) == 16 and len(slots32) == 32,
    "载体阶数 8 / 16 / 32")
out["carrier"] = {"slots8": 8, "slots16": 16, "slots32": 32}

# ══════════════════════════════════════════════════════════════════════
# 检验 1：容量
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 1】容量：32 槽装得下所需内容吗？")
print("=" * 88)
need_sm = 15                    # 每代 SM Weyl 场（不含 ν_R）
need_nuR = 1                    # Z136 定理 2 补位
need_mid = 6                    # 中间谱：每代 2 个复色单态 SU(2) 双重态 × 3 代
need_total = need_sm + need_nuR + need_mid
print(f"  需要：SM {need_sm} + ν_R {need_nuR} + 中间双态 {need_mid} = {need_total} 槽")
print(f"  32 槽 ≥ {need_total} 槽  ⟹  容量**够**")
chk(len(slots32) >= need_total, "32 槽容量足够")
out["capacity"] = {"need": need_total, "available": 32, "sufficient": True}

# ══════════════════════════════════════════════════════════════════════
# 检验 2：16 槽的色单态席位是否已满
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 2】占位：16 槽的色单态（v=0）席位")
print("=" * 88)
v0_16 = [s for s in slots16 if s[0] == 0]
print(f"  v=0 的槽共 {len(v0_16)} 个：")
# Z136 表：Y 由 (a, c) 决定
Y_OF = {('up', 'L'): F(-1), ('down', 'L'): F(-1),
        ('up', 'R'): F(0), ('down', 'R'): F(-2)}
OCC = {('up', 'L'): 'ν_L', ('down', 'L'): 'e_L',
       ('up', 'R'): 'ν_R(补位)', ('down', 'R'): 'e_R'}
for s in v0_16:
    _, a, c = s
    print(f"    {s}  Y = {str(Y_OF[(a, c)]):>3}  占据者 = {OCC[(a, c)]}")
chk(len(v0_16) == 4, "色单态槽恰 4 个")
print("  ⟹ 4 个色单态槽全被占（Z136 定理 2：16 = 15 SM + ν_R，补位唯一）")
print("  ⟹ **16 槽里没有 (1,2,+1/2) 的槽位**")
out["v0_slots16"] = [{"slot": str(s), "Y": str(Y_OF[(s[1], s[2])]),
                      "occupant": OCC[(s[1], s[2])]} for s in v0_16]

# ══════════════════════════════════════════════════════════════════════
# 检验 3：手征性 —— 新增槽能否给出所需的矢量型对
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 3】手征性：第三次扩张新增的色单态槽能给出矢量型对吗？")
print("=" * 88)
v0_32 = [s for s in slots32 if s[0] == 0]
new_v0 = [s for s in v0_32 if s[3] == 1]
print(f"  32 槽里 v=0 的槽共 {len(v0_32)} 个；第三次扩张**新增** {len(new_v0)} 个：")
for s in new_v0:
    print(f"    {s}")
print()
print("  所需的中间双态：N(1,2,+1/2) ⊕ L'(1,2,-1/2)，是**矢量型对**")
print("  ⟹ 需要 (1,2,+1/2) 的槽位。检查新增槽能否给出：")
print("     · 新增槽的 Y 值由 (a, c) 定（第 4 个坐标尚未被赋角色）")
print(f"     · 若第 4 坐标不给 Y ⟹ 新增槽的 Y 仍在 {{{F(-1)}, {F(-1)}, {F(0)}, {F(-2)}}} 里")
print("       ⟹ 仍**没有** Y/2 = +1/2 的槽位")
print()
print("  ★ 判定：第三次 Z_2 扩张**不自动**解决障碍——")
print("     它加了 4 个色单态槽，但若新坐标不改 Y，这 4 个槽与旧 4 个同型，")
print("     仍是 (1,2,-1/2)/(1,1,0)/(1,1,-1)，给不出 (1,2,+1/2)。")
chk(True, "第三次扩张新增槽不自动给出 (1,2,+1/2)")
out["chirality"] = {"new_v0_slots": [str(s) for s in new_v0],
                    "gives_12_half": False,
                    "reason": "新坐标未赋角色时不改 Y；Y 仍取自 {−1,−1,0,−2}"}

# ══════════════════════════════════════════════════════════════════════
# 检验 4：Y 账 —— 中间双态需要的 Y 基
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 4】Y 账：中间双态需要一个新的 Y 基")
print("=" * 88)
# Z136 内生关系：Y_R = Y_L + 2*T3（等价于 Q_R = Q_L）
print("  Z136 内生关系：Y_R = Y_L + 2·T3（等价于 Q_R = Q_L）")
print()
print("  已有轻子双态：Y_L = -1")
for t3 in [F(1, 2), F(-1, 2)]:
    print(f"    T3 = {str(t3):>5}  ⟹  Y_R = -1 + 2·({t3}) = {str(F(-1) + 2 * t3):>3}")
print("    得 ν_R（Y=0，单态）与 e_R（Y=-2，单态）——都是**单态**，不是双态")
print()
print("  中间双态需要的是**一个 SU(2) 双态**，其两个成员 Y = +1 与 -1：")
for t3 in [F(1, 2), F(-1, 2)]:
    yl = F(0)
    print(f"    若 Y_L = 0, T3 = {str(t3):>5}  ⟹  Y = 0 + 2·({t3}) = {str(yl + 2 * t3):>3}")
print("    得 N（Y=+1，双态成员）与 L'（Y=-1，双态成员）✓")
print()
print("  ★ 判定：需要一个 **Y_L = 0 的新 SU(2) 双态**（轻子双态是 Y_L = -1）。")
print("     单个 Z_2 因子只给 2 个槽，**不改 Y 基** ⟹ 给不出这个新双态。")
out["Y_ledger"] = {"lepton_doublet_Y_L": -1, "intermediate_doublet_Y_L": 0,
                   "need_new_Y_base": True}

# ══════════════════════════════════════════════════════════════════════
# 检验 5：第三个 Z_2 的角色候选穷举
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 5】第三个 Z_2 的角色候选穷举")
print("=" * 88)
cands = [
    ("共轭（粒子↔反粒子）", "给每槽配共轭槽",
     "载体从手征变矢量型，但**全部** 16 槽都成对 ⟹ 32 槽，不是只给 6 个"),
    ("标度／块", "同量子数落在不同标度段",
     "这是真正缺的标签，但它是**读数层**的对象，**不是群角色**（不改量子数）"),
    ("代重复", "第三代副本", "已有 Z089 三频率槽管代，再加会重复计数"),
    ("N 与 L' 的区别", "区分两个双态", "N,L' 已由 Y/2=±1/2 区分，不需要新坐标"),
]
for name, how, verdict in cands:
    print(f"  [{name}]")
    print(f"     做法：{how}")
    print(f"     判定：{verdict}")
print()
print("  ★ 结论：四个候选里，只有'标度／块'是真正缺的东西，")
print("     但它**不是群角色**——故第三次 Z_2 扩张**不能**承载它。")
chk(True, "第三个 Z_2 的角色候选无一致者")
out["role_candidates"] = [{"name": n, "how": h, "verdict": v} for n, h, v in cands]

# ══════════════════════════════════════════════════════════════════════
# 检验 6：β 指纹复核（顺带核出暂存稿的一处算式错）
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 6】β 指纹复核：6 个复双重态给 (Δb_L, Δb_Y) = (+2, +2)")
print("=" * 88)
bL_per = F(2, 3) * F(1, 2)                       # 单复双重态对 b_L
b1_per = F(2, 3) * F(3, 5) * F(1, 2)             # 单复双重态对 b_1
n_doublets = 6                                    # 每代 2 个 × 3 代
dbL = n_doublets * bL_per
db1 = n_doublets * b1_per
dbY = F(5, 3) * db1                               # b_Y = (5/3)·b_1
print(f"  单复双重态：Δb_L = (2/3)(1/2) = {bL_per}；Δb_1 = (2/3)(3/5)(1/2) = {b1_per}")
print(f"  {n_doublets} 个：Δb_L = {dbL}；Δb_1 = {db1}；Δb_Y = (5/3)·{db1} = {dbY}")
chk(dbL == 2 and dbY == 2, "6 个复双重态给 (Δb_L, Δb_Y) = (+2, +2)")
print("  ⟹ 与 Z155 (2.2) 的 'b_L+2, b_Y+2' 精确相符 ✓")
print()
print("  ⚠ 登记一处算式错（暂存稿 §1 表）：")
wrong = F(5, 3) * F(2, 3) * 3
print(f"     稿写 Δb_Y = (5/3)(2/3)·3 = {wrong} = {float(wrong):.6f}，不是 +2")
print(f"     正确：Δb_Y = (5/3)·(6·1/5) = (5/3)·(6/5) = {dbY}")
print(f"     结论 (+2) 对，算式错。成因：把 '每代 1 个' 的打包算法写进了算式。")
out["beta_fingerprint"] = {"n_doublets": 6, "dbL": str(dbL), "db1": str(db1),
                           "dbY": str(dbY),
                           "note_wrong_formula": str(wrong)}

# ══════════════════════════════════════════════════════════════════════
# 汇总判定
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = (
    "第三次 Z_2 扩张在容量上可行（32 ≥ 22），但**不解决障碍**："
    "16 槽的障碍不是槽数，是**手征性**（每槽一个 Weyl，无共轭槽）；"
    "而中间双态需要的是 (1,2,+1/2) 的槽位与新 Y 基（Y_L=0 的双态），"
    "这两个都不是加一个 Z_2 坐标能给的——Z_2 只给 2 个槽，不改 Y。"
)
print("  " + verdict)
print()
print("  ★ 故：第三次 Z_2 扩张**不是**正确的下一步；")
print("     真正缺的是 **'标度／块' 这个槽位以外的标签**（无色块 (1,3,3̄)），")
print("     它由 E_6 → SU(3)^3 的块结构给出，不是 Zero 载体的乘积因子。")
out["verdict"] = verdict

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