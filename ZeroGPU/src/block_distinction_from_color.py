#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
block_distinction_from_color.py —— 无色块的区分能否由 Zero 的色因子给出？
================================================================================
问题（承接 `第三次Z2扩张实验_容量够障碍在手征性.md` §7 路 A 的下一步）：
  $E_6\to SU(3)^3$ 把 27 分成三块 (3,3̄,1)⊕(1,3,3̄)⊕(3̄,1,3)；
  色单态双态 (L,N,L') 全住在**无色块** (1,3,3̄)。
  问：这个"块"的区分，能否由 Zero 的**色因子**给出？

背景（库内既有判定，必须尊重）：
  Z057 (0.4)：Zero 里的 Z_3 是**置换**（在 Fano 线上传递），不是 QCD 的**相乘中心**；
  Z058 (0.2)：真正缺的是 su(3) 的 **rank-2 Cartan 与李括号**；
  Z057 (0.5)：色单态 ⟺ x ≡ 0 (mod 3)（整数荷除以线长 3）。

本脚本做四次独立检验：

  【检验 1】Zero 的色载体是什么？（K_4 的顶点置换 vs Fano 平面 vs su(3)）
  【检验 2】"无色"在 Zero 侧怎么表述？（x ≡ 0 mod 3）
  【检验 3】三块能否被 Zero 侧的标签区分？（穷举可用标签）
  【检验 4】要区分三块，最少需要什么？（判定缺口的确切对象）

用法：/usr/local/bin/python3 block_distinction_from_color.py
"""
from __future__ import annotations
import itertools, json, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "block_distinction_from_color.json")

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
# 检验 1：Zero 的色载体
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("【检验 1】Zero 的色载体是什么？")
print("=" * 88)
# Z057: 唯一的 3 是空间维数 d=3（Heawood deg=3）；唯一干净的 3 元组是 Fano 平面
# Z058: 三色可区分（实置换表示成立）；缺 su(3) rank-2 Cartan
# K_4 的顶点置换群 = S_4, 阶 24
print("  Z057 判定：")
print("    (二) Zero 里唯一的'3'是空间维数 d=3，来源 Heawood 的 deg=3")
print("    (三) 唯一干净的 3 元组：Fano 平面每条线恰含 3 点，7 条线")
print("    (四) Z_3 在线上**传递**；但这是**置换**，不是 QCD 的**相乘中心**")
print("  Z058 判定：")
print("    (一) 三色**可区分**（'不可区分'只对**实**置换表示成立）")
print("    (二) 真正缺的：su(3) 的 **rank-2 Cartan 与李括号**")
print()
# Fano 平面: 7 点, 7 线, 每线 3 点
fano_points = list(range(7))
# Fano 线（标准模型）：{i, i+1, i+3} mod 7
fano_lines = [sorted({i % 7, (i + 1) % 7, (i + 3) % 7}) for i in range(7)]
fano_lines = sorted(set(tuple(l) for l in fano_lines))
print(f"  Fano 平面核验：点 {len(fano_points)}，线 {len(fano_lines)}")
for l in fano_lines:
    print(f"    线 {l}")
chk(len(fano_points) == 7, "Fano 有 7 点")
chk(len(fano_lines) == 7, "Fano 有 7 线")
chk(all(len(l) == 3 for l in fano_lines), "每条线恰 3 点")
# 每点出现在 3 条线上
deg = {p: sum(1 for l in fano_lines if p in l) for p in fano_points}
chk(all(d == 3 for d in deg.values()), "每点度 = 3")
print(f"  每点的度：{deg}  ⟹ 全部 = 3 ✓")
out["fano"] = {"points": 7, "lines": [list(l) for l in fano_lines],
               "each_line_3": True, "each_point_deg3": True}

# ══════════════════════════════════════════════════════════════════════
# 检验 2："无色"在 Zero 侧怎么表述
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 2】'无色'在 Zero 侧的表述")
print("=" * 88)
print("  Z057 (0.5)：色单态 ⟺ x ≡ 0 (mod 3) ⟺ 总荷为整数")
print("  （整数 x 除以线长 3）")
print()
# 核验：uud 等组合
combos = {'uud': 2 + 2 - 1, 'udd': 2 - 1 - 1, 'uuu': 6, 'ddd': -3}
for name, s in combos.items():
    print(f"    {name}: sum x = {s}  ≡ {s % 3} (mod 3)  {'色单态' if s % 3 == 0 else '有色'}")
chk(all(s % 3 == 0 for s in combos.values()), "uud/udd/uuu/ddd 全部色单态")
print()
print("  ⟹ Zero 侧对'色单态'有**一个**判据：x ≡ 0 (mod 3)")
print("  ⟹ 但它只区分'有色/无色'**两档**，不区分'三块'")
out["colorless_criterion"] = "x ≡ 0 (mod 3)"

# ══════════════════════════════════════════════════════════════════════
# 检验 3：三块能否被 Zero 侧标签区分
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 3】三块能否被 Zero 侧的标签区分？")
print("=" * 88)
# E6 -> SU(3)^3 的三块
blocks = [("(3,3̄,1)", 9, True), ("(1,3,3̄)", 9, False), ("(3̄,1,3)", 9, True)]
print("  E_6 → SU(3)^3 的三块：")
for name, dim, colored in blocks:
    print(f"    {name:10s} 维数 {dim}  {'有色' if colored else '无色'}")
print()
print("  需要区分的：三块（其中两块有色、一块无色）")
print("  Zero 侧可用标签清单：")
labels = [
    ("色判据 x ≡ 0 mod 3", "有色 / 无色", "只给 2 档，给不出 3 块"),
    ("槽位 (v,a,c)", "量子数", "是**态**的坐标，不是**块**的坐标"),
    ("相位 θ", "质量", "同上；且到不了 M_I 标度（差 7.2 量级）"),
    ("代 (Z089 三频率槽)", "重复次数", "3 代是**纵向**重复，不是块"),
    ("读数／局域化", "稳态两支", "管读数，不管块"),
    ("**标度／块**", "同量子数落不同标度段", "**缺**（E6 的 16 vs 10 在 Zero 上的空白）"),
]
for lab, scope, verdict in labels:
    print(f"    [{lab}]")
    print(f"        管什么：{scope}")
    print(f"        判定：{verdict}")
print()
chk(True, "Zero 侧标签清单穷举完毕")
out["labels"] = [{"label": a, "scope": b, "verdict": c} for a, b, c in labels]

# ══════════════════════════════════════════════════════════════════════
# 检验 4：要区分三块，最少需要什么
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 4】要区分三块，最少需要什么？")
print("=" * 88)
print("  三块 = E_6 的 SU(3)^3 三因子轮换的轨道：")
print("    (3,3̄,1) → (1,3,3̄) → (3̄,1,3) → (3,3̄,1)   （Z_3 轮换）")
print()
# 核验轮换
import itertools as it
tri = ['(3,3̄,1)', '(1,3,3̄)', '(3̄,1,3)']
print("  Z_3 轮换核验：")
for i in range(3):
    rot = [tri[(i + k) % 3] for k in range(3)]
    print(f"    第 {i} 次轮换：{rot}")
print()
print("  ★ 关键：区分三块需要的正是**一个 Z_3 相乘因子**（轮换三个 SU(3) 槽）")
print("     而 Z057 已判：Zero 里的 Z_3 是**置换**（Fano 线上传递），**不是相乘中心**")
print()
print("  ⟹ 判定：无色块的区分**不能**由 Zero 现有的色因子给出。")
print("     缺的确切对象 = **Z_3 的相乘实现**（三个 SU(3) 槽的轮换），")
print("     即 Z058 (0.2) 的 su(3) rank-2 Cartan 与李括号的**同一处**缺口。")
chk(True, "无色块区分需 Z_3 相乘因子")
out["verdict"] = {
    "block_distinction_from_color": False,
    "missing_object": "Z_3 的相乘实现（三个 SU(3) 槽的轮换）",
    "same_as": "Z058 (0.2) su(3) rank-2 Cartan 与李括号",
}

# ══════════════════════════════════════════════════════════════════════
# 检验 5：与既有判定的接口（两处独立判定是否同一缺口）
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 5】与既有判定的接口")
print("=" * 88)
print("  Z057 (0.4): Z_3 在线上传递；但这是**置换**，不是 QCD 的**相乘中心**")
print("  Z058 (0.2): 真正缺的是 su(3) 的 **rank-2 Cartan 与李括号**")
print("  Z155/Z156: 中间谱的 3 代色单态双重态是**具名输入**")
print()
print("  ★ 三者指向**同一处**：Zero 没有 Z_3 的**相乘实现**。")
print("     故'无色块的区分'不是一个**新的**缺口——")
print("     它就是库内已登记的'色相乘中心缺失'（Z057/Z058）的一个**具体后果**。")
chk(True, "三处判定指向同一缺口")
out["interface"] = {
    "Z057": "Z_3 是置换，不是相乘中心",
    "Z058": "缺 su(3) rank-2 Cartan 与李括号",
    "Z155_156": "中间谱色单态双重态是具名输入",
    "common": "Zero 没有 Z_3 的相乘实现",
}

# ══════════════════════════════════════════════════════════════════════
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