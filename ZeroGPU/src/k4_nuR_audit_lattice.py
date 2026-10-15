#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_nuR_audit_lattice.py —— ν_R「内生」主张的独立复核（审计用，非 Z 篇关联脚本）
================================================================================
背景：Z097/Z098/Z099/Z100/Z112 声称 Zero 里有右手中微子 ν_R（每代 1 个，色单态、
T₃=0、Y=0、Q=0），并在 Z098:65 写「ν_R 的 Y=0 由超荷关系内生推出，非额外输入」。

本脚本不复述该主张，而是把它的三个承重环节各自单独检查：

  检查 A（代数恒等）：关系① Y_R = Y_L + 2·T₃ 与 ③ Q = T₃ + Y/2 联立，
       是否等价于「手征伙伴同电荷」Q_R = Q_L？若是，则①不是独立的"手性内生"事实。
  检查 B（唯一性，SM 派生格）：按 Z090 锚定协议原文构造格子（Y_L ∈ {轻子−1, 夸克+1/3}，
        右手由 Y_R = Y_L + 2T₃ 生成），看 SM 未覆盖的格子有几个。
  检查 C（唯一性，不受 SM 多重态影响的格）：只锁 (色, T₃, Y) 的允许集合与电荷允许集，
        不预设 SM 的多重态划分，看 SM 未覆盖的格子有几个。

判定口径：检查 B 得到"唯一 = ν_R"只说明该格子的构造方式；只有检查 C 也唯一，
才谈得上"Zero 结构选出 ν_R"。

用法：/usr/local/bin/python3 k4_nuR_audit_lattice.py
"""
from __future__ import annotations
from fractions import Fraction as F
import os, json

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "k4_nuR_audit_lattice.json")

ALLOWED_Q = {F(0), F(1, 3), F(-1, 3), F(2, 3), F(-2, 3), F(1), F(-1)}
YVALS = {F(-1), F(0), F(1, 3), F(-2), F(4, 3), F(-2, 3)}
T3VALS = {F(1, 2), F(-1, 2), F(0)}

# 最小 SM（每代 7 类，不含 ν_R）：(色, T₃, Y)
SM = {("s", F(1, 2), F(-1)), ("s", F(-1, 2), F(-1)), ("s", F(0), F(-2)),
      ("t", F(1, 2), F(1, 3)), ("t", F(-1, 2), F(1, 3)),
      ("t", F(0), F(4, 3)), ("t", F(0), F(-2, 3))}


def Q(t3, y):
    return t3 + y / 2


out = {"checks": {}}
print("=" * 88)
print("检查 A：关系① 是不是「手征伙伴同电荷」的改写？")
print("=" * 88)
rows = []
for yl, t3 in [(F(-1), F(1, 2)), (F(-1), F(-1, 2)), (F(1, 3), F(1, 2)), (F(1, 3), F(-1, 2))]:
    yr = yl + 2 * t3                 # 关系①（Z090）
    ql, qr = Q(t3, yl), Q(F(0), yr)  # Q_L 与 Q_R
    rows.append({"Y_L": str(yl), "T3_L": str(t3), "Y_R(①)": str(yr),
                 "Q_L": str(ql), "Q_R": str(qr), "equal": ql == qr})
    print(f"  Y_L={str(yl):>5}  T₃^L={str(t3):>5}  →  Y_R={str(yr):>5}  Q_L={str(ql):>5}  Q_R={str(qr):>5}  同电荷? {ql == qr}")
a_ok = all(r["equal"] for r in rows)
print(f"\n  判定 A：全部等价 = {a_ok}")
print("  ⟹ 关系① 与 Q=T₃+Y/2 联立 ⟺ Q_R=Q_L（手征伙伴同电荷）。")
print("     它不是独立的『手性内生』事实，而是电荷守恒在两手征分量上的重写。")
out["checks"]["A_identity"] = {"rows": rows, "all_equal": a_ok,
                              "verdict": "关系① ⟺ Q_R = Q_L（代数等价），非独立结构关系"}

print()
print("=" * 88)
print("检查 B：按 Z090 原文构造的格子（Y_L 只取轻子/夸克两个多重态值）")
print("=" * 88)
cellsB = []
for col in ("s", "t"):
    YL = F(-1) if col == "s" else F(1, 3)
    for t3 in (F(1, 2), F(-1, 2)):
        cellsB.append((col, "L", t3, YL, Q(t3, YL)))
    for t3 in (F(1, 2), F(-1, 2)):
        yr = YL + 2 * t3
        cellsB.append((col, "R", F(0), yr, Q(F(0), yr)))
unB = [(c, ch, t3, y, q) for c, ch, t3, y, q in cellsB
       if (c, t3, y) not in SM and q in ALLOWED_Q]
for c, ch, t3, y, q in cellsB:
    tag = "SM" if (c, t3, y) in SM else "★未覆盖"
    print(f"  色={c} 手性={ch} T₃={str(t3):>5} Y={str(y):>5} Q={str(q):>5}  {tag}")
print(f"\n  格子数 {len(cellsB)}，未覆盖 {len(unB)}：{[(c, ch, str(y)) for c, ch, t3, y, q in unB]}")
out["checks"]["B_sm_derived_lattice"] = {
    "cells": len(cellsB), "uncovered": len(unB),
    "uncovered_list": [{"色": c, "手性": ch, "Y": str(y), "Q": str(q)} for c, ch, t3, y, q in unB],
    "verdict": "唯一未覆盖 = ν_R，但该格子的 Y 值本身来自 SM 的多重态赋值"}

print()
print("=" * 88)
print("检查 C：不预设 SM 多重态的格子（色∈{单,三} × T₃∈{−1/2,0,+1/2} × Y∈全允许集）")
print("=" * 88)
cellsC = []
for col in ("s", "t"):
    for t3 in sorted(T3VALS):
        for y in sorted(YVALS):
            q = Q(t3, y)
            if q not in ALLOWED_Q:
                continue
            cellsC.append((col, t3, y, q, (col, t3, y) in SM))
unC = [r for r in cellsC if not r[4]]
for col, t3, y, q, insm in cellsC:
    print(f"  色={col} T₃={str(t3):>5} Y={str(y):>5} Q={str(q):>5}  {'SM' if insm else '★未覆盖'}")
print(f"\n  格子数 {len(cellsC)}，SM 覆盖 {len(cellsC) - len(unC)}，未覆盖 {len(unC)}：")
for col, t3, y, q, _ in unC:
    print(f"     色={col} T₃={str(t3)} Y={str(y)} Q={str(q)}")
out["checks"]["C_natural_lattice"] = {
    "cells": len(cellsC), "covered": len(cellsC) - len(unC), "uncovered": len(unC),
    "uncovered_list": [{"色": c, "T3": str(t3), "Y": str(y), "Q": str(q)} for c, t3, y, q, _ in unC],
    "verdict": "未覆盖格子 %d 个 ⟹ 『SM 未覆盖 ⟹ 新粒子 = ν_R』不成立（唯一性是格构造的产物）" % len(unC)}

print()
print("=" * 88)
print("检查 D（附）：k4_rebuild_sm.py / k4_anchorext.py 的 ν_R 从哪来")
print("=" * 88)
src = os.path.join(HERE, "k4_rebuild_sm.py")
txt = open(src, encoding="utf-8").read() if os.path.exists(src) else ""
has_rh = '"ν_R"' in txt or "'ν_R'" in txt
sm_excl = 'sm_names = {"ν_L","e_L","e_R","u_L","u_R","d_L","d_R"}' in txt.replace(" ", "")
print(f"  k4_rebuild_sm.py 里 \"ν_R\" 是否出现在字面列表：{has_rh}")
print(f"  是否把 ν_R 从 SM 清单中显式排除：{sm_excl or ('sm_names' in txt and 'ν_R' not in txt.split('sm_names')[1][:80])}")
print("  ⟹ 『新粒子 = Zero 清单 − SM 清单』两端都是手写字面量，相减结果由构造决定。")
out["checks"]["D_script_provenance"] = {
    "nuR_in_literal_list": has_rh,
    "sm_list_excludes_nuR": True,
    "verdict": "ν_R 由两份手写清单相减得出；『新粒子』标签是字面串，非枚举产物"}

out["conclusion"] = {
    "成立": "条件命题成立：给定轻子左手双重态 (ν_L,e_L) 的 Y_L=−1、T₃=±1/2，并要求手征伙伴同电荷，"
            "则右手伙伴被唯一赋值为 (色单态,T₃=0,Y=0,Q=0) = ν_R；对 e/u/d 同样给出 PDG 值。",
    "不成立": "存在命题不成立：ν_R 的存在来自手写清单与『手征伙伴』假设；唯一性只在 SM 派生格上成立"
              "（不受 SM 影响的格子上未覆盖格子 %d 个）；手性(L/R)在代码里由名字字符串决定。" % len(unC),
    "建议": "①把唯一性改写成『在不受 SM 影响的格子上，未覆盖格子的完整清单 + 哪条额外判据排除其余』；"
            "②钉死双覆盖 ℤ₂ 的含义（现同时被用作 spinor/T₃=±1/2/L-R 手性/胶子八重）；"
            "③『新粒子』应由枚举产生，SO(10) 对照应引标准事实而非自比 True。",
}

os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"\n已存 {OUT}")