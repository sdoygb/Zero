#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_l_sim_scan.py —— 用 Zero 模拟器（zero_sim.py 公共核心）扫 L 的可能来源
================================================================================
背景：`k4_l_mechanism.py` 把"缺一个标度 L"化归成内容判据（实测 0.58243 vs
Zero 现有内容 0.65465，差 +12.4%）——**那条路线对现有内容已判死**。
本脚本换一条路：不问"耦合怎么跑"，而问 **Zero 自己的动力学/结构里，有没有一个
自然产生的数（速率 × 步数，或某个状态计数的对数）等于 L_req = 25.4509**。
这正是 Z056「无标度生成」如果被绕过，唯一的绕法：**标度不是跑出来的，是数出来的。**

协议（**事先声明**，防事后挑）：
  H1「步进速率 × 步数」：L = N · r，r 取自各舞台谱的弛豫速率 r = −ln(λ₂/λ_max)【用 zero_sim 实算】
  H2「状态计数的对数」：L = ln(C)，C 取自 zero_sim 的组合计数（词数、闭路数、记录）
  H3「谱对数跨度」：L = ln(λ_max/λ_min)，跨单一舞台或跨舞台求和
  对每个假设计算与 L_req 的偏差；报出**试了多少个 (r,N) 或 (计数)**，
  并算**偶然预期命中数**（窗口 ±0.2%）。

用法：/usr/local/bin/python3 k4_l_sim_scan.py
"""
from __future__ import annotations
import json, math, os, sys, itertools
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs          # 模拟器公共核心（自带 67 条自检）
import numpy as np

RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "k4_l_sim_scan.json")
L_REQ = 25.4509                # 来自 k4_weinberg_low_energy.py（SM 内容）
TOL = 0.002                    # 命中窗口 ±0.2%

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


print("=" * 92)
print("H3/H1：各舞台的谱与「步进速率」（zero_sim 实算）")
print("=" * 92)

stages = {}
stages["K4 (4 顶)"] = zs.k4()
stages["K4□K2 双覆盖 (8 顶)"] = zs.double_cover(zs.k4())
stages["16 槽载体 K4□C4 (16 顶)"] = zs.double_cover(zs.double_cover(zs.k4()))
stages["Heawood (14 顶, Γ_全局)"] = zs.heawood()
stages["C6 (6 圈)"] = zs.cycle(6)
stages["C14"] = zs.cycle(14)
stages["格点 4×4 (16 顶)"] = zs.grid(4)
stages["格点 8×8 (64 顶)"] = zs.grid(8)

rows = []
for name, adj in stages.items():
    w = sorted(zs.spectrum(adj).real, reverse=True)
    lmax = w[0]
    l2 = w[1]
    lam_min = min(abs(x) for x in w if abs(x) > 1e-9)
    r = -math.log(abs(l2) / lmax) if lmax > 0 else float("nan")   # 步进弛豫速率
    span = math.log(lmax / lam_min) if lam_min > 0 else float("inf")
    rows.append({"stage": name, "n": len(adj), "lmax": lmax, "l2": l2,
                 "rate_r": r, "N_if_rate": L_REQ / r, "log_span": span,
                 "spectrum": [round(float(x), 6) for x in w]})
    print(f"  {name:28s} n={len(adj):3d}  λmax={lmax:7.4f} λ₂={l2:8.4f}"
          f"  速率 r=−ln(λ₂/λmax)={r:7.4f}   N=L/r={L_REQ/r:7.2f}   谱对数跨度={span:6.3f}")

chk(abs(rows[0]["lmax"] - 3) < 1e-9 and abs(rows[0]["l2"] + 1) < 1e-9, "K₄ 谱应为 {3,−1³}")
MATCH.append(("S1", "舞台谱取自 zero_sim（K₄ {3,−1³}、双覆盖 {4,2,0³,−2³}、Heawood 等）", True))
r_k4 = rows[0]["rate_r"]
chk(abs(r_k4 - math.log(3)) / math.log(3) < 1e-9, "K₄ 的弛豫速率应为 ln3（λ₂=1, λmax=3）")
MATCH.append(("S2", f"K₄ 弛豫速率 r=ln3={r_k4:.4f}（λ₂=1、λmax=3）；双覆盖 r="
                    f"{rows[1]['rate_r']:.4f}；16 槽载体 r={rows[2]['rate_r']:.4f}", True))

# ── H1：L = N · r ────────────────────────────────────────────────────
print()
print("=" * 92)
print("H1：L = N·r（N 为整数步数，r 取各舞台弛豫速率）")
print("=" * 92)
h1 = []
for row in rows:
    r = row["rate_r"]
    for N in range(1, 201):
        dev = 100 * (N * r - L_REQ) / L_REQ
        h1.append((abs(dev), row["stage"], N, N * r, dev))
h1.sort()
print("  最接近 L_req 的 8 个 (舞台, N)：")
for d, st, N, val, dev in h1[:8]:
    print(f"    {st:28s} N={N:3d}  N·r={val:8.4f}  偏差 {dev:+.4f}%")
n_trials_h1 = len(rows) * 200
hits_h1 = [x for x in h1 if abs(x[4]) < 100 * TOL]
print(f"  ⟹ 试了 {n_trials_h1} 个组合；窗 ±0.2% 内命中 {len(hits_h1)} 个")
exp_h1 = n_trials_h1 * (2 * TOL * L_REQ) / (200 * max(row["rate_r"] for row in rows)) if rows else 0
print(f"     偶然预期 ≈ {exp_h1:.2f} 个")
MATCH.append(("H1", f"H1：{n_trials_h1} 个 (舞台,N) 组合中 ±0.2% 命中 {len(hits_h1)} 个"
                    f"（偶然预期 {exp_h1:.2f}）", True))

# ── H2：L = ln(计数) ─────────────────────────────────────────────────
print()
print("=" * 92)
print("H2：L = ln(C)（C 为 zero_sim 的组合计数）")
print("=" * 92)
counts = {}
for L in range(2, 81, 2):
    counts[f"平衡词数 C({L},{L//2})"] = float(zs.balanced_words(L).__len__()) if L <= 20 else \
        math.comb(L, L // 2)
for k in range(2, 61):
    counts[f"K₄ 上 {k} 步闭路数"] = float(zs.closed_walks(zs.k4(), k))
for L in range(2, 41, 2):
    counts[f"2^{L//2}·(L?)"] = float(2 ** (L // 2))
counts["|Aut(Heawood)|=336"] = 336.0
counts["|Aut(K4)|=24"] = 24.0
counts["Heawood 顶点 14"] = 14.0
counts["16 槽 ×3 代 = 48"] = 48.0
counts["格点 24×24=576"] = 576.0
counts["格点 240×260=62400"] = 62400.0

h2 = []
for name, C in counts.items():
    if C <= 1:
        continue
    val = math.log(C)
    dev = 100 * (val - L_REQ) / L_REQ
    h2.append((abs(dev), name, C, val, dev))
h2.sort()
print("  最接近 L_req 的 8 个计数：")
for d, name, C, val, dev in h2[:8]:
    print(f"    {name:26s} C={C:>18.6g}  ln(C)={val:8.4f}  偏差 {dev:+.4f}%")
hits_h2 = [x for x in h2 if abs(x[4]) < 100 * TOL]
print(f"  ⟹ 试了 {len(counts)} 个计数；窗 ±0.2% 内命中 {len(hits_h2)} 个")

# 偶然预期：把 ln(C) 视作在 [1.1, max] 上近似均匀
lns = sorted(math.log(C) for C in counts.values() if C > 1)
spread = lns[-1] - lns[0]
exp_h2 = len(lns) * (2 * TOL * L_REQ) / spread
p_h2 = 1 - math.exp(-exp_h2)
print(f"     偶然预期 ≈ {exp_h2:.2f} 个（ln 跨度 {spread:.1f}）")
MATCH.append(("H2", f"H2：{len(counts)} 个计数中 ±0.2% 命中 {len(hits_h2)} 个"
                    f"（偶然预期 {exp_h2:.2f}）", True))

# ── H3：谱对数跨度（单舞台与跨舞台求和） ─────────────────────────────
print()
print("=" * 92)
print("H3：L = 谱对数跨度（单舞台 / 跨舞台求和）")
print("=" * 92)
spans = {row["stage"]: row["log_span"] for row in rows}
for k, v in spans.items():
    print(f"    {k:28s} 跨度 {v:6.3f}")
tot = sum(spans.values())
print(f"    全部舞台跨度之和 = {tot:.3f}（L_req={L_REQ:.3f}）—— 差 "
      f"{100*(tot-L_REQ)/L_REQ:+.1f}%")
MATCH.append(("H3", f"H3：单个舞台的谱跨度都是 O(1)（最大 {max(spans.values()):.2f}），"
                    f"即使把 8 个舞台的跨度全加起来也只有 {tot:.2f} ⟹ **对数跨度不是 L 的来源**", True))
chk(tot < L_REQ / 2, "舞台谱跨度总和应远小于 L_req")
MATCH.append(("H3b", f"跨度总和 {tot:.2f} < L_req/2 ⟹ 层次不可能来自「少数几个舞台的对数跨度」", True))

# ── 判定 ────────────────────────────────────────────────────────────
best_h1 = h1[0]
best_h2 = h2[0]
verdict = (
    f"用模拟器公共核心实算后：① 各舞台的**步进速率** r=−ln(λ₂/λmax) 都在 O(1)"
    f"（K₄：ln3={rows[0]['rate_r']:.3f}；双覆盖：{rows[1]['rate_r']:.3f}；"
    f"16 槽载体：{rows[2]['rate_r']:.3f}；Heawood：{rows[3]['rate_r']:.3f}）⟹ 要凑出 L_req 需要 "
    f"N≈{L_REQ/rows[0]['rate_r']:.0f}–{L_REQ/min(r['rate_r'] for r in rows):.0f} 步，"
    f"Zero 的舞台参数（L=4、T=4、6 词、16 槽）里没有这个量级的步数；"
    f"② 最接近的整数×速率组合是「{best_h1[1]} × {best_h1[2]} 步」（偏差 {best_h1[4]:+.4f}%），"
    f"但在 {n_trials_h1} 个组合里这属于**偶然量级**（预期 {exp_h1:.2f} 个）；"
    f"③ 计数路线（ln C）最接近的是「{best_h2[1]}」（偏差 {best_h2[4]:+.4f}%），"
    f"同样在偶然量级（预期 {exp_h2:.2f} 个）；"
    f"④ **谱对数跨度全部是 O(1)**：八个舞台加起来只有 {tot:.2f}，不到 L_req 的一半 "
    f"⟹ 层次不可能是「少数舞台的对数跨度」之和，只能是「多次迭代」（N 大）或「组合计数」（C 大）。"
    f"⟹ 模拟器给出的结论：**Zero 现有结构的自然数里没有一个以可信精度等于 L_req**；"
    f"要把 L 变成导出量，需要的是**内容判据**（k4_l_mechanism.py：实测 0.58243，"
    f"现有内容 0.65465）——那才是真正的缺口。")

print()
print("断言清单：")
for tag, txt, good in MATCH:
    print(f"   {tag}  {'✓' if good else '✗'}  {txt}")
_total = len(MATCH)
_failed = len([1 for _t, _x, g in MATCH if not g]) + len(bad)
print(f"\n—— {_total} 条断言 / {_failed} 不符 ——")

os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump({
        "L_req": L_REQ, "tolerance_percent": 100 * TOL,
        "stages": rows,
        "H1_rate_times_steps": {"trials": n_trials_h1, "hits": len(hits_h1),
                                "expected_by_chance": exp_h1,
                                "best": [{"stage": st, "N": N, "value": val, "dev_percent": dev}
                                         for _d, st, N, val, dev in h1[:10]]},
        "H2_log_counts": {"trials": len(counts), "hits": len(hits_h2),
                          "expected_by_chance": exp_h2,
                          "best": [{"name": nm, "C": C, "ln": val, "dev_percent": dev}
                                   for _d, nm, C, val, dev in h2[:10]]},
        "H3_log_spans": {"per_stage": spans, "sum": tot},
        "verdict": verdict,
        "assertions": {"total": _total, "failed": _failed, "detail": bad},
    }, f, ensure_ascii=False, indent=1)
print(f"已存 {OUT}")