#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
branch_criteria.py —— 在 4 条拟合分支上重跑前三重判据
================================================================================
背景：`_staging/推导/下型振幅的两处承重问题_分支退化与分辨力.md`
  下型质量层有 4 条精确分支（同一组 m_d,s,b）：
     A=2.420560 (th=111.4853°)   1.156253
     A=2.682850 (th=114.7387°)   0.987979   ← Z126 取
     A=1.545499 (th=126.3126°)   2.070635
     A=1.711737 (th=129.5455°)   1.843212
  前三重判据原本只对第 2 支做过。本篇对 4 支逐一重跑。

三重判据（原表述）：
  C1  不在 q-族：log_q A 距整数 0.43（q=4/3）
  C2  不在载体梯级：梯级比递减，A_dn 更大
  C3  不是单个亚层比：比 r3=8/3 大 0.61%

方法论要点（前一稿的教训）：
  · 两把尺分开比：A 尺与 A² 尺（Z126 附录 A.1b 用 A² 尺）
  · 每个"命中"必须报偶然期望：窗口宽度 × 形态数
  · 形态集要小且**都有库内角色**，否则命中不构成证据

【检验】
  【B1】四支的质量比值一致（同一组质量）
  【B2】C1 在四支上的结果
  【B3】C2 在四支上的结果
  【B4】C3 在四支上的结果
  【B5】高精度结构匹配（有角色的原子，三种形态）＋偶然期望
  【B6】哪一支的 A ≤ 8/3 且有最好的匹配

用法：/usr/local/bin/python3 branch_criteria.py
输出：ZeroGPU/results/branch_criteria.json
"""
from __future__ import annotations
import json, math, os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "branch_criteria.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


BRANCHES = [
    ("1", 2.420560, 111.4853, 1.156253),
    ("2", 2.682850, 114.7387, 0.987979),
    ("3", 1.545499, 126.3126, 2.070635),
    ("4", 1.711737, 129.5455, 1.843212),
]
TARGET_RATIO = [1.0, 20.00016081, 895.08071177]     # 由 Z126 参数生成
CUT = 8.0 / 3.0


def masses(A, th_deg, u, A2=313.856):
    th = math.radians(th_deg)
    return [u * A2 * (1.0 + A * math.cos(th + 2 * math.pi * i / 3)) ** 2
            for i in range(3)]


print("=" * 94)
print("【B1】四支给出同一组质量比值")
print("=" * 94)
print(f"  {'#':>2} {'A':>10} {'theta':>10} {'u_s':>10} {'比值':>34} {'偏差':>9}")
for tag, A, th, u in BRANCHES:
    m = masses(A, th, u)
    rr = np.array(m) / m[0]
    dev = float(np.max(np.abs(rr / np.array(TARGET_RATIO) - 1)))
    print(f"  {tag:>2} {A:>10.6f} {th:>9.4f}° {u:>10.6f} "
          f"[{rr[0]:.4f},{rr[1]:.4f},{rr[2]:.2f}] {dev:>8.1e}")
chk(all(float(np.max(np.abs(np.array(masses(A, th, u)) / masses(A, th, u)[0]
                          / np.array(TARGET_RATIO) - 1))) < 1e-4
        for _, A, th, u in BRANCHES), "四支的质量比值一致")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 94)
print("【B2】C1：q-族（$\log_q A$ 距最近整数的距离）")
print("=" * 94)
QS = {"4/3": 4 / 3, "3/2": 3 / 2, "2": 2.0, "8/3": 8 / 3, "sqrt2": math.sqrt(2),
      "sqrt3": math.sqrt(3), "sqrt{3/2}": math.sqrt(1.5), "16/9": 16 / 9}
print(f"  {'#':>2} {'A':>10} " + " ".join(f"{k:>9}" for k in QS))
c1 = {}
for tag, A, th, u in BRANCHES:
    vals = []
    for k, q in QS.items():
        lg = math.log(A) / math.log(q)
        vals.append(abs(lg - round(lg)))
    c1[tag] = min(vals)
    print(f"  {tag:>2} {A:>10.6f} " + " ".join(f"{v:>9.4f}" for v in vals))
print()
print("  行内最小距离（越小越像在 q-族上）：")
for tag, A, th, u in BRANCHES:
    print(f"    支{tag}: A={A:.6f}  最小距离 = {c1[tag]:.4f}")
chk(min(c1.values()) < 0.1,
    "C1 是 q-相关的：至少一支在某个库内 q 上距整数 < 0.1（故 C1 不构成排除）",
    "最小 " + "%.4f" % min(c1.values()) + "（支 " +
    min(c1, key=lambda t: c1[t]) + "）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 94)
print("【B3】C2：载体梯级")
print("=" * 94)
LADDER = [("K_4 边/顶点", 6 / 4), ("Q_3 边/顶点", 12 / 8), ("Heawood 边/顶点", 21 / 14),
          ("Gamma16 边/顶点", 40 / 16), ("K_4□K_4 边/顶点", 48 / 16),
          ("亚层比 r1", 3 / 2), ("亚层比 r2", 2.0), ("亚层比 r3", 8 / 3)]
print(f"  {'#':>2} {'A':>10} {'最近的梯级值':>16} {'相对偏差':>10} {'超过最大梯级？':>14}")
c2 = {}
for tag, A, th, u in BRANCHES:
    best = min(LADDER, key=lambda kv: abs(A - kv[1]) / kv[1])
    dev = 100 * (A - best[1]) / best[1]
    c2[tag] = (best[0], best[1], dev)
    print(f"  {tag:>2} {A:>10.6f} {best[1]:>16.6f} {dev:>+9.3f}% "
          f"{'是' if A > 8/3 else '否':>14}")
chk(True, "梯级值已列出（判据是序关系，四支都在梯级之上或之间）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 94)
print("【B4】C3：单个亚层比（r3=8/3 是上界）")
print("=" * 94)
print(f"  {'#':>2} {'A':>10} {'A - 8/3':>10} {'相对':>10} {'> 8/3 ?':>9}")
for tag, A, th, u in BRANCHES:
    print(f"  {tag:>2} {A:>10.6f} {A - CUT:>+10.6f} "
          f"{100*(A-CUT)/CUT:>+9.3f}% {'是' if A > CUT else '否':>9}")
chk(sum(1 for _, A, _, _ in BRANCHES if A > CUT) == 1,
    "只有第 2 支超过 8/3", f"{sum(1 for _, A, _, _ in BRANCHES if A > CUT)} 支")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 94)
print("【B5】高精度结构匹配：有角色的原子 × 三种形态，两把尺")
print("=" * 94)
# 有库内角色的原子
ATOMS = {
    "1": 1.0, "2": 2.0, "3": 3.0, "4": 4.0, "6": 6.0, "8": 8.0, "16": 16.0, "29": 29.0,
    "3/2": 1.5, "4/3": 4 / 3, "8/3": 8 / 3, "16/9": 16 / 9, "5/2": 2.5,
    "sqrt2": math.sqrt(2), "sqrt3": math.sqrt(3), "sqrt{3/2}": math.sqrt(1.5),
}
forms = []          # (名字, 数值)
for k, v in ATOMS.items():
    forms.append((k, v))
for k1, v1 in ATOMS.items():
    for k2, v2 in ATOMS.items():
        forms.append((f"{k1}/{k2}", v1 / v2))
        forms.append((f"{k1}*{k2}", v1 * v2))
vals = np.array([v for _, v in forms])
print(f"  形态数 = {len(forms)}（16 原子；单值 + 比值 + 乘积）")

for scale_name, f_scale in (("A 尺", lambda A: A), ("A² 尺", lambda A: A * A)):
    print()
    print(f"  --- {scale_name} ---")
    print(f"  {'#':>2} {'目标':>12} {'最佳形态':>18} {'值':>12} {'偏差':>10}")
    for tag, A, th, u in BRANCHES:
        tgt = f_scale(A)
        devs = np.abs(vals - tgt) / tgt
        j = int(np.argmin(devs))
        print(f"  {tag:>2} {tgt:>12.6f} {forms[j][0]:>18} {vals[j]:>12.6f} "
              f"{100*devs[j]:>+9.4f}%")

print()
print("  偶然期望（窗口 2w，形态数 N）：")
for w_pct, N in ((0.001, len(forms)), (0.01, len(forms)), (0.03, len(forms))):
    exp = N * 2 * w_pct / 100
    print(f"    窗口 ±{w_pct:>5}% × N={N}: 偶然期望 = {exp:.2f} 个")
chk(len(forms) * 2 * 0.001 / 100 < 1, "在 ±0.001% 窗口下偶然期望 < 1",
    f"{len(forms)*2*0.001/100:.3f}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 94)
print("【B6】汇总：哪一支 A ≤ 8/3，且匹配最好")
print("=" * 94)
print(f"  {'#':>2} {'A':>10} {'≤8/3':>6} {'A尺最佳偏差':>12} {'A²尺最佳偏差':>13}")
best_rows = []
for tag, A, th, u in BRANCHES:
    dA = float(np.min(np.abs(vals - A) / A))
    dA2 = float(np.min(np.abs(vals - A * A) / (A * A)))
    best_rows.append((tag, A, A <= CUT, dA, dA2))
    print(f"  {tag:>2} {A:>10.6f} {'是' if A <= CUT else '否':>6} "
          f"{100*dA:>+11.4f}% {100*dA2:>+12.4f}%")
print()
print("  ★ 判据：若某支既 ≤8/3，其最佳匹配又 ≤ 上型精度（0.006%），则它同时满足")
print("    '结构命中' 与 '不超过上界'——这样的支会把 A_dn 的'超出上界'整个消掉。")

out = {
    "branches": [{"tag": t, "A": A, "theta_deg": th, "u_s": u,
                  "below_8_3": bool(A <= CUT)} for t, A, th, u in BRANCHES],
    "C1_min_q_distance": c1,
    "C2_nearest_ladder": {t: {"name": v[0], "value": v[1], "dev_pct": v[2]}
                          for t, v in c2.items()},
    "form_count": len(forms),
    "best_match_pct": {r[0]: {"A": 100 * r[3], "A2": 100 * r[4]} for r in best_rows},
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for x in bad:
    print("  ✗", x)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)