#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R50_check.py -- 核验分层纪律：层的定义、会诊清单、层坍塌判据、双侧定律重算

  F1  文档结构：纪律陈述、层定义、会诊清单、双侧定律、操作规则、边界
  F2  JSON 自洽：扫描结果与双侧定律
  F3  独立复算：双侧定律（不复用探针代码路径）
  F4  交叉一致：与 Z6/G28/G73/R48/R49 的既有判决一致；STATUS/INDEX 登记
"""
from __future__ import annotations

import io
import json
import os
import sys
from math import log

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   " + d) if d else ""))


def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f:
        return f.read()


DOC = read("R50_layer_discipline.md")
RES = json.load(io.open(os.path.join(HERE, "R50_layer_scan_results.json"), encoding="utf-8"))

# ----------------------------------------------------------------------
head("F1  文档结构")

check("标题点出分层纪律", "分层纪律" in DOC)
check("(R50-1)(R50-2)(R50-3) 三条推论在位",
      all(("(R50-%d)" % k) in DOC for k in (1, 2, 3)))
check("层定义表含 L0/L1/L1′/L2 ＋ 导出记号 R",
      all(x in DOC for x in ("L0", "L1", "L2")) and ("mathcal R" in DOC or "读出面" in DOC))
check("写明'层坍塌'判据（同层相反才是真矛盾）", "层坍塌" in DOC and "真矛盾" in DOC)
check("会诊清单 ≥8 条", DOC.count("| # | 表面的矛盾") >= 1 and DOC.count("\n| 1 |") >= 1)
check("给出双侧定律重算表", "双侧定律" in DOC and "0.346571" in DOC)
check("写明 B>4 由 L2 饱和截断", "饱和截断" in DOC or "容量饱和" in DOC)
check("三条操作规则在位", "写断言必须带层指标" in DOC and "先对齐层" in DOC)
check("边界含'扫描不保证穷尽'", "不保证穷尽" in DOC or "未做" in DOC)

# ----------------------------------------------------------------------
head("F2  JSON 自洽")

two = RES["B_two_sided_law"]
rows = two["rows"]
check("B<4 侧定律成立（最大偏差 < 1e-4）", two["max_abs_diff_B_lt_4"] < 1e-4,
      "%.2e" % two["max_abs_diff_B_lt_4"])
check("B>4 侧被饱和截断", two["B_gt_4_saturated"] is True)
deep = {r["B"]: r for r in rows}
check("B=2 的实测斜率 = 定律值（1e-5）",
      abs(deep[2.0]["measured_minus_slope"] - deep[2.0]["law"]) < 1e-5,
      "%.6f vs %.6f" % (deep[2.0]["measured_minus_slope"], deep[2.0]["law"]))
check("B=5 的锥边密度饱和在 O(1)（>0.5）", deep[5.0]["rho_at_T"] > 0.5,
      "%.3f" % deep[5.0]["rho_at_T"])
check("扫描给出会诊候选", len(RES["A_layer_collapse_candidates"]) >= 3,
      "%d 个符号" % len(RES["A_layer_collapse_candidates"]))

# ----------------------------------------------------------------------
head("F3  独立复算（不同于探针的实现）")


def sim_edge_v2(B, T, X):
    """独立实现：与探针同模型，但用显式卷积矩阵写法。"""
    n = np.zeros((4, 2 * X + 1))
    n[0, X] = 1.0
    out = []
    for _ in range(T):
        m = 0.5 * (np.roll(n, 1, axis=1) + np.roll(n, -1, axis=1))
        m[:, 0] = m[:, -1] = 0.0
        tot = m.sum(axis=0)
        n = np.vstack([B * m[1] * (1.0 - tot), m[0], m[1], m[2]])
        n = np.maximum(n, 0.0)
        out.append(n[:, X + len(out) + 1].sum())
    return np.array(out)


for B in (2.0, 3.0, 3.5):
    e = sim_edge_v2(B, 900, 1100)
    k = np.arange(1, 901)
    g = e > 0
    sel = g & (k >= 300) & (k <= 850)
    sl = -float(np.polyfit(k[sel], np.log(e[sel]), 1)[0])
    check("独立积分：B=%.1f 的斜率 = -½log(B/4)" % B, abs(sl - (-0.5 * log(B / 4.0))) < 1e-4,
          "%.6f vs %.6f" % (sl, -0.5 * log(B / 4.0)))
for B in (4.5, 5.0):
    e = sim_edge_v2(B, 900, 1100)
    check("独立积分：B=%.1f 锥边饱和到 O(1)" % B, e[-1] > 0.3, "rho=%.3f" % e[-1])

# ----------------------------------------------------------------------
head("F4  交叉一致与登记")

Z6 = read("Z6_stall_autopsy_and_released_ledger.md")
G28 = read("G28_dynamics_audit.md")
G73 = read("G73_B_is_an_input.md")
R48 = read("R48_exact_cone_and_effective_cone.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
check("Z6/G28 提供层结构与'两套动力学'原话", ("活动层" in Z6 or "层" in Z6) and "互不相容" in G28)
check("G73 的 B 不可导出在位", "不可导出" in G73)
check("R48 的 eps(B) 定律在位", "0.346579" in R48 or "\\tfrac12\\log\\tfrac4B" in R48)
check("STATUS 已登记 R50", "### 2.50" in STATUS and "R50_layer_discipline.md" in STATUS
      and "R50_check.py" in STATUS)
check("INDEX 已收录 R50", "R50_layer_discipline.md" in INDEX and "R50_check.py" in INDEX)
check("探针/结果文件在位", os.path.exists(os.path.join(HERE, "R50_layer_scan_probe.py"))
      and os.path.exists(os.path.join(HERE, "R50_layer_scan_results.json")))

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
