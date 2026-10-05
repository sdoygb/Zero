#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R35_check.py -- 核验 P2″ 判据与轮廓分类结果。

独立复算：差集、循环残差 rho、六种轮廓的类型判定，以及"几何因子不足以强制 III_lambda"。
"""

from __future__ import annotations

import io
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(title: str) -> None:
    print("")
    print("=" * 72)
    print(title)
    print("=" * 72)


def check(name: str, cond: bool, detail: str = "") -> None:
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def read(name: str) -> str:
    with io.open(os.path.join(HERE, name), encoding="utf-8") as handle:
        return handle.read()


DOC = read("R35_type_iii_classification.md")
R34 = read("R34_finite_dimensional_boost_obstruction.md")
R12 = read("R12_zero_native_gap_filling.md")
G29 = read("G29_probability_as_derived_not_postulated.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R35_type_iii_results.json"), encoding="utf-8"))


def classify(w: np.ndarray, window: float = 30.0):
    w = np.asarray(w, dtype=float)
    w = w / w.max()
    logw = np.log(w)
    logw = logw - logw[0]
    if np.ptp(logw) < 1e-9:
        return "trivial", 0.0, None
    vals = sorted({round(abs(logw[i] - logw[j]), 12)
                   for i in range(len(logw)) for j in range(i + 1, len(logw))
                   if 1e-12 < abs(logw[i] - logw[j]) <= window})
    if not vals:
        return "undecided", None, None
    D = np.array(vals)
    c = float(D.min())
    rho = float(np.abs(D / c - np.round(D / c)).max())
    return ("III_lambda" if rho < 1e-3 else "III_1"), rho, c


T = 48
cases = {
    "uniform": np.ones(T + 1),
    "native_2^a": 2.0 ** np.arange(T + 1),
    "geometric_0.99": 0.99 ** np.arange(T + 1),
    "polynomial": (np.arange(T + 1) + 1.0) ** 2,
    "factorial": np.array([math.factorial(a + 1) for a in range(T + 1)], dtype=float),
    "geometric_x_poly": (2.0 ** np.arange(T + 1)) * (np.arange(T + 1) + 1.0),
}

# ======================================================================
head("F1  文档结构：判据、读数、新约束")

check("标题与性质为 P2″ 结果＋对 π 的新约束",
      "P2″ 结果" in DOC
      and "几何模流反过来**约束 π**" in DOC
      and "对缺失输入的新约束" in DOC)
check("判据 (R35-1)(R35-2) 与关键读数 (R35-3) 在位",
      "(R35-1)" in DOC and "(R35-2)" in DOC and "(R35-3)" in DOC
      and "只要轮廓}\\textbf{不是恰好等差}" in DOC)
check("新约束 (R35-5) 在位",
      "(R35-5)" in DOC and "非等差" in DOC)
check("决定性否证判据在位（III_{1/2} ⇒ 非 QFT 模论 ⇒ 无 BW）",
      "III$_{1/2}$" in DOC and "BW 不适用" in DOC)
check("没有声称算出原生轮廓",
      "没有算出**原生**轮廓" in DOC and "没有证明极限含 III$_1$ 因子" in DOC)

# ======================================================================
head("F2  独立复算：轮廓分类")

res = {name: classify(w) for name, w in cases.items()}
check("均匀轮廓给平凡（G={0}）", res["uniform"][0] == "trivial")
check("原生满分支 2^a 给 III_lambda 且 λ=1/2",
      res["native_2^a"][0] == "III_lambda"
      and abs(math.exp(-res["native_2^a"][2]) - 0.5) < 1e-9,
      "rho=%.2e lambda=%.6f" % (res["native_2^a"][1], math.exp(-res["native_2^a"][2])))
check("幂律给 III_1（rho≈0.5）",
      res["polynomial"][0] == "III_1" and res["polynomial"][1] > 0.4,
      "rho=%.4f" % res["polynomial"][1])
check("阶乘给 III_1", res["factorial"][0] == "III_1", "rho=%.4f" % res["factorial"][1])
check("几何×多项式仍给 III_1（指数因子不足以强制 III_lambda）",
      res["geometric_x_poly"][0] == "III_1" and res["geometric_x_poly"][1] > 0.4,
      "rho=%.4f" % res["geometric_x_poly"][1])
check("几何 0.99^a 给 III_0.99",
      res["geometric_0.99"][0] == "III_lambda"
      and abs(math.exp(-res["geometric_0.99"][2]) - 0.99) < 1e-6)

# ======================================================================
head("F3  与结果 JSON 一致")

by_name = {p["name"]: p for p in RES["profiles"]}
check("JSON 记录原生分支为 III_lambda，λ=1/2",
      by_name["native_branching_2^a"]["type"] == "III_lambda"
      and abs(by_name["native_branching_2^a"]["lambda"] - 0.5) < 1e-9)
check("JSON 记录多项式／阶乘／几何×多项式为 III_1",
      by_name["polynomial_(a+1)^2"]["type"].startswith("III_1")
      and by_name["factorial_(a+1)!"]["type"].startswith("III_1")
      and by_name["stairs_2^a_times_(a+1)"]["type"].startswith("III_1"))
check("JSON 的 verdict 写明对 π 的约束",
      "非指数" in RES["verdict"]["consequence"] or "非等差" in RES["verdict"]["consequence"])
check("JSON 的 BW 需求写明 III_1", RES["verdict"]["bw_requires"] == "III_1")

# ======================================================================
head("F4  与既有账本一致")

check("G29 的 4 块例在位且比值递增（非等差）",
      "0.787402" in G29 and "50.0" in G29)
check("R12 的常数剖面机制被引为 {0} 情形",
      "常数剖面" in R12 and "常数轮廓" in DOC)
check("R34 的 P2″ 纲领被正确回指",
      "P2″（可测指纹，新增）" in R34 and "P2″" in DOC)

# ======================================================================
head("F5  账本登记")

check("STATUS 已登记 R35",
      "### 2.35" in STATUS
      and "R35_type_iii_classification.md" in STATUS
      and "R35_check.py" in STATUS
      and "非等差" in STATUS)
check("INDEX 已收录 R35",
      "R35_type_iii_classification.md" in INDEX
      and "R35_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
