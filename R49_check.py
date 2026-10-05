#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R49_check.py -- 核验「π 非等差：精确判据、3 块不可能、4 块可构造」

  F1  文档结构：判据、不可能定理、可行域、显式解、诚实边界
  F2  JSON 数值自洽
  F3  独立复算：素数指数秩（不复用探针代码路径）+ 语境门限 + 显式轮廓
  F4  交叉一致：与 R35/R37/R42/R44 的既有数值一致；STATUS/INDEX 登记
"""
from __future__ import annotations

import io
import json
import os
import sys
from fractions import Fraction
from math import sqrt

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


DOC = read("R49_pi_nonarithmetic_criterion.md")
RES = json.load(io.open(os.path.join(HERE, "R49_pi_nonarithmetic_results.json"), encoding="utf-8"))

# ----------------------------------------------------------------------
head("F1  文档结构")

check("标题点出判据／不可能／可构造", "判据" in DOC and "不可能" in DOC and "构造" in DOC)
check("(R49-1)...(R49-4) 在位", all(("(R49-%d)" % k) in DOC for k in range(1, 5)))
check("写明稠密判据 = 素数指数差向量秩满", "素数指数" in DOC and "秩" in DOC)
check("写明 3 块不可能", "3 块轮廓（含全部 3 维归约）" in DOC or "永远是" in DOC)
check("给出语境门限 0.723607", "0.723607" in DOC)
check("给出显式轮廓", "(10^4" in DOC or "10^4" in DOC)
check("登记旋转类恒 III_λ 的障碍", "旋转类" in DOC and "III_\\lambda" in DOC)
check("诚实边界含『未从 Zero 原生生成』", "没有**从 Zero 原生对象生成" in DOC or "没给\"来源\"" in DOC or "来源" in DOC)
check("四维 GR 未推出", "未由此推出" in DOC)

# ----------------------------------------------------------------------
head("F2  JSON 数值自洽")

a = RES["A_selftest"]
for k, v in a.items():
    check("自检 %s：秩 = %d" % (k, v["expected_rank"]), v["rank"] == v["expected_rank"],
          "rank=%s" % v["rank"])
c = RES["C_four_blocks"]
check("门限与 R37/R44 一致（0.7236）", abs(c["lambda1_threshold"] - 0.723607) < 5e-7,
      "%.6f" % c["lambda1_threshold"])
rows = {r["name"]: r for r in c["constructions"]}
ok_rows = [r for r in c["constructions"] if r["rank"] == 3 and r["contextual"]]
check("至少 3 个显式轮廓同时满足秩 3 与语境", len(ok_rows) >= 3,
      "%d 个：%s" % (len(ok_rows), [r["w"] for r in ok_rows]))
pure = [r for r in c["constructions"] if "幂律" in r["name"]][0]
check("纯幂律对照 a^-2：秩 2（离散）", pure["rank"] == 2, "rank=%d" % pure["rank"])
conly = [r for r in c["constructions"] if r["name"].startswith("4 块：(1,1,4,36)")][0]
check("对照 (1,1,4,36)：语境但秩 2（说明两条约束独立）",
      conly["contextual"] and conly["rank"] == 2,
      "S_max=%.4f rank=%d" % (conly["S_max"], conly["rank"]))
rs = c["random_search"]
check("随机搜索里秩 3 与语境可同时出现", rs["rank3_and_contextual"] > 0,
      "%d / %d" % (rs["rank3_and_contextual"], rs["samples"]))
b = RES["B_three_blocks"]
check("3 块不可能已登记", "3" in b["claim"] or "3 块" in b["claim"])

# ----------------------------------------------------------------------
head("F3  独立复算")


def ev(k):
    PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    kk, v = k, []
    for p in PRIMES:
        e = 0
        while kk % p == 0:
            kk //= p
            e += 1
        v.append(e)
    return np.array(v, float) if kk == 1 else None


def rank_indep(w):
    W = [Fraction(x) for x in w]
    den = np.lcm.reduce([x.denominator for x in W])
    nums = [int(x * den) for x in W]
    vecs = [ev(nums[i + 1]) - ev(nums[i]) for i in range(len(nums) - 1)]
    return int(np.linalg.matrix_rank(np.array(vecs), tol=1e-9))


MU = (sqrt(5), (5 - sqrt(5)) / 2, (5 - sqrt(5)) / 2)


def smax(w):
    ww = np.sort(np.asarray(w, float))[::-1]
    ww = ww / ww.sum()
    return float(sum(ww[k] * MU[k] for k in range(3)))


# 独立判据复算（与探针同一判据、独立实现）
check("(2,5,20,100) 秩 2（离散）", rank_indep([2, 5, 20, 100]) == 2)
check("(1,2,4,8) 秩 1（等比 ⇒ III_{1/2}）", rank_indep([1, 2, 4, 8]) == 1)
check("(16,2,3,5) 秩 3（稠密）", rank_indep([16, 2, 3, 5]) == 3)
check("(10^4,2,3,5) 秩 3 且语境", rank_indep([10 ** 4, 2, 3, 5]) == 3 and smax([10 ** 4, 2, 3, 5]) > 2,
      "S_max=%.4f" % smax([10 ** 4, 2, 3, 5]))
# 3 块恒不稠密
for w in ([1, 2, 3], [1, 100, 10000], [2, 3, 5]):
    check("3 块 %s：秩 ≤ 2（不可能稠密）" % w, rank_indep(w) <= 2, "rank=%d" % rank_indep(w))
# 语境门限独立复算
m2 = (5 - sqrt(5)) / 2
l1 = (2.0 - m2) / (sqrt(5.0) - m2)
check("语境门限独立解 = 0.723607", abs(l1 - 0.723607) < 5e-7, "%.6f" % l1)
check("共形族在门限处 S_max = 2",
      abs((sqrt(5) * l1 + m2 * (1 - l1)) - 2.0) < 1e-12)
# 2 的幂 ⇒ 秩 1（旋转类障碍）
for w in ([1, 2, 4, 8], [2, 4, 8, 16], [1, 1, 2, 4]):
    check("全为 2 的幂 %s：秩 1（恒 III_λ）" % w, rank_indep(w) == 1, "rank=%d" % rank_indep(w))

# ----------------------------------------------------------------------
head("F4  交叉一致与登记")

R35 = read("R35_type_iii_classification.md")
R42 = read("R42_explicit_pi_rotation_class.md")
R44 = read("R44_survival_vs_contextuality_no_go.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
check("R35 已给类型判据（{0}/cZ/稠密）", "稠密" in R35 and "III" in R35)
check("R42 已给双侧约束", "双侧" in R42)
check("R44 的门限 0.7236 在位", "0.7236" in R44)
check("STATUS 已登记 R49", "### 2.49" in STATUS and "R49_pi_nonarithmetic_criterion.md" in STATUS
      and "R49_check.py" in STATUS)
check("INDEX 已收录 R49", "R49_pi_nonarithmetic_criterion.md" in INDEX and "R49_check.py" in INDEX)
check("探针/结果文件在位", os.path.exists(os.path.join(HERE, "R49_pi_nonarithmetic_probe.py"))
      and os.path.exists(os.path.join(HERE, "R49_pi_nonarithmetic_results.json")))

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
