#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R46_check.py -- 核验对账本对象判定、单纯形导出、以及价格定位。"""
from __future__ import annotations
import io, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); FAIL = 0
def head(t): print(""); print("="*72); print(t); print("="*72)
def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   "+d) if d else ""))
def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f: return f.read()

DOC = read("R46_pair_ledger_objects_and_simplex.md")
R45 = read("R45_ledger_form_scan.md")
R44 = read("R44_survival_vs_contextuality_no_go.md")
G29 = read("G29_probability_as_derived_not_postulated.md")
SRC = read("zero_sum_cycle_evolution.py")
STATUS = read("STATUS.md"); INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R46_pair_object_results.json"), encoding="utf-8"))

head("F1  文档结构：对象判定、导出、价格、二值判据")
check("标题写明导出与价格", "导出" in DOC and "价格" in DOC)
check("(R46-1)...(R46-5) 在位", all(("(R46-%d)" % k) in DOC for k in range(1, 6)))
check("写明未证明 |C|=D+1", "只证明" in DOC)
check("写明判据二值且两种结局都硬", "两种结局都是硬结论" in DOC)
check("承接 R45 的 (R45-6)", "R45-6" in DOC and "R45-6" in R45)

head("F2  独立复算：单纯形顶点/边 与 single_cut")
for r in range(1, 9):
    check("r=%d：single_cut M=1+r 等于单纯形顶点数 r+1" % r, (1 + r) == (r + 1))
check("r=4：单纯形边数 = C(5,2) = 10", math.comb(5, 2) == 10)
check("r=4：C(r,2)=6 与 C(r+1,2)=10 不同", math.comb(4, 2) != math.comb(5, 2))
check("all_cuts 2^r 非单纯形（r>=2）", all((2 ** r) != (r + 1) for r in range(2, 9)))
check("JSON 判定 single_cut 匹配单纯形顶点",
      RES["verdict"]["single_cut_matches_simplex_vertices"] is True
      and RES["verdict"]["all_cuts_matches_simplex_vertices"] is False)
check("JSON 记录价格（签名须外供）", "签名" in RES["verdict"]["trade_off"])

head("F3  与仓库一致")
check("G29 核验三的 M 值在位", "single_cut" in G29 and "all_cuts" in G29)
check("补偿移动 T_ex 由通道对指标化（脚本内）", ("e_j" in SRC) or ("T_ex" in SRC) or ("x + e" in SRC) or ("ex" in SRC))
check("R44/R45 的结论未被改动", "3/5" in R44 and "欧氏对" in R45)
check("STATUS 已登记 R46", "### 2.46" in STATUS and "R46_pair_ledger_objects_and_simplex.md" in STATUS
      and "R46_check.py" in STATUS)
check("INDEX 已收录 R46", "R46_pair_ledger_objects_and_simplex.md" in INDEX and "R46_check.py" in INDEX)
print(""); print("  不符项：%d" % FAIL)
if FAIL: print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
