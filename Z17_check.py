#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z17_check.py -- 核验 Z17（底层条款（Z0 条款 ＋ Z1–Z5 定理）退场后的空缺账本）：定理 Z17.1 与 §9 的本轮更新。
"""
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

Z17 = read("Z17_A0_A5_retirement_vacancy_ledger.md")
G89 = read("G89_dimension_no_go_and_the_balance_condition.md") if os.path.exists(
    os.path.join(HERE, "G89_dimension_no_go_and_the_balance_condition.md")) else ""
R44 = read("R44_survival_vs_contextuality_no_go.md")
R46 = read("R46_pair_ledger_objects_and_simplex.md")
R47 = read("R47_signature_from_causal_cone.md")
G29 = read("G29_probability_as_derived_not_postulated.md")
STATUS = read("STATUS.md"); INDEX = read("INDEX.md")

head("F1  定理 Z17.1 与逐条款核验表在位")
check("标题与空缺账本定位", "退场后的空缺账本" in Z17)
check("定理 Z17.1 在位且标【已证】", "定理 Z17.1" in Z17 and "（no-go 的 Z0 化）【已证】" in Z17)
check("Z0①/②/③ 与 Z1–Z5、识别 U 的逐条款表在位",
      all(k in Z17 for k in ("**Z0①**", "**Z0②**", "**Z0③**", "**Z1**", "**Z2**", "**Z3**", "**Z4**", "**Z5**")))
check("m-一致性条款在位（Z0 动力学条款救不了选维）",
      "Z0 的动力学条款" in Z17 and "排除（no-go）" in Z17)
check("§4 的单一空缺结论在位", "仅一项" in Z17 and "不被固定" in Z17)

head("F2  §9 的本轮更新在位")
check("(Z17-3) 双侧约束与 no-go", "(Z17-3)" in Z17 and "互斥" in Z17)
check("(Z17-4) 对象＝通道", "(Z17-4)" in Z17 and "通道" in Z17)
check("(Z17-5) 签名由因果锥供出", "(Z17-5)" in Z17 and "离散不变量" in Z17)
check("(Z17-6) 净更新：|C| 收窄为二值，空缺仍在",
      "(Z17-6)" in Z17 and "二值" in Z17 and "仍" in Z17)
check("仓库证据 single_cut 单纯形在位", "single_cut" in Z17 and "单纯形" in Z17)

head("F3  与 R 系列一致（数值可复核）")
check("R44 的 q=3/5 与互斥结论", "3/5" in R44)
check("R46 的 single_cut 匹配单纯形顶点（G29 核验三）", "single_cut" in G29 and "1+r" in G29)
check("R47 的签名表 D=2..6", "1, 3, 0" in R47 or "(1, 3, 0)" in R47 or "1,3,0" in R47.replace(" ", ""))
mu1, mu2 = math.sqrt(5), (5-math.sqrt(5))/2
q_star = 0.6; lam1 = (1+math.sqrt(2*q_star-1))/2
S = mu2 + (mu1-mu2)*lam1
check("独立复算：S_max(3/5)=2", abs(S-2) < 1e-12, "S=%.12f" % S)
check("独立复算：C(5,2)-C(4,2)=4", math.comb(5,2)-math.comb(4,2) == 4)

head("F4  账本登记")
check("STATUS 已登记 Z17", "Z17" in STATUS and "Z17_A0_A5_retirement_vacancy_ledger.md" in STATUS)
check("INDEX 已收录 Z17", "Z17_A0_A5_retirement_vacancy_ledger.md" in INDEX and "Z17_check.py" in INDEX)
print(""); print("  不符项：%d" % FAIL)
if FAIL: print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
