#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R44_check.py -- 核验"选维（生存）与单体语境性互斥"的解析 no-go。"""
from __future__ import annotations
import io, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); FAIL = 0
def head(t): print(""); print("="*72); print(t); print("="*72)
def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   "+d) if d else ""))
def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f: return f.read()

DOC = read("R44_survival_vs_contextuality_no_go.md")
R31 = read("R31_phase_ledger_and_lifetime_selection.md")
R32 = read("R32_ledger_readout_selection_and_L8_resolution.md")
R37 = read("R37_kcbs_contextuality.md")
STATUS = read("STATUS.md"); INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R44_survival_vs_contextuality_results.json"), encoding="utf-8"))
MU1 = math.sqrt(5.0); MU2 = (5.0 - math.sqrt(5.0)) / 2.0

def peak_D(q, Dmax=40):
    F = [math.comb(D, 2)*q**D for D in range(2, Dmax+1)]
    return 2 + int(np.argmax(F))
def sb(q):
    return None if 2*q-1 < 0 else MU2 + (MU1-MU2)*(1+math.sqrt(2*q-1))/2
def s_spec(lam):
    lam = np.sort(np.asarray(lam, float))[::-1]
    if lam.size < 3: lam = np.concatenate([lam, np.zeros(3-lam.size)])
    return float(np.dot(lam[:3], [MU1, MU2, MU2]))

head("F1  文档结构：两条要求、解析临界、no-go、逃逸出口")
check("标题与性质为解析 no-go", "互斥" in DOC and "解析 no-go" in DOC)
check("(R44-1)...(R44-7) 在位", all(("(R44-%d)" % k) in DOC for k in range(1, 8)))
check("写明 no-go 只在 R31/R32 账本框架内", "只在 `R31`/`R32` 的账本框架内成立" in DOC)
check("写明没有推翻 R31/R32/R37 各自结论",
      "没有**推翻" in DOC or "没有推翻" in DOC)
check("给出逃逸出口（改账本形式）", "改账本形式" in DOC and "R44-7" in DOC)

head("F2  独立复算：生存窗口、S 上界、临界点")
check("q=0.5556（路线A L=4）峰在 D=4", peak_D(5/9) == 4)
check("q=0.50 峰不在 D=4（临界下沿）", peak_D(0.50) != 4)
check("q=0.60 峰仍在 D=4 且 F5/F4=1.0000",
      peak_D(0.60) == 4 and abs(math.comb(5,2)*0.6**5/(math.comb(4,2)*0.6**4) - 1.0) < 1e-12)
check("q=0.65 峰在 D=5", peak_D(0.65) == 5)
check("S_max(3/5)=2 恰好", abs(sb(0.6) - 2.0) < 1e-12, "S=%.12f" % sb(0.6))
check("lambda1*(3/5)=(2-mu2)/(mu1-mu2)=0.723607",
      abs((1+math.sqrt(0.2))/2 - (2-MU2)/(MU1-MU2)) < 1e-12)
check("S_max 单调增：q<3/5 => S<2；q>3/5 => S>2",
      sb(0.55) < 2 < sb(0.65))
check("窗口互补：q<3/5 峰在 4；q>3/5 峰>=5",
      peak_D(0.58) == 4 and peak_D(0.62) == 5)

head("F3  两条既有结论落在互斥两侧")
check("路线 A L=4：生存 ✅ / 语境 ✗",
      peak_D(5/9) == 4 and s_spec([4/6, 2/6]) < 2, "S=%.4f" % s_spec([4/6, 2/6]))
docw = [0.0157, 0.0394, 0.1575, 0.7874]; qd = sum(w*w for w in docw)
check("文档例：语境 ✅ / 生存 ✗",
      s_spec(docw) > 2 and peak_D(qd) >= 5, "q=%.4f peak=%d S=%.4f" % (qd, peak_D(qd), s_spec(docw)))
check("JSON 记录互补性", RES["verdict"]["windows_are_complementary"] is True)
check("JSON 记录逃逸出口", "账本形式" in RES["verdict"]["escape"])

head("F4  账本一致")
check("R31/R32 的账本形式在被引文档中", "q_L" in R31 or "q^D" in R31 or "q_L^D" in R31)
check("R32 的 D=4 唯一存活结论在位", "D=4" in R32)
check("R37 的语境性结论在位", "2.0327" in R37 or "2.0146" in R37)
check("STATUS 已登记 R44", "### 2.44" in STATUS and "R44_survival_vs_contextuality_no_go.md" in STATUS
      and "R44_check.py" in STATUS)
check("INDEX 已收录 R44", "R44_survival_vs_contextuality_no_go.md" in INDEX and "R44_check.py" in INDEX)
print(""); print("  不符项：%d" % FAIL)
if FAIL: print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
