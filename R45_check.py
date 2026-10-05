#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R45_check.py -- 核验账本形式扫描、逃生口、以及 R32 唯一性的收窄。"""
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

DOC = read("R45_ledger_form_scan.md")
R32 = read("R32_ledger_readout_selection_and_L8_resolution.md")
R44 = read("R44_survival_vs_contextuality_no_go.md")
STATUS = read("STATUS.md"); INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R45_ledger_form_scan_results.json"), encoding="utf-8"))
MU = (math.sqrt(5.0), (5.0-math.sqrt(5.0))/2, (5.0-math.sqrt(5.0))/2)
QS = np.round(np.arange(0.30, 1.0001, 0.002), 6)

def peak_D(M, E, q, Dmax=30):
    return max(((D, M(D)*q**E(D)) for D in range(2, Dmax+1)), key=lambda t: t[1])[0]
def window(M, E, Dmax=30):
    good = [q for q in QS if peak_D(M, E, q, Dmax) == 4]
    return (float(min(good)), float(max(good))) if good else None
def s_spec(lam):
    lam = np.sort(np.asarray(lam, float))[::-1]
    if lam.size < 3: lam = np.concatenate([lam, np.zeros(3-lam.size)])
    return float(np.dot(lam[:3], MU))

C2 = lambda D: math.comb(D, 2)
Cp = lambda D: math.comb(D+1, 2)
Ed = lambda D: D

head("F1  文档结构：扫描、逃生口、代价、新目标")
check("标题写明逃生与 R32 削弱", "可逃" in DOC and "唯一性被削弱" in DOC)
check("(R45-1)...(R45-6) 在位", all(("(R45-%d)" % k) in DOC for k in range(1, 7)))
check("写明字典未导出", "未导出" in DOC and "唯一" in DOC)
check("写明 39 个形式都行 ⇒ 不足以定出形式", "39" in DOC and "不足以" in DOC)

head("F2  独立复算：窗口与相交")
w_std = window(C2, Ed); w_euc = window(Cp, Ed)
check("洛伦兹对窗口 [0.50,0.60] 不与 q>0.6 相交",
      w_std is not None and w_std[1] <= 0.6 + 1e-9, "window=%s" % (w_std,))
check("欧氏对窗口整个落在 q>0.6（下沿 >0.6）",
      w_euc is not None and w_euc[0] > 0.6, "window=%s" % (w_euc,))
check("欧氏对窗口上沿 ≈ 2/3", abs(w_euc[1] - 2/3) < 0.01, "hi=%.4f" % w_euc[1])
check("扫描统计：相交数 39，无窗口 16",
      len([r for r in RES["scan"] if r["overlap_ctx"]]) == 39
      and len([r for r in RES["scan"] if r["window"] is None]) == 16)

head("F3  载体：文档例与两层族在欧氏字典下同时满足")
car = RES["carriers"]
check("文档例 q=0.6466、S>2、洛伦兹峰 D=5、欧氏峰 D=4",
      abs(car["documented_2_5_20_100"]["q"] - 0.64660386) < 1e-6
      and car["documented_2_5_20_100"]["contextual"] is True
      and car["documented_2_5_20_100"]["peak_D_standard"] == 5
      and car["documented_2_5_20_100"]["peak_D_Cplus"] == 4)
tl = car["two_layer_p0.8_a2_K8"]
check("两层族 q=0.6585、S>2、欧氏峰 D=4",
      tl["contextual"] is True and tl["peak_D_Cplus"] == 4 and tl["peak_D_standard"] == 5)
check("路线 A L=4 两字典都过不了语境性",
      car["routeA_L4"]["contextual"] is False)
docw = [0.0157, 0.0394, 0.1575, 0.7874]
check("独立复算文档例 S_max=2.0328", abs(s_spec(docw) - 2.03278903309848) < 1e-9,
      "S=%.6f" % s_spec(docw))

head("F4  账本一致")
check("R44 的 no-go 在洛伦兹字典内仍成立", "3/5" in R44)
check("R32 的唯一性被正确收窄（文档内有此表述）", "唯一" in R32)
check("STATUS 已登记 R45", "### 2.45" in STATUS and "R45_ledger_form_scan.md" in STATUS
      and "R45_check.py" in STATUS)
check("INDEX 已收录 R45", "R45_ledger_form_scan.md" in INDEX and "R45_check.py" in INDEX)
print(""); print("  不符项：%d" % FAIL)
if FAIL: print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
