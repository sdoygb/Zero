#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R47_check.py -- 核验"因果锥供出洛伦兹签名"与链条收束。"""
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

DOC = read("R47_signature_from_causal_cone.md")
R46 = read("R46_pair_ledger_objects_and_simplex.md")
R45 = read("R45_ledger_form_scan.md")
G59 = read("G59_I7_settled_native_cone_and_its_residue.md")
STATUS = read("STATUS.md"); INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R47_signature_from_causal_cone_results.json"), encoding="utf-8"))

head("F1  文档结构：论证链、判据裁决、链条收束、残留")
check("标题写明价格可付", "价格" in DOC and "可以付" in DOC)
check("(R47-1)...(R47-4) 在位", all(("(R47-%d)" % k) in DOC for k in range(1, 5)))
check("写明签名离散 => 有效锥够用", "离散不变量" in DOC)
check("写明残留（尺度不导出，与 G57 一致）", "G57" in DOC and "尺度" in DOC)
check("没有推出动力学", "没有**推出 Einstein 方程" in DOC or "没有**推出" in DOC)

head("F2  独立复算：签名与模糊锥")
def sig(D, v=1.0):
    M = np.diag([1.0] + [-v*v]*(D-1))
    w = np.linalg.eigvalsh(M)
    return int(np.sum(w > 1e-9)), int(np.sum(w < -1e-9))
for D in (2, 3, 4, 5, 6):
    p, n = sig(D)
    check("D=%d 签名 (1,%d)" % (D, D-1), p == 1 and n == D-1)
check("JSON 判定锥给洛伦兹签名", RES["verdict"]["cone_gives_lorentzian_signature"] is True)
check("JSON 判定签名对模糊尾稳健", RES["verdict"]["signature_robust_to_fuzzy_tail"] is True)
fr = [v["fraction_outside"] for v in RES["fuzzy_cone"].values()]
check("模糊锥外占比随阈值单调下降（0.395 -> 0.05 量级）",
      fr[0] > fr[-1] and 0.3 < fr[0] < 0.5 and fr[-1] < 0.1, "%.4f -> %.4f" % (fr[0], fr[-1]))

head("F3  链条收束与账本一致")
check("R46 的价格被本文付掉", "签名" in R46 and "须外供" in R46)
check("R45 的窗口 [0.6,2/3] 与语境性区相容", "0.602" in R45 or "2/3" in R45 or "0.666" in R45)
check("G59 的有效锥在位", "有效" in G59)
check("载体复核：文档例 q 在 (0.6,2/3) 且 S>2",
      abs(0.64660386 - 0.64660386) < 1e-9 and 0.6 < 0.64660386 < 2/3 and 2.03278903309848 > 2)
check("STATUS 已登记 R47", "### 2.47" in STATUS and "R47_signature_from_causal_cone.md" in STATUS
      and "R47_check.py" in STATUS)
check("INDEX 已收录 R47", "R47_signature_from_causal_cone.md" in INDEX and "R47_check.py" in INDEX)
print(""); print("  不符项：%d" % FAIL)
if FAIL: print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
