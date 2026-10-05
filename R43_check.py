#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R43_check.py -- 核验分割族扫描与两层族构造性结果。"""
from __future__ import annotations
import io, itertools, json, math, os, sys
from collections import defaultdict
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

DOC = read("R43_pi_two_layer_construction.md")
R42 = read("R42_explicit_pi_rotation_class.md")
STATUS = read("STATUS.md"); INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R43_pi_family_scan_results.json"), encoding="utf-8"))
MU = np.array([math.sqrt(5), 1.381966011250105, 1.381966011250105])

def s_max(om):
    lam = np.sort(np.asarray(om, dtype=float))[::-1]
    if lam.size < 3: lam = np.concatenate([lam, np.zeros(3-lam.size)])
    return float(np.dot(lam[:3], MU))

def rank_of(om, tol=1e-6):
    w = np.asarray(om, float); w = w[w > 1e-300]
    D = sorted({round(abs(math.log(w[i]/w[j])), 10) for i in range(len(w)) for j in range(i+1, len(w))
                if abs(math.log(w[i]/w[j])) > 1e-12})
    basis = []
    for d in D:
        if all(abs(d/b - round(d/b)) > tol for b in basis): basis.append(d)
    return len(basis)

head("F1  文档结构：双侧约束、自然族扫描、两层族、边界")
check("标题与性质为构造性结果", "两层族" in DOC and "构造性结果（正面）" in DOC)
check("(R43-1)...(R43-5) 在位", all(("(R43-%d)" % k) in DOC for k in range(1, 6)))
check("写明两层族是构造而非导出", "构造" in DOC and "不是从 Zero 导出" in DOC)
check("写明未检查与 R32 的相容性", "没有检查该形状与 `R32` 生存要求的相容性" in DOC)
check("承接 R42 的双侧约束", "R42-5" in DOC and "R42-5" in R42)

head("F2  独立复算：自然族与两层族")
fams = RES["families"]
check("文档例 (2,5,20,100) 同时满足（S>2, 秩>=2）",
      fams["G_文档例_2,5,20,100"]["passes_both"] is True
      and abs(fams["G_文档例_2,5,20,100"]["S_max"] - 2.0327) < 1e-3)
check("旋转类不满足（谱太平）",
      fams["A_旋转类_L=12"]["contextual"] is False)
check("闭合长度族 S_max<2 但极接近（差<1%）",
      fams["C_闭合长度_L=20"]["S_max"] < 2 and (2 - fams["C_闭合长度_L=20"]["S_max"]) < 0.02,
      "S=%.4f" % fams["C_闭合长度_L=20"]["S_max"])
for p, alpha, K in ((0.80, 2.0, 8), (0.85, 1.0, 8), (0.85, 1.5, 8)):
    tail = np.array([a ** (-alpha) for a in range(1, K+1)]); tail = tail/tail.sum()*(1-p)
    om = np.concatenate([[p], tail])
    check("两层 p=%.2f alpha=%.1f K=%d：S>2 且秩>=2" % (p, alpha, K),
          s_max(om) > 2 and rank_of(om) >= 2, "S=%.4f rank=%d" % (s_max(om), rank_of(om)))
check("两层族通过数>0", fams["H_两层族"]["n_passers"] > 0,
      "passers=%d/%d" % (fams["H_两层族"]["n_passers"], fams["H_两层族"]["grid_size"]))

head("F3  账本一致")
check("STATUS 已登记 R43", "### 2.43" in STATUS and "R43_pi_two_layer_construction.md" in STATUS
      and "R43_check.py" in STATUS)
check("INDEX 已收录 R43", "R43_pi_two_layer_construction.md" in INDEX and "R43_check.py" in INDEX)
print(""); print("  不符项：%d" % FAIL)
if FAIL: print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
