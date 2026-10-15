#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E07_dynamics_lookup_check.py -- 核验《资料查证：动力学在旧语料库里早已建好》

纪律：只读 / 纯标准库 / 不 import 项目脚本 / 不写文件

分组：
  A  D_L 的不可约表示计数与 sum d^2 = 2L
  B  群代数的 M_2 数 vs A_L 的 M_2 数（数目对不上）
  C  旧库资料的存在性与关键数字
  D  概率推前的非均匀性（复算）
"""
from __future__ import annotations
import math
import os
import re

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

OLD="/Users/oygb/Downloads/Zero-261010"

def find(name):
    for root,dirs,files in os.walk(OLD):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git")]
        if name in files: return os.path.join(root,name)
    return None

def has(path,pat):
    try: return re.search(pat,open(path,encoding="utf-8").read()) is not None
    except Exception: return False

print("=== E07 核验：动力学资料查证 ===")

# ---------------- A
print("\n[A] D_L 的不可约表示计数")
for L in [3,4,5,6,7,8,12]:
    if L%2==0: n1,n2=4,(L-2)//2
    else:      n1,n2=2,(L-1)//2
    s=n1*1+n2*4
    ck(s==2*L, "A: L=%-3d n1=%d n2=%d sum d^2=%d = 2L=%d"%(L,n1,n2,s,2*L))
    print("     L=%-3d 1维 %-2d 2维 %-2d sum d^2 = %-3d = 2L ✓"%(L,n1,n2,s))
ck(4%2==0, "A8: 框架自定 L=4 为偶")

# ---------------- B
print("\n[B] 群代数的 M_2 数 vs A_L 的 M_2 数")
for L in [4,6,8,12]:
    ga=(L-2)//2; aL=L//2
    ck(ga!=aL, "B: L=%-3d 群代数 %d 个 vs A_L %d 个（不等）"%(L,ga,aL))
    print("     L=%-3d 群代数的 M_2 数 = %-3d ; A_L 的 M_2 数 = %-3d ; 差 %d"%(
        L,ga,aL,aL-ga))
ck((4-2)//2==1 and 4//2==2, "B5: L=4 时 1 vs 2，确实不等")

# ---------------- C
print("\n[C] 旧库资料的存在性与关键数字")
checks=[
 ("G28_dynamics_audit.md", r"微观动力学", "C1: G28 存在且有'微观动力学'"),
 ("G28_dynamics_audit.md", r"建好并运行", "C2: G28 判'建好并运行'"),
 ("G28_dynamics_audit.md", r"max_total_charge_error\s*=\s*0", "C3: 零约束每步精确"),
 ("G28_dynamics_audit.md", r"确定性全分支", "C4: G28 有'确定性全分支'"),
 ("G28_dynamics_audit.md", r"违反 Z0", "C5: G28 记概率型违反公理③"),
 ("G29_probability_as_derived_not_postulated.md", r"推前", "C6: G29 有'推前'"),
 ("G29_probability_as_derived_not_postulated.md", r"K=-\\log\\omega|K\s*=\s*-\\log", "C7: G29 有 K = -log omega"),
 ("G29_probability_as_derived_not_postulated.md", r"非平凡", "C8: G29 判非平凡"),
 ("G27_purification_attempt.md", r"二面体", "C9: G27 有二面体群"),
 ("G27_purification_attempt.md", r"M_2\(\\mathbb C\)", "C10: G27 有 M_2(C)"),
 ("G43_dynamics_is_metric_free.md", r"metric-free", "C11: G43 metric-free"),
 ("G55_dynamics_line_degeneration_to_GR.md", r"记忆核", "C12: G55 记忆核"),
]
for fn,pat,msg in checks:
    p=find(fn)
    ck(p is not None, msg+" —— 文件存在")
    if p: ck(has(p,pat), msg+" —— 措辞匹配")

print("\n[C'] R34 的 K_omega 谱（JSON）")
p=find("R34_mixing_generator_results.json")
ck(p is not None, "C13: R34 结果 JSON 存在")
if p:
    import json
    d=json.load(open(p,encoding="utf-8"))
    sp=d.get("Q3_K_omega_span", None)
    ck(sp is not None and abs(sp-3.912023005428146)<1e-9,
       "C14: K_omega 跨度 = %.15f（非零）"%sp if sp else "C14: 缺 Q3_K_omega_span")
    print("     K_omega 跨度 = %s（非零 => 模 Hamiltonian 非平凡）"%sp)

# ---------------- D
print("\n[D] 概率推前的非均匀性（G29 §2 的复算）")
blocks=[2,5,20,100]; tot=sum(blocks)
ws=[b/tot for b in blocks]
ck(abs(ws[0]-0.015748)<1e-5, "D1: 最小权重 ≈ 0.015748（实得 %.6f）"%ws[0])
ck(abs(ws[3]-0.787402)<1e-5, "D2: 最大权重 ≈ 0.787402（实得 %.6f）"%ws[3])
ratio=ws[3]/ws[0]
ck(abs(ratio-50.0)<1e-6, "D3: w_max/w_min = %.1f"%ratio)
span=math.log(ratio)
ck(abs(span-3.912023005428146)<1e-9, "D4: log 跨度 = %.15f"%span)
ck(span>0, "D5: 跨度 > 0 => K 非平凡")
print("     权重 %s ; 比 %.1f ; log 跨度 %.15f" % (
    ["%.6f"%w for w in ws], ratio, span))

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
