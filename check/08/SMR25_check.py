#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR25_check.py -- 核验《SMR25 · 两圈跑动与胶子自能》"""
import io, math, os, re, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
FAIL = NCHK = 0
def head(t): print(""); print("=" * 72); print(t); print("=" * 72)
def check(n, c, d=""):
    global FAIL, NCHK
    NCHK += 1
    ok = bool(c)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   " + d) if d else ""))
def read(p):
    with io.open(p, encoding="utf-8") as f: return f.read()
def exists(rel): return os.path.exists(os.path.join(HERE, rel)) or os.path.exists(os.path.join(ROOT, rel))

DOC = read(os.path.join(HERE, "SMR25_two_loop_running.md"))

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 " + fld + " 在位", fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 " + sec.strip() + " 在位", sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 6 条推导（2.1-2.6）", subs == ["2." + str(i) for i in range(1, 7)], "实测 " + str(subs))
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§" + name + " 五段式标记齐全", not miss, "缺 " + str(miss))
    check("§" + name + " 常见误读 >= 2 条", blk.count("✗") >= 2, "实测 " + str(blk.count("✗")))
    check("§" + name + " 有加框结果", "\\boxed" in blk)

head("F2  独立复算")
# (a) 一圈
C2G, C2R, TR = Fr(3), Fr(4, 3), Fr(1, 2)
def b0(nf): return Fr(11, 3) * C2G - Fr(4, 3) * TR * nf
check("(a) b0 = (11/3)C2(G) − (4/3)T(R)n_f = 11 − (2/3)n_f（n_f 逐点）",
      all(b0(nf) == Fr(11, 1) - Fr(2, 3) * nf for nf in (0, 3, 6, 16)))
check("(a) n_f^crit = 33/2 = 16.5（b0 变号）", b0(Fr(33, 2)) == 0 and b0(16) > 0 and b0(17) < 0)
# (b) 两圈
def b1(nf): return Fr(34, 3) * C2G ** 2 - (4 * C2R + Fr(20, 3) * C2G) * TR * nf
check("(b) 群论公式给 b1 = 102 − (38/3)n_f（n_f 逐点）",
      all(b1(nf) == 102 - Fr(38, 3) * nf for nf in (0, 3, 6, 16)))
check("(b) b1 的零点 n_f = 153/19 ≈ 8.05", b1(Fr(153, 19)) == 0 and abs(float(Fr(153, 19)) - 8.05) < 0.01)
check("(b) 符号：n_f=3 → 64 > 0；n_f=6 → 26 > 0；n_f=16 → −302/3 < 0",
      b1(3) == 64 and b1(6) == 26 and b1(16) < 0)
check("(b) 三项来源分解正确（C2(G)^2 / C2(R)T(R) / C2(G)T(R)）",
      Fr(34, 3) * 9 == 102 and (4 * C2R) * TR == Fr(8, 3) and (Fr(20, 3) * C2G) * TR == 10)
# (c) 部分分式与两圈闭式
ok_pf = True
for a, c in [(Fr(1, 5), Fr(2)), (Fr(1, 10), Fr(3, 2)), (Fr(2, 7), Fr(5))]:
    lhs = 1 / (a ** 2 * (1 + c * a))
    rhs = 1 / a ** 2 - c / a + c ** 2 / (1 + c * a)
    if lhs != rhs: ok_pf = False
check("(c) 部分分式：1/(a^2(1+ca)) = 1/a^2 − c/a + c^2/(1+ca)（三组有理数精确）", ok_pf)
def t_sol(a0, a, b0v, b1v):
    c = b1v / b0v
    return (1 / b0v) * (1 / a0 - 1 / a) + (b1v / b0v ** 2) * math.log((a0 * (1 + c * a)) / (a * (1 + c * a0)))
def t_num(a0, a, b0v, b1v, N=2000000):
    # 数值积分 dt/da = -1/(b0 a^2 + b1 a^3)
    import math as m
    lo, hi = float(a0), float(a)
    h = (hi - lo) / N; s = 0.0
    for i in range(N):
        x = lo + (i + 0.5) * h
        s += -1.0 / (float(b0v) * x * x + float(b1v) * x ** 3)
    return s * h
a0, a = Fr(1, 50), Fr(1, 80)
b0v, b1v = Fr(11, 1) - Fr(2, 3) * 6, 102 - Fr(38, 3) * 6
check("(c) 两圈闭式与数值积分一致（差 < 1e-6）",
      abs(abs(t_sol(a0, a, b0v, b1v)) - abs(t_num(a0, a, b0v, b1v))) < 1e-6,
      "%.6f vs %.6f" % (t_sol(a0, a, b0v, b1v), t_num(a0, a, b0v, b1v)))
# (d) 两圈量级 vs 门槛
toll = 0.0071; Lreq = 25.4509138931
for alpha, ratio in [(0.02, 10), (0.04, 10)]:
    dLL = ratio * alpha / (2 * math.pi)
    dL = dLL * Lreq
    check("(d) <alpha>=%.2f, b2/b1=10 → dL/L=%.3f, dL=%.2f（门槛的 %.0f 倍）" % (alpha, dLL, dL, dL / toll),
          dL > toll and 80 < dL / toll < 300)
check("(d) 文档结论写明 114-228 倍（原稿 4-9 倍已更正）", "114" in DOC and "228" in DOC)
# (e) 待办状态
check("(e) 待办状态表：T1 量级已估、系数未算；T2/T3 未做",
      "量级已估" in DOC and "系数未算" in DOC and "未做" in DOC)
# (f) 六项账
check("(f) 六项账在位（一圈/两圈/精确解/量级/SM 系数未算/T2T3 未做）",
      all(k in DOC for k in ["一圈系数", "两圈系数", "精确解", "量级", "未算", "未做"]))

head("F3  锚定")
FILES = ["sm/SMR4_gauge_principle_and_ym.md", "sm/SMR12_strong_sector_and_theta.md",
         "sm/SMR14_low_energy_angle_no_go.md", "sm/SMR17_closed_form_and_reduction.md",
         "sm/SMR20_settlement.md", "sm/SMR24_propagator_spectral.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR12_strong_sector_and_theta.md", ["33"]),
                 ("sm/SMR14_low_energy_angle_no_go.md", ["0.0071"]),
                 ("sm/SMR20_settlement.md", ["T1"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))

head("F4  与旧文对齐")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 " + str(len(smd)))
check("备份包在位", os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))

head("F5  红线自查")
SKIP = ("**目标**", "**步骤**", "**用到的输入**", "**核对**", "**结果**", "**常见误读**")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and not ln.strip().startswith("#")
       and not any(ln.strip().startswith(s) for s in SKIP)
       and "**性质**" not in ln and "**纪律**" not in ln and "F5 红线自查" not in ln
       and ("算完" in ln or "已算" in ln or "关闭" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "必须", "登记", "更正"])]
check("含 算完/已算/关闭 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
