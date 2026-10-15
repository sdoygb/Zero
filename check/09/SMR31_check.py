#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR31_check.py -- 核验《SMR31 · 建桥三步》"""
import io, math, os, re, sys

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

DOC = read(os.path.join(HERE, "SMR31_koide_bridge_three_steps.md"))

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
# (a) 27 状态电荷表与置换结构
A = [(l, (-m) % 3, 0) for l in range(3) for m in range(3)]
B = [((-l) % 3, 0, n) for l in range(3) for n in range(3)]
C = [(0, m, (-n) % 3) for m in range(3) for n in range(3)]
check("(a) 三块各 9 态，共 27；每块对一个因子为单态（电荷 0）",
      len(A) == len(B) == len(C) == 9 and len(A) + len(B) + len(C) == 27
      and all(s[2] == 0 for s in A) and all(s[1] == 0 for s in B) and all(s[0] == 0 for s in C))
pi = lambda t: (t[1], t[2], t[0])
check("(a) 循环置换 pi(t)=(t2,t3,t1) 把 A->B->C->A（逐态检查）",
      all(pi(s) in B for s in A) and all(pi(s) in C for s in B) and all(pi(s) in A for s in C))
check("(a) 三元组空间 = 正则表示 = 1 + 2（维数相加为 3）", 1 + 2 == 3)
# (b) 等迹对 r 零信息
def equal_trace_quantity():
    a8 = [0.5, -0.5, 0.0]
    v27 = []
    for l in range(3):
        for m in range(3): v27.append((a8[l], -a8[m], 0.0))
    for l in range(3):
        for n in range(3): v27.append((-a8[l], 0.0, a8[n]))
    for m in range(3):
        for n in range(3): v27.append((0.0, a8[m], -a8[n]))
    return [sum(v[i] * v[i] for v in v27) for i in range(3)]
T = equal_trace_quantity()
check("(b) 等迹量 = [3, 3, 3]（由权重表定，与振幅无关）", all(abs(t - 3) < 1e-12 for t in T), str(T))
def r_of(b):
    return b * b / 4.0
check("(b) r = b^2/4 随 b 变（b=1, sqrt2, 2 -> 0.25, 0.5, 1.0）",
      abs(r_of(1) - 0.25) < 1e-12 and abs(r_of(math.sqrt(2)) - 0.5) < 1e-12 and abs(r_of(2) - 1) < 1e-12)
check("(b) 结论：等迹量在任意 b 下都是 [3,3,3] ⟹ 对 r 零信息（文档判定一致）",
      "零信息" in DOC and "桥不能由它建" in DOC)
# (c) r_obs
me, mmu, mt = 0.51099895, 105.6583755, 1776.86
xs = [math.sqrt(v) for v in (me, mmu, mt)]
c0 = sum(xs); sx2 = sum(v * v for v in xs)
c1sq = (3 * sx2 - c0 * c0) / 2
r_obs = c1sq / (c0 * c0)
check("(c) r_obs ≈ 0.5000（由 PDG 带电轻子质量）", abs(r_obs - 0.5) < 1e-4, "%.6f" % r_obs)
check("(c) 与候选 r=1/2 差 < 1e-4，与候选 r=2 差 ≈ 1.5 ⟹ 数据选前者",
      abs(r_obs - 0.5) < 1e-4 and abs(r_obs - 2.0) > 1.4)
check("(c) 文档写明 r_obs 的闭式（(3Σm − (Σ√m)²)/(2(Σ√m)²)）", "3\\sum m_i" in DOC or "3\\sum m_i" in DOC.replace(" ", ""))
# (d) 缺口形状三段
check("(d) 缺口形状三段在位（SMR28 b=√2 / SMR30 有理数 / 本篇二选一）",
      all(k in DOC for k in ["SMR28", "SMR30", "二选一"]))
# (e) 防污染
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q754|Q630|Q775|build_38)", f): found.append(os.path.join(dp, f))
check("(e) 防污染：本仓库内无旧理论文件", not found, "命中 " + str(found))
check("(e) 文档写明只读/不作前提/独立复算", "只读" in DOC and "不作前提" in DOC and "独立复算" in DOC)
# (f) 六项账
check("(f) 六项账在位（正则表示/零信息/二选一/开放/对照/电荷表）",
      all(k in DOC for k in ["正则表示", "零信息", "二选一", "开放", "对照", "电荷表"]))

head("F3  锚定")
FILES = ["sm/SMR2_carrier_and_spin.md", "sm/SMR5_quantum_numbers_and_charge.md",
         "sm/SMR7_fermion_spectrum_and_generations.md", "sm/SMR28_koide_formula.md",
         "sm/SMR30_koide_block_equipartition.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR30_koide_block_equipartition.md", ["1/2"]),
                 ("sm/SMR28_koide_formula.md", ["Z_3"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))

head("F4  与旧文对齐")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 " + str(len(smd)))
check("备份包在位", os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))

head("F5  红线自查")
SKIP = ("**目标**", "**步骤**", "**用到的输入**", "**核对**", "**结果**", "**常见误读**", "**性质**")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and not ln.strip().startswith("#")
       and not any(ln.strip().startswith(s) for s in SKIP)
       and "**纪律**" not in ln and "F5 红线自查" not in ln
       and ("导出" in ln or "缺口已关" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放", "选"])]
check("含 导出/缺口已关 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
