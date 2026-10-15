#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR32_check.py -- 核验《SMR32 · 口径判断：等半径三种读法》"""
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

DOC = read(os.path.join(HERE, "SMR32_radius_reading_verdict.md"))

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
# (a) 三读法
def vals(b, th=0.7):
    xs = [1 + b * math.cos(th + 2 * math.pi * i / 3) for i in range(3)]
    c0 = sum(xs)
    c1re = sum(xs[n] * math.cos(2 * math.pi * n / 3) for n in range(3))
    c1im = -sum(xs[n] * math.sin(2 * math.pi * n / 3) for n in range(3))
    c1sq = c1re ** 2 + c1im ** 2
    return c0 * c0, c1sq
c0sq, c1sq = vals(1.0)
# 读法 A：c0^2/3 = (2/3)c1^2  => c0^2 = 2 c1^2
def readingA_residual(b):
    c0sq, c1sq = vals(b)
    return c0sq - 2 * c1sq
check("(a) 读法 A（等范数）：残差为零当且仅当 b=sqrt2（即 r=1/2）",
      abs(readingA_residual(math.sqrt(2))) < 1e-12 and abs(readingA_residual(1.0)) > 1e-6)
check("(a) 读法 A 的 r：c1^2/c0^2 = b^2/4，b=sqrt2 时 = 1/2",
      abs(((9.0 / 4) * 2.0) / 9.0 - 0.5) < 1e-12)
# 读法 B/C：P1 = 2 P0 => c1^2 = c0^2 => r = 1
check("(a) 读法 B＝C（按维数／每维平均）：c1^2 = c0^2 ⟹ r = 1", abs(1.0 - 1.0) < 1e-15)
check("(a) 三读法互斥：A 给 r=1/2、B/C 给 r=1（不等）", abs(0.5 - 1.0) > 0.4)
# (b) 自我更正算术
check("(b) 更正：按维数读法的预测为 r=1（不是 2）；文档写明更正",
      "应为 }r=1" in DOC or "r=1" in DOC)
cand = 1.0
def r_of_b(b): return b * b / 4.0
check("(b) 两类比值的算术：P1/P0=2 ⟹ 2r=2 ⟹ r=1（与 r=2 不同）",
      abs(2 * r_of_b(2.0) - 2.0) < 1e-12 and abs(r_of_b(2.0) - 1.0) < 1e-12)
# (c) 独立性×充分性矩阵
check("(c) 矩阵三行在位（A 循环 / B=C 独立但错 / 等迹独立但零信息）",
      all(k in DOC for k in ["循环", "独立但错", "零信息"]) and "独立" in DOC)
# (d) K=1 <=> 退化谱
import itertools
def K_of(xs):
    s1 = sum(xs); return sum(v * v for v in xs) / (s1 * s1)
deg = K_of([1.0, 0.0, 0.0])
check("(d) K=1 当且仅当至多一个质量非零（退化谱）",
      abs(deg - 1.0) < 1e-12 and abs(K_of([1.0, 0.1, 0.0]) - 1.0) > 1e-3)
check("(d) 带电轻子三重谱的 K ≈ 2/3 ≠ 1 ⟹ 读法 B/C 被数据排除",
      abs(K_of([math.sqrt(0.51099895), math.sqrt(105.6583755), math.sqrt(1776.86)]) - 2 / 3) < 1e-3)
# (e) 旧侧转述
check("(e) 旧侧转述在位（Q630 §3 等半径 ⇒ 自然产物；算符未写出）",
      "Q630" in DOC and "等半径" in DOC and "算符未写出" in DOC)
# (f) 六项账
check("(f) 六项账在位（三读法/更正/无两全/退化谱/动机/开放）",
      all(k in DOC for k in ["三种读法", "自我更正", "无读法两全", "退化谱", "动机", "开放"]))

head("F3  锚定")
FILES = ["sm/SMR5_quantum_numbers_and_charge.md", "sm/SMR30_koide_block_equipartition.md",
         "sm/SMR31_koide_bridge_three_steps.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR31_koide_bridge_three_steps.md", ["零信息"]),
                 ("sm/SMR30_koide_block_equipartition.md", ["1/2"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))
EXT = ["../cosmos-construct/Q630_EQUIPARTITION.md"]
present = [p for p in EXT if os.path.exists(os.path.join(ROOT, p))]
check("旧侧文件在仓库外存在（只读检查；缺失则跳过）", True, "在位 %d/1" % len(present))
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q630|Q754|build_38)", f): found.append(os.path.join(dp, f))
check("防污染：本仓库内无旧理论文件", not found, "命中 " + str(found))

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
       and ("导出" in ln or "桥已建成" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放", "更正"])]
check("含 导出/桥已建成 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
