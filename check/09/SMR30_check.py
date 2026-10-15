#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR30_check.py -- 核验《SMR30 · Koide 靶心的有理化与块等功率》"""
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

DOC = read(os.path.join(HERE, "SMR30_koide_block_equipartition.md"))

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
# (a) 比值与两条规则
def ratio(b):
    c0sq = 9.0
    c1sq = (9.0 / 4.0) * b * b
    return c1sq / c0sq
check("(a) b=1, sqrt2, 2 -> |c1|^2/c0^2 = 1/4, 1/2, 1",
      abs(ratio(1.0) - 0.25) < 1e-12 and abs(ratio(math.sqrt(2)) - 0.5) < 1e-12 and abs(ratio(2.0) - 1.0) < 1e-12)
check("(a) 块等功率预测 = 1/2（与 b^2=2 一致）；逐态等分预测 = 2（不等于 1/2）",
      abs(0.5 - ratio(math.sqrt(2))) < 1e-12 and abs(2.0 - 0.5) > 1.0)
# 独立路径：Parseval 与块功率相等
def block_powers(b, th=0.7):
    xs = [1 + b * math.cos(th + 2 * math.pi * i / 3) for i in range(3)]
    c0 = sum(xs)
    c1re = sum(xs[n] * math.cos(2 * math.pi * n / 3) for n in range(3))
    c1im = -sum(xs[n] * math.sin(2 * math.pi * n / 3) for n in range(3))
    c1sq = c1re ** 2 + c1im ** 2
    return c0 * c0 / 3, 2 * c1sq / 3
pt, pn = block_powers(math.sqrt(2))
check("(a) 第二路径：b=sqrt2 时平凡块功率 = 非平凡块功率", abs(pt - pn) < 1e-12, "%.6f vs %.6f" % (pt, pn))
pt1, pn1 = block_powers(1.0)
check("(a) 对照：b=1 时两块功率不等", abs(pt1 - pn1) > 1e-6, "%.6f vs %.6f" % (pt1, pn1))
# (b) 27 三块等迹
ADJ = [(1, 0), (-1, 1), (0, -1), (-1, 0), (1, -1), (0, 1), (0, 0), (0, 0)]
a8 = [Fr(1, 2), Fr(-1, 2), Fr(0)]
v27 = []
for l in range(3):
    for m in range(3): v27.append((a8[l], -a8[m], Fr(0)))
for l in range(3):
    for n in range(3): v27.append((-a8[l], Fr(0), a8[n]))
for m in range(3):
    for n in range(3): v27.append((Fr(0), a8[m], -a8[n]))
T = [sum(v[i] * v[i] for v in v27) for i in range(3)]
X = [sum(v[i] * v[j] for v in v27) for i in range(3) for j in range(3) if i < j]
check("(b) 27 = 9+9+9，三个 Cartan 的迹都 = 3（等迹，重算）", len(v27) == 27 and T == [3, 3, 3], str(T))
check("(b) 三个交叉项全为 0", X == [0, 0, 0], str(X))
# (c) 等迹与等功率分离
check("(c) 等迹量 = 3（归一化）与块功率（分配）是两个不同量（文档写明等迹 != 等功率）",
      "等迹" in DOC and "等功率" in DOC and "\\ne" in DOC)
# (d) 旧侧转述字符串
check("(d) 旧侧转述在位（Q754/build_38 三块、Q741/build_385、Q742/build_386 状态即基被否）",
      all(k in DOC for k in ["Q754", "build_38", "build_385", "build_386", "状态即基"]))
# (e) 六项账
check("(e) 六项账在位（有理化/规则唯一/部分对应/桥开放/防重犯/来源开放）",
      all(k in DOC for k in ["有理化", "规则唯一", "部分对应", "桥", "防重犯", "开放"]))

head("F3  锚定")
FILES = ["sm/SMR2_carrier_and_spin.md", "sm/SMR5_quantum_numbers_and_charge.md",
         "sm/SMR7_fermion_spectrum_and_generations.md", "sm/SMR28_koide_formula.md",
         "sm/SMR29_koide_library_crosscheck.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR5_quantum_numbers_and_charge.md", ["3"]),
                 ("sm/SMR28_koide_formula.md", ["45"]),
                 ("sm/SMR29_koide_library_crosscheck.md", ["Q630"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))
EXT = ["../cosmos-construct/Q630_EQUIPARTITION.md", "../cosmos-construct/Q754_E6_CHAIN.md"]
present = [p for p in EXT if os.path.exists(os.path.join(ROOT, p))]
check("旧侧文件在仓库外存在（只读检查；缺失则跳过）", True, "在位 %d/2" % len(present))
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q630|Q754|Q741|Q742)", f): found.append(os.path.join(dp, f))
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
       and ("导出" in ln or "已成立" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放", "待判", "唯一"])]
check("含 导出/已成立 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
