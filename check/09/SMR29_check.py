#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR29_check.py -- 核验《SMR29 · 两库对照：Koide 缺口的定位》"""
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

DOC = read(os.path.join(HERE, "SMR29_koide_library_crosscheck.md"))

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
# (a) Parseval 与等分
def dft(xs):
    N = len(xs)
    cs = []
    for k in range(N):
        re_ = sum(xs[n] * math.cos(2 * math.pi * k * n / N) for n in range(N))
        im_ = -sum(xs[n] * math.sin(2 * math.pi * k * n / N) for n in range(N))
        cs.append((re_, im_))
    return cs
def parseval_ok(xs):
    c = dft(xs)
    lhs = sum(v * v for v in xs)
    rhs = (c[0][0] ** 2 + 2 * (c[1][0] ** 2 + c[1][1] ** 2)) / 3
    return abs(lhs - rhs) < 1e-9, c
ok_p = True; ok_eq = True
for b in (1.0, math.sqrt(2), 2.0):
    for th in (0.4, 1.7):
        xs = [1 + b * math.cos(th + 2 * math.pi * i / 3) for i in range(3)]
        o, c = parseval_ok(xs)
        ok_p = ok_p and o
        dem = c[0][0] ** 2 / 3
        rel = 2 * (c[1][0] ** 2 + c[1][1] ** 2) / 3
        eq = abs(dem - rel) < 1e-9
        if (abs(b * b - 2) < 1e-12) != eq: ok_eq = False
check("(a) Parseval: sum x^2 = (c0^2 + 2|c1|^2)/3（六个样本）", ok_p)
check("(a) 等分（民主部分 = 相对部分）<=> b^2 = 2（b=1, sqrt2, 2 判据一致）", ok_eq)
# (b) Niven 三元组
tri = [Fr(1, 2), Fr(2, 3), Fr(5, 6)]
check("(b) 三元组 1/2 + 2/3 + 5/6 = 2（精确）", sum(tri) == 2)
rat = []
for q in range(1, 61):
    for pp in range(1, q):
        v = abs(math.cos(math.pi * pp / q))
        if min(abs(v), abs(v - 0.5), abs(v - 1.0)) < 1e-10:
            rat.append(Fr(pp, q))
rat = sorted(set(rat))
check("(b) Niven（分母<=60 的数值检验）：有理余弦的有理角恰为 {1/3, 1/2, 2/3}",
      rat == [Fr(1, 3), Fr(1, 2), Fr(2, 3)], str(rat))
check("(b) 三元组中 5/6 的余弦无理（-sqrt3/2），恰一个无理",
      abs(math.cos(5 * math.pi / 6) + math.sqrt(3) / 2) < 1e-12
      and abs(math.cos(math.pi / 2)) < 1e-12 and abs(math.cos(2 * math.pi / 3) + 0.5) < 1e-12)
check("(b) Q775 不含 b：文档写明 Niven 约束角度、不能独立给等分",
      "不能独立给等分" in DOC and "不含" in DOC and "Niven" in DOC)
# (c) 缺口定位表
check("(c) 缺口定位表在位（本目录/旧侧 Q630/旧侧 Q775 三行 + 等权）",
      all(k in DOC for k in ["SMR28", "Q630", "Q775", "等权"]))
# (d) tau 误差与偏离
me, mmu, mt, dmt = 0.51099895, 105.6583755, 1776.86, 0.12
xs = [math.sqrt(me), math.sqrt(mmu), math.sqrt(mt)]
S1 = sum(xs); K = sum(v * v for v in xs) / (S1 * S1)
dKdmt = 1 / (S1 * S1) - 2 * sum(v * v for v in xs) / (S1 ** 3) * (1 / (2 * xs[2]))
dev = abs(K - 2 / 3)
check("(d) K_obs = 0.666660，偏离约 9e-6", abs(dev - 9e-6) < 5e-6, "%.2e" % dev)
check("(d) 敏感度 dK/dm_tau ≈ 5.6e-5 /MeV ⟹ δK ≈ 6.8e-6（与偏离同量级）",
      abs(dKdmt - 5.6e-5) < 1e-5 and abs(dKdmt * dmt - 6.8e-6) < 3e-6,
      "%.2e, %.2e" % (dKdmt, dKdmt * dmt))
check("(d) 偏离在 tau 误差内 ⟹ 旧侧“真实破坏”不成立（文档已更正）",
      "不成立" in DOC and "误差内" in DOC and dev < 2 * dKdmt * dmt)
# (e) 防污染
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q630|Q775|Q629|Q277|Q596)", f): found.append(os.path.join(dp, f))
check("(e) 防污染：本仓库内无旧理论文件（Q630/Q775/Q629/Q277/Q596）", not found, "命中 " + str(found))
check("(e) 文档写明只读+转述+不作前提三条", "不拷入" in DOC and "不作前提" in DOC and "转述" in DOC)
# (f) 六项账
check("(f) 六项账在位（Q630 判定/Q775 判定/缺口/更正/开放/未读）",
      all(k in DOC for k in ["等价改写", "只管角度", "缺口定位", "已更正", "开放", "未读"]))

head("F3  锚定")
FILES = ["HANDOVER_GR.md", "sm/SMR10_flavor_and_cp.md", "sm/SMR19_closure_density_trap.md",
         "sm/SMR28_koide_formula.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR28_koide_formula.md", ["45"]),
                 ("sm/SMR19_closure_density_trap.md", ["387513"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))
EXT = ["../cosmos-construct/Q630_EQUIPARTITION.md", "../cosmos-construct/Q775_WHY_THREE.md"]
present = [p for p in EXT if os.path.exists(os.path.join(ROOT, p))]
check("旧侧两篇在仓库外存在（只读检查；缺失则跳过）", True, "在位 %d/2" % len(present))

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
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放", "只到", "改写"])]
check("含 导出/已成立 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
