#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR15_check.py -- 核验《SMR15 · 专攻 L：原生池零命中，L_req 只依赖 b1-b2》"""
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

DOC = read(os.path.join(HERE, "SMR15_hierarchy_attack_on_L.md"))
Lreq, TOL = 25.4509138931, 0.0071

def L_from(b1, b2, target=0.23122, aem=127.95):
    lo, hi = 0.0, 60.0
    for _ in range(200):
        mid = (lo + hi) / 2
        A = Fr(3, 8) * (Fr(aem) + (b1 - b2) * Fr(mid) / (2 * math.pi)); B = Fr(aem) - Fr(5, 3) * A
        if float(1 / (1 + Fr(5, 3) * (A / B))) > target: lo = mid
        else: hi = mid
    return (lo + hi) / 2

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 4 条推导（2.1-2.4）", subs == ["2.%d" % i for i in range(1, 5)], "实测 %s" % subs)
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

head("F2  独立复算")
pool = {}
for n in range(2, 61, 2): pool["C(%d,%d)" % (n, n // 2)] = math.log(math.comb(n, n // 2))
for n, v in {2:1,4:2,6:4,8:10,10:26,12:80,14:246,16:810,18:2704}.items(): pool["K(%d)" % n] = math.log(v)
for n in range(3, 13): pool["%d!" % n] = math.log(math.factorial(n))
for nm, c in [("pi",math.pi),("2pi",2*math.pi),("4pi",4*math.pi),("8pi",8*math.pi),("e^3",math.e**3),("1/alpha",137.036)]:
    pool["ln(%s)" % nm] = math.log(c)
for nm, v in [("5/3",5/3),("3/8",3/8),("41/10",41/10),("109/15",109/15),("12",12),("24",24),("48",48)]:
    pool[nm] = math.log(v)
vals = list(pool.values()); spread = max(vals) - min(vals)
hits = [k for k, v in pool.items() if abs(v - Lreq) <= TOL]
exp_hits = len(pool) * 2 * TOL / spread
best = min(pool.items(), key=lambda kv: abs(kv[1] - Lreq))
check("(a) 池大小 62、跨 40.3 个数量级", len(pool) == 62 and abs(spread - 40.3) < 0.5, "%.2f" % spread)
check("(a) 门槛 0.028%% 下命中 0 个（期望 %.2f）" % exp_hits, hits == [] and exp_hits < 0.1)
check("(a) 最近者 ln C(40,20) = 25.649407，偏离门槛的 28 倍",
      abs(best[1] - 25.649407) < 1e-5 and abs(abs(best[1] - Lreq) / TOL - 28.0) < 0.5,
      "%s=%.6f（%.1f 倍）" % (best[0], best[1], abs(best[1] - Lreq) / TOL))

base = L_from(Fr(41, 10), Fr(-19, 6))
check("(b) 基准 L_req = 25.450914", abs(base - Lreq) < 1e-5, "%.6f" % base)
ok = all(abs(L_from(Fr(41,10)+Fr(d), Fr(-19,6)+Fr(d)) - base) < 1e-9 for d in (Fr(1), Fr(4,3), Fr(-7,2)))
check("(b) b1、b2 同加 δ ⟹ L_req 逐位不变（只依赖差）", ok)
check("(b) 加一整代（Δb1=Δb2=4/3）⟹ L_req 不变",
      abs(L_from(Fr(41,10)+Fr(4,3), Fr(-19,6)+Fr(4,3)) - base) < 1e-9)
d10 = L_from(Fr(41,10)*Fr(11,10), Fr(-19,6)*Fr(11,10))
check("(c) b 差 +10% ⟹ L_req = 23.137194（变化 −2.3137）",
      abs(d10 - 23.137194) < 1e-5 and abs((d10 - base) + 2.3137) < 1e-4, "%.6f" % d10)
sec24 = DOC.split("### 2.4")[1].split("### ")[0] if "### 2.4" in DOC else DOC
check("(d) R1–R5 五条要求齐（无量纲／只经差／精度／不依赖代数／给等耦合能标）",
      all(("R%d" % i) in DOC for i in range(1, 6)) and "0.028" in DOC and "整代" in DOC)

head("F3  锚定")
FILES = ["G57_unreachability_of_absolute_normalization.md", "G60_dimensionless_ledger_and_one_free_unit.md",
         "sm/SMR9_scale_gauge_unification.md", "sm/SMR13_weak_angle_derivation.md",
         "sm/SMR14_low_energy_angle_no_go.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
for rel, kws in [("sm/SMR14_low_energy_angle_no_go.md", ["0.0071"]),
                 ("sm/SMR9_scale_gauge_unification.md", ["109"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

head("F4  与旧文对齐")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 %d" % len(smd))
check("备份包在位", os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))

head("F5  红线自查")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and "**性质**" not in ln and "**纪律**" not in ln
       and "F5 红线自查" not in ln and "L" in ln and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "不可导出", "式肯定句"])]
check("含 L + 导出的行都带否定限定", not bad, "命中 %s" % bad)
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all(("%d. **不许**" % i) in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
