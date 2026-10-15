#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR19_check.py -- 核验《SMR19 · 表达闭包的密度陷阱：命中不是证据》"""
import io, math, os, re, sys, time

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

DOC = read(os.path.join(HERE, "SMR19_closure_density_trap.md"))
TGTL, TOLL = 25.4509138931, 0.0071
TGTA = 42.411622
TOLA = (11.0 / (16 * math.pi)) * 0.0071

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
S0 = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0, 26.0, 27.0, 80.0,
      float(math.comb(32, 16)), float(math.comb(36, 18)), float(math.comb(40, 20)),
      float(math.comb(44, 22)), float(math.comb(48, 24)),
      1 / 3, 2 / 3, 1 / 2, 5 / 3, 3 / 8, 41 / 10, 109 / 15, 5 / 9, 7 / 25, 19 / 175,
      23 / 50, 197 / 1225, 313 / 7938,
      math.pi, 2 * math.pi, 4 * math.pi, 8 * math.pi, 16 * math.pi, math.e, math.e ** 3,
      math.pi ** 2 / 6, math.sqrt(2), math.sqrt(3), math.sqrt(5)]
UN = [math.log, math.sqrt, lambda x: 1.0 / x]
BIN = [lambda a, b: a + b, lambda a, b: a - b, lambda a, b: a * b,
       lambda a, b: (a / b if b else None)]
def safe(f, *a):
    try:
        v = f(*a)
        if v != v or v in (float("inf"), float("-inf")) or abs(v) < 1e-8 or abs(v) > 1e8:
            return None
        return v
    except Exception:
        return None
K = lambda v: round(v, 9)
t0 = time.time()
L0 = {K(v): v for v in S0}
L1 = dict(L0)
for v in L0.values():
    for f in UN:
        w = safe(f, v)
        if w is not None: L1[K(w)] = w
for a in L0.values():
    for b in L0.values():
        for f in BIN:
            w = safe(f, a, b)
            if w is not None: L1[K(w)] = w
L2 = dict(L1)
for v in L1.values():
    for f in UN:
        w = safe(f, v)
        if w is not None: L2.setdefault(K(w), w)
for a in list(L0.values()):
    for b in list(L1.values()):
        for f in BIN:
            for w in (safe(f, a, b), safe(f, b, a)):
                if w is not None: L2.setdefault(K(w), w)
el = time.time() - t0
check("(a) 层大小 42 / 3330 / 387513（±2%）",
      abs(len(L0) - 42) <= 1 and abs(len(L1) - 3330) <= 70 and abs(len(L2) - 387513) <= 8000,
      "%d / %d / %d，%.1fs" % (len(L0), len(L1), len(L2), el))
vals = list(L2.values())
hL = [v for v in vals if abs(v - TGTL) <= TOLL]
hA = [v for v in vals if abs(v - TGTA) <= TOLA]
bL = min(vals, key=lambda v: abs(v - TGTL)); bA = min(vals, key=lambda v: abs(v - TGTA))
check("(b) L 命中 28 个（±3）", abs(len(hL) - 28) <= 3, "实测 " + str(len(hL)))
check("(b) alpha_U^-1 命中 1 个", len(hA) == 1, "实测 " + str(len(hA)))
check("(b) L 最接近者偏离 0.02 倍门槛（即比门槛准 50 倍）",
      abs(abs(bL - TGTL) / TOLL - 0.02) < 0.01, "偏离 %.5f = %.3f 倍" % (abs(bL - TGTL), abs(bL - TGTL) / TOLL))
check("(b) alpha_U^-1 最接近者偏离 0.08 倍门槛",
      abs(abs(bA - TGTA) / TOLA - 0.08) < 0.02, "偏离 %.6f = %.3f 倍" % (abs(bA - TGTA), abs(bA - TGTA) / TOLA))

def local_exp(tgt, tol):
    lo, hi = tgt / 2, tgt * 2
    n = sum(1 for v in vals if lo <= v <= hi)
    dens = n / (hi - lo)
    return dens * 2 * tol, n
eL, nL = local_exp(TGTL, TOLL); eA, nA = local_exp(TGTA, TOLA)
spread = max(vals) - min(vals)
naiveL = len(vals) * 2 * TOLL / spread
check("(c) 局部密度期望命中 O(1)-O(10)（L：%.1f，alpha_U^-1：%.2f）" % (eL, eA),
      1.0 <= eL <= 40.0 and 0.05 <= eA <= 20.0, "L 期望 %.1f（区间内 %d 个值）" % (eL, nL))
check("(c) 跨度法错值 << 1（说明跨度法会误判为惊人命中）", naiveL < 1e-3, "%.2e" % naiveL)
check("(c) 实测命中与局部密度期望同量级（28 vs %.1f；1 vs %.2f）" % (eL, eA),
      len(hL) <= 4 * eL + 5 and len(hA) <= 4 * eA + 5)

check("(d) 文档登记受限判死不可用", "受限判死" in DOC and "不可用" in DOC and "自我否决" in DOC)
check("(e) V1 / V2 / V3 三条新判据在位", all(("V" + str(i)) in DOC for i in (1, 2, 3)) and "结构性多对一" in DOC)
check("(f) 战略调整表在位（池扫描更正 / 受限判死 / 维度论 / 外部旧链）",
      all(k in DOC for k in ["池扫描", "受限判死", "维度论", "外部旧链"]))

head("F3  锚定")
FILES = ["sm/SMR14_low_energy_angle_no_go.md", "sm/SMR15_hierarchy_attack_on_L.md",
         "sm/SMR16_lead_audit_scale_unification.md", "sm/SMR17_closed_form_and_reduction.md",
         "sm/SMR18_three_attacks_alpha_coupling_structural.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR15_hierarchy_attack_on_L.md", ["R1"]),
                 ("sm/SMR18_three_attacks_alpha_coupling_structural.md", ["42.4116"]),
                 ("sm/SMR16_lead_audit_scale_unification.md", ["J1"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))

head("F4  与旧文对齐")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 " + str(len(smd)))
check("备份包在位", os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))

head("F5  红线自查")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and not ln.strip().startswith("#") and "**性质**" not in ln and "**纪律**" not in ln
       and "F5 红线自查" not in ln and ("命中" in ln or "判死" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "廉价", "更正",
                                     "登记", "否决", "不可用", "密度", "期望", "适用", "是否", "命中表"])]
check("含 命中/判死 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
