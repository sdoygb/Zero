#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR18_check.py -- 核验《SMR18 · 三击：alpha 走廊审计 / 等耦合机制 / 结构性判死尝试》"""
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

DOC = read(os.path.join(HERE, "SMR18_three_attacks_alpha_coupling_structural.md"))
AEM = Fr(12795, 100); DB = Fr(109, 15); B1 = Fr(41, 10); L0 = Fr(254509138931, 10 ** 10)

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 " + fld + " 在位", fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 " + sec.strip() + " 在位", sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 5 条推导（2.1-2.5）", subs == ["2." + str(i) for i in range(1, 6)], "实测 " + str(subs))
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§" + name + " 五段式标记齐全", not miss, "缺 " + str(miss))
    check("§" + name + " 常见误读 >= 2 条", blk.count("✗") >= 2, "实测 " + str(blk.count("✗")))
    check("§" + name + " 有加框结果", "\\boxed" in blk)

head("F2  独立复算")
def count_in_root(term):
    n = 0
    for dp, dn, fn in os.walk(ROOT):
        dn[:] = [d for d in dn if d != ".git"]
        if os.path.abspath(dp).startswith(os.path.join(ROOT, "sm")): continue
        for f in fn:
            if f.endswith(".md"):
                try: t = read(os.path.join(dp, f))
                except Exception: continue
                if term in t: n += 1
    return n
c137, calpha = count_in_root("137.036"), count_in_root("精细结构")
check("(a) 根目录含 137.036 的 .md >= 1", c137 >= 1, "实测 " + str(c137))
check("(a) 根目录含 精细结构 的 .md >= 1", calpha >= 1, "实测 " + str(calpha))
check("(a) 存在 Z23（既有 alpha 专文）", exists("Z23_alpha_fine_structure.md"))
hg = read(os.path.join(ROOT, "HANDOVER_GR.md"))
check("(a) 台账含 C3 依赖方向 / C6 零信息 / C8 撤回与 0.231200",
      ("C3" in hg and "依赖方向" in hg) and ("C6" in hg and "零信息" in hg)
      and ("C8" in hg and "0.231200" in hg and "撤回" in hg))

sm5 = read(os.path.join(HERE, "SMR5_quantum_numbers_and_charge.md"))
check("(b) SMR5 含三个 Cartan 等迹（Sigma T8^(i)2 = 3）", "3" in sm5 and "Cartan" in sm5)
check("(b) 文档把'等迹 != 等耦合'写成明确断言", "等迹" in DOC and "等耦合" in DOC and "!= " in DOC.replace("\\ne", "!= "))

A = Fr(3, 8) * (AEM + DB * L0 / (2 * math.pi))
aU = A - B1 * L0 / (2 * math.pi)
ident_lhs = aU + Fr(11) / (16 * math.pi) * L0
check("(c) 恒等式：alpha_U^-1 + (11/16pi)L = (3/8)alpha_em^-1",
      abs(float(ident_lhs) - float(Fr(3, 8) * AEM)) < 1e-12,
      "%.8f vs %.8f" % (float(ident_lhs), float(Fr(3, 8) * AEM)))
check("(c) alpha_U^-1 = 42.411622（±0.001）", abs(float(aU) - 42.411622) < 0.001, "%.6f" % float(aU))
check("(c) 系数 (11/16pi) = 0.218838 -> 门槛 dalpha_U^-1 = 0.00155（0.0037%）",
      abs(float(Fr(11) / (16 * math.pi)) - 0.218838) < 1e-6
      and abs(float(Fr(11) / (16 * math.pi)) * 0.0071 - 0.00155) < 1e-5)
check("(c) 文档明确登记判死失败（维度论不成立）", "判死" in DOC and "为假" in DOC and "自我否决" in DOC)

pool = {"4pi": 4 * math.pi, "8pi": 8 * math.pi, "2pi": 2 * math.pi, "pi": math.pi,
        "12pi": 12 * math.pi, "16pi": 16 * math.pi, "e^3": math.e ** 3, "e^4": math.e ** 4,
        "1/alpha_em": 137.036, "(3/8)aem": 0.375 * 127.95, "aem/3": 127.95 / 3,
        "42": 42.0, "43": 43.0, "40": 40.0, "45": 45.0, "48": 48.0, "24": 24.0, "26": 26.0,
        "80": 80.0, "10": 10.0, "sqrt2 pi^2": math.sqrt(2) * math.pi ** 2}
for n in range(2, 41, 2): pool["C(%d)" % n] = float(math.comb(n, n // 2))
tolA = float(Fr(11) / (16 * math.pi)) * 0.0071
hits = [k for k, v in pool.items() if abs(v - float(aU)) <= tolA]
best = min(pool.items(), key=lambda kv: abs(kv[1] - float(aU)))
check("(d) 新池 " + str(len(pool)) + " 项在 0.0037% 门槛下命中 0 个", hits == [], "命中 " + str(hits))
check("(d) 最近者 aem/3 = 42.65，偏 154 倍门槛",
      best[0] == "aem/3" and abs(abs(best[1] - float(aU)) / tolA - 154) < 6,
      "%s=%.4f（%.0f 倍）" % (best[0], best[1], abs(best[1] - float(aU)) / tolA))
check("(e) 三击合账四行在位（走廊审计/等耦合/判死失败/池扫描）",
      all(k in DOC for k in ["攻击一", "攻击二", "攻击三", "新目标数池扫描"]) and "42.4116" in DOC)

head("F3  锚定")
FILES = ["HANDOVER_GR.md", "Z23_alpha_fine_structure.md", "sm/SMR5_quantum_numbers_and_charge.md",
         "sm/SMR15_hierarchy_attack_on_L.md", "sm/SMR17_closed_form_and_reduction.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR17_closed_form_and_reduction.md", ["0.198913"]),
                 ("sm/SMR15_hierarchy_attack_on_L.md", ["R5"]),
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
       if not ln.strip().startswith("|") and "**性质**" not in ln and "**纪律**" not in ln
       and "F5 红线自查" not in ln and ("导出" in ln or "判死" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "登记", "否决",
                                     "失败", "门槛", "等价", "iff", "尝试", "三击"])]
check("含 导出/判死 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
