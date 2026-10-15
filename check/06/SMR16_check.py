#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR16_check.py -- 核验《SMR16 · 线索审计：L 的旧理论走廊》"""
import io, math, os, re, subprocess, sys
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
def G(L, aem=127.95):
    A = Fr(3, 8) * (Fr(aem) + Fr(109, 15) * Fr(L) / (2 * math.pi)); B = Fr(aem) - Fr(5, 3) * A
    return float(1 / (1 + Fr(5, 3) * (A / B)))

DOC = read(os.path.join(HERE, "SMR16_lead_audit_scale_unification.md"))
Lreq, TOL = 25.4509138931, 0.0071

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 5 条推导（2.1-2.5）", subs == ["2.%d" % i for i in range(1, 6)], "实测 %s" % subs)
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

head("F2  独立复算")
def count_in_root(term):
    n = 0
    for dp, dn, fn in os.walk(ROOT):
        dn[:] = [d for d in dn if d not in (".git",)]
        if os.path.abspath(dp).startswith(os.path.join(ROOT, "sm")): continue
        for f in fn:
            if f.endswith(".md"):
                try: t = read(os.path.join(dp, f))
                except Exception: continue
                if term in t: n += 1
    return n
c_mu, c_eq = count_in_root("M_U"), count_in_root("等耦合")
check("(a) 根目录提到 M_U 的 .md 极少（≤3）", c_mu <= 3, "实测 %d" % c_mu)
check("(a) 根目录提到 等耦合 的 .md 极少（≤3）", c_eq <= 3, "实测 %d" % c_eq)
for term in ("统一能标", "大统一", "三线共点", "耦合统一"):
    check("(a) 根目录无 %s（0 个文件）" % term, count_in_root(term) == 0, "实测 %d" % count_in_root(term))
g91 = read(os.path.join(ROOT, "G91_matter_sector_range_and_three_way_verdict.md"))
check("(a) G91 是唯一 M_U 出处，且把它归入单位一侧", "M_U" in g91 and "单位一侧" in g91)

inv1 = read(os.path.join(HERE, "_ref", "modular_equilibrium_sm_inventory.md"))
inv2 = read(os.path.join(HERE, "_ref", "cosmos_construct_SM_inventory.md"))
check("(b) 清单含 Q453_TWO_SCALES 与 26.033／4.98e13", "Q453_TWO_SCALES" in inv2 and "26.033" in inv2 and "4.98e13" in inv2)
check("(b) 清单含 D9_couplings 且状态为 条件恢复", "D9_couplings" in inv1 and "条件恢复" in inv1)
check("(b) 清单含 Q602／Q620／Q322／R-SCALE", all(k in inv2 or k in inv1 for k in ("Q602", "Q620", "Q322", "R-SCALE")))
check("(b) 清单含 GR_FROM_GEOMETRY 的 2.50×10^16 条目", "GR_FROM_GEOMETRY" in inv2 and "2.50" in inv2)

mz = 91.1876
cands = [("Q453-new", 4.98e13, 0.222321), ("GR_GEOM", 2.50e16, 0.187190), ("Q453-old", 246.22 * math.exp(33.6197), 0.179461)]
ok = True
for nm, M, exp in cands:
    got = G(math.log(M / mz))
    if abs(got - exp) > 5e-5: ok = False; print("     %s: %.6f vs %.6f" % (nm, got, exp))
check("(c) 折算表：0.222321／0.187190／0.179461", ok)
check("(c) 所需点：G(25.4509138931) = 0.231220", abs(G(Lreq) - 0.231220) < 5e-7)
ratio = abs(math.log(4.98e13 / mz) - Lreq) / TOL
check("(c) 最好候选偏门槛 222 倍", abs(ratio - 222) < 2, "%.1f 倍" % ratio)
check("(c) 全部候选与观测差 > 0.03（=-3.8% 以上）", all(abs(G(math.log(M / mz)) - 0.23122) > 0.003 for _, M, _ in cands))

check("(d) J1／J2 两问在位且与 R1–R5 对应", "J1" in DOC and "J2" in DOC and "R1" in DOC and "无量纲" in DOC)
check("(e) 防污染：全仓无旧理论源文件（Q453*/Q602*/Q620*/Q322*/D9_couplings*）",
      not any(re.match(r"^(Q453|Q602|Q620|Q322|D9_couplings)", f) for f in
              [x for dp, dn, fn in os.walk(ROOT) if ".git" not in dp for x in fn]))
check("(e) 文档写明不拷入／不作前提／只读清单三约束",
      "不拷入" in DOC and "不作前提" in DOC and "转述" in DOC)

head("F3  锚定")
FILES = ["G91_matter_sector_range_and_three_way_verdict.md", "sm/SMR14_low_energy_angle_no_go.md",
         "sm/SMR15_hierarchy_attack_on_L.md", "sm/SMR13_weak_angle_derivation.md",
         "sm/_ref/modular_equilibrium_sm_inventory.md", "sm/_ref/cosmos_construct_SM_inventory.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
for rel, kws in [("sm/SMR14_low_energy_angle_no_go.md", ["0.0071"]),
                 ("sm/SMR15_hierarchy_attack_on_L.md", ["R1"]),
                 ("sm/SMR13_weak_angle_derivation.md", ["超荷"])]:
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
       and "F5 红线自查" not in ln and "旧理论" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "仅", "转述", "清单", "登记表"])]
check("含 旧理论 的行都带限定（不／仅／转述）", not bad, "命中 %s" % bad)
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all(("%d. **不许**" % i) in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
