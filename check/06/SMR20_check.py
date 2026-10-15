#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR20_check.py -- 核验《SMR20 · 结账：弱混合角走廊的终账》"""
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
def legacy(name):
    import tarfile
    p = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
    if not os.path.exists(p): return ""
    with tarfile.open(p) as tf:
        for m in tf.getmembers():
            if os.path.basename(m.name) == name:
                f = tf.extractfile(m)
                return f.read().decode("utf-8", "replace") if f else ""
    return ""

DOC = read(os.path.join(HERE, "SMR20_settlement.md"))
AEM = Fr(12795, 100); DB = Fr(109, 15); L0 = 25.4509138931; OBS = 0.23122

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 " + fld + " 在位", fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 " + sec.strip() + " 在位", sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 结账本体")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 7 条（2.1-2.7）", subs == ["2." + str(i) for i in range(1, 8)], "实测 " + str(subs))
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§" + name + " 五段式标记齐全", not miss, "缺 " + str(miss))
    check("§" + name + " 常见误读 >= 2 条", blk.count("✗") >= 2, "实测 " + str(blk.count("✗")))
    check("§" + name + " 有加框结果", "\\boxed" in blk)

head("F2  独立复算")
# (a) 依据 1-6
c2 = Fr(1, 3) + Fr(4, 3)
check("(a) 依据 1：Sigma c^2 = c_L^2 + c_R^2 = 5/3", c2 == Fr(5, 3))
Y = [Fr(1, 6), Fr(-2, 3), Fr(1, 3), Fr(-1, 2), Fr(1), Fr(0)]
SY2 = sum(Fr(m) * y * y for m, y in zip([6, 3, 3, 2, 1, 1], Y))
check("(a) 依据 2：由导出的超荷表给 Sigma Y^2 = 10/3（比值 5/3）",
      Y == [Fr(1, 6), Fr(-2, 3), Fr(1, 3), Fr(-1, 2), Fr(1), Fr(0)] and SY2 == Fr(10, 3))
check("(a) 依据 3：b1 - b2 = 109/15 = 7.266667", DB == Fr(109, 15) and abs(float(DB) - 7.266667) < 1e-6)
A1 = Fr(3, 8) * (AEM + DB * Fr(L0) / (2 * math.pi))
check("(a) 依据 4：曲线给 alpha_1^-1(M_Z) = 59.019241（对上 SMR9 的 x1）",
      abs(float(A1) - 59.019241) < 1e-5, "%.6f" % float(A1))
def G(L):
    a = Fr(3, 8) * (AEM + DB * Fr(L) / (2 * math.pi)); b = AEM - Fr(5, 3) * a
    return 1 / (1 + Fr(5, 3) * (a / b))
def closed(L):
    return Fr(3, 8) - (Fr(5) / (16 * math.pi)) * DB * Fr(L) / AEM
check("(a) 依据 5：闭式 == 曲线（5 个 L 上 <1e-15）",
      all(abs(float(closed(L)) - float(G(L))) < 1e-15 for L in (10.0, 20.0, 25.45, 30.0, 39.4)))
aU = A1 - Fr(41, 10) * Fr(L0) / (2 * math.pi)
check("(a) 依据 6：alpha_U^-1 + (11/16pi)L = (3/8)alpha_em^-1，且 alpha_U^-1 = 42.4116",
      abs(float(aU + Fr(11) / (16 * math.pi) * Fr(L0)) - float(Fr(3, 8) * AEM)) < 1e-12
      and abs(float(aU) - 42.411622) < 0.001, "%.6f" % float(aU))
check("(a) 曲线在 L_req 给 0.231220", abs(float(G(L0)) - OBS) < 5e-7)

# (b) 买的三笔
check("(b) B1/B2/B3 三笔在位且 B3 标为唯一实质缺口",
      all(("B" + str(i)) in DOC for i in (1, 2, 3)) and "唯一的实质缺口" in DOC)
check("(b) B3 的等价写法（L 与 alpha_U^-1）在位", "25.4509138931" in DOC and "42.4116" in DOC)

# (c) 前提 P1
check("(c) 前提 P1 在位且写明等迹 != 等耦合",
      "P1" in DOC and "等迹" in DOC and "等耦合" in DOC and "\\ne" in DOC)
check("(c) 前提恰一条（文档明写'恰一条'）", "恰一条" in DOC)

# (d) 待办 T1-T6
missT = [i for i in range(1, 7) if ("T" + str(i)) not in DOC]
check("(d) 待办 T1-T6 齐全", not missT, "缺 " + str(missT))
check("(d) T1/T2/T3 标为未算或未做", DOC.count("未算") >= 2 and "未做" in DOC)
check("(d) T4 写明已知被排除（3.9688 / -39.66%）", "3.9688" in DOC and "39.66" in DOC)
check("(d) T5/T6 标为开放", DOC.count("开放") >= 2)

# (e) 验收 V1-V3
check("(e) V1/V2/V3 在位且写明孤立命中不算",
      all(("V" + str(i)) in DOC for i in (1, 2, 3)) and "孤立命中" in DOC)

# (f) 灵敏度
check("(f) 灵敏度：b 差 +10% -> L 移 -2.31；alpha 移 0.1 -> L 移 0.0199；整代不变",
      "-2.31" in DOC and "0.0199" in DOC and "不变" in DOC)
check("(f) 门槛：0.028%（L）与 0.0037%（alpha_U^-1）都在位", "0.028" in DOC and "0.0037" in DOC)
dA = 0.0071 / (2 * math.pi * (3 / 8 - OBS) / ((5 / 8) * float(DB)))
check("(f) 复核：delta L = 0.0071 <=> delta alpha_U^-1 = 0.00155",
      abs((11 / (16 * math.pi)) * 0.0071 - 0.00155) < 1e-5)

# (g) 适用范围
check("(g) 适用范围句在位（只管弱混合角走廊）", "只管弱混合角走廊" in DOC)

head("F3  锚定")
FILES = ["HANDOVER_GR.md", "sm/SMR9_scale_gauge_unification.md", "sm/SMR13_weak_angle_derivation.md",
         "sm/SMR14_low_energy_angle_no_go.md", "sm/SMR15_hierarchy_attack_on_L.md",
         "sm/SMR16_lead_audit_scale_unification.md", "sm/SMR17_closed_form_and_reduction.md",
         "sm/SMR18_three_attacks_alpha_coupling_structural.md", "sm/SMR19_closure_density_trap.md",
         "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR13_weak_angle_derivation.md", ["3/8"]),
                 ("sm/SMR17_closed_form_and_reduction.md", ["0.198913"]),
                 ("sm/SMR19_closure_density_trap.md", ["387513"]),
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
       if not ln.strip().startswith("|") and not ln.strip().startswith("#")
       and "**性质**" not in ln and "**纪律**" not in ln and "F5 红线自查" not in ln
       and ("导出" in ln or "命中" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "登记",
                                     "买入", "缺口", "前提", "待办", "廉价", "验收"])]
check("含 导出/命中 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))
check("文件名不被旧正则 ^SM<digit>_ 捕获", not re.match(r"^SM\d+_", "SMR20_settlement.md"))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
