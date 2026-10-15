#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR17_check.py -- 核验《SMR17 · 精确闭式与归约：L 与 α_em⁻¹ 同源》"""
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

DOC = read(os.path.join(HERE, "SMR17_closed_form_and_reduction.md"))
AEM = Fr(12795, 100); DB = Fr(109, 15); L0 = 25.4509138931; OBS = 0.23122

def G(L):
    A = Fr(3, 8) * (AEM + DB * Fr(L) / (2 * math.pi)); B = AEM - Fr(5, 3) * A
    return 1 / (1 + Fr(5, 3) * (A / B))
def closed(L):
    return Fr(3, 8) - (Fr(5) / (16 * math.pi)) * DB * Fr(L) / AEM

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

head("F2  独立复算")
check("(a) 闭式与 G(L) 在 L_req 处一致（<1e-16）", abs(float(closed(L0)) - float(G(L0))) < 1e-16,
      "%.10f vs %.10f" % (float(closed(L0)), float(G(L0))))
ok = all(abs(float(closed(L)) - float(G(L))) < 1e-15 for L in (10.0, 20.0, 25.45, 30.0, 39.4))
check("(a) 闭式与 G(L) 在 5 个 L 上一致（<1e-15）", ok)
A = Fr(3, 8) * (AEM + DB * Fr(L0) / (2 * math.pi))
check("(a) 恒等式 sin²θ_W = 1 − (5/3)α₁⁻¹/α_em⁻¹",
      abs(float(1 - Fr(5, 3) * A / AEM) - float(G(L0))) < 1e-15)

c = 2 * math.pi * (3 / 8 - OBS) / ((5 / 8) * float(DB))
check("(b) c = 2π(3/8−0.23122)/((5/8)(109/15)) = 0.198913", abs(c - 0.198913) < 5e-7, "%.6f" % c)
check("(b) L = c·α_em⁻¹ 回代：0.198913×127.95 = 25.4509", abs(c * 127.95 - L0) < 5e-4,
      "%.6f" % (c * 127.95))
check("(b) L/α_em⁻¹ = 1/5.0273", abs(L0 / 127.95 - 1 / 5.0273) < 1e-6, "%.6f" % (L0 / 127.95))

dA = 0.0071 / c
check("(c) 门槛转移：δL=0.0071 ⟹ δα_em⁻¹ = 0.0357（0.028%）",
      abs(dA - 0.0357) < 5e-4 and abs(100 * dA / 127.95 - 0.028) < 5e-3, "%.5f" % dA)
check("(c) 灵敏度：α 移 0.1 ⟹ L 移 0.0199", abs(c * 0.1 - 0.0198913) < 1e-6)

hg = read(os.path.join(ROOT, "HANDOVER_GR.md"))
check("(d) 台账含 C3（R(r) 依赖方向反）", "C3" in hg and "依赖方向" in hg)
check("(d) 台账含 C6（零信息）", "C6" in hg and "零信息" in hg)
check("(d) 台账含 C8（0.231200 已被 build_83 撤回）", "C8" in hg and "0.231200" in hg and "撤回" in hg)
check("(d) 0.231200 与观测值的偏差 = −0.0088%（落在 0.028% 门槛内）",
      abs((0.231200 - OBS) / OBS + 0.000088) < 5e-6,
      "%.6f%%" % (100 * (0.231200 - OBS) / OBS))
check("(d) 台账写明撤回原因是 β 归一化错（b_Y 须为 (5/3)b_1）", "归一化" in hg and "b_Y" in hg)
check("(e) R1–R5 在归约后的新形式在位（R1..R5 全出现）",
      all(("R%d" % i) in DOC for i in range(1, 6)) and "0.0357" in DOC)

head("F3  锚定")
FILES = ["HANDOVER_GR.md", "sm/SMR13_weak_angle_derivation.md", "sm/SMR14_low_energy_angle_no_go.md",
         "sm/SMR15_hierarchy_attack_on_L.md", "sm/SMR16_lead_audit_scale_unification.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
for rel, kws in [("sm/SMR14_low_energy_angle_no_go.md", ["0.0071"]),
                 ("sm/SMR16_lead_audit_scale_unification.md", ["J1"]),
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
       and "F5 红线自查" not in ln and ("导出" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "归约", "转述", "登记", "门槛", "iff", "\\iff", "等价"])]
check("含 导出 的行都带限定（不／未／归约／登记／门槛）", not bad, "命中 %s" % bad)
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all(("%d. **不许**" % i) in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
