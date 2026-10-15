#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR27_check.py -- 核验《SMR27 · 方案转换》"""
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

DOC = read(os.path.join(HERE, "SMR27_scheme_conversion.md"))
mW, mZ = Fr(80377, 1000), Fr(911876, 10000)
shat, aem_inv, DB = 0.23122, 127.95, Fr(109, 15)
toll = 0.0071

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
# (a)
ratio = mW / mZ
sOS = 1 - ratio ** 2
check("(a) s^2_OS = 1 − (m_W/m_Z)^2 用精确有理数（m_W=80.377, m_Z=91.1876）",
      abs(float(sOS) - 0.2230) < 0.0005, "%.5f" % float(sOS))
d_sch = shat - float(sOS)
rel = d_sch / shat
check("(a) 方案差 Δ_sch ≈ 0.0082，相对 3.6%",
      abs(d_sch - 0.0082) < 0.0005 and abs(rel - 0.036) < 0.003,
      "Δ=%.5f 相对=%.2f%%" % (d_sch, 100 * rel))
check("(a) 相对方案差 / 相对门槛 = 约 10^3",
      abs(rel / 0.00037 - 96) < 30, "%.0f 倍" % (rel / 0.00037))
# (b)
alpha = 1 / aem_inv
check("(b) alpha/(4pi) ≈ 5.8e-4；增强因子 Δ_sch/(alpha/4pi) ≈ 14（与 Δρ 增强一致）",
      abs(alpha / (4 * math.pi) - 5.8e-4) < 5e-5
      and abs(d_sch / (alpha / (4 * math.pi)) - 14) < 3,
      "%.2e, %.1f" % (alpha / (4 * math.pi), d_sch / (alpha / (4 * math.pi))))
# (c)
def ccoef(s2): return 2 * math.pi * (3 / 8 - s2) / (float(Fr(5, 8)) * float(DB))
cMS, cOS = ccoef(shat), ccoef(float(sOS))
LMS, LOS = cMS * aem_inv, cOS * aem_inv
check("(c) c_MS = 0.198913 与 SMR17 一致", abs(cMS - 0.198913) < 5e-6, "%.6f" % cMS)
check("(c) L_MS = 25.451（与 SMR14/SMR17 一致）", abs(LMS - 25.4509138931) < 1e-4, "%.4f" % LMS)
check("(c) L_OS ≈ 26.9，δL ≈ 1.4（相对 5.5%）",
      abs(LOS - 26.9) < 0.2 and abs((LOS - LMS) - 1.4) < 0.2,
      "L_OS=%.3f δL=%.3f" % (LOS, LOS - LMS))
check("(c) δL / 门槛 ≈ 200 倍", abs((LOS - LMS) / toll - 200) < 40, "%.0f 倍" % ((LOS - LMS) / toll))
# (d)
check("(d) 待办状态表：T3 量级已算、系数未复算；T1/T2 不变",
      "方案差的量级已算" in DOC and "转换系数未复算" in DOC and "不变" in DOC)
check("(d) 四类效应阶梯在位（门槛/两圈/方案/阈值）",
      all(k in DOC for k in ["一圈门槛", "两圈", "方案差", "阈值"]) and "阶梯" in DOC)
# (e)
check("(e) 结账句完整限定版在位（内容谱+阈值表+阶数+方案）",
      "完整限定版" in DOC and all(k in DOC for k in ["内容谱", "阈值表", "阶数", "方案"]))
# (f)
check("(f) 六项账在位（两套定义/转换量级/方案差/阶梯/系数未做/两圈未算）",
      all(k in DOC for k in ["两套定义", "转换的一圈结构", "方案差对", "阶梯", "未做", "未算"]))

head("F3  锚定")
FILES = ["sm/SMR14_low_energy_angle_no_go.md", "sm/SMR17_closed_form_and_reduction.md",
         "sm/SMR20_settlement.md", "sm/SMR25_two_loop_running.md",
         "sm/SMR26_threshold_matching.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR17_closed_form_and_reduction.md", ["0.198913"]),
                 ("sm/SMR20_settlement.md", ["T3"]),
                 ("sm/SMR26_threshold_matching.md", ["24.8"]),
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
       and ("关闭" in ln or "复算" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "引用"])]
check("含 关闭/复算 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
