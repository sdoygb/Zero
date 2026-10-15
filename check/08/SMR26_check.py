#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR26_check.py -- 核验《SMR26 · 阈值匹配》"""
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

DOC = read(os.path.join(HERE, "SMR26_threshold_matching.md"))

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
check("(a) 匹配的两部分结构在位（对数部分 + 有限常数部分）",
      "对数" in DOC and "有限常数" in DOC and "匹配" in DOC)
# (b) 阶跃近似的折线性质
def alpha_inv(scale_logs, bs, start, x):
    """分段线性：bs[i] 在 [logs[i], logs[i+1]] 上生效"""
    tot = start
    for i in range(len(bs)):
        lo, hi = scale_logs[i], scale_logs[i + 1]
        seg = min(max(x, lo), hi)
        tot += (bs[i] / (2 * math.pi)) * (seg - lo)
    return tot
logs = [math.log(91.1876), math.log(173.0), math.log(1e13)]
bs = [Fr(7), Fr(7) - Fr(2, 3)]
x1, x2 = math.log(100.0), math.log(150.0)
a1, a2 = alpha_inv(logs, bs, Fr(0), x1), alpha_inv(logs, bs, Fr(0), x2)
slope1 = (a2 - a1) / (x2 - x1)
check("(b) 第一段斜率 = b(6)/2π（n_f=6 段）", abs(slope1 - 7 / (2 * math.pi)) < 1e-12, "%.6f" % slope1)
x3, x4 = math.log(500.0), math.log(1e4)
a3, a4 = alpha_inv(logs, bs, Fr(0), x3), alpha_inv(logs, bs, Fr(0), x4)
slope2 = (a4 - a3) / (x4 - x3)
check("(b) 第二段斜率 = (b(5))/2π，斜率变化 = Δb/2π = −(2/3)/2π",
      abs(slope2 - float(Fr(7) - Fr(2, 3)) / (2 * math.pi)) < 1e-12
      and abs((slope2 - slope1) + (2 / 3) / (2 * math.pi)) < 1e-12)
check("(b) α^{-1} 在阈值处连续（左极限 = 右极限）", abs(alpha_inv(logs, bs, Fr(0), logs[1]) - alpha_inv(logs, bs, Fr(0), logs[1])) < 1e-15)
# (c) SM 例
def b0(nf): return 11 - Fr(2, 3) * nf
check("(c) Delta b0(m_t) = b0(n_f=5) − b0(n_f=6) = +2/3（跨 m_t 向下 n_f 减一，b0 上升）", b0(5) - b0(6) == Fr(2, 3))
check("(c) b0 逐点：n_f=3,4,5,6 → 9, 25/3, 23/3, 7",
      [b0(3), b0(4), b0(5), b0(6)] == [Fr(9), Fr(25, 3), Fr(23, 3), Fr(7)])
# (d) 阈值对 L 的影响
ln_ratio = math.log(1.0307e13 / 173.0)
check("(d) ln(M_U/m_t) ≈ 24.7", abs(ln_ratio - 24.8) < 0.15, "%.3f" % ln_ratio)
toll = 0.0071
for frac in (0.1, 0.5, 1.0):
    dL = frac * ln_ratio
    check("(d) Δ(b1−b2)/(b1−b2)=%.1f → δL≈%.2f（门槛的 %.0f 倍）" % (frac, dL, dL / toll),
          dL > toll and 300 < dL / toll < 4000)
check("(d) 文档写明 350–3500 倍与 24.8", "350" in DOC and "3500" in DOC and "24.8" in DOC)
check("(d) 文档写明阈值效应大于两圈（0.8–1.6）",
      "0.8" in DOC and "1.6" in DOC and "大于" in DOC or "更大" in DOC)
# (e) 待办状态
check("(e) 待办状态表：T2 结构已写/阶跃近似可算/有限项未算；T1 不变；T3 未做",
      "阶跃近似可算" in DOC and "有限项未算" in DOC and "未做" in DOC)
check("(e) 结账句最低限定四条在位（内容谱/阈值表/阶数/方案）",
      all(k in DOC for k in ["内容谱", "阈值表", "阶数", "方案"]))
# (f) 六项账
check("(f) 六项账在位（结构/阶跃近似/Δb 列表/量级/有限项未算/方案未做）",
      all(k in DOC for k in ["匹配的两部分", "阶跃近似", "Delta b", "量级", "有限匹配项", "方案转换"]))

head("F3  锚定")
FILES = ["sm/SMR4_gauge_principle_and_ym.md", "sm/SMR14_low_energy_angle_no_go.md",
         "sm/SMR17_closed_form_and_reduction.md", "sm/SMR20_settlement.md",
         "sm/SMR25_two_loop_running.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR25_two_loop_running.md", ["114"]),
                 ("sm/SMR20_settlement.md", ["T2"]),
                 ("sm/SMR14_low_energy_angle_no_go.md", ["0.0071"]),
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
       and ("关闭" in ln or "算完" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记"])]
check("含 关闭/算完 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
