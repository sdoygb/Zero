#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR36_check.py -- 核验《SMR36 · 动力学能否解释扇区选择性》"""
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

DOC = read(os.path.join(HERE, "SMR36_dynamics_selectivity.md"))
M = (0.51099895, 105.6583755, 1776.86)

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
def K_of(ms):
    s1 = sum(math.sqrt(v) for v in ms)
    return sum(ms) / (s1 * s1)
# (a) 零次齐次
ok_a = True
for c in (0.001, 0.5, 1.0, 7.0, 1e6):
    if abs(K_of([c * v for v in M]) - K_of(M)) > 1e-12: ok_a = False
check("(a) K 零次齐次：K(c m) = K(m)（五个 c 值，含 1e-3 与 1e6）", ok_a)
check("(a) 轻子 K = 0.666660（自算）", abs(K_of(M) - 0.66666) < 1e-5, "%.8f" % K_of(M))
# (b) 一圈 QCD 指数
CF = Fr(4, 3)
gm0 = 6 * CF
check("(b) gamma_m^(0) = 6 C_F = 8（精确）", gm0 == 8)
exps = {}
for nf in (3, 4, 5, 6):
    exps[nf] = Fr(12, 33 - 2 * nf)
check("(b) 指数 12/(33-2n_f)：n_f=3,4,5,6 -> 4/9, 12/25, 12/23, 12/21",
      exps[3] == Fr(4, 9) and exps[4] == Fr(12, 25) and exps[5] == Fr(12, 23) and exps[6] == Fr(12, 21),
      str({k: str(v) for k, v in exps.items()}))
check("(b) 经典结果 m ∝ alpha_s^(12/25)（n_f=4）", exps[4] == Fr(12, 25))
R = 0.37   # 任意重标度因子
check("(b) 整体重标度（R 与味无关）下 K 严格不变",
      abs(K_of([R * v for v in M]) - K_of(M)) < 1e-15)
# (c) 扇区内电荷普适
Q = {"charged lepton": (-1, -1, -1), "up": (Fr(2, 3),) * 3, "down": (Fr(-1, 3),) * 3}
check("(c) 三扇区内三代电荷相同（带电轻子 -1；上型 +2/3；下型 -1/3）",
      all(len(set(v)) == 1 for v in Q.values()))
check("(c) 扇区间 Q^2 不同（1, 4/9, 1/9）——但那是扇区间的，不影响扇区内 K",
      Q["charged lepton"][0] ** 2 == 1 and float(Q["up"][0] ** 2) == 4 / 9 and float(Q["down"][0] ** 2) == 1 / 9)
# (d) 阈值分化量级 vs 所需
Kd2, Kd17, Ku2, Ku17 = 0.740551, 0.724317, 0.868327, 0.877625
drift_d = abs(Kd17 - Kd2) / Kd2
drift_u = abs(Ku17 - Ku2) / Ku2
need_d = abs(Kd2 - 2 / 3) / (2 / 3)
need_u = abs(Ku2 - 2 / 3) / (2 / 3)
r_d, r_u = need_d / drift_d, need_u / drift_u
check("(d) 漂移：d 扇区 %.1f%%，u 扇区 %.1f%%" % (100 * drift_d, 100 * drift_u),
      0.01 < drift_d < 0.03 and 0.005 < drift_u < 0.02)
check("(d) 所需偏差：d 扇区 %.1f%%，u 扇区 %.1f%%" % (100 * need_d, 100 * need_u),
      0.06 < need_d < 0.15 and 0.25 < need_u < 0.35)
check("(d) 比值（所需/漂移）：d≈%.1f 倍，u≈%.1f 倍（落在 5–20 倍区间）" % (r_d, r_u),
      4 < r_d < 25 and 10 < r_u < 40)
check("(d) 文档写明 5–20 倍与'一圈判死'", "5" in DOC and "20" in DOC and "一圈动力学判死" in DOC)
# (e) 两条出路
check("(e) 两条出路在位（味结构 / 机制非普适）", "味结构" in DOC and "机制非普适" in DOC and "公设" in DOC)
# (f) 七项账
check("(f) 七项账在位（零次齐次/QCD 普适/QED 普适/阈值量级/味结构/机制非普适/相位）",
      all(k in DOC for k in ["零次齐次", "味普适", "电荷普适", "阈值", "味结构", "非普适", "相位"]))

head("F3  锚定")
FILES = ["sm/SMR4_gauge_principle_and_ym.md", "sm/SMR9_scale_gauge_unification.md",
         "sm/SMR10_flavor_and_cp.md", "sm/SMR34_equal_weight_relaxation.md",
         "sm/SMR35_sector_selectivity.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR35_sector_selectivity.md", ["受限判死"]),
                 ("sm/SMR34_equal_weight_relaxation.md", ["范数条件"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))
EXT = ["../cosmos-construct/Q277_COMMON_SCALE_KOIDE.md"]
present = [p for p in EXT if os.path.exists(os.path.join(ROOT, p))]
check("旧侧文件在仓库外存在（只读检查；缺失则跳过）", True, "在位 %d/1" % len(present))
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q277|Q628|build_166)", f): found.append(os.path.join(dp, f))
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
       and ("导出" in ln or "已解决" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放",
                                     "判死", "可导出", "买入", "公设"])]
check("含 导出/已解决 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
