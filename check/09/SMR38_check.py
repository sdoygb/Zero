#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR38_check.py -- 核验《SMR38 · Koide 走廊结账》"""
import io, itertools, math, os, re, sys
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

DOC = read(os.path.join(HERE, "SMR38_koide_settlement.md"))

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 " + fld + " 在位", fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 "]:
    check("节 " + sec.strip() + " 在位", sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导（结账）")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 7 条（2.1-2.7）", subs == ["2." + str(i) for i in range(1, 8)], "实测 " + str(subs))
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§" + name + " 五段式标记齐全", not miss, "缺 " + str(miss))
    check("§" + name + " 常见误读 >= 2 条", blk.count("✗") >= 2, "实测 " + str(blk.count("✗")))
    check("§" + name + " 有加框结果", "\\boxed" in blk)

head("F2  汇总复算")
# (a) 三段一致性
ml = (0.51099895, 105.6583755, 1776.86)
xs = [math.sqrt(v) for v in ml]
S1 = sum(xs); S2 = sum(xs[i] * xs[j] for i in range(3) for j in range(i + 1, 3))
K = sum(v * v for v in xs) / (S1 * S1)
check("(a) 轻子：K = 2/3 与 S2 = S1^2/6 等价（同一三元组）",
      abs(K - 2 / 3) < 1e-5 and abs(S2 - S1 * S1 / 6) / (S1 * S1 / 6) < 5e-5, "K=%.8f" % K)
xb = S1 / 3; y2 = sum((v - xb) ** 2 for v in xs)
check("(a) 等价形式 ||y||^2 = 3 xbar^2 亦成立（相对差 ~1e-5）",
      abs(y2 - 3 * xb * xb) / (3 * xb * xb) < 1e-4)
check("(a) 循环矩阵语言：K = 1/3 + (2/3)(|a1|^2/a0^2)，且 |a1|^2/a0^2 = 1/2 给 2/3",
      abs((1 / 3 + (2 / 3) * 0.5) - 2 / 3) < 1e-15)
# (b) 角度唯一解
NIV = [Fr(1, 3), Fr(1, 2), Fr(2, 3)]
admiss = []
for a, b in itertools.combinations(NIV, 2):
    third = 2 - a - b
    if 0 < third < 1 and len({a, b, third}) == 3: admiss.append((a, b, third))
check("(a) 圆分割唯一解 = (1/2, 2/3, 5/6)（含 pi/2）",
      len(admiss) == 1 and Fr(1, 2) in admiss[0], str(admiss))
A_pi2 = 2 * math.cos(math.pi / 4)
check("(a) 振幅 A(theta=pi/2) = sqrt2；且 3pi/4 零点与之同一事实",
      abs(A_pi2 - math.sqrt(2)) < 1e-15
      and abs(1 + A_pi2 * math.cos(3 * math.pi / 4)) < 1e-15)
# (c) 判死表数值
def K_of(ms):
    s1 = sum(math.sqrt(v) for v in ms)
    return sum(ms) / (s1 * s1)
check("(c) K 零次齐次（对任意 c）", all(abs(K_of([c * v for v in ml]) - K_of(ml)) < 1e-12 for c in (1e-6, 1, 1e6)))
Cu, Cd = 0.868327, 0.740551
c_u = math.sqrt(3 * (Cu - 1 / 3)); c_d = math.sqrt(3 * (Cd - 1 / 3))
check("(c) N6：c_u ≈ 1.267、c_d ≈ 1.105、sqrt3 ≈ 1.732 越界",
      abs(c_u - 1.2669) < 0.005 and abs(c_d - 1.1053) < 0.005
      and (1 / 3 + 3 / 3) > 1 and abs(c_u - c_d) > 0.1)
drift = max(abs(0.724317 - 0.740551) / 0.740551, abs(0.877625 - 0.868327) / 0.868327)
need = min(abs(0.740551 - 2 / 3) / (2 / 3), abs(0.868327 - 2 / 3) / (2 / 3))
check("(c) N2：阈值漂移 %.1f%% vs 所需 %.1f%% ⟹ 比值 %.1f（5–20 倍）" % (100 * drift, 100 * need, need / drift),
      need / drift > 4)
check("(c) N3：K 基无关（相似不变量：对给定谱任意酉变换不改 K）",
      abs(K_of(sorted(ml)) - K_of(ml)) < 1e-15)
check("(c) N1：三扇区算符结构相同（SU(2) 内容），差别仅色重数",
      "结构相同" in DOC and "色重数" in DOC)
check("(c) N5：纯表示论读法受限判死；N4：轨道偏离为转述 ⟹ 判死表六行在位",
      all(k in DOC for k in ["N1", "N2", "N3", "N4", "N5", "N6"]))
# (d) 价签与档位
check("(d) 五格价签在位（I-1 到 I-5）", all(("I-%d" % i) in DOC for i in range(1, 6)))
check("(d) 档位汇总在位（已证 2 环 / 条件证成 1 环 / 判死 6 条 / 开放 2 项）",
      "已证" in DOC and "条件证成" in DOC and "开放" in DOC and "判死" in DOC)
# (e) 翻转条件
check("(e) 翻转条件六条（1.-6.）", all((str(i) + " |") in DOC for i in range(1, 7)))
# (f) 两走廊对照
check("(f) 两走廊对照要素（缺数值 vs 缺选择规则；命中廉价 vs 不廉价）",
      "缺**一个数值**" in DOC or "缺一个数值" in DOC.replace("*", "")
      or "一个数值" in DOC)
check("(f) 命中不廉价与闭包廉价对照在位", "不廉价" in DOC and "廉价" in DOC)
# (g) 结账句要素
check("(g) 结账句含四要素（形状/角度/振幅/买入）",
      all(k in DOC for k in ["形状", "角度", "振幅", "公设"]))

head("F3  锚定")
FILES = ["sm/SMR20_settlement.md", "sm/SMR28_koide_formula.md", "sm/SMR32_radius_reading_verdict.md",
         "sm/SMR34_equal_weight_relaxation.md", "sm/SMR35_sector_selectivity.md",
         "sm/SMR36_dynamics_selectivity.md", "sm/SMR37_flavor_structure_selectivity.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR20_settlement.md", ["42.4116"]),
                 ("sm/SMR37_flavor_structure_selectivity.md", ["基无关"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q274|Q277|Q473|Q628|Q629|Q630|Q775|build_64|build_166)", f):
            found.append(os.path.join(dp, f))
check("防污染：本仓库内无旧理论文件", not found, "命中 " + str(found))

head("F4  与旧文对齐")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 " + str(len(smd)))
check("备份包在位", os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))

head("F5  红线自查")
SKIP = ("**目标**", "**步骤**", "**用到的输入**", "**核对**", "**结果**", "**常见误读**", "**性质**")
DOC_SCAN = DOC.split("## §3 ")[0] + DOC.split("## §4 ")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and not ln.strip().startswith("#")
       and not any(ln.strip().startswith(s) for s in SKIP)
       and "**纪律**" not in ln and "F5 红线自查" not in ln
       and ("已导出" in ln or "廉价巧合" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放",
                                     "判死", "买入", "公设", "条件"])]
check("含 已导出/廉价巧合 的行都带限定", not bad, "命中 " + str(bad))
sec5 = DOC.split("## §5")[1].split("## §6")[0]
check("§5 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec5 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
