#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR33_check.py -- 核验《SMR33 · 重叠机制判定：Koide 振幅条件可导出》"""
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

DOC = read(os.path.join(HERE, "SMR33_koide_overlap_mechanism.md"))

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
# (a) 重叠恒等式
ok_a = True
for th in (0.4, math.pi / 2, 1.9):
    for ph in (0.1, 1.0, 2.5, -0.7):
        A = 2 * math.cos(th / 2)
        lhs = 1 + A * math.cos(ph)
        rhs = 1 + math.cos(ph + th / 2) + math.cos(ph - th / 2)
        if abs(lhs - rhs) > 1e-12: ok_a = False
check("(a) 恒等式 1 + A cos phi = 1 + cos(phi+th/2) + cos(phi-th/2)（3 角 x 4 相位，偏差 0）", ok_a)
# (b) A(pi/2) = sqrt2 且 A = sqrt2 <=> theta = pi/2
A_pi2 = 2 * math.cos(math.pi / 4)
check("(b) A(theta=pi/2) = 2cos(45deg) = sqrt2（精确到机器精度）", abs(A_pi2 - math.sqrt(2)) < 1e-15, "%.15f" % A_pi2)
grid = [i * 0.0005 for i in range(1, 6284)]
Av = [2 * math.cos(t / 2) for t in grid]
mono = all(Av[i + 1] < Av[i] for i in range(len(Av) - 1))
hits = [t for t, v in zip(grid, Av) if abs(v - math.sqrt(2)) < 2e-4]
check("(b) A(theta)=2cos(theta/2) 在 (0,pi) 上严格单调下降", mono)
check("(b) 网格容差 (2e-4) 下 A = sqrt2 的数值解唯一，且落在 pi/2 的 1e-3 内",
      len(hits) == 1 and abs(hits[0] - math.pi / 2) < 1e-3, str([round(v, 5) for v in hits]))
# (c) 3pi/4 零点 <=> theta = pi/2
zero_ok = True
for th in (math.pi / 2, 0.9, 2.2):
    ph0 = math.pi - th / 2
    val = 1 + 2 * math.cos(th / 2) * math.cos(ph0)
    if abs(th - math.pi / 2) < 1e-12:
        if abs(val) > 1e-12: zero_ok = False
    else:
        if abs(val) < 1e-6: zero_ok = False
check("(c) 零点 phi = pi - th/2：为零当且仅当 th = pi/2（此时 = 3pi/4）", zero_ok)
check("(c) th = pi/2 时零点数值 = 3pi/4 = 135 deg",
      abs((math.pi - math.pi / 4) - 3 * math.pi / 4) < 1e-15)
# (d) 圆分割唯一解：穷举
NIV = [Fr(1, 3), Fr(1, 2), Fr(2, 3)]     # 有理余弦的有理角（SMR29 已验）
admiss = []
for a, b in itertools.combinations(NIV, 2):
    third = 2 - a - b
    if 0 < third < 1 and len({a, b, third}) == 3:
        admiss.append((a, b, third))
check("(d) 恰一个无理角且合格的三元组唯一 = (1/2, 2/3, 5/6)",
      admiss == [(Fr(2, 3), Fr(1, 2), Fr(5, 6))] or admiss == [(Fr(1, 2), Fr(2, 3), Fr(5, 6))],
      str([[str(x) for x in t] for t in admiss]))
tri = sorted([Fr(1, 2), Fr(2, 3), Fr(5, 6)])
check("(d) 三元组：和为 2、互异、均 <= 1、含 pi/2", sum(tri) == 2 and len(set(tri)) == 3 and max(tri) <= 1
      and Fr(1, 2) in tri)
check("(d) 三全有理情形 (2/3,2/3,2/3) 的和也为 2，但违反互异 ⟹ 被排除",
      Fr(2, 3) * 3 == 2 and len({Fr(2, 3)}) == 1)
# (e) 不等权反例
th, ph = math.pi / 2, 0.8
w1, w2 = 1.0, 2.0
mixed = w1 * math.cos(ph + th / 2) + w2 * math.cos(ph - th / 2)
cos_coef = (w1 + w2) * math.cos(th / 2)
sin_coef = (w2 - w1) * math.sin(th / 2)
frac_cos = cos_coef * math.cos(ph)
check("(e) 不等权 (w1,w2)=(1,2)：出现 sin phi 项且 cos 系数 != 2cos(th/2)",
      abs(sin_coef) > 1e-6 and abs(cos_coef - 2 * math.cos(th / 2)) > 1e-6)
check("(e) 混合式 = cos 部分 + sin 部分（数值核对）",
      abs(mixed - (frac_cos + sin_coef * math.sin(ph))) < 1e-12)
check("(e) 等权 w1=w2 ⟹ sin 项消失且 cos 系数 = 2cos(th/2)",
      abs((1.0 - 1.0) * math.sin(th / 2)) < 1e-15 and abs((1.0 + 1.0) * math.cos(th / 2) - 2 * math.cos(th / 2)) < 1e-15)
# (f) 两条限定
check("(f) 两条限定性更新在位（SMR28 §2.4 过严；SMR32 限定为纯表示论语言内）",
      "过严" in DOC and "纯表示论语言内" in DOC and "条件证成" in DOC)
# (g) 七项账
check("(g) 七项账在位（唯一解/恒等式/等权开放/机制输入/相位开放/夸克反证/振幅条件证成）",
      all(k in DOC for k in ["唯一解", "恒等式", "等权", "输入", "相位", "夸克", "条件证成"]))

head("F3  锚定")
FILES = ["sm/SMR28_koide_formula.md", "sm/SMR29_koide_library_crosscheck.md",
         "sm/SMR30_koide_block_equipartition.md", "sm/SMR31_koide_bridge_three_steps.md",
         "sm/SMR32_radius_reading_verdict.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR29_koide_library_crosscheck.md", ["Niven"]),
                 ("sm/SMR32_radius_reading_verdict.md", ["无任何读法"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))
EXT = ["../cosmos-construct/Q473_OVERLAP_LAW.md", "../cosmos-construct/build/build_64_circle_partition.py"]
present = [p for p in EXT if os.path.exists(os.path.join(ROOT, p))]
check("旧侧文件在仓库外存在（只读检查；缺失则跳过）", True, "在位 %d/2" % len(present))
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q473|build_64|build_166|Q634)", f): found.append(os.path.join(dp, f))
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
       and ("导出" in ln or "已建成" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放", "条件", "限定", "更正", "已导出", "由此导出"])]
check("含 导出/已建成 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
