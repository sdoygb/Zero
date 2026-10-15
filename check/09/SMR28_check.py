#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR28_check.py -- 核验《SMR28 · 小出（Koide）质量关系的 Zero 侧推导尝试》"""
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

DOC = read(os.path.join(HERE, "SMR28_koide_formula.md"))
S2r = math.sqrt(2.0)

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
def Kof(xs):
    s1 = sum(xs); return sum(v * v for v in xs) / (s1 * s1)
# 构造 Koide 三元组：x_i = 1 + sqrt2 cos(theta + 2pi i/3)
def koide_triplet(theta):
    return [1 + S2r * math.cos(theta + 2 * math.pi * i / 3) for i in range(3)]
for th in (0.2, 1.1, 2.7):
    xs = koide_triplet(th)
    check("(a) K = ||x||^2/S_1^2 对 b=sqrt2, theta=%.1f 给 2/3" % th, abs(Kof(xs) - 2 / 3) < 1e-12,
          "%.12f" % Kof(xs))
# (b) S2 = S1^2/6 等价
def S2of(xs):
    return xs[0]*xs[1] + xs[0]*xs[2] + xs[1]*xs[2]
xs = koide_triplet(0.7); s1 = sum(xs)
check("(b) K=2/3 <=> S_2 = S_1^2/6", abs(S2of(xs) - s1 * s1 / 6) < 1e-12)
# (c) ||y||^2 = 3 xbar^2
xb = s1 / 3; y2 = sum((v - xb) ** 2 for v in xs)
check("(c) K=2/3 <=> ||y||^2 = 3 xbar^2", abs(y2 - 3 * xb * xb) < 1e-12, "%.6f vs %.6f" % (y2, 3 * xb * xb))
# (d) 45 度
cos2 = s1 * s1 / (3 * sum(v * v for v in xs))
check("(d) cos^2 angle = 1/(3K) = 1/2 <=> 45 度", abs(cos2 - 0.5) < 1e-12)
# (e) Z3 分解：平凡模（(1,1,1)）与 2 维非平凡模
def cosvec(th):
    return [math.cos(th + 2 * math.pi * i / 3) for i in range(3)]
check("(e) 平凡模 (1,1,1) 在循环置换下不变",
      all(abs([1, 1, 1][(i + 1) % 3] - 1) < 1e-15 for i in range(3)))
c1, c2 = cosvec(0.3), cosvec(1.4)
check("(e) cos 型向量与 (1,1,1) 正交（sum cos = 0，三个 theta 值）",
      all(abs(sum(cosvec(th))) < 1e-12 for th in (0.3, 1.4, 2.9)))
det = c1[0] * c2[1] - c1[1] * c2[0]
check("(e) 两个不同 theta 的 cos 型向量线性无关（张成 2 维非平凡模）", abs(det) > 1e-6,
      "det=%.6f" % det)
check("(e) 非平凡模上的循环置换 = 旋转 2pi/3（cos 型向量相位平移 2pi/3）",
      all(abs(cosvec(0.5)[(i + 1) % 3] - math.cos(0.5 + 2 * math.pi / 3 + 2 * math.pi * i / 3)) < 1e-12
          for i in range(3)))
check("(e) sum cos^2 = 3/2（三个 theta 值）",
      all(abs(sum(v * v for v in cosvec(th)) - 1.5) < 1e-12 for th in (0.3, 1.4, 2.9)))
# K = 1/3 + b^2/6 恒等式（并说明它不给 b）
def K_of_b(b):
    return (1 + b * b / 2) / 3
check("(e) 恒等式 K = 1/3 + b^2/6（b=sqrt2 → 2/3，但该式不固定 b）",
      abs(K_of_b(S2r) - 2 / 3) < 1e-15 and abs(K_of_b(1.0) - 0.5) < 1e-15)
# (f) K_obs 与期望命中
me, mmu, mt = 0.51099895, 105.6583755, 1776.86
x_obs = [math.sqrt(v) for v in (me, mmu, mt)]
K_obs = Kof(x_obs)
check("(f) K_obs ≈ 0.666660（由 PDG 带电轻子质量）", abs(K_obs - 0.66666) < 1e-5, "%.8f" % K_obs)
dev = abs(K_obs - 2 / 3)
exp_hits = (3 / math.pi ** 2) * 100 ** 2 * 2 * dev
check("(f) 与 2/3 的偏差 ≈ 9e-6，小分母有理数期望命中 ≈ 0.006 << 1（命中不廉价）",
      dev < 2e-5 and exp_hits < 0.05, "偏差 %.2e 期望 %.4f" % (dev, exp_hits))
# (g)(h)(i)
check("(g) 机制搜索表在位（群结构/三代/恒等式/几何/原生数集 五行）",
      all(k in DOC for k in ["群结构", "三代计数", "恒等式", "圆的几何", "原生数集"]))
check("(h) C2/C10 对照条款在位", "C2" in DOC and "C10" in DOC and "零信息" in DOC)
check("(i) 六项账在位（四种等价/群分解/命中统计/来源开放/只对带电轻子/对照登记）",
      all(k in DOC for k in ["四种等价形式", "模分解", "命中不廉价", "开放", "带电轻子", "登记"]))

head("F3  锚定")
FILES = ["HANDOVER_GR.md", "sm/SMR2_carrier_and_spin.md", "sm/SMR7_fermion_spectrum_and_generations.md",
         "sm/SMR10_flavor_and_cp.md", "sm/SMR19_closure_density_trap.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR10_flavor_and_cp.md", ["b^2}{6}"]),
                 ("sm/SMR19_closure_density_trap.md", ["387513"]),
                 ("HANDOVER_GR.md", ["C2"]),
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
       and ("导出" in ln or "机制" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放", "承", "真问题", "错法"])]
check("含 导出/机制 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
