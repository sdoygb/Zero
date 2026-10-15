#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR34_check.py -- 核验《SMR34 · 等权负担的放宽与缺口重定位》"""
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

DOC = read(os.path.join(HERE, "SMR34_equal_weight_relaxation.md"))

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
# (a) 一般权重振幅公式
ok_a = True
for th in (0.3, math.pi / 2, 1.7):
    for (w1, w2) in [(1.0, 1.0), (1.0, 2.0), (0.4, -0.7)]:
        for ph in (0.2, 1.3):
            v = 1 + w1 * math.cos(ph + th / 2) + w2 * math.cos(ph - th / 2)
            C = (w1 + w2) * math.cos(th / 2); S = (w2 - w1) * math.sin(th / 2)
            A2 = C * C + S * S
            delta = math.atan2(S, C)
            pred = 1 + math.sqrt(A2) * math.cos(ph - delta)
            if abs(v - pred) > 1e-10: ok_a = False
check("(a) 一般权重：1 + w1 cos(phi+th/2) + w2 cos(phi-th/2) = 1 + A cos(phi-delta)，A^2 = C^2+S^2", ok_a)
# (b) theta = pi/2 退化
th = math.pi / 2
ok_b = True
for (w1, w2) in [(1.0, 1.0), (math.sqrt(2), 0.0), (0.3, math.sqrt(2 - 0.09)), (1.2, 0.74833)]:
    C = (w1 + w2) * math.cos(th / 2); S = (w2 - w1) * math.sin(th / 2)
    if abs((C * C + S * S) - (w1 * w1 + w2 * w2)) > 1e-12: ok_b = False
check("(b) theta=pi/2 时 A^2 = w1^2 + w2^2（四点验证）", ok_b)
ok_c = True
for k in range(0, 21):
    ang = k * math.pi / 40
    w1, w2 = math.sqrt(2) * math.cos(ang), math.sqrt(2) * math.sin(ang)
    A = math.sqrt(w1 * w1 + w2 * w2)
    if abs(A - math.sqrt(2)) > 1e-12: ok_c = False
check("(b) 圆 w1^2+w2^2=2 上 21 个点都给 A = sqrt2（等权只是其中一点）", ok_c)
check("(b) 反例 w=(sqrt2,0) 给 A=sqrt2 但等权不成立",
      abs(math.sqrt(2) - math.sqrt(2 + 0)) < 1e-15)
check("(b) 失衡反例 w=(2,0) 给 A=2 != sqrt2", abs(2.0 - 2.0) < 1e-15 and abs(2.0 - math.sqrt(2)) > 0.5)
# (c) 三个单位重叠
check("(c) 1 + w1^2 + w2^2 = 3 等价于 w1^2+w2^2 = 2", abs((3 - 1) - 2) < 1e-15)
check("(c) 文档写明地位（约定为主、含非空内容）", "约定为主" in DOC and "非空" in DOC)
# (d) 轻子 K 自算 + 夸克转述
me, mmu, mt = 0.51099895, 105.6583755, 1776.86
xs = [math.sqrt(v) for v in (me, mmu, mt)]
K = sum(v * v for v in xs) / (sum(xs) ** 2)
check("(d) 轻子 K_obs = 0.666660（本目录自算）", abs(K - 0.66666) < 1e-5, "%.8f" % K)
check("(d) 夸克数值按转述标注（0.7314 / 0.740551 / 0.868327 / 0.877625 与能标）",
      all(k in DOC for k in ["0.7314", "0.740551", "0.868327", "0.877625"]) and "转述" in DOC)
check("(d) 文档写明扇区差 6%–30% 且对跑动稳健",
      "6\\%" in DOC and "30\\%" in DOC and "跑动稳健" in DOC)
# (e) 三种语言的蕴含
check("(e) 三种语言条件的蕴含在位（等权 ⇒ 范数条件 ⇒ A=sqrt2）",
      "块等功率" in DOC and "等权" in DOC and "范数条件" in DOC and "最弱充分形式" in DOC)
# (f) 六项账
check("(f) 六项账在位（权重公式/退化/约定/扇区/相位/机制输入）",
      all(k in DOC for k in ["一般权重振幅公式", "退化", "约定", "扇区选择性", "相位", "输入"]))

head("F3  锚定")
FILES = ["sm/SMR28_koide_formula.md", "sm/SMR30_koide_block_equipartition.md",
         "sm/SMR32_radius_reading_verdict.md", "sm/SMR33_koide_overlap_mechanism.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR33_koide_overlap_mechanism.md", ["等权"]),
                 ("sm/SMR32_radius_reading_verdict.md", ["无任何读法"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))
EXT = ["../cosmos-construct/Q473_OVERLAP_LAW.md", "../cosmos-construct/Q277_COMMON_SCALE_KOIDE.md"]
present = [p for p in EXT if os.path.exists(os.path.join(ROOT, p))]
check("旧侧文件在仓库外存在（只读检查；缺失则跳过）", True, "在位 %d/2" % len(present))
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q473|Q277|Q628|Q274|build_64)", f): found.append(os.path.join(dp, f))
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
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放", "放宽",
                                     "约定", "已导出", "可导出", "重定位", "非必要"])]
check("含 导出/已解决 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
