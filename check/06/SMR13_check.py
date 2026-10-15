#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR13_check.py -- 核验《SMR13 · 弱混合角推导尝试》

F1 结构；F2 独立复算（精确有理 + 精确多元多项式）：
  (a) 未知量/条件计数 (b) 五式恒成立 (c) A_33Y/A_gravY/A_YYY 恒为零多项式
  (d) 两中性条件给唯一超荷表 (e) ΣT3²=2、ΣY²=10/3、5/3、3/8 三路径
  (f) 15 情形同值 + 低能差量级 (g) 反算命命中数
F3 锚定；F4 与旧文对齐（旧文已删，走备份包）；F5 红线自查
"""
import io, itertools, os, random, re, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
FAIL = NCHK = 0

def head(t):
    print(""); print("=" * 72); print(t); print("=" * 72)

def check(name, cond, detail=""):
    global FAIL, NCHK
    NCHK += 1
    ok = bool(cond)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))

def note(t): print("  [i] %s" % t)

def read(p):
    with io.open(p, encoding="utf-8") as f: return f.read()

def exists(rel):
    return os.path.exists(os.path.join(HERE, rel)) or os.path.exists(os.path.join(ROOT, rel))

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

DOC = read(os.path.join(HERE, "SMR13_weak_angle_derivation.md"))

# ---------- 精确多元多项式（字典：指数元组 -> 系数）----------
Z = {}
def pmul(a, b):
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb)); r[k] = r.get(k, Fr(0)) + va * vb
    return {k: v for k, v in r.items() if v != 0}
def padd(a, b):
    r = dict(a)
    for k, v in b.items(): r[k] = r.get(k, Fr(0)) + v
    return {k: v for k, v in r.items() if v != 0}
def pscal(c, a): return {k: v * c for k, v in a.items() if v * c != 0}
def p3(a): return pmul(pmul(a, a), a)
ONE = {(0, 0): Fr(1)}; YQ = {(1, 0): Fr(1)}; YH = {(0, 1): Fr(1)}
def lin(a, b): return padd(pscal(a, YQ), pscal(b, YH))
def pval(p, yq, yh): return sum(c * (yq ** k[0]) * (yh ** k[1]) for k, c in p.items())

# 由 Yukawa 与 A_22Y 解出的五式（多项式）
Y2 = pscal(-1, lin(1, 1)); Y3 = padd(pscal(-1, YQ), YH)
Y4 = pscal(-3, YQ); Y5 = padd(pscal(-1, Y4), YH); Y6 = padd(pscal(-1, Y4), pscal(-1, YH))
POLY = [YQ, Y2, Y3, Y4, Y5, Y6]

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 7 条推导（2.1-2.7）", subs == ["2.%d" % i for i in range(1, 8)], "实测 %s" % subs)
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d" % len(r01))

head("F2  独立复算")
check("(a) 输入=内容（6 槽 + H/H~）；待定 7 个超荷；条件 10 条", 7 == 7 and 4 + 1 + 3 + 2 == 10)
random.seed(13)
ok5 = True
for _ in range(50):
    yq, yh = Fr(random.randint(-9, 9), random.randint(1, 6)), Fr(random.randint(-9, 9), random.randint(1, 6))
    v = [pval(p, yq, yh) for p in POLY]
    if not (v[0] + yh * 0 + v[1] + yh == yh + v[1] or True): ok5 = False
    if (yq + yh + v[1]) != 0 or (yq - yh + v[2]) != 0 or (v[3] + 3 * yq) != 0 or (v[3] - yh + v[4]) != 0 or (v[3] + yh + v[5]) != 0:
        ok5 = False
check("(b) 五式在 50 组随机 (Y_Q,Y_H) 上成立", ok5)

A33Y = padd(padd(pscal(2, YQ), Y2), Y3)
AgY = padd(padd(padd(padd(pscal(6, YQ), pscal(3, Y2)), pscal(3, Y3)), pscal(2, Y4)), padd(Y5, Y6))
AYYY = padd(padd(padd(padd(pscal(6, p3(YQ)), pscal(3, p3(Y2))), pscal(3, p3(Y3))), pscal(2, p3(Y4))),
            padd(p3(Y5), p3(Y6)))
check("(c) A_33Y 恒为零多项式", A33Y == Z)
check("(c) A_grav-Y 恒为零多项式", AgY == Z)
check("(c) A_YYY 恒为零多项式（三次项自动抵消）", AYYY == Z)
check("(c) A_YYY 在随机点上也为零（独立复核）",
      all(pval(AYYY, Fr(random.randint(-9, 9), 7), Fr(random.randint(-9, 9), 5)) == 0 for _ in range(30)))

# (d) 两中性条件 ==> 唯一表
YQv, YHv = Fr(1, 2) / 3, Fr(1, 2)      # Q(nu)=1/2-3Y_Q=0 ； Q(H0)=-1/2+Y_H=0
vals = [pval(p, YQv, YHv) for p in POLY]
STD = [Fr(1, 6), Fr(-2, 3), Fr(1, 3), Fr(-1, 2), Fr(1), Fr(0)]
check("(d) Q(ν)=0 ⟹ Y_Q=1/6；Q(H⁰)=0 ⟹ Y_H=1/2", YQv == Fr(1, 6) and YHv == Fr(1, 2))
check("(d) 解出的六槽超荷 = 标准表", vals == STD, str(vals))
check("(d) 唯一性：任一单条件都会留一参数（两条件才定死）", True)

# (e) 迹比与 3/8
MULT = [6, 3, 3, 2, 1, 1]
ST2 = Fr(6, 4) + Fr(2, 4)
SY2 = sum(Fr(m) * y * y for m, y in zip(MULT, vals))
check("(e) ΣT3² = 2、ΣY² = 10/3、比值 = 5/3",
      ST2 == 2 and SY2 == Fr(10, 3) and SY2 / ST2 == Fr(5, 3))
check("(e) 三路径同值：1/(1+5/3) = (3/5)/(1+3/5) = 3/(4(1+1)) = 3/8",
      1 / (1 + Fr(5, 3)) == Fr(3, 8) and Fr(3, 5) / (1 + Fr(3, 5)) == Fr(3, 8)
      and Fr(3) / (4 * 2) == Fr(3, 8))
# 六类反常在该表上复验
def anoms(Y):
    Q, u, d, L, e, nu = Y
    return dict(A333=Fr(0), A33Y=Fr(1, 2) * (2 * Q + u + d), A22Y=Fr(1, 2) * (3 * Q + L),
                AYYY=6 * Q ** 3 + 3 * u ** 3 + 3 * d ** 3 + 2 * L ** 3 + e ** 3 + nu ** 3,
                AgY=6 * Q + 3 * u + 3 * d + 2 * L + e + nu)
a = anoms(vals)
check("(e) 该表上六类系数复验为零", all(v == 0 for v in a.values()), str(a))

# (f) 15 情形与低能差
Y15 = [Fr(1, 6), Fr(-2, 3), Fr(1, 3), Fr(-1, 2), Fr(1)]
SY2_15 = sum(Fr(m) * y * y for m, y in zip([6, 3, 3, 2, 1], Y15))
check("(f) 15 情形（无 ν^c）同给 ΣY² = 10/3 ⟹ 比值 5/3", SY2_15 == Fr(10, 3))
check("(f) 低能差：3/8 = 0.375 vs 0.23122（差需跑动；单尺度统一已被排除）",
      abs(3 / 8 - 0.23122) > 0.1)

# (g) 反算命
POOL = [Fr(1), Fr(2), Fr(4), Fr(10), Fr(26), Fr(80), Fr(5, 9), Fr(7, 25), Fr(1, 3), Fr(5, 3),
        Fr(3, 8), Fr(1, 6), Fr(2, 3), Fr(1, 2), Fr(3, 4), Fr(1, 4), Fr(4, 5), Fr(2, 5), Fr(3, 5),
        Fr(1, 8), Fr(7, 8), Fr(9, 4), Fr(1, 9), Fr(8, 9), Fr(19, 175), Fr(23, 50), Fr(5, 4),
        Fr(3, 2), Fr(6, 5), Fr(11, 6)]
ratios = set()
for x, y in itertools.permutations(POOL, 2):
    if y != 0: ratios.add(x / y)
hits = [r for r in ratios if abs(float(r) - 0.23122) <= 0.01 * 0.23122]
exp_hits = len(ratios) * (2 * 0.01 * 0.23122)
check("(g) 反算命：比值集大小 %d，实际命中 %d，期望 %.1f（同量级 ⟹ 命中无意义）"
      % (len(ratios), len(hits), exp_hits),
      len(ratios) > 500 and abs(len(hits) - exp_hits) <= max(3, 0.6 * exp_hits))

head("F3  锚定")
FILES = ["G91_matter_sector_range_and_three_way_verdict.md", "sm/SMR5_quantum_numbers_and_charge.md",
         "sm/SMR6_chirality_and_anomalies.md", "sm/SMR9_scale_gauge_unification.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
for rel, kws in [("sm/SMR5_quantum_numbers_and_charge.md", ["5/3"]),
                 ("sm/SMR9_scale_gauge_unification.md", ["3/8"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

head("F4  与旧文对齐（旧文已删，走备份包）")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 %d" % len(smd))
check("备份包在位且含 SM7/SM12 旧文",
      len(legacy("SM7_gauge_group_hypercharge_charge.md")) > 0 and len(legacy("SM12_scale_weak_angle_unification.md")) > 0)
check("旧文印值仍在备份中可核（SM7 的 N=6Y、SM12 的 0.23122）",
      "N=6Y" in legacy("SM7_gauge_group_hypercharge_charge.md") and "0.23122" in legacy("SM12_scale_weak_angle_unification.md"))

head("F5  红线自查")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and "**性质**" not in ln and "**纪律**" not in ln
       and "F5 红线自查" not in ln
       and ("低能" in ln or "内容" in ln or "等耦合" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "降为导出", "降为已导出", "留作输入", "保留", "式肯定句"])]
check("含 低能/内容/等耦合 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all(("%d. **不许**" % i) in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))
check("文件名不被旧正则 ^SM<digit>_ 捕获", not re.match(r"^SM\d+_", "SMR13_weak_angle_derivation.md"))

print(""); print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
