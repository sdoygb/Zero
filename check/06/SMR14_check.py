#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR14_check.py -- 核验《SMR14 · 低能弱混合角：0.23122 等价于一个数（L）》

F1 结构；F2 独立复算：(a) G(L) 与 SMR9 的 x1 对上 (b) 二分反解 L_req (c) 灵敏度/δL
(d) L 的尺度比属性 (e) 候选尺度表 (f) 0.028% 门槛下的命中统计与最佳有理逼近
F3 锚定；F4 与旧文对齐（走备份包）；F5 红线自查
"""
import io, itertools, math, os, re, sys
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

DOC = read(os.path.join(HERE, "SMR14_low_energy_angle_no_go.md"))
B1, B2 = Fr(41, 10), Fr(-19, 6)
DB = B1 - B2                       # 109/15
AEM_INV = 127.95
OBS, OBS_ERR = 0.23122, 4e-5

def G(L):
    A = Fr(3, 8) * (Fr(AEM_INV) + DB * Fr(L) / (2 * math.pi))
    B = Fr(AEM_INV) - Fr(5, 3) * A
    return 1 / (1 + Fr(5, 3) * (A / B)), float(A), float(B)

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 6 条推导（2.1-2.6）", subs == ["2.%d" % i for i in range(1, 7)], "实测 %s" % subs)
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

head("F2  独立复算")
L0 = 25.4509138931
s0, A0, B0 = G(L0)
check("(a) b1-b2 = 109/15 = 7.266667", DB == Fr(109, 15) and abs(float(DB) - 7.266667) < 1e-6)
check("(a) G(25.4509138931) 给 α1⁻¹ = 59.019241（逐位对上 SMR9 的 x1）",
      abs(A0 - 59.019241) < 1e-5, "%.6f" % A0)
check("(a) 同一 L 给 α2⁻¹ = 29.584599 与 sin²θ_W = 0.231220",
      abs(B0 - 29.584599) < 1e-5 and abs(float(s0) - OBS) < 5e-7, "%.6f %.6f" % (B0, float(s0)))

lo, hi = 20.0, 32.0
for _ in range(200):
    mid = (lo + hi) / 2
    if float(G(mid)[0]) > OBS: lo = mid
    else: hi = mid
Lreq = (lo + hi) / 2
check("(b) 二分反解：0.23122 ⟺ L_req = 25.4509138931（10 位）",
      abs(Lreq - L0) < 5e-10, "%.10f" % Lreq)
check("(b) 回代 G(L_req) = 0.23122（4×10⁻⁵ 内）", abs(float(G(Lreq)[0]) - OBS) < 1e-9)
check("(b) M_U,req = 91.1876·e^L = 1.030688e13 GeV",
      abs(91.1876 * math.exp(Lreq) - 1.0306880687e13) / 1.0306880687e13 < 1e-6)

d = float(G(Lreq + 1)[0]) - float(G(Lreq)[0])
dL = OBS_ERR / abs(d)
check("(c) dsin²/dL ≈ −0.005649", abs(d + 0.005649) < 5e-5, "%.6f" % d)
check("(c) 实验精度 4e-5 ⟹ δL ≈ 0.0071（0.028%）",
      abs(dL - 0.0071) < 5e-4 and abs(100 * dL / Lreq - 0.028) < 5e-3, "%.5f" % dL)

ls = [ln for ln in DOC.splitlines() if "ln(M_U/M_Z)" in ln or "尺度比" in ln]
check("(d) 文档把 L 写成尺度比并归入绝对尺度账本", len(ls) >= 2 and "绝对尺度" in DOC)

CAND = [("M_Pl(1.22e19)", math.log(1.22e19 / 91.1876), 0.152219),
        ("10^16", math.log(1e16 / 91.1876), 0.192367),
        ("10^6", math.log(1e6 / 91.1876), 0.322447)]
ok = True
for nm, L, exp in CAND:
    got = float(G(L)[0])
    if abs(got - exp) > 5e-5: ok = False; print("     %s: 实测 %.6f 期望 %.6f" % (nm, got, exp))
check("(e) 候选尺度表：M_Pl→0.152219、10^16→0.192367、10^6→0.322447", ok)
check("(e) 三个候选与观测值的差都 > 0.03（远离）",
      all(abs(float(G(L)[0]) - OBS) > 0.03 for _, L, _ in CAND))

POOL = [Fr(1), Fr(2), Fr(4), Fr(10), Fr(26), Fr(80), Fr(5, 9), Fr(7, 25), Fr(1, 3), Fr(5, 3),
        Fr(3, 8), Fr(1, 6), Fr(2, 3), Fr(1, 2), Fr(3, 4), Fr(1, 4), Fr(4, 5), Fr(2, 5), Fr(3, 5),
        Fr(1, 8), Fr(7, 8), Fr(9, 4), Fr(1, 9), Fr(8, 9), Fr(19, 175), Fr(23, 50), Fr(5, 4),
        Fr(3, 2), Fr(6, 5), Fr(11, 6)]
ratios = set()
for x, y in itertools.permutations(POOL, 2):
    if y != 0: ratios.add(x / y)
tol = 0.00071 * OBS                     # 0.028% 门槛（对应 δL≈0.0071）
hits = [r for r in ratios if abs(float(r) - OBS) <= tol]
exp_hits = len(ratios) * 2 * tol
check("(f) 0.028%% 门槛下：比值集 %d 个，命中 %d 个（期望 %.2f ⟹ <1 才说明门槛有区分力）"
      % (len(ratios), len(hits), exp_hits),
      exp_hits < 1.0, "%.3f" % exp_hits)
rhits = []
for q in range(1, 1001):
    pp = round(OBS * q)
    if pp > 0 and abs(pp / q - OBS) <= OBS_ERR: rhits.append((pp, q))
best = min((abs(pp / q - OBS), pp, q) for pp, q in rhits)
check("(f) 更正：分母 ≤1000 落在 ±4e-5 内的有理数有 %d 个（期望≈24）⟹ 有理命中廉价、不算证据"
      % len(rhits), len(rhits) >= 5,
      "最佳 %d/%d 差 %.2e（%.3fσ）" % (best[1], best[2], best[0], best[0] / OBS_ERR))
check("(f) 3/13 距观测值 4.5×10⁻⁴ ⟹ 约 11σ，被排除",
      abs(3 / 13 - OBS) > 1e-4, "%.2e" % abs(3 / 13 - OBS))

head("F3  锚定")
FILES = ["G57_unreachability_of_absolute_normalization.md",
         "G60_dimensionless_ledger_and_one_free_unit.md",
         "sm/SMR4_gauge_principle_and_ym.md", "sm/SMR9_scale_gauge_unification.md",
         "sm/SMR13_weak_angle_derivation.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
for rel, kws in [("sm/SMR9_scale_gauge_unification.md", ["25.4509138931", "127.95"]),
                 ("sm/SMR13_weak_angle_derivation.md", ["3/8"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

head("F4  与旧文对齐（旧文已删，走备份包）")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 %d" % len(smd))
sm12 = legacy("SM12_scale_weak_angle_unification.md")
check("备份包内 SM12 可读，且含 0.23122／25.4509138931（可核对）",
      len(sm12) > 0 and "0.23122" in sm12 and "25.4509138931" in sm12)

head("F5  红线自查")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and "**性质**" not in ln and "**纪律**" not in ln
       and "F5 红线自查" not in ln
       and ("0.23122" in ln or "预言" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "反解", "不可导出", "式肯定句", "iff", "\\iff"])]
check("含 0.23122/预言 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all(("%d. **不许**" % i) in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))
check("文件名不被旧正则 ^SM<digit>_ 捕获", not re.match(r"^SM\d+_", "SMR14_low_energy_angle_no_go.md"))

print(""); print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
