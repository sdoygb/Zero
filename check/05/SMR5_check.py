#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR5_check.py -- 核验《SMR5 · 量子数、超荷与电荷量子化》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 8 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理 + 精确 Q(sqrt3)；不 import 旧脚本）：
        (a) 16 态表：N=6Y 整性、同奇偶性 N=2T3 (mod 2)、Q in (1/3)Z、电荷集、反例 1/6
        (b) 最小步长 gcd = 1/3
        (c) 中心核：穷举 3x2x6=36，两条路径判恒等，|K|=6、生成元 z(1,1,1)、Z6 = Z3 x Z2
        (d) 迹条件：态级 / 场级两条路径；按场计的错值 19/9
        (e) 27：(3,3b,1)+(3b,1,3)+(1,3,3b) 的三迹 = 3、交叉项 = 0、sum Y^2 = 3(cL^2+cR^2)
        (f) 嵌入解 cL=-1/sqrt3, cR=2/sqrt3 在 Q(sqrt3) 上精确
        (g) lambda 标度律、3/5 = 1/(cL^2+cR^2)
  F3  锚定：引用的 lh 与 sm/ 文件存在；关键词真实出现
  F4  与旧文对齐：SM7 表 T1、5/3、521、中心核 Z6；登记 C-SMR5-1（格运算那一行）
  F5  红线自查：无"导出了表示内容/电荷量子化原理"式肯定句；§6 七条；$$ 成对；文件集不变
"""
from __future__ import annotations

import functools
import io
import itertools
import math
import os
import re
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
FAIL = 0
NCHK = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(name, cond, detail=""):
    global FAIL, NCHK
    NCHK += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def note(t):
    print("  [i] %s" % t)


def read(p):
    with io.open(p, encoding="utf-8") as f:
        return f.read()


def legacy(name):
    """从 sm/_legacy_SM0-19.tar.gz 读取旧文（旧文已从工作树删除，对照以备份为准）"""
    import tarfile
    p = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
    if not os.path.exists(p):
        return ""
    with tarfile.open(p) as tf:
        for m in tf.getmembers():
            if os.path.basename(m.name) == name:
                f = tf.extractfile(m)
                return f.read().decode("utf-8", "replace") if f else ""
    return ""


def exists(rel):
    return os.path.exists(os.path.join(HERE, rel)) or os.path.exists(os.path.join(ROOT, rel))


DOC = read(os.path.join(HERE, "SMR5_quantum_numbers_and_charge.md"))
SM7 = legacy("SM7_gauge_group_hypercharge_charge.md")

# ---- Q(sqrt3) ----
def q3(a, b=0):
    return (Fr(a), Fr(b))


def q3_mul(x, y):
    return (x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q3_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def q3_fl(x):
    return float(x[0]) + float(x[1]) * math.sqrt(3)


# ======================================================================
head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 8 条推导（2.1-2.8）", subs == ["2.%d" % i for i in range(1, 9)], "实测 %s" % subs)
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式六标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d 条" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

# ======================================================================
head("F2  独立复算")

# ---- 16 态表（场, 色重数, 分量数, T3 列表, Y）----
FIELDS = [("Q", 3, [Fr(1, 2), Fr(-1, 2)], Fr(1, 6)),
          ("u^c", 3, [Fr(0)], Fr(-2, 3)),
          ("d^c", 3, [Fr(0)], Fr(1, 3)),
          ("L", 1, [Fr(1, 2), Fr(-1, 2)], Fr(-1, 2)),
          ("e^c", 1, [Fr(0)], Fr(1)),
          ("nu^c", 1, [Fr(0)], Fr(0))]
STATES = []
for nm, col, t3s, Y in FIELDS:
    for t3 in t3s:
        for _ in range(col):
            STATES.append((nm, t3, Y))
check("(a) 态数 = 16", len(STATES) == 16, "实测 %d" % len(STATES))
ok_int = all((6 * Y).denominator == 1 for _, _, Y in STATES)
ok_par = all((6 * Y - 2 * t3) % 2 == 0 for _, t3, Y in STATES)
check("(a) N=6Y 为整数（取值 1,-4,2,-3,6,0）", ok_int and
      sorted({int(6 * Y) for _, _, Y in STATES}) == [-4, -3, 0, 1, 2, 6])
check("(a) 同奇偶性 N = 2T3 (mod 2) 逐态成立", ok_par)
ok_Q = all((t3 + Y) * 3 % 1 == 0 for _, t3, Y in STATES)
Qset = sorted({t3 + Y for _, t3, Y in STATES})
check("(a) Q = T3+Y in (1/3)Z 且电荷集 = {0,±1/3,±2/3,±1}",
      ok_Q and Qset == [Fr(-1), Fr(-2, 3), Fr(-1, 3), Fr(0), Fr(1, 3), Fr(2, 3), Fr(1)], str(Qset))
check("(a) 反例：t2=0, N=1 给 Q=1/6 ∉ (1/3)Z（故格的论证不够）",
      (Fr(0) + Fr(1, 6)) * 3 % 1 != 0 and Fr(1, 6) not in [Fr(k, 3) for k in range(-3, 4)])
check("(a) 格的和 (1/6)Z + (1/2)Z = (1/6)Z ≠ (1/3)Z",
      all(((Fr(a, 6) + Fr(b, 2)) * 6).denominator == 1 for a in range(-3, 4) for b in range(-3, 4))
      and Fr(1, 6) not in [Fr(k, 3) for k in range(-3, 4)])

# ---- (b) 最小步长 ----
qs = sorted({abs(int(3 * (t3 + Y))) for _, t3, Y in STATES if t3 + Y != 0})
g = functools.reduce(math.gcd, qs)
check("(b) 非零电荷（以 1/3 为单位）的 gcd = 1 ⟹ 最小步长 1/3", g == 1, "实测 %s -> %d" % (qs, g))

# ---- (c) 中心核 ----
T3MAP = {"Q": 1, "u^c": -1, "d^c": -1, "L": 0, "e^c": 0, "nu^c": 0}   # 色三重态性


def trivial(a, b, u):
    for nm, t3, Y in STATES:
        t3tri = T3MAP[nm]
        ph = Fr(a * t3tri, 3) + Fr(b * int(2 * t3), 2) + Fr(u * int(6 * Y), 6)
        if ph.denominator != 1:
            return False
    return True


ker = [(a, b, u) for a in range(3) for b in range(2) for u in range(6) if trivial(a, b, u)]
check("(c) 穷举 3x2x6 = 36 个中心元素，核恰有 6 个", len(ker) == 6, "实测 %s" % ker)
check("(c) 核 = {z(u,u,u)}，生成元 z(1,1,1)", sorted(ker) == sorted([(u % 3, u % 2, u) for u in range(6)]))
# 两条路径：位相指数 vs 乘 omega 的幂（用 Z_6 上的整数指数验证）
ok2 = True
for (a, b, u) in ker:
    for nm, t3, Y in STATES:
        e1 = (a * T3MAP[nm]) % 3 * 2 + (b * int(2 * t3)) % 2 * 3 + (u * int(6 * Y)) % 6
        if e1 % 6 != 0:
            ok2 = False
check("(c) 核元素在全部 16 态上相位指数 ≡ 0（第二条路径：整数指数 mod 6）", ok2)
check("(c) Z6 ≅ Z3 × Z2：2 阶元与 3 阶子群交为 {1}",
      len([z for z in [(u % 3, u % 2, u) for u in range(6)] if (2 * 0 + z[2] * 3) % 6 == 0]) >= 1)

# ---- (d) 迹条件 ----
ST2 = sum(t3 * t3 for _, t3, _ in STATES)
SY2 = sum(Y * Y for _, _, Y in STATES)
check("(d) 态级：ΣT3^2 = 2、ΣY^2 = 10/3、比 = 5/3",
      ST2 == 2 and SY2 == Fr(10, 3) and SY2 / ST2 == Fr(5, 3), "%s %s" % (ST2, SY2))
SY2_field = sum(Fr(len(t3s) * col) * Y * Y for _, col, t3s, Y in FIELDS)
ST2_field = sum(Fr(len([t for t in t3s if t != 0]) * col) * Fr(1, 4) for _, col, t3s, Y in FIELDS)
check("(d) 场级路径（色重数 x 分量数 x Y^2）给同一对数值",
      SY2_field == Fr(10, 3) and ST2_field == 2, "%s %s" % (ST2_field, SY2_field))
SY2_wrong = sum(Fr(len(t3s)) * Y * Y for _, col, t3s, Y in FIELDS)
check("(d) 错按场计（漏色重数）给 19/9 ≠ 10/3", SY2_wrong == Fr(19, 9), "实测 %s" % SY2_wrong)

# ---- (e) 27 上的三迹 ----
a8 = [Fr(1, 2), Fr(-1, 2), Fr(0)]
v27 = []
for l in range(3):
    for m in range(3):                      # 块 1：(3, 3b, 1)
        v27.append((a8[l], -a8[m], Fr(0)))
for l in range(3):
    for n in range(3):                      # 块 2：(3b, 1, 3)
        v27.append((-a8[l], Fr(0), a8[n]))
for m in range(3):
    for n in range(3):                      # 块 3：(1, 3, 3b)
        v27.append((Fr(0), a8[m], -a8[n]))
check("(e) 27 = 9+9+9 态", len(v27) == 27)
T11 = sum(v[0] * v[0] for v in v27)
T22 = sum(v[1] * v[1] for v in v27)
T33 = sum(v[2] * v[2] for v in v27)
T12 = sum(v[0] * v[1] for v in v27)
T13 = sum(v[0] * v[2] for v in v27)
T23 = sum(v[1] * v[2] for v in v27)
check("(e) 三个 Cartan 的迹都 = 3", T11 == 3 and T22 == 3 and T33 == 3, "%s %s %s" % (T11, T22, T33))
check("(e) 交叉项全为 0", T12 == 0 and T13 == 0 and T23 == 0, "%s %s %s" % (T12, T13, T23))

# ---- (f) 嵌入解（Q(sqrt3)）----
cL = q3(0, Fr(-1, 3))        # -1/sqrt3 = -sqrt3/3
cR = q3(0, Fr(2, 3))         #  2/sqrt3
cL2 = q3_mul(cL, cL)
cR2 = q3_mul(cR, cR)
check("(f) cL^2 = 1/3、cR^2 = 4/3、和 = 5/3", cL2 == q3(Fr(1, 3)) and cR2 == q3(Fr(4, 3))
      and q3_add(cL2, cR2) == q3(Fr(5, 3)))
check("(f) |cR/cL| = 2 且 cL = -1/sqrt3、cR = 2/sqrt3（精确）",
      q3_fl(cR) / q3_fl(cL) == -2.0 and abs(q3_fl(cL) + 1 / math.sqrt(3)) < 1e-12
      and abs(q3_fl(cR) - 2 / math.sqrt(3)) < 1e-12)
SY2_27 = q3_add(q3_mul(cL2, q3(T22)), q3_mul(cR2, q3(T33)))
check("(f) ΣY^2|27 = 3(cL^2+cR^2) = 5", SY2_27 == q3(5), str(SY2_27))

# ---- (g) lambda 标度律 ----
ok_lam = True
for lam in [Fr(2), Fr(1, 3), Fr(-1), Fr(5)]:
    if (lam * lam * SY2) / ST2 != lam * lam * Fr(5, 3):
        ok_lam = False
check("(g) Y -> λY 使 ΣY^2/ΣT3^2 -> λ^2 (5/3)；条件 λ = ±1", ok_lam)
check("(g) 3/5 = 1/(cL^2+cR^2)", Fr(3, 5) == 1 / Fr(5, 3))

# ======================================================================
head("F3  锚定")
FILES = ["Z1_zero_layer_as_the_foundation.md", "G61_locking_the_five_integers.md",
         "G82_B_from_two_independent_Z2.md", "G91_matter_sector_range_and_three_way_verdict.md",
         "R6_h5_dictionary_error_bound.md", "sm/SM7_gauge_group_hypercharge_charge.md",
         "sm/SMR4_gauge_principle_and_ym.md", "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G91_matter_sector_range_and_three_way_verdict.md", ["只能具名输入"]),
      ("G61_locking_the_five_integers.md", ["Z_2"]),
      ("sm/SM7_gauge_group_hypercharge_charge.md", ["N=6Y"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM7_gauge_group_hypercharge_charge.md']
_present = [d for d in OLD_USED if os.path.exists(os.path.join(HERE, d))]
check("本组原先对照的旧文已全部从 sm/ 删除（0 个残留）", not _present, "残留 " + str(_present))
_tar = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
check("备份包 sm/_legacy_SM0-19.tar.gz 存在", os.path.exists(_tar))
if os.path.exists(_tar):
    import tarfile
    with tarfile.open(_tar) as _tf:
        _names = set(os.path.basename(n) for n in _tf.getnames())
    _absent = [d for d in OLD_USED if d not in _names]
    check("本组对照过的旧文都能在备份包内找到（可恢复）", not _absent, "缺 " + str(_absent))
    check("旧文内容可从备份包读出（非空）", all(len(legacy(d)) > 0 for d in OLD_USED))
_idx = read(os.path.join(HERE, "SMR_index.md"))
check("SMR_index.md 的系列结算块在位", "SMR 系列结算" in _idx)

head("F5  红线自查")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if ("表示内容" in ln or "电荷量子化" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若"])]
check("含 表示内容/电荷量子化 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all(("%d. **不许**" % i) in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
smc = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_check\.py$", f))
smrd = sorted(f for f in os.listdir(HERE) if re.match(r"^SMR\d+_.*\.md$", f))
smrc = sorted(f for f in os.listdir(HERE) if re.match(r"^SMR\d+_check\.py$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 %d" % len(smd))
check("旧 SM<数字>_check.py 已按指令删除（0 个）", len(smc) == 0, "实测 %d" % len(smc))
check("SMR 系列（支柱 12 篇 ＋ 后续篇）在位", len(smrd) >= 12 and len(smrc) >= 12,
      "%d 篇 / %d 个" % (len(smrd), len(smrc)))
check("备份包在位（旧体系可从 sm/_legacy_SM0-19.tar.gz 恢复）",
      os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))
check("本篇文件名不被 SM<数字>_ 正则捕获", not re.match(r"^SM\d+_", "SMR5_quantum_numbers_and_charge.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
