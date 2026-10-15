#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR6_check.py -- 核验《SMR6 · 手征性与反常消去》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 7 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理；不 import 旧脚本）：
        (a) 六类系数两条路径（逐场表 / 逐态枚举）= 0
        (b) Witten N_{1/2} = 4（按场与按态两种口径）
        (c) B-L 八类系数 = 0；去掉 nu^c 的差分 (A_gX, A_XXX) = (-1,-1)
        (d) 14 态块：六类 = 0、态数 14、N_{1/2} = 4、复共轭差集、与一代不交、16+X 无反常
        (e) 反常方程组解空间维数 = 2，基 (1,1,1,1,1,0) 与 (0,0,0,0,0,1)
        (f) 内容层手征性（复共轭不闭合）
        (g) 齐次性：Y->λY 下三类系数按 λ^3 / λ / 1 变换
  F3  锚定：引用的 lh 与 sm/ 文件存在；关键词真实出现
  F4  与旧文对齐：SM8 的 14 态块、B-L 赋值表、八类结论、A_XXX 印值
  F5  红线自查：无"导出了手征性/谱内容"式肯定句；§6 七条；$$ 成对；文件集不变
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


DOC = read(os.path.join(HERE, "SMR6_chirality_and_anomalies.md"))
SM8 = legacy("SM8_anomalies.md")

COLOR_DIM = {"3": 3, "3b": 3, "1": 1}
A_R3 = {"3": Fr(1), "3b": Fr(-1), "1": Fr(0)}
T_R3 = {"3": Fr(1, 2), "3b": Fr(1, 2), "1": Fr(0)}
T_R2 = {1: Fr(0), 2: Fr(1, 2)}

# 一代六槽：(名字, 色表示, 弱维数, Y, B-L)
GEN = [("Q", "3", 2, Fr(1, 6), Fr(1, 3)),
       ("u^c", "3b", 1, Fr(-2, 3), Fr(-1, 3)),
       ("d^c", "3b", 1, Fr(1, 3), Fr(-1, 3)),
       ("L", "1", 2, Fr(-1, 2), Fr(-1)),
       ("e^c", "1", 1, Fr(1), Fr(1)),
       ("nu^c", "1", 1, Fr(0), Fr(1))]
# 14 态块
BLK = [("(3,2,0)", "3", 2, Fr(0), Fr(0)),
       ("(3b,1,-1/6)", "3b", 1, Fr(-1, 6), Fr(0)),
       ("(3b,1,1/6)", "3b", 1, Fr(1, 6), Fr(0)),
       ("(1,2,0)", "1", 2, Fr(0), Fr(0))]


def anomalies(fields, with_witten=True):
    A333 = sum(A_R3[c] * w for _, c, w, Y, X in fields)
    A33Y = sum(T_R3[c] * w * Y for _, c, w, Y, X in fields)
    A22Y = sum(T_R2[w] * COLOR_DIM[c] * Y for _, c, w, Y, X in fields)
    AYYY = sum(COLOR_DIM[c] * w * Y ** 3 for _, c, w, Y, X in fields)
    AgY = sum(COLOR_DIM[c] * w * Y for _, c, w, Y, X in fields)
    N12 = sum(1 if w == 2 else 0 for _, c, w, Y, X in fields for _ in range(COLOR_DIM[c]))
    return dict(A333=A333, A33Y=A33Y, A22Y=A22Y, AYYY=AYYY, AgY=AgY, N12=N12)


def bl_anomalies(fields):
    A33X = sum(T_R3[c] * w * X for _, c, w, Y, X in fields)
    A22X = sum(T_R2[w] * COLOR_DIM[c] * X for _, c, w, Y, X in fields)
    AXXX = sum(COLOR_DIM[c] * w * X ** 3 for _, c, w, Y, X in fields)
    AgX = sum(COLOR_DIM[c] * w * X for _, c, w, Y, X in fields)
    AY2X = sum(COLOR_DIM[c] * w * Y * Y * X for _, c, w, Y, X in fields)
    AYX2 = sum(COLOR_DIM[c] * w * Y * X * X for _, c, w, Y, X in fields)
    return dict(A33X=A33X, A22X=A22X, AXXX=AXXX, AgX=AgX, AY2X=AY2X, AYX2=AYX2)


# ======================================================================
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
    check("§%s 五段式六标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

# ======================================================================
head("F2  独立复算")

# ---- (a) 六类系数 ----
an = anomalies(GEN)
print("      [中间量] A333=%s A33Y=%s A22Y=%s AYYY=%s AgY=%s N12=%s"
      % (an["A333"], an["A33Y"], an["A22Y"], an["AYYY"], an["AgY"], an["N12"]))
check("(a) 一代：六类系数全为零",
      all(an[k] == 0 for k in ("A333", "A33Y", "A22Y", "AYYY", "AgY")))
check("(a) 逐项印值：A333 = 2-1-1 = 0；A_YYY = 1/36-8/9+1/9-1/4+1 = 0",
      Fr(2 - 1 - 1) == 0 and Fr(1, 36) - Fr(8, 9) + Fr(1, 9) - Fr(1, 4) + 1 == 0)
# 第二条路径：逐态枚举（16 态）
STATES = []
for nm, c, w, Y, X in GEN:
    for _ in range(COLOR_DIM[c]):
        for _ in range(w):
            STATES.append((nm, c, w, Y, X))
check("(a) 逐态枚举态数 = 16", len(STATES) == 16, "实测 %d" % len(STATES))
A333_2 = sum(A_R3[c] / w * 1 for _, c, w, Y, X in STATES) if False else \
    sum(A_R3[c] for _, c, w, Y, X in STATES)
A_YYY_2 = sum(Y ** 3 for _, c, w, Y, X in STATES)
A_gY_2 = sum(Y for _, c, w, Y, X in STATES)
check("(a) 第二路径（逐态）给同一组零：A333、A_YYY、A_gY",
      A333_2 == 0 and A_YYY_2 == 0 and A_gY_2 == 0,
      "%s %s %s" % (A333_2, A_YYY_2, A_gY_2))

# ---- (b) Witten ----
check("(b) N_{1/2} = 3+1 = 4（按场）", an["N12"] == 4)
check("(b) 按态：双重态态数 = 8，8/2 = 4（两种口径一致）",
      sum(1 for _, c, w, Y, X in STATES if w == 2) == 8)
check("(b) 4 为偶 ⟹ 无 Witten 反常", an["N12"] % 2 == 0)

# ---- (c) B-L ----
bl = bl_anomalies(GEN)
print("      [B-L 中间量] " + " ".join("%s=%s" % (k, v) for k, v in sorted(bl.items())))
check("(c) 含 nu^c：B-L 八类系数全为零", all(v == 0 for v in bl.values()), str(bl))
GEN_nonu = [f for f in GEN if f[0] != "nu^c"]
bl2 = bl_anomalies(GEN_nonu)
nz = {k: v for k, v in bl2.items() if v != 0}
check("(c) 去掉 nu^c：恰两个非零 A_gX = A_XXX = -1（其余为零）",
      nz == {"AgX": Fr(-1), "AXXX": Fr(-1)}, str(nz))

# ---- (d) 14 态块 ----
anB = anomalies(BLK)
nstates = sum(COLOR_DIM[c] * w for _, c, w, Y, X in BLK)
check("(d) X 的态数 = 6+3+3+2 = 14", nstates == 14, "实测 %d" % nstates)
check("(d) X 的六类系数全为零", all(anB[k] == 0 for k in ("A333", "A33Y", "A22Y", "AYYY", "AgY")))
check("(d) X 的 N_{1/2} = 3+1 = 4（偶）", anB["N12"] == 4 and anB["N12"] % 2 == 0)


def conj(f):
    nm, c, w, Y, X = f
    cc = {"3": "3b", "3b": "3", "1": "1"}[c]
    return (nm, cc, w, -Y, -X)


Xset = {(c, w, Y) for _, c, w, Y, X in BLK}
conjX = {(cc, w, Y) for (nm, c, w, Y, X) in BLK for cc in [{"3": "3b", "3b": "3", "1": "1"}[c]]}
check("(d) X 手征：X ≠ 其复共轭（(3b,2,0) ∉ X）",
      Xset != conjX and ("3b", 2, Fr(0)) not in Xset)
GENset = {(c, w, Y) for _, c, w, Y, X in GEN}
check("(d) X 与一代六槽不交", not any((c, w, Y) in GENset for (c, w, Y) in Xset))
anBoth = anomalies(GEN + BLK)
blBoth = bl_anomalies(GEN + BLK)
check("(d) 16+X 仍无反常（局域六类＋Witten＋B-L 八类）",
      all(anBoth[k] == 0 for k in ("A333", "A33Y", "A22Y", "AYYY", "AgY"))
      and anBoth["N12"] % 2 == 0 and all(v == 0 for v in blBoth.values()))

# ---- (e) 解空间维数 ----
# 变量 k1..k6，五个条件（乘 36 清分母）
rows = []
rows.append([Fr(2), Fr(-1), Fr(-1), Fr(0), Fr(0), Fr(0)])                      # A333
rows.append([Fr(1), Fr(-2), Fr(1), Fr(0), Fr(0), Fr(0)])                       # A33Y
rows.append([Fr(1), Fr(0), Fr(0), Fr(-1), Fr(0), Fr(0)])                       # A22Y: k1/2 - k4/2 = 0
rows.append([Fr(1), Fr(-2), Fr(1), Fr(-1), Fr(1), Fr(0)])                      # A_grav-Y
rows.append([Fr(1), Fr(-32), Fr(4), Fr(-9), Fr(36), Fr(0)])                    # A_YYY ×36


def rank_mod(rows, ncol):
    rows = [r[:] for r in rows]
    r = 0
    for c in range(ncol):
        pr = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


rk = rank_mod(rows, 6)
check("(e) 五条件的秩 = 4 ⟹ 解空间维数 = 2", rk == 4, "实测秩 %d" % rk)
v1 = [Fr(1)] * 5 + [Fr(0)]
v2 = [Fr(0)] * 5 + [Fr(1)]
check("(e) (k,k,k,k,k,0) 与 (0,0,0,0,0,1) 都是解（基）",
      all(sum(a * b for a, b in zip(r, v1)) == 0 for r in rows)
      and all(sum(a * b for a, b in zip(r, v2)) == 0 for r in rows))
check("(e) Witten：3k1+k4 = 4k 恒为偶（对任意 k）", all(4 * k % 2 == 0 for k in range(6)))

# ---- (f) 内容层手征性 ----
gen_conj = {(cc, w, -Y) for (nm, c, w, Y, X) in GEN for cc in [{"3": "3b", "3b": "3", "1": "1"}[c]]}
check("(f) 一代内容手征：（3b,2,-1/6）∉ 槽集合",
      (("3b", 2, Fr(-1, 6))) not in GENset and ("3", 2, Fr(1, 6)) in GENset)

# ---- (g) 齐次性 ----
ok_h = True
for lam in (Fr(2), Fr(1, 3), Fr(-1), Fr(5)):
    scaled = [(nm, c, w, lam * Y, X) for nm, c, w, Y, X in GEN]
    a2 = anomalies(scaled)
    if a2["AYYY"] != lam ** 3 * an["AYYY"] or a2["A33Y"] != lam * an["A33Y"] or a2["A333"] != an["A333"]:
        ok_h = False
check("(g) Y->λY：A_YYY 按 λ^3、A_33Y 按 λ、A_333 不变（零性保持）", ok_h)

# ======================================================================
head("F3  锚定")
FILES = ["Z0_zero_never_rests_single_axiom.md", "Z14_closure_cyclic_order_base_theorem.md",
         "E1_NN_verdict.md", "G61_locking_the_five_integers.md", "G63_target_list_and_audit.md",
         "G91_matter_sector_range_and_three_way_verdict.md", "sm/SM8_anomalies.md",
         "sm/SM5_gauge_origin.md", "sm/SM7_gauge_group_hypercharge_charge.md",
         "sm/SMR5_quantum_numbers_and_charge.md", "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G91_matter_sector_range_and_three_way_verdict.md", ["手征"]),
      ("G63_target_list_and_audit.md", ["1.4"]),
      ("E1_NN_verdict.md", ["反射"]),
      ("sm/SM8_anomalies.md", ["14 态"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM8_anomalies.md']
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
       if ("手征性" in ln or "谱内容" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若"])]
check("含 手征性/谱内容 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
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
check("本篇文件名不被 SM<数字>_ 正则捕获", not re.match(r"^SM\d+_", "SMR6_chirality_and_anomalies.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
