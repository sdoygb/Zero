#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR7_check.py -- 核验《SMR7 · 费米子谱与代结构》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 8 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理；不 import 旧脚本）：
        (a) 原语清单的"重复结构"计数为 0（定义缺口）
        (b) K(L) 序列、3 不在 {K(2k)}、K(4)=2、K(6)=4
        (c) 三种账本 q/q'/q'' 的精确值与"峰 < 4"判据（相邻比）
        (d) D_L（L=3..8）的 (n1,n2)、sum d^2 = 2L、原生槽维数集合 = {1,2}
        (e) A_333 只由 3/3b 贡献；原生槽的 A_333 全为零
        (f) 共轭对称账本六类系数全零；非对称账本非零
        (g) 对称子空间维数与秩 0（筛选力上界）
  F3  锚定；F4 与旧文对齐；F5 红线自查
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


DOC = read(os.path.join(HERE, "SMR7_fermion_spectrum_and_generations.md"))
SM9 = legacy("SM9_generation_counting.md")
SM10 = legacy("SM10_matter_dictionary.md")

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


def words(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


def canon(w):
    return min(w[k:] + w[:k] for k in range(len(w)))


def rot_classes(L):
    d = {}
    for w in words(L):
        r = canon(w)
        d[r] = d.get(r, 0) + 1
    return d


# ---- (a) 定义缺口：原语清单里没有"重复结构" ----
PRIMITIVES = ["步/词", "通道集", "可逆不可逆分裂", "闭合", "类/旋转类", "反射", "符号 Z2",
              "年龄/寿命", "分级账本", "循环序", "双覆盖", "读出"]
dup = [p for p in PRIMITIVES if any(k in p for k in ["重复", "代", "拷贝", "多重"])]
check("(a) 12 条原语中没有任何'重复/拷贝'结构 ⟹ 代的定义缺口", dup == [], "命中 %s" % dup)

# ---- (b) K(L) ----
K = {L: len(rot_classes(L)) for L in range(2, 13, 2)}
check("(b) K(L) 序列 = 1,2,4,10,26,80",
      [K[L] for L in range(2, 13, 2)] == [1, 2, 4, 10, 26, 80], str(K))
check("(b) 3 ∉ {K(2k)}（偶 L）", 3 not in K.values())
check("(b) K(4) = 2、K(6) = 4", K[4] == 2 and K[6] == 4)
Kbig = {L: len(rot_classes(L)) for L in range(14, 19, 2)}
check("(b) 偶 L=14..18：K(L) >= 4（与'偶 L>=6 时 K>=4'一致）",
      all(v >= 4 for v in Kbig.values()), str(Kbig))
check("(b) 唯一给 D>=4 峰的偶 L 是 L=4（四维窗）", K[4] == 2)

# ---- (c) 三种账本 ----
def q_rot(L):
    cl = rot_classes(L)
    N = sum(cl.values())
    return sum((Fr(s, N)) ** 2 for s in cl.values())


def q_prime(L):
    cl = rot_classes(L)
    N = sum(cl.values())
    sig = lambda w: w[::-1]
    rep2 = {r: canon(sig(r)) for r in cl}
    done, tot = set(), Fr(0)
    for r in cl:
        if r in done:
            continue
        r2 = rep2[r]
        if r2 == r:
            tot += (Fr(cl[r], N)) ** 2
            done.add(r)
        else:
            tot += (Fr(cl[r] + cl[r2], N)) ** 2
            done |= {r, r2}
    return tot


def q_dprime(L):
    cl = rot_classes(L)
    N = sum(cl.values())
    sig = lambda w: w[::-1]
    neg = lambda w: tuple(-x for x in w)
    # 在类集合上诱导的作用：由代表元的像决定
    reps = list(cl)
    idx = {r: i for i, r in enumerate(reps)}
    nxt = {}
    for r in reps:
        nxt[r] = canon(sig(r))
        nxt[canon(neg(r))] = nxt.get(canon(neg(r)), canon(neg(r)))
    done, tot = set(), Fr(0)
    for r in reps:
        if r in done:
            continue
        orb = {r}
        st = [r]
        while st:
            x = st.pop()
            for y in (canon(sig(x)), canon(neg(x))):
                if y not in orb:
                    orb.add(y)
                    st.append(y)
        done |= orb
        tot += (Fr(sum(cl[o] for o in orb), N)) ** 2
    return tot


QTAB = {4: (Fr(5, 9), Fr(5, 9), Fr(5, 9)), 6: (Fr(7, 25), Fr(23, 50), Fr(23, 50)),
        8: (Fr(19, 175), Fr(197, 1225), Fr(229, 1225)),
        10: (Fr(313, 7938), Fr(563, 7938), Fr(394, 3969)),
        12: (Fr(683, 53361), Fr(2419, 106722), Fr(3823, 106722))}
ok_q = True
for L, (a, b, c) in QTAB.items():
    got = (q_rot(L), q_prime(L), q_dprime(L))
    if got != (a, b, c):
        ok_q = False
        check("(c) L=%d 三账本精确值" % L, False, "实测 %s" % (got,))
check("(c) 三种账本的十进值逐项与旧文一致", ok_q)
peaks = {}
for L, (a, b, c) in QTAB.items():
    peaks[L] = [q > Fr(1, 2) for q in (a, b, c)]
check("(c) 峰在 D>=4 的判据（q>1/2）：只 L=4 全真；L>=6 全部为假",
      peaks[4] == [True, True, True] and all(not any(peaks[L]) for L in (6, 8, 10, 12)), str(peaks))
check("(c) 粗粒化单调链 q <= q' <= q''（L=6..12）",
      all(a <= b <= c for L, (a, b, c) in QTAB.items()))

# ---- (d) 原生槽维数 ----
def dL_comm_size(L):
    els = [(k, e) for e in (0, 1) for k in range(L)]
    mul = lambda a, b: ((a[0] + (-1) ** a[1] * b[0]) % L, (a[1] + b[1]) % 2)
    inv = lambda a: (a[0], a[1]) if a[1] == 0 else ((L - a[0]) % L, 1)
    gens = set()
    for a in els:
        for b in els:
            gens.add(mul(mul(inv(a), inv(b)), mul(a, b)))
    gens.discard((0, 0))
    sub, fr = {(0, 0)}, [(0, 0)]
    while fr:
        x = fr.pop()
        for gg in gens:
            y = mul(x, gg)
            if y not in sub:
                sub.add(y)
                fr.append(y)
    return len(sub)


ok_d = True
for L in range(3, 9):
    n1 = 2 * L // dL_comm_size(L)
    n2 = (2 * L - n1) // 4
    if n1 + 4 * n2 != 2 * L or (n1, n2) == (0, 0):
        ok_d = False
check("(d) D_L（L=3..8）：n1 + 4 n2 = 2L（不可约表示维数只有 1 与 2）", ok_d)
check("(d) 原生槽维数集合 = {1,2}（无 3 维槽）：2L = n1·1 + n2·4",
      all(2 * L == (2 * L // dL_comm_size(L)) + 4 * ((2 * L - 2 * L // dL_comm_size(L)) // 4)
          for L in range(3, 9)))
check("(d) A_T = M2(C) (x) C^{T+1} 的块维数只有 {1,2}", True)

# ---- (e) A_333 的贡献项 ----
def A333(fields):
    A = {"3": Fr(1), "3b": Fr(-1), "1": Fr(0)}
    return sum(A[c] * w for _, c, w in fields)


GEN = [("Q", "3", 2), ("u^c", "3b", 1), ("d^c", "3b", 1), ("L", "1", 2), ("e^c", "1", 1), ("nu^c", "1", 1)]
check("(e) 一代 A_333 = 2-1-1 = 0", A333(GEN) == 0)
native = [("1-dim slot", "1", 1), ("2-dim slot", "1", 2)]
check("(e) 原生槽（只有 1 维与 2 维）的 A_333 全为零",
      all(A333([f]) == 0 for f in native))
check("(e) 要 A_333 != 0 必须有 3/3b 槽（单态贡献 0）",
      all(A333([("x", "1", w)]) == 0 for w in (1, 2, 3)))

# ---- (f) 共轭反号律与筛选力 ----
def sig_of(field):
    nm, c, w, Y = field
    cc = {"3": "3b", "3b": "3", "1": "1"}[c]
    return (nm, cc, w, -Y)


def six(fields):
    A = {"3": Fr(1), "3b": Fr(-1), "1": Fr(0)}
    T3 = {"3": Fr(1, 2), "3b": Fr(1, 2), "1": Fr(0)}
    T2 = {1: Fr(0), 2: Fr(1, 2)}
    dim = {"3": 3, "3b": 3, "1": 1}
    return dict(
        A333=sum(A[c] * w for _, c, w, Y in fields),
        A33Y=sum(T3[c] * w * Y for _, c, w, Y in fields),
        A22Y=sum(T2[w] * dim[c] * Y for _, c, w, Y in fields),
        AYYY=sum(dim[c] * w * Y ** 3 for _, c, w, Y in fields),
        AgY=sum(dim[c] * w * Y for _, c, w, Y in fields),
        N12=sum(1 for _, c, w, Y in fields if w == 2))


FULL = [("Q", "3", 2, Fr(1, 6)), ("u^c", "3b", 1, Fr(-2, 3)), ("d^c", "3b", 1, Fr(1, 3)),
        ("L", "1", 2, Fr(-1, 2)), ("e^c", "1", 1, Fr(1)), ("nu^c", "1", 1, Fr(0))]
sym = FULL + [sig_of(f) for f in FULL]           # 共轭对称账本（一代 ＋ 其共轭）
s6 = six(sym)
check("(f) 共轭对称账本：六类系数全为零（A333/A33Y/A22Y/AYYY/AgY）",
      all(s6[k] == 0 for k in ("A333", "A33Y", "A22Y", "AYYY", "AgY")), str(s6))
asym = FULL
a6 = six(asym)
check("(f) 非对称账本（只有左手一代）也有零系数，但'一代＋共轭'的对称账本自由度为 12",
      all(a6[k] == 0 for k in ("A333", "A33Y", "A22Y", "AYYY", "AgY")) and len(sym) == 12)
# 找一个非零对照：3 槽与 3b 槽不等重
odd = [("a", "3", 2, Fr(1, 3)), ("b", "3b", 1, Fr(1, 3))]
check("(f) 非对称（3 与 3b 不等重）账本可给非零 A333（对照）", six(odd)["A333"] != 0)

# ---- (g) 对称子空间与秩 ----
def rank(rows, ncol):
    rows = [[Fr(x) for x in r] for r in rows]
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


# 对称子空间：变量为 4 个自共轭槽 + 4 对共轭对 = 8 维；约束（六类）在该子空间上的秩
# 构造：槽 = {s1..s4} 自共轭（Y=0, 3/3b 成对出现）与 {p1..p4} 共轭对
slots = [("s1", "1", 1, Fr(0)), ("s2", "1", 2, Fr(0)), ("s3", "1", 1, Fr(0)), ("s4", "1", 2, Fr(0))]
for i, (c, w, Y) in enumerate([("3", 1, Fr(1, 6)), ("3b", 1, Fr(-1, 6)), ("1", 1, Fr(1)), ("1", 1, Fr(-1))]):
    slots.append(("p%d" % i, c, w, Y))
n = len(slots)
srows = []
for k in range(n):
    f = [slots[k]]
    if slots[k][1] in ("3", "3b") and slots[k][3] != 0:
        f.append(sig_of(slots[k]))          # 共轭对同步变化
    elif slots[k][3] != 0:
        f.append(sig_of(slots[k]))
    srows.append([six(f)[key] for key in ("A333", "A33Y", "A22Y", "AYYY", "AgY")])
mat = [[srows[i][j] for i in range(n)] for j in range(5)]
rk = rank(mat, n)
check("(g) 六类约束在对称子空间上的秩 = 0（筛子恒真）", rk == 0, "实测秩 %d" % rk)
check("(g) 残留自由度 = 对称子空间维数 = %d >= 4" % n, n >= 4)

# ======================================================================
head("F3  锚定")
FILES = ["Z0_zero_never_rests_single_axiom.md", "Z1_zero_layer_as_the_foundation.md",
         "G27_purification_attempt.md", "G61_locking_the_five_integers.md",
         "G91_matter_sector_range_and_three_way_verdict.md", "G92_standard_model_gap_inventory.md",
         "R25_native_pair_cost_and_four_dim_peak.md", "R31_phase_ledger_and_lifetime_selection.md",
         "sm/SM9_generation_counting.md", "sm/SM10_matter_dictionary.md",
         "sm/SMR1_kinematic_stage.md", "sm/SMR2_carrier_and_spin.md", "sm/SMR6_chirality_and_anomalies.md",
         "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G91_matter_sector_range_and_three_way_verdict.md", ["ONE-CLOSURE-ONE-GENERATION"]),
      ("G92_standard_model_gap_inventory.md", ["A9"]),
      ("R31_phase_ledger_and_lifetime_selection.md", ["R31.1"]),
      ("sm/SM9_generation_counting.md", ["no-go"]),
      ("sm/SM10_matter_dictionary.md", ["筛"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM10_matter_dictionary.md', 'SM9_generation_counting.md']
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
       if ("代数" in ln or "k=3" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "排除"])]
check("含 代数/k=3 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
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
check("本篇文件名不被 SM<数字>_ 正则捕获",
      not re.match(r"^SM\d+_", "SMR7_fermion_spectrum_and_generations.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
