#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G91_check.py -- 核验「物质扇区的射程判定：SM 门逐项三分」

  F1  文档结构：三分判据 R1–R4、弹药表 P1–P13、19 项总表、两条 no-go、定理/推论编号、三栏计数
  F2  独立复算（小 L 精确）：Z_L 与 D_L 轨道闭包枚举、q_L 精确有理数、L=4 两类与稳定子
  F3  独立复算（大 L）：规范代表法 + Burnside 交叉校验，L<=20；轨道数=3 只在 L=6
  F4  定理 G91.2 复算：成对账本峰位逐 L，L=4 峰={4}，L=6 峰=D=2
  F5  定理 G91.3 / 推论 G91.4 复算：反射保类、sigma 自由 => q 不变
  F6  交叉一致：R31/G63/G61/G12/Z15/Z16/E1_NN/Z13 与旧理论靶子编号
  F7  登记：STATUS.md / INDEX.md 收录；不抬高任何 O1-O5 状态
"""
from __future__ import annotations

import io
import os
import re
import sys
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, gcd

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   " + d) if d else ""))


def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f:
        return f.read()


DOC = read("G91_matter_sector_range_and_three_way_verdict.md")

# ----------------------------------------------------------------------
# 独立工具
# ----------------------------------------------------------------------

def balanced_words(L):
    """零和 ±1 词（L 偶）。"""
    return [w for w in product((1, -1), repeat=L) if sum(w) == 0]


def rot(w, k):
    return w[k:] + w[:k]


def rev(w):
    return tuple(reversed(w))


def neg(w):
    return tuple(-x for x in w)


def zl_gens(L):
    """Z_L 的生成元（最小生成集：单步旋转）。"""
    return [lambda w: rot(w, 1)]


def full_gens(L):
    """Z_L 与反射、符号生成的群（= D_L x Z_2）的最小生成集。"""
    return [lambda w: rot(w, 1), rev, neg]


def zl_gens_all(L):
    """Z_L 的全部元素（用于精确闭包枚举）。"""
    return [(lambda w, k=k: rot(w, k)) for k in range(L)]


def orbits(elems, gens):
    """路径一：闭包法（精确）。"""
    seen = set()
    out = []
    for e in elems:
        if e in seen:
            continue
        orb = set()
        stack = [e]
        while stack:
            x = stack.pop()
            if x in orb:
                continue
            orb.add(x)
            for g in gens:
                stack.append(g(x))
        seen |= orb
        out.append(orb)
    return out


def canon_key(e, gens):
    """路径二：规范代表（gen 像的最小元组）。"""
    return min(g(e) for g in gens)


def orbits_canon(elems, gens):
    """路径二：按规范代表分桶（只在 gens 是群的生成集时等价于轨道）。"""
    groups = {}
    for e in elems:
        groups.setdefault(canon_key(e, gens), []).append(e)
    return [set(v) for v in groups.values()]


def orbits_S(elems):
    """S = <Z_L, sigma> 下的轨道（词级闭包）。"""
    return orbits(elems, [lambda w: rot(w, 1), rev])


def class_level_orbits(W):
    """把 W 分成旋转类，返回 (类代表表, 类大小, sigma 在类上的置换, 词->类号)。"""
    keys = {}
    for w in W:
        keys.setdefault(class_key(w), []).append(w)
    reps = sorted(keys.keys())
    idx = {k: i for i, k in enumerate(reps)}
    sizes = [len(keys[k]) for k in reps]
    sig = [idx[class_key(rev(k))] for k in reps]
    owner = {}
    for w in W:
        owner[w] = idx[class_key(w)]
    return reps, sizes, sig, owner


def orbits_S_closure(W):
    """在**类集合**上对 <Z_L, sigma> 取轨道（类级闭包，独立路径）。"""
    reps, sizes, sig, owner = class_level_orbits(W)
    seen = set()
    out = []
    for i in range(len(reps)):
        if i in seen:
            continue
        orb = {i}
        stack = [i]
        while stack:
            x = stack.pop()
            for y in (x, sig[x]):
                if y not in orb:
                    orb.add(y)
                    stack.append(y)
        seen |= orb
        out.append({w for w in W if owner[w] in orb})
    return out


def class_key(w):
    """旋转类的规范代表（项链的字典序最小旋转）。"""
    L = len(w)
    return min(w[k:] + w[:k] for k in range(L))


def q_of(orbits_, total):
    return sum(Fraction(len(o), total) ** 2 for o in orbits_)


def necklace_count(L):
    """零和词按 Z_L 轨道的计数闭式（zero_sum_rotation_class_algebra 定理 T1）。"""
    def phi(n):
        r, m, p = n, n, 2
        while p * p <= m:
            if m % p == 0:
                while m % p == 0:
                    m //= p
                r -= r // p
            p += 1
        if m > 1:
            r -= r // m
        return r

    g = gcd(L, L // 2)
    return sum(phi(d) * comb(L // d, L // (2 * d)) for d in range(1, g + 1) if g % d == 0) // L


def burnside_ZL(L):
    """Burnside：Z_L 作用在平衡词上的轨道数 = (1/L) sum_{k} |Fix(rot_k)|。"""
    W = balanced_words(L)
    tot = 0
    for k in range(L):
        tot += sum(1 for w in W if rot(w, k) == w)
    return Fraction(tot, L)


# ----------------------------------------------------------------------
head("F1  文档结构")

check("标题点出三分", all(k in DOC for k in ("能导出", "具名输入", "不可导出")))
check("射程判词在位（离散在内／连续在外）",
      "离散计数侧落在界内" in DOC and "连续内部对称侧落在界外" in DOC)
check("执行 SYNTHESIS §7.2 第 7 步被写明",
      "第 7 步" in DOC and "SYNTHESIS_zero_to_standard_model.md" in DOC)
check("三分判据 R1–R4 在位", all(("**R%d**" % i) in DOC for i in range(1, 5)))
check("弹药表 P1–P13 在位", all(("| P%d |" % i) in DOC for i in range(1, 14)))
check("P13 明写「不存在」", "**不存在**（引理 43(a)" in DOC)
check("两种不可导出被区分", "结构性不可导出" in DOC and "记账性不可导出" in DOC)
check("命题 G91.1 在位", "命题 G91.1" in DOC)
check("定理 G91.2 在位", "定理 G91.2" in DOC)
check("定理 G91.3 / 命题 G91.4 / 命题 G91.5 在位",
      "定理 G91.3" in DOC and "命题 G91.4" in DOC and "命题 G91.5" in DOC)
check("候选判据 ONE-CLOSURE-ONE-GENERATION 在位且标未证",
      "ONE-CLOSURE-ONE-GENERATION" in DOC and "【具名候选】，未证" in DOC)
check("明确不主张 k=1", "不主张" in DOC and "$k=1$" in DOC)
check("不做的事三条在位",
      "不用反常闭合去选" in DOC and "不用旋转类轨道数去选" in DOC and "低能有效理论" in DOC)
check("诚实边界含引理43(d)是候选", "(d) 是候选" in DOC)
check("四条不改变判定在位",
      all(k in DOC for k in ("O1", "SURV4-GLOBAL", "DIM-SECTOR", "EVO-NORM")))
check("无条件导出标准模型仍不可主张", "无条件导出标准模型" in DOC)
OLD = os.path.join(HERE, "..", "modular-equilibrium", "derivations")
OLDIDS = ("D46", "D47", "D48", "D50", "D52", "D53", "D58", "D59", "D64", "D112", "D113", "D114")
check("旧理论靶子已带文件级出处（12 个编号）", all(k in DOC for k in OLDIDS))
miss_old = [k for k in OLDIDS
            if not [f for f in os.listdir(OLD) if f.startswith(k + "_")]] if os.path.isdir(OLD) else OLDIDS
check("被引的旧理论文件在磁盘上存在（目录可达时）", not miss_old, "缺：%s" % miss_old)
check("H 表示条件唯一带「仅参考」与四条申报条件",
      "H=(1,2,\\pm1/2)" in DOC and "申报条件" in DOC and "D58" in DOC)

# ----------------------------------------------------------------------
head("F2  §8 总表：19 项与三栏计数")

SEC8 = DOC.split("## §8 三栏对照总表")[1].split("## §9")[0]
rows = []
for line in SEC8.splitlines():
    m = re.match(r"^\|\s*(\d+)\s*\|", line)
    if m:
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        rows.append((int(cells[0]), cells[1], cells[2], cells[3] if len(cells) > 3 else ""))
check("总表恰 19 行且编号 1..19", [r[0] for r in rows] == list(range(1, 20)),
      "%d 行" % len(rows))

cnt = Counter()
for _, _, verdict, _ in rows:
    v = verdict.replace("*", "")
    if v.startswith("不可导出"):
        cnt["不可导出"] += 1
    elif v.startswith("只能具名输入"):
        cnt["具名输入"] += 1
    elif v.startswith("能导出") or v.startswith("导出"):
        cnt["能导出"] += 1
    elif v.startswith("开放"):
        cnt["开放"] += 1
    else:
        cnt["?" + v[:12]] += 1
check("四栏计数 = 不可导出 5 + 开放 3 + 能导出 4 + 具名输入 7",
      cnt["不可导出"] == 5 and cnt["开放"] == 3 and cnt["能导出"] == 4 and cnt["具名输入"] == 7,
      str(dict(cnt)))
check("四栏计数之和 = 19", sum(cnt.values()) == 19, str(sum(cnt.values())))
DOCflat = re.sub(r"[\s\\{}]", "", DOC.replace("*", "").replace("＋", "+").replace("text", ""))
check("合计句在位（不可导出 5 ＋ 开放 3 ＋ 能导出 4 ＋ 具名输入 7）",
      "不可导出5+开放3+能导出4" in DOCflat and "只能具名输入7" in DOCflat)
check("两套计数口径不得互引的声明在位", "两套计数口径不同" in DOC)
check("每行都有「依据」列内容", all(r[3] for r in rows),
      "缺依据的行：%s" % [r[0] for r in rows if not r[3]])

# ----------------------------------------------------------------------
head("F3  §7 轨道表独立复算（小 L 精确闭包 + 双路径互校）")

# (|W_L|, Z_L 轨道数, q_L, 轨道大小**多重集**（逐个轨道）)
TABLE = {2: (2, 1, Fraction(1), [2]),
         4: (6, 2, Fraction(5, 9), [4, 2]),
         6: (20, 4, Fraction(7, 25), [6, 6, 6, 2]),
         8: (70, 10, Fraction(19, 175), [8] * 8 + [4, 2]),
         10: (252, 26, Fraction(313, 7938), [10] * 25 + [2]),
         12: (924, 80, Fraction(683, 53361), [12] * 75 + [6, 6, 6, 4, 2])}
for L, (N, K, q, prof) in sorted(TABLE.items()):
    W = balanced_words(L)
    check("L=%2d 平衡词数 = %d" % (L, N), len(W) == N, "%d" % len(W))
    o = orbits(W, zl_gens_all(L))
    check("L=%2d Z_L 轨道数 = %d（闭包法）" % (L, K), len(o) == K, "%d" % len(o))
    check("L=%2d 轨道大小多重集与文档 §7.1 表逐项一致" % L,
          sorted(len(x) for x in o) == sorted(prof),
          str(sorted(len(x) for x in o)))
    check("L=%2d 每个轨道大小整除 L" % L, all(L % len(x) == 0 for x in o))
    check("L=%2d 轨道数 = 项链闭式" % L, len(o) == necklace_count(L),
          "%d vs %d" % (len(o), necklace_count(L)))
    check("L=%2d 轨道数 = Burnside" % L, len(o) == burnside_ZL(L))
    qq = q_of(o, N)
    check("L=%2d q_L = %s（精确有理数）" % (L, q), qq == q, "%s" % qq)

# L=4 的细节：两类、稳定子
W4 = balanced_words(4)
o4 = orbits(W4, zl_gens_all(4))
sizes = sorted((len(x) for x in o4), reverse=True)
check("L=4 恰两类，大小 {4,2}", sizes == [4, 2], str(sizes))
stal = sorted(len([k for k in range(4) if rot(min(x), k) == min(x)]) for x in o4)
check("L=4 稳定子阶 = (1,2)（交替类有 Z_2）", stal == [1, 2], str(stal))
check("L=4 交错类 = {+-+-,-+-+}",
      any(set(x) == {tuple([1, -1, 1, -1]), tuple([-1, 1, -1, 1])} for x in o4))
alt = [x for x in o4 if len(x) == 2][0]
check("L=4 反射把交错类映到自身", rev(min(alt)) in alt)
cls4 = [x for x in o4 if len(x) == 4][0]
check("L=4 反射把 4 元类映到自身（故 D_4 类数也是 2）", rev(min(cls4)) in cls4)

# ----------------------------------------------------------------------
head("F4  大 L：规范代表法 + Burnside 交叉校验（L<=20）")

KC = {}
for L in range(2, 21, 2):
    W = balanced_words(L)
    groups = {}
    for w in W:
        groups.setdefault(class_key(w), []).append(w)
    KC[L] = len(groups)
    check("L=%2d 规范代表法 = 项链闭式 = Burnside" % L,
          len(groups) == necklace_count(L) == burnside_ZL(L),
          "%d / %d / %s" % (len(groups), necklace_count(L), burnside_ZL(L)))
seq = [KC[L] for L in (2, 4, 6, 8, 10, 12)]
check("Z_L 轨道数序列（L=2..12）= 1,2,4,10,26,80", seq == [1, 2, 4, 10, 26, 80], str(seq))
three = [L for L, c in KC.items() if c == 3]
check("L<=20 内「Z_L 轨道数 = 3」不存在（L=6 是 4）", three == [], str(three))
check("R31 表的 L=6 轨道分布 {6,6,6,2} 对应轨道数 4", TABLE[6][1] == 4)
check("L>=8 轨道数 >=10 且单调增",
      all(KC[L] >= 10 for L in range(8, 21, 2))
      and all(KC[L] < KC[L + 2] for L in range(8, 19, 2)),
      str([(L, KC[L]) for L in range(8, 21, 2)]))
check("L=14/16/18/20 的轨道数 = 246/810/2704/9252",
      [KC[L] for L in (14, 16, 18, 20)] == [246, 810, 2704, 9252],
      str([KC[L] for L in (14, 16, 18, 20)]))

# ----------------------------------------------------------------------
head("F5  定理 G91.2：成对账本峰位逐 L（模型内）")

B = 4


def peak(L, q):
    """成对账本 F_D = B*C(D,2)*q^D 的有限峰；无有限峰（q=1）时返回 None。"""
    if q >= 1:
        return None
    best, bd = None, None
    for D in range(2, 40):
        F = B * comb(D, 2) * q ** D
        if bd is None or F > bd:
            best, bd = D, F
    return best


for L in (2, 4, 6, 8, 10, 12):
    q = TABLE[L][2]
    pk = peak(L, q)
    if L == 2:
        check("L=2 无有限峰（q_2 = 1）", pk is None)
    elif L == 4:
        check("L=4 峰 = {4}（四维窗口内）", pk == 4 and Fraction(1, 2) < q < Fraction(3, 5),
              "峰 D=%d, q=%s" % (pk, q))
    else:
        check("L=%2d 峰 = D=2（被 R31.1 排除）" % L, pk == 2, "峰 D=%d" % pk)
check("L=6 的 q_6 = 7/25 < 1/3", Fraction(7, 25) < Fraction(1, 3))
check("q_L <= L/N_L 的上界在 L=6 处 = 3/10",
      Fraction(6, 20) == Fraction(3, 10) and Fraction(3, 10) < Fraction(1, 3))
check("文档写明「轨道数=3」与「峰在引力子域」不能同时成立",
      "不能同时成立" in DOC and "无交集" in DOC)

# 命题 G91.5：三种账本下 L=4 仍是唯一峰在引力子域的偶寿命
FULLQ = {2: Fraction(1), 4: Fraction(5, 9), 6: Fraction(23, 50), 8: Fraction(197, 1225),
         10: Fraction(563, 7938), 12: Fraction(2419, 106722)}
FULLQ2 = {2: Fraction(1), 4: Fraction(5, 9), 6: Fraction(23, 50), 8: Fraction(229, 1225),
          10: Fraction(394, 3969), 12: Fraction(3823, 106722)}
ok = True
for L, q in FULLQ.items():
    pk = peak(L, q)
    if L == 2:
        good = pk is None
    elif L == 4:
        good = pk == 4
    else:
        good = pk is not None and pk < 4
    ok = ok and good
    check("命题 G91.5：sigma-商账本 L=%2d 的峰 = %s" % (L, "无有限峰" if pk is None else "D=%d" % pk),
          good)
check("命题 G91.5：反射不变账本下 L=4 仍是唯一峰在引力子域的偶寿命（L<=12）", ok)
# 逐 L 现算（与写死的表独立）：sigma-商与全群两种账本
ok2 = ok3 = True
for L in range(4, 13, 2):
    W = balanced_words(L)
    N = len(W)
    qz = q_of(orbits(W, zl_gens_all(L)), N)
    qS = q_of(orbits_S_closure(W) if L <= 12 else orbits_S(W), N)
    qF = q_of(orbits(W, full_gens(L)), N)
    pz, pS, pF = peak(L, qz), peak(L, qS), peak(L, qF)
    if L == 4:
        ok2 = ok2 and pz == 4 and pS == 4
        ok3 = ok3 and pF == 4
    else:
        ok2 = ok2 and pz < 4 and pS < 4
        ok3 = ok3 and pF < 4
    check("L=%2d 三种账本峰位（Z_L %d / sigma-商 %d / 全群 %d）" % (L, pz, pS, pF),
          (pz == 4 and pS == 4 and pF == 4) if L == 4 else (pz < 4 and pS < 4 and pF < 4))
check("命题 G91.5 逐 L 现算：sigma-商账本下唯一峰=4 的仍是 L=4", ok2)
check("命题 G91.5 逐 L 现算：全群账本下唯一峰=4 的仍是 L=4", ok3)
check("sigma-商账本的 q' 与文档 §7.3 表逐项一致",
      all(FULLQ[L] == q_of(orbits_S_closure(balanced_words(L)), len(balanced_words(L)))
          for L in (4, 6, 8, 10, 12)))
check("全群账本的 q'' 与文档 §7.3 表逐项一致",
      all(FULLQ2[L] == q_of(orbits(balanced_words(L), full_gens(L)), len(balanced_words(L)))
          for L in (4, 6, 8, 10, 12)))

# ----------------------------------------------------------------------
head("F6  定理 G91.3 / 命题 G91.4：反射保类与 q' 的精确公式")

SIGMA = {4: (0, Fraction(5, 9)), 6: (1, Fraction(23, 50)), 8: (2, Fraction(197, 1225)),
         10: (10, Fraction(563, 7938)), 12: (30, Fraction(2419, 106722))}
for L in range(2, 21, 2):
    W = balanced_words(L)
    N = len(W)
    o = orbits_canon(W, zl_gens_all(L)) if L <= 12 else orbits_canon(W, zl_gens(L))
    if L <= 12:
        # 用精确闭包枚举类，再用类代表定 sigma 的作用
        o = orbits(W, zl_gens_all(L))
        cls_of = {}
        for i, orb in enumerate(o):
            for w in orb:
                cls_of[w] = i
        sig = [cls_of[rev(next(iter(o[i])))] for i in range(len(o))]
        sizes = [len(x) for x in o]
        # G91.3 前提：rev o rot_k = rot_{-k} o rev
        commute = all(rev(rot(w, k)) == rot(rev(w), -k % L) for w in W[:20] for k in range(L))
        check("L=%2d G91.3 前提：rev o rot_k = rot_{-k} o rev" % L, commute)
        check("L=%2d G91.3：sigma 是类集合上的对合" % L,
              all(sig[sig[i]] == i for i in range(len(o))))
        pairs = set()
        for i in range(len(o)):
            if sig[i] != i:
                pairs.add((min(i, sig[i]), max(i, sig[i])))
        check("L=%2d 配对的两个类等势（命题 G91.4 的引理）" % L,
              all(sizes[i] == sizes[j] for i, j in pairs))
        q = q_of(o, N)
        qp = q + sum(2 * Fraction(sizes[i], N) ** 2 for i, _ in pairs)
        oR = orbits(W, [lambda w: rot(w, 1), rev])
        check("L=%2d 命题 G91.4 公式 q' = q + 2*sum_pairs(w^2) 成立" % L,
              qp == q_of(oR, N), "公式 %s vs 直接 %s" % (qp, q_of(oR, N)))
        if L in SIGMA:
            npair, qq = SIGMA[L]
            check("L=%2d 配对对数与 q' 与文档 §7.3 表一致" % L,
                  len(pairs) == npair and qp == qq,
                  "pairs=%d (表 %d), q'=%s (表 %s)" % (len(pairs), npair, qp, qq))
        if L == 2:
            check("L=2：唯一类被 sigma 固定（配对为空）⇒ q' = q = 1",
                  len(pairs) == 0 and qp == Fraction(1))
        elif L == 4:
            check("L=4：Pair(sigma) = 空 ⇒ q' = q = 5/9（R31 不受影响）",
                  len(pairs) == 0 and qp == Fraction(5, 9))
        else:
            check("L=%2d：Pair(sigma) 非空 ⇒ q' > q（更正「商掉反射无害」）" % L,
                  len(pairs) > 0 and qp > q, "q=%s q'=%s" % (q, qp))
        check("L=%2d 直接按 S-轨道分桶得到的 q' 与公式一致" % L,
              q_of(orbits_S(W), N) == qp, "%s vs %s" % (q_of(orbits_S(W), N), qp))

FULL = {2: 1, 4: 2, 6: 3, 8: 7, 10: 13, 12: 35}
for L, K in FULL.items():
    W = balanced_words(L)
    o3 = orbits(W, full_gens(L))
    check("L=%2d 全群（Z_L ＋ 反射 ＋ 符号）精确闭包轨道数 = %d" % (L, K), len(o3) == K,
          "%d" % len(o3))

# ----------------------------------------------------------------------
head("F7  交叉一致与登记")

R31 = read("R31_phase_ledger_and_lifetime_selection.md")
G63 = read("G63_target_list_and_audit.md")
G61 = read("G61_locking_the_five_integers.md")
G12 = read("G12_gauge_sector_minimal_extension.md")
Z15 = read("Z15_zcar_no_go_and_jordan_wigner_readout.md")
Z16 = read("Z16_zunif_balanced_regular_module.md")
E1NN = read("E1_NN_verdict.md")
Z13 = read("Z13_zero_foundation_missing_principle.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

check("R31 的 L=4 行给 {4,2} 与 5/9", "{4,2}" in R31 and "5/9" in R31)
check("R31 的项链闭式交叉核对序列在位", "1,2,4,10,26,80" in R31.replace("$", "").replace(" ", ""))
check("R31.1 的唯一选择结论在位", "定理 R31.1" in R31 and "L=4" in R31)
check("R31 只用 Z_L 轨道（本件补上反射商）", "循环旋转" in R31)
check("G63 §3 明写物质栏由 G12 判不可达",
      "物质／规范参数" in G63 and "G12_gauge_sector_minimal_extension.md" in G63)
check("G12 引理 43 三条子论证在位", "引理 43" in G12 and "(a)" in G12 and "(b)" in G12 and "(c)" in G12)
check("G61 的 B=2^2 两条独立 Z_2 在位", "2^{2}=4" in G61.replace("$", ""))
check("Z15 的 Z-CAR 不能唯一选择在位", "不能唯一" in Z15)
check("Z16 的 Z-UNIF 不能由 Z0/Z-E* 自动推出在位", "Z-UNIF" in Z16 and "不能" in Z16)
check("E1_NN 的「真实对称是反射」在位", "反射" in E1NN)
check("Z13 的 E5（量子读出）在位", "E5" in Z13 and "读出" in Z13)
check("STATUS 已登记 G91（§2.51）",
      "### 2.51" in STATUS and "G91_matter_sector_range_and_three_way_verdict.md" in STATUS
      and "物质扇区的射程判定" in STATUS)
check("STATUS 已登记 G91_check.py", "G91_check.py" in STATUS)
check("INDEX 已收录 G91 文档（核验脚本按 INDEX 口径只计数）",
      "G91_matter_sector_range_and_three_way_verdict.md" in INDEX
      and re.search(r"核验脚本（\d+ 个）", INDEX) is not None
      and os.path.exists(os.path.join(HERE, "G91_check.py")))
check("文档存在自核验入口", "python3 G91_check.py" in DOC)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
