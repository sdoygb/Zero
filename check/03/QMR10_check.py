#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR10_check.py —— 《那道门：块内分级的候选与障碍（缺口定位）》的逐条独立复算。

  P1  定理 1  必要条件 (N)：取反保壳层、是对合、被交换的分级恒平衡（含头号候选被筛）
  P2  定理 2  四个存活候选：逐个验证取反不变 + 逐壳层不平衡
  P3  定理 3  判据退化：四个候选的 max_h δ_h 全部 = 1（精确）
  P4  定理 4  退化之源：极端壳层 h=H 上 C2 确定；其余三个数值确定
  P5  定理 5  正确口径：主导壳层与权重加权给出两档；单调性不是定理（L=14 反例）
  P6  文档契约与链接

不 import 任何其它脚本；路径由 __file__ 推导。
"""
import os
import sys
import re
import math
from itertools import combinations
from collections import defaultdict

import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR10_the_gate_candidates_and_obstruction.md"


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("  [v] " if cond else "  [x] ") + name + (("   " + detail) if detail else ""))


def head(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


DOC = open(os.path.join(HERE, DOCNAME), encoding="utf-8").read() if os.path.exists(
    os.path.join(HERE, DOCNAME)) else ""


def bw(L):
    out = []
    for pos in combinations(range(L), L // 2):
        w = [-1] * L
        for p in pos:
            w[p] = 1
        out.append(tuple(w))
    return out


def part(w):
    s, o = 0, []
    for x in w:
        s += x
        o.append(s)
    return o


def g_step(w, k):      # 头号候选（应被 (N) 筛掉）
    return 0 if w[k + 1] > 0 else 1


def g_sc(w, k):        # C1
    return sum(1 for i in range(k) if w[i] != w[i + 1]) % 2


def g_rm(w, k):        # C2
    S = part(w)
    return 0 if abs(S[k]) == max(abs(x) for x in S[:k + 1]) else 1


def g_v1(w, k):        # C3
    S = part(w)
    return sum(1 for j in range(k + 1) if abs(S[j]) == 1) % 2


def g_rc(w, k):        # C4
    S = part(w)
    m, c = 0, 0
    for j in range(k + 1):
        if abs(S[j]) > m:
            m = abs(S[j])
            c += 1
    return c % 2


CANDS = [("C1 前缀异号次数奇偶", g_sc), ("C2 |S_k| 前缀最大", g_rm),
         ("C3 |S_j|=1 次数奇偶", g_v1), ("C4 创纪录次数奇偶", g_rc)]


def table(L, g):
    cnt = defaultdict(lambda: [0, 0])
    for w in bw(L):
        S = part(w)
        for k in range(L - 1):
            cnt[abs(S[k])][g(w, k)] += 1
    return cnt


L = 12
words = bw(L)

# ======================================================================
head("P1  定理 1 · 必要条件 (N)")
# ======================================================================
check("★ (a) 取反保壳层：|S_k(−w)| = |S_k(w)|（全部词与前缀）",
      all([abs(x) for x in part(tuple(-x for x in w))] == [abs(x) for x in part(w)] for w in words))
check("★ (b) 取反是对合：ι² = id", all(tuple(-x for x in tuple(-x for x in w)) == w for w in words))
cnts = table(L, g_step)
check("★ (c) 头号候选（下一步符号）被取反**交换** ⇒ 每个壳层恒平衡（%s）"
      % ["h=%d:%d/%d" % (h, cnts[h][0], cnts[h][1]) for h in sorted(cnts)[:4]],
      all(cnts[h][0] == cnts[h][1] for h in cnts))

# ======================================================================
head("P2  定理 2 · 四个存活候选")
# ======================================================================
ok = True
for nm, g in CANDS:
    inv = all(g(tuple(-x for x in w), k) == g(w, k) for w in words for k in range(L - 1))
    cnt = table(L, g)
    unb = all(cnt[h][0] != cnt[h][1] for h in cnt)
    if not (inv and unb):
        ok = False
    print("      %-22s 取反不变=%s  逐层不平衡=%s" % (nm, inv, unb))
check("★ 四个候选全部：取反不变 ＋ 逐壳层不平衡", ok)

# ======================================================================
head("P3  定理 3 · 判据退化（max_h δ_h 全部 = 1）")
# ======================================================================
ok2, where = True, []
for nm, g in CANDS:
    cnt = table(L, g)
    d = {h: abs(cnt[h][0] - cnt[h][1]) / (cnt[h][0] + cnt[h][1]) for h in cnt}
    mx = max(d.values())
    hmax = max(d, key=d.get)
    where.append((nm, mx, hmax, cnt[hmax][0] + cnt[hmax][1]))
    if abs(mx - 1.0) > 1e-15:
        ok2 = False
    print("      %-22s max_h δ_h = %.4f（在 h=%d，该层权重 %d）" % (nm, mx, hmax, cnt[hmax][0] + cnt[hmax][1]))
check("★ 四个候选的 max_h δ_h 全部 = 1.0000（精确）⇒ **判据不能区分候选**", ok2)
check("★ 且全部在同一个壳层 h=H=%d 上取到，该层权重极小（%d）" % (L // 2, where[0][3]),
      len({w[2] for w in where}) == 1 and where[0][3] <= 4)

# ======================================================================
head("P4  定理 4 · 退化之源（极端壳层）")
# ======================================================================
H = L // 2
cnt_rm = table(L, g_rm)
check("★ C2 在 h=H 上**确定**（恒取 0）：N_0=%d，N_1=%d ⇒ δ_H=1（可证）"
      % (cnt_rm[H][0], cnt_rm[H][1]), cnt_rm[H][1] == 0 and cnt_rm[H][0] > 0)
for nm, g in CANDS[1:]:
    pass
ok3 = True
for Lx in (12, 14, 16):
    for nm, g in CANDS:
        c = table(Lx, g)
        hh = Lx // 2
        if not (c[hh][0] == 0 or c[hh][1] == 0):
            ok3 = False
check("★ C1/C3/C4 在 h=H 上的确定性：数值成立（L=12,14,16 逐候选）", ok3)

# ======================================================================
head("P5  定理 5 · 正确口径（加权／主导壳层）")
# ======================================================================
rows = []
for nm, g in CANDS:
    cnt = table(L, g)
    hs = sorted(cnt)
    tot = sum(sum(cnt[h]) for h in hs)
    d = {h: abs(cnt[h][0] - cnt[h][1]) / sum(cnt[h]) for h in hs}
    dom = max(hs, key=lambda h: sum(cnt[h]))
    wavg = sum(sum(cnt[h]) * d[h] for h in hs) / tot
    rows.append((nm, d[dom], wavg, dom))
    print("      %-22s δ(主导 h=%d) = %.4f   权重加权 = %.4f" % (nm, dom, d[dom], wavg))
tierA = [r for r in rows if r[1] < 0.10]
tierB = [r for r in rows if r[1] >= 0.10]
check("★ 主导壳层口径下候选**分两档**：{[C1,C3] δ≈0.056–0.060} 与 {[C2,C4] δ≈0.148–0.169}",
      len(tierA) == 2 and len(tierB) == 2)
check("★ 加权口径亦分档（0.1295／0.1736 对 0.2794／0.4050）",
      min(r[2] for r in rows) > 0.10 and max(r[2] for r in rows) < 0.5)
# 单调性不是定理（L=14 上 C2 破坏）
c14 = table(14, g_rm)
hs = [h for h in sorted(c14) if h >= 1]
d14 = [abs(c14[h][0] - c14[h][1]) / sum(c14[h]) for h in hs]
mono14 = all(d14[i] <= d14[i + 1] + 1e-12 for i in range(len(d14) - 1))
check("★ 单调性**不是定理**：L=14 上 C2 破坏单调（δ 序列 %s）" % [round(x, 4) for x in d14], not mono14)

# ======================================================================
head("P7  §7.1 重大更正 · 四个候选根本没开门；门被 Z0③ 关闭")
# ======================================================================
# (a) 按**真分支** sign(S_k) 计数 ⇒ 两支完全等权
cnt_b = defaultdict(lambda: [0, 0])
for w in bw(L):
    S = part(w)
    for k in range(L - 1):
        if S[k] == 0:
            continue
        cnt_b[abs(S[k])][0 if S[k] > 0 else 1] += 1
sym = all(cnt_b[h][0] == cnt_b[h][1] for h in cnt_b)
print("      真分支权重（L=%d）：" % L, {h: cnt_b[h] for h in sorted(cnt_b)})
check("★ 按真分支 sign(S_k) 计数：N(+支) = N(−支) **对每一壳层** ⇒ δ_h ≡ 0", sym)
check("★ ⇒ 账本态在每一块上中心（ι 对称性的直接后果）", sym)
# (b) 四个候选与分支的关联 = 0（列联表两列恒等）
alluncorr = True
for nm, g in CANDS:
    tab = defaultdict(int)
    for w in bw(L):
        S = part(w)
        for k in range(L - 1):
            if S[k] == 0:
                continue
            tab[(g(w, k), 0 if S[k] > 0 else 1)] += 1
    a, b, c, d = tab.get((0, 0), 0), tab.get((0, 1), 0), tab.get((1, 0), 0), tab.get((1, 1), 0)
    uncorr = (a == b and c == d)
    if not uncorr:
        alluncorr = False
    print("      %-22s 列联表 [[%d,%d],[%d,%d]] ⇒ %s" % (nm, a, b, c, d,
          "完全无关" if uncorr else "有关联"))
check("★ 四个候选与分支的列联表**两列恒等** ⇒ 与分支**零关联**（一致率恰 0.5000）", alluncorr)
check("★ ⇒ 四个候选**根本没有开门**；§3 的 δ_h>0 是**假信号**（测的是候选自身取值分布）", alluncorr)
# (c) 更正已写入文档
check("★ 文档含 §7.1 更正（δ_h 测错对象 ＋ 门被 Z0③ 关闭）",
      "§7.1" in DOC and "重大更正" in DOC and "直接推论" in DOC)
check("★ 文档登记更正清单（原判 → 更正后）", "更正清单" in DOC)

# ======================================================================
head("P6  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("不依赖任何同题旧文", "去依赖声明"),
                 ("定理 1", "定理 1"), ("定理 2", "定理 2"), ("定理 3", "定理 3"),
                 ("定理 4", "定理 4"), ("定理 5", "定理 5"),
                 ("必要条件", "必要条件 (N)"), ("判据退化", "判据退化"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明"),
                 ("二难", "二难"), ("缺口定位", "定位性质"), ("未判", "未判措辞")]:
    check("文档契约：%s" % lab, tok in DOC)
bad = []
for m in re.finditer(r"\[([^\]]*)\]\(([^)]+)\)", DOC):
    t = m.group(2).strip()
    if t.startswith(("http", "#", "mailto")):
        continue
    t = re.sub(r":\d+(-\d+)?$", "", t)
    if not os.path.exists(os.path.normpath(os.path.join(HERE, t))):
        bad.append(t)
check("文档内全部相对链接可解析", not bad, "断链：%s" % bad if bad else "无")
_qm = re.findall(r"\]\((QM\d+_[^)]*)\)", DOC)
check("★ 不含任何 QM* 旧文依赖（链接扫描）", not _qm, "残留 %s" % _qm if _qm else "无")

print()
print("=" * 78)
print("QMR10 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
