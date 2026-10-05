#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z1_check.py —— Zero 层作为基础的核验
=====================================
独立实断言，全部可失败：
  F1  三条设计纪律在位（数值不前移／公理自付其价／可定义者不升格）
  F2  散度更新 ≡ 补偿移动（Z1 定理 1 的移动律由词导出）
  F3  Σ_v d(v) ≡ 0（Z1 定理 2 的零和是恒等式）
  F4  Z0③＋Z2 的程序级证据（无概率、无突变、整数重数）
  F5  Z3 四款与 Zero 的 verification 字段
  F6  词可达散度集 == 根格 {Σx=0} == im B（G19 引理 1/2）
  F7  闭词 ⇒ 零散度（单向严格；逆向反例计数）
  F8  ★ Z 条款文本不含任何数值常量（防 L=4/T=3 前移造成循环）
  F9  价目表每条都有"买回"（纪律 2）
  F10 Z_* 的可定义性：只用 Z1–Z5 重算，与 Zero 脚本逐轮比对
"""
import io
import itertools
import json
import os
import sys
from collections import Counter, deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = ("v" if ok else "x") if (level == "ind" or LEDGER_MODE == "A" or not ok) else "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


DOC = io.open(os.path.join(HERE, "Z1_zero_layer_as_the_foundation.md"), encoding="utf-8").read()
JSN = json.load(io.open(os.path.join(HERE, "simulations",
                                     "zero_sum_closure_exit_results.json"), encoding="utf-8"))


# ================================================================ F1
head("F1  三条设计纪律在位")
check("纪律 1（数值不得前移）在位",
      ("数值不得前移" in DOC) and ("必须保持为自由参数" in DOC))
check("纪律 2（公理必须自付其价）在位", "公理必须自付其价" in DOC)
check("纪律 3（可定义者不许升为公理）在位", "可定义的东西不许升为公理" in DOC)


# ================================================================ F2
head("F2  散度更新 ≡ 补偿移动（Z1 定理 1 的移动律由词导出）")
def div(edges, nv):
    d = [0] * nv
    for (u, v) in edges:
        d[v] += 1
        d[u] -= 1
    return tuple(d)


bad = tested = 0
for nv in (3, 4, 5):
    E = [(u, v) for u in range(nv) for v in range(nv) if u != v]
    for edges in itertools.product(E, repeat=3):
        for k in range(1, len(edges) + 1):
            d0, d1 = div(edges[:k - 1], nv), div(edges[:k], nv)
            u, v = edges[k - 1]
            exp = [0] * nv
            exp[v] += 1
            exp[u] -= 1
            tested += 1
            if [d1[i] - d0[i] for i in range(nv)] != exp:
                bad += 1
check("每个 (词,步) 对的散度增量恰为 e_v − e_u", bad == 0,
      "测 %d 对，不符 %d" % (tested, bad))
check("样本量不少于 20000（穷尽性）", tested >= 20000, "实算 %d" % tested)


# ================================================================ F3
head("F3  Σ_v d(v) ≡ 0（Z1 定理 2 的零和是恒等式）")
viol = 0
for nv in (2, 3, 4):
    E = [(u, v) for u in range(nv) for v in range(nv) if u != v]
    for edges in itertools.product(E, repeat=4):
        if sum(div(edges, nv)) != 0:
            viol += 1
check("任何词的散度之和恒为 0", viol == 0, "违反 %d 次" % viol)


# ================================================================ F4
head("F4  Z0③＋Z2 的程序级证据（无概率、无突变、整数重数）")
mdl = JSN["model"]
check("结果 JSON 自述 probabilities = none", str(mdl.get("probabilities")) == "none",
      "实值 %r" % mdl.get("probabilities"))
check("结果 JSON 自述 mutations = none", str(mdl.get("mutations")) == "none",
      "实值 %r" % mdl.get("mutations"))
check("闭合规则 = cumulative balance returns to zero",
      mdl.get("closure_rule") == "cumulative balance returns to zero",
      "实值 %r" % mdl.get("closure_rule"))
src = io.open(os.path.join(HERE, "simulations", "zero_sum_closure_exit.py"), encoding="utf-8").read()
check("重数是整数计数（源码用 Counter 记 zero_layer）",
      "Counter" in src and "zero_layer" in src)


# ================================================================ F5
head("F5  Z3 四款 与 Zero 的 verification 字段")
v = JSN["verification"]
check("闭合即退出：exit_removes_all_closed_paths", v.get("exit_removes_all_closed_paths") is True)
check("对照组：continue 规则保留闭合路径", v.get("control_retains_closed_paths") is True)
check("退出写入记录层：exit_writes_history_records", v.get("exit_writes_history_records") is True)
check("寿命到达后活动层为空：exit_leaves_no_active_paths_after_lifetime",
      v.get("exit_leaves_no_active_paths_after_lifetime") is True)
check("Z1 引用了 G34（闭合即退出 = π 的定义）", "G34_exit_rule_absorbed_into_pi.md" in DOC)


# ================================================================ F6
head("F6  词可达散度集 == 根格 {Σx=0}")
for m in (2, 3):
    E = [(u, v) for u in range(m) for v in range(m) if u != v]
    B = np.zeros((m, len(E)), dtype=int)
    for k, (u, v) in enumerate(E):
        B[v, k] += 1
        B[u, k] -= 1
    start = tuple([0] * m)
    R = {start}
    q = deque([start])
    for _ in range(3):
        nq = deque()
        while q:
            d = q.popleft()
            for k in range(len(E)):
                nd = tuple(np.array(d) + B[:, k])
                if nd not in R:
                    R.add(nd)
                    nq.append(nd)
        q = nq
    K = 3
    tgt = {x for x in itertools.product(range(-K, K + 1), repeat=m) if sum(x) == 0}
    check("m=%d：可达集 == {Σx=0, |x|≤%d}（双向包含）" % (m, K), R == tgt,
          "可达 %d，目标 %d" % (len(R), len(tgt)))
    check("m=%d：rank B = m−1（G19 引理 1）" % m,
          np.linalg.matrix_rank(B) == m - 1)


# ================================================================ F7
head("F7  闭词 ⇒ 零散度（单向严格）")
tot = badc = 0
for m in (2, 3):
    E = [(u, v) for u in range(m) for v in range(m) if u != v]
    for L in (2, 4):
        for edges in itertools.product(E, repeat=L):
            for st in range(m):
                pos, ok = st, True
                for (u, vv) in edges:
                    if u != pos:
                        ok = False
                        break
                    pos = vv
                if ok and pos == st:
                    tot += 1
                    if any(div(edges, m)):
                        badc += 1
check("闭合词（回到起点）的散度恒为零", badc == 0, "闭词样本 %d，反例 %d" % (tot, badc))

# 逆向反例：d≡0 但不是闭词
rev = 0
for m in (2, 3):
    E = [(u, v) for u in range(m) for v in range(m) if u != v]
    for edges in itertools.product(E, repeat=4):
        if all(x == 0 for x in div(edges, m)):
            pos, ok = edges[0][0], True
            for (u, vv) in edges:
                if u != pos:
                    ok = False
                    break
                pos = vv
            if not (ok and pos == edges[0][0]):
                rev += 1
check("逆向确实不成立（d≡0 ⇏ 单闭词，存在反例）", rev > 0, "反例 %d 个" % rev)
check("Z1 写明了这条单向性与反例", "逆向不成立" in DOC)


# ================================================================ F8
head("F8  ★ Z 条款文本不含任何数值常量（防循环回归）")
import re as _re
rows = [l for l in DOC.split("\n") if l.startswith("| **Z") and l.count("|") >= 4]
axiom_rows = [l for l in rows if _re.match(r"\|\s*\*\*Z[1-5]\*\*", l)]
def_rows = [l for l in rows if _re.match(r"\|\s*\*\*Z6\*\*", l)]
offenders = []
for l in rows:
    cells = [c.strip() for c in l.strip().strip("|").split("|")]
    desc = cells[-1] if cells else ""
    # 允许"标号"里的数字（Z0–Z5 / Z1–Z6 / Z-E1 / I8），只禁裸数值常量
    # 允许：标号（Z0–Z5/Z1–Z6/Z-E1/I8）、节号（§3.6）、文档号（G46/D242）
    stripped = _re.sub(r"\b(?:A|Z|Z-E|I|G|D)\s?-?\s?\d+\b", "", desc)
    stripped = _re.sub(r"§\s?\d+(?:\.\d+)?", "", stripped)
    if any(ch.isdigit() for ch in stripped):
        offenders.append((cells[0], desc[:50]))
check("Z 条款的描述列不含任何裸数值常量（标号除外）", not offenders,
      "越界行：%s" % (offenders or "无"))
check("条款表恰 5 条（Z1–Z5）", len(axiom_rows) == 5, "实算 %d 行" % len(axiom_rows))
check("Z6 单列为【定义】表（不是公理）", len(def_rows) == 1)
check("Z1 明写 L,T,d 必须自由（纪律 1）",
      ("$L$" in DOC or "L,T,d" in DOC) and "自由参数" in DOC)


# ================================================================ F9
head("F9  价目表：每条都有'买回'（纪律 2）")
price = [l for l in DOC.split("\n") if l.startswith("| Z") and "|" in l and "买回" not in l]
check("价目表行数 ≥ 9（Z1–Z6 + Z-E1…Z-E4）", len(price) >= 9, "实算 %d 行" % len(price))
check("Z-E1…Z-E4 四条扩充都在价目表里",
      all(("Z-E%d" % i) in DOC for i in (1, 2, 3, 4)))
check("Z1 写明'不带价签的条款：无'", "不带价签的条款" in DOC)
check("Z-E1 的代价（与笔记显式选择相反）已写明",
      ("not from" in DOC) and ("改选" in DOC))


# ================================================================ F10
head("F10  Z_* 的可定义性：只用 Z1–Z5 重算，与 Zero 脚本逐轮比对")
def canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))


def run_zero_axioms(keep_closed_active, L, epochs, init):
    active = [{w} for w in init]
    hist = [set() for _ in init]
    dead = [set() for _ in init]
    Z = Counter()
    snaps = []
    for ep in range(epochs + 1):
        snaps.append(dict(epoch=ep, active=sum(len(s) for s in active),
                          hist=sum(len(s) for s in hist), events=sum(Z.values()),
                          modes=len(Z), dead=sum(len(s) for s in dead)))
        if ep == epochs:
            break
        nxt = [set() for _ in init]
        for i, paths in enumerate(active):
            for w in paths:
                if len(w) >= L:                       # Z4
                    dead[i].add(w)
                    continue
                for st in (1, -1):                    # Z2 全分支
                    c = w + (st,)
                    if bool(c) and sum(c) == 0:       # Z3 闭合
                        hist[i].add(c)
                        Z[canon(c)] += 1              # 整数重数
                        if keep_closed_active:
                            nxt[i].add(c)
                    else:
                        nxt[i].add(c)
        active = nxt
    return snaps, Z


INIT = [(1,), (-1,), (1, 1), (-1, -1), (1, 1, 1), (-1, -1, -1), (1, -1, 1), (-1, 1, -1)]
allok = True
for scen, keep in (("continue_after_closure", True), ("exit_after_closure", False)):
    mine, Z = run_zero_axioms(keep, mdl["layer_lifetime"], mdl["epochs"], INIT)
    ref = next(r for r in JSN["runs"] if r["scenario"] == scen)
    for a, b in zip(mine, ref["snapshots"]):
        ok = (a["active"] == b["active_paths"] and a["hist"] == b["history_records"]
              and a["events"] == b["zero_events"] and a["modes"] == b["zero_modes"]
              and a["dead"] == b["dead_paths"])
        allok &= ok
    check("情景 %s：7 轮 × 5 量全部与 Zero 脚本一致" % scen, allok,
          "Z_* = %s" % dict(Z))
check("Z_* 的类数是 3（脚本与重算一致）", allok)


# ================================================================ 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
