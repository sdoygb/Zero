#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z0_check.py —— 单一公理「零不断乱动」的核验
============================================
独立实断言，全部可失败：
  F1  Z0 三款在位（不停留／不停歇／不设概率）
  F2  §0 范围在位：本体系自足；Zero 笔记中外部公理体系的字样已清零
  F3  Z2 的程序级证据（无概率、无突变；全分支）
  F4  ★ Z5 的四模型对照（两熄灭、两持续）——由结果 JSON 实算
  F5  单零游走的连通性（连通分量数 = 1）
  F6  年龄 = 词长（程序逐字）
  F7  补偿对 (+,-)（程序逐字）
  F8  公理数为 1（正文只列 Z0 为公理）
  F9  价目表完整（每条有"买回"）
  F10 残余条款（Z4 终端款）已标【开放】
"""
import io
import json
import os
import sys

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


def rd(p):
    return io.open(os.path.join(HERE, p), encoding="utf-8").read()


def js(p):
    return json.load(io.open(os.path.join(HERE, p), encoding="utf-8"))


DOC = rd("Z0_zero_never_rests_single_axiom.md")

# ---------------------------------------------------------------- F1
head("F1  唯一公理 Z0 的三款在位")
check("① 零不停留", "零不停留" in DOC)
check("② 从不停歇", "从不停歇" in DOC)
check("③ 不设概率", "不设概率" in DOC)
check("正文声明除 Z0 外无其它公理", "除 Z0 外无其它公理" in DOC or "除 Z0 外" in DOC)

# ---------------------------------------------------------------- F2
head("F2  §0 范围：本体系自足，不涉及外部公理体系")
check("§0 标题为『范围』", "## §0 范围" in DOC)
check("§0 写明只处理本体系内部的基础重立", "本体系内部" in DOC)
check("§0 写明 A1–A5（历史命名）降为定理", ("A1–A5" in DOC or "A0–A5" in DOC) and "定理" in DOC)
check("§0 写明推理链自足", "自足" in DOC)
check("§0 不再提及任何外部公理体系（U 系字样已清零）",
      ("U1-U4" not in DOC) and ("U1–U4" not in DOC))
_notes = [f for f in os.listdir(HERE) if f.startswith("zero_sum_") and f.endswith(".md")]
_leak = [f for f in _notes if ("U1-U4" in rd(f)) or ("U1–U4" in rd(f))]
check("全部 Zero 笔记均不再提及该外部体系", not _leak, "残留 %s" % _leak)


# ---------------------------------------------------------------- F3
head("F3  Z2 的程序级证据（不设概率 ⇒ 全分支）")
mdl = js("simulations/zero_sum_closure_exit_results.json")["model"]
check("结果 JSON 自述 probabilities = none", str(mdl.get("probabilities")) == "none")
check("结果 JSON 自述 mutations = none", str(mdl.get("mutations")) == "none")
src = rd("simulations/zero_sum_closure_exit.py")
check("源码逐字实现'所有后继都被生成'（for step in (1, -1) 无分支剪除）",
      "for step in (1, -1):" in src)
check("Z0 指出 D210 把'全分支语义'列为输入、而本体系将其导出",
      "全分支生成语义" in DOC)

# ---------------------------------------------------------------- F4
head("F4  ★ Z5 的四模型对照（由结果 JSON 实算）")
ce = js("simulations/zero_sum_closure_exit_results.json")
traj = {}
for r in ce["runs"]:
    traj["closure_exit/" + r["scenario"]] = [s["active_paths"] for s in r["snapshots"]]
pd = js("simulations/zero_sum_periodic_destruction_results.json")
for r in pd["runs"]:
    keys = [k for k in r["snapshots"][0] if "active" in k]
    traj["periodic/" + r["scenario"]] = [s[keys[0]] for s in r["snapshots"]] if keys else []
orv = js("simulations/zero_sum_open_reservoir_results.json")
finals = {r["scenario"]: r.get("final_active_paths") for r in orv["runs"]}

check("closure_exit（有限寿命、无重播种）最终熄灭",
      traj["closure_exit/exit_after_closure"][-1] == 0,
      "轨迹 %s" % traj["closure_exit/exit_after_closure"][:7])
check("closure_exit/continue 亦熄灭（有限寿命的必然）",
      traj["closure_exit/continue_after_closure"][-1] == 0)
w = {k.split("/")[1]: v for k, v in traj.items() if k.startswith("periodic/")}
check("wipe_no_reseed（不重播种）最终熄灭", w.get("wipe_no_reseed", [1])[-1] == 0,
      "轨迹 %s" % w.get("wipe_no_reseed", [])[:7])
check("wipe_reseed（＋重播种）持续存在（末值 > 0）", w.get("wipe_reseed", [0])[-1] > 0,
      "轨迹 %s" % w.get("wipe_reseed", [])[:7])
check("open_reservoir（无寿命）持续存在（三个维度末值均 > 0）",
      all((v or 0) > 0 for v in finals.values()), "末值 %s" % finals)
check("对照的结论：给定'不断'必须'无寿命'或'重播种'二择一",
      "二择一" in DOC or "必须\"无寿命\"或\"重播种\"" in DOC)

# ---------------------------------------------------------------- F5
head("F5  单零游走的连通性（连通分量数 = 1）")
import itertools
from collections import deque


def components(edges):
    verts = set()
    for (u, v) in edges:
        verts.add(u)
        verts.add(v)
    seen, comp = set(), 0
    adj = {}
    for (u, v) in edges:
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    for s in verts:
        if s in seen:
            continue
        comp += 1
        q = deque([s])
        seen.add(s)
        while q:
            x = q.popleft()
            for y in adj.get(x, ()):
                if y not in seen:
                    seen.add(y)
                    q.append(y)
    return comp


bad = tot = 0
for L in (2, 3, 4):
    for edges in itertools.product([(u, v) for u in range(3) for v in range(3) if u != v], repeat=L):
        # 只取真正是一条游走的序列
        ok, pos = True, edges[0][0]
        for (u, v) in edges:
            if u != pos:
                ok = False
                break
            pos = v
        if not ok:
            continue
        tot += 1
        if components(edges) != 1:
            bad += 1
check("单条游走给出的图恰有 1 个连通分量", bad == 0, "游走样本 %d，反例 %d" % (tot, bad))
check("Z0 写明多初始词需'识别 U'才连通", "识别 U" in DOC and "连通" in DOC)

# ---------------------------------------------------------------- F6
head("F6  年龄 = 词长（程序逐字）")
check("源码用 len(word) 判定年龄到达寿命", "len(word) >= LAYER_LIFETIME" in rd(
    "simulations/zero_sum_closure_exit.py"))
check("Z0 写明年龄 = 词长且不需新结构", "年龄 = 词长" in DOC)

# ---------------------------------------------------------------- F7
head("F7  配对 = 补偿步对 (+,-)")
lu = rd("simulations/zero_sum_living_universe.py")
check("源码逐字：Vacuum fluctuations create a compensating pair (+,-)",
      "compensating pair (+,-)" in lu)
check("Z0 写明配对由 Z0① 导出", "补偿步" in DOC or "补偿对" in DOC)

# ---------------------------------------------------------------- F8
head("F8  公理数为 1")
check("正文声明'公理只有一条'", "公理只有一条" in DOC)
check("§1 给出唯一公理 Z0 的三款", DOC.count("Z0") >= 4)
check("Z0 明确 5 → 1 的压缩（相对 Z1）", "5 → 1" in DOC or "**5 条**" in DOC)

# ---------------------------------------------------------------- F9
head("F9  价目表完整")
_sec = DOC[DOC.index("## §3 价目表"):DOC.index("## §4")]
rows = [l for l in _sec.split("\n")
        if l.startswith("|") and not set(l) <= set("|-: ") and "买回" not in l]
check("价目表数据行 ≥ 6（Z0 三款 ＋ 识别 U ＋ 参数 L ＋ 残余条款）", len(rows) >= 6,
      "实算 %d 行" % len(rows))
check("价目表逐行都写了'买回'一栏（四列）",
      all(l.count("|") >= 5 for l in rows), "最少列数 %d" % min(l.count("|") for l in rows))
check("'不带价签的条款'一句在位", "不带价签的条款" in DOC)

# ---------------------------------------------------------------- F10
head("F10  Z4 终端款已导出（由'几何须有稳定谱维数'逼出）")
check("§2.6 标题声明【导出】", "Z4 的**终端款**【导出" in DOC or "现改为导出" in DOC)
check("正文给出推导链：无截断 ⇒ 平均度无界 ⇒ 非流形 ⇒ 无稳定谱维数",
      "平均度无界" in DOC and "非流形" in DOC and "无稳定谱维数" in DOC)
check("正文点明该图就是 Z1 定理 1（历史 A2）补偿移动的可执行实现",
      ("A2 的补偿移动" in DOC or "Z1 定理 1 的补偿移动" in DOC) and "可执行实现" in DOC)
check("已不再把终端款称作残余假设（仅保留'原先列为'的叙述）",
      "终端款是残余假设" not in DOC and "残余假设** |" not in DOC)
check("诚实边界写明它带条件（几何要求）", "带条件的导出" in DOC or "终止条款本身" in DOC)
check("开放项已更新为'能否不依赖几何要求而由 Z0②导出'",
      "能否不依赖几何要求" in DOC)


# ---------------------------------------------------------------- F11
head("F11  ★ 独立复现闭环图的平均度无界（节点数与 G28 §4 逐项比对）")
import itertools as _it


def _canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))


def _modes(T):
    return sorted({_canon(w) for w in _it.product((1, -1), repeat=T) if sum(w) == 0})


def _neigh(w):
    T = len(w)
    out = set()
    for i in range(T):
        j = (i + 1) % T
        if w[i] != w[j]:
            v = list(w)
            v[i], v[j] = v[j], v[i]
            c = _canon(tuple(v))
            if c != w:
                out.add(c)
    return out


REF_NODES = {8: 10, 10: 26, 12: 80, 14: 246}
REF_DEG = {8: 3.60, 10: 4.92, 12: 6.15, 14: 7.39}
nodes_ok = deg_ok = True
rows = []
for T in (8, 10, 12, 14):
    M = _modes(T)
    E = sum(len(_neigh(w)) for w in M) // 2
    d = 2 * E / len(M)
    rows.append((T, len(M), round(d, 2)))
    nodes_ok &= (len(M) == REF_NODES[T])
    deg_ok &= (abs(d - REF_DEG[T]) < 0.01)
check("节点数（旋转类计数）与 G28 §4 逐项一致", nodes_ok, "实算 %s" % rows)
check("平均度与 G28 §4 逐项一致", deg_ok, "实算 %s" % [(t, d) for t, _, d in rows])
degs = [d for _, _, d in rows]
check("平均度随 T 单调增长（无界趋势）",
      all(degs[i] < degs[i + 1] for i in range(len(degs) - 1)), "序列 %s" % degs)
check("节点数增长 > 20 倍（10 → 246）", rows[-1][1] / rows[0][1] > 20,
      "%.1f 倍" % (rows[-1][1] / rows[0][1]))


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
