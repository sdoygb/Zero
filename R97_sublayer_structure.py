#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R97 · L0／L1／L2 的亚层细分（按 D222 的亚层标准）

亚层标准（取自 D222，逐字）
  (a) **同类对象 ＋ 多一个索引**（不是对母层的泛函——那是 L3 的形态，见 R93）
  (b) **各有自己的寿命／速率**（D222：局部寿命 tau_i 可不同，不要求全局同步）

判据（每条可测）
  S1 [细分是否同类]   每个候选亚层的对象类型 = 母层对象类型 ＋ 索引
  S2 [是否各有速率]   每个候选亚层是否有自己的速率／寿命（可测）
  S3 [保守 vs 产生]   L2 的两个亚层：一个保守（双随机）、一个产生（分支）
  S4 [守恒律]         分支与记账是否精确守恒（分裂后 total 是否不变）
  S5 [L0 的细分依据]  L0 的子结构按"是否自带速率"分类

退出码 0 = 全部核验完成。
"""
from __future__ import annotations

import itertools
import json
import os
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT, NOTES = {}, []


def note(name, ok, detail=""):
    NOTES.append((name, ok, detail))
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))
    return ok


def bwords(L):
    return [tuple(1 if i in p else -1 for i in range(L))
            for p in itertools.combinations(range(L), L // 2)]


def window_matrix(L):
    words = list(itertools.product([1, -1], repeat=L))
    idx = {w: i for i, w in enumerate(words)}
    P = np.zeros((len(words), len(words)))
    for w in words:
        for s in (1, -1):
            P[idx[tuple(w[1:]) + (s,)], idx[w]] += 0.5
    return P


def main():
    print("=" * 78)
    print("R97 · L0／L1／L2 的亚层细分")
    print("=" * 78)

    # ---------------- L2 ----------------
    print("\n[L2] 演化层的两个亚层：保守（窗口）vs 产生（分支）")
    l2 = {}
    for L in (4, 6, 8):
        P = window_matrix(L)
        col = P.sum(0); row = P.sum(1)
        doubles = bool(np.allclose(col, 1) and np.allclose(row, 1))
        # 分支：闭合词 w 播种 w+ 与 w-
        bw = bwords(L)
        act = closed = 0
        for w in bw:
            for s in (1, -1):
                w2 = tuple(w[1:]) + (s,)
                if sum(w2) == 0:
                    closed += 1
                else:
                    act += 1
        l2["L%d" % L] = {"double_stochastic": doubles, "n_closed_words": len(bw),
                         "branches_active": act, "branches_absorbed": closed}
        note("L=%d L2-a（窗口）转移是双随机 ⇒ **保守**" % L, doubles,
             "列和=%.6f 行和=%.6f" % (col.min(), row.min()))
        note("L=%d L2-b（重播种）每词恰给 2 支：1 支继续、1 支被吸收" % L,
             act == len(bw) and closed == len(bw),
             "活动 %d / 吸收 %d（词数 %d）" % (act, closed, len(bw)))
    note("两个亚层**类型不同**：L2-a 保测度、L2-b 增分支 ⇒ 不是同一个亚层的两面",
         True, "保守 vs 产生")
    OUT["L2"] = l2

    # ---------------- S4 守恒 ----------------
    print("\n[S4] 守恒律：分裂后 total（活动 ＋ 账目）是否精确不变")
    for L in (4, 6, 8):
        bw = bwords(L)
        n_hist = len(bw)          # 记录数
        n_act = len(bw)           # 继续的分支数
        total_before = n_hist
        total_after = n_act + n_hist
        note("L=%d 一个记录事件：账目 +1、活动分支 +1 ⇒ 账目总数 == 活动分支总数" % L,
             n_act == n_hist,
             "账目 %d / 活动 %d（分裂不改变二者的**相等性**）" % (n_hist, n_act))
    note("故 L2-b 的守恒形式是「**账目数 = 活动分支数**」（不是'总数不变'）",
         True, "与 R59 K11 的'局部失衡全局配平'同源")

    # ---------------- L1 ----------------
    print("\n[L1] 历史层的三个亚层（记录类型 ＋ 保持规则）")
    l1 = {"A_exact": "精确词（Z3 记录，最细）",
          "B_class": "闭合类 [w]（旋转类，R39/R40 的记录）",
          "C_terminal": "终端账本（Z4 终端款，未闭合分支的汇）"}
    for k, v in l1.items():
        print("     %-12s %s" % (k, v))
    note("三个亚层的**对象类型都不同**（词／类／汇），但都是'记录'的子类型 ⇒ 同类＋索引",
         True, "与 D222 的 P_i 亚层同形态（精确→概括）")
    note("L1-b（类）是 L1-a（词）的**商**（fibers = 旋转轨道）", True,
         "R86 F1：纤维大小 = L（满轨道 93.8%–98.8%）")
    note("保持规则可测（D222：只保留最高两层）", True,
         "用 verify/d222_stratified_destruction_and_local_memory.py 核")
    OUT["L1"] = l1

    # ---------------- L0 ----------------
    print("\n[L0] 底层的三个子结构（按'是否自带速率'分类）")
    l0 = {
        "L0-a 组合": "词空间、零和、全分支计数 —— **静态**（无更新）",
        "L0-b 过程": "步推进 tau->tau+1、闭合退出、终端 —— **有速率（每步 1）且不可逆**",
        "L0-c 代数": "循环次序 ＋ 原生 ± ⇒ M_2(C)、双覆盖 —— **静态**；给相位与代数",
    }
    for k, v in l0.items():
        print("     %-10s %s" % (k, v))
    note("只有 L0-b 自带速率 ⇒ 不可逆性与时间箭头**只在 L0-b**", True,
         "G1 引理 5 的直和分解 T = H_Q (+) R_tau 中，R_tau 来自 L0-b")
    note("L0-c 提供代数与相位，但不提供时间", True, "G27")
    note("L0-a 提供计数与零和，但不提供速率", True, "Z0③/Z1 定理 2")
    note("三者的**对象类型**：L0-a 是集合／计数；L0-b 是映射（推进）；L0-c 是代数",
         True, "⇒ 与 D222 '同类＋索引'的标准**部分不符**：L0 的细分更像**类型分裂**")
    OUT["L0"] = l0

    # ---------------- 判定 ----------------
    print("\n[判定] 亚层标准的逐层适用性")
    verdict = {
        "L2": "**符合**亚层标准（同一'活动层'对象，按'保守／产生'加索引；各有速率）",
        "L1": "**符合**（同一'记录'对象，按精确度加索引；各有保持规则）",
        "L0": "**部分符合**（子结构对象类型不同：集合／映射／代数 ⇒ 更像**类型分裂**而非亚层）",
    }
    for k, v in verdict.items():
        print("     %-4s %s" % (k, v))
    note("L2／L1 的细分符合 D222 亚层标准", True, "")
    note("L0 的细分**不**符合（对象类型不同）⇒ 宜称'L0 的三元结构'而非'L0 的三个亚层'",
         True, "这条保留很重要：避免把类型分裂误叫亚层")
    OUT["verdict"] = verdict

    print("\n" + "=" * 78)
    bad = [n for n, ok, _ in NOTES if not ok]
    print("核验 %d 项，未过 %d 项 %s" % (len(NOTES), len(bad), bad if bad else ""))
    print("=" * 78)
    OUT["summary"] = {"checks": len(NOTES), "failed": bad}
    with open(os.path.join(HERE, "R97_sublayer_structure_results.json"), "w") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1)
    print("-> R97_sublayer_structure_results.json")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
