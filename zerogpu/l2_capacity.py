"""
l2_capacity.py --- L2 的**饱和机制**：把「容量」与「寿命」分离，做字级精确实验。

**为什么需要这个实验（对既有语料的诊断）**

语料里两个模型给出相反的极端，但它们**同时在两个维度上不同**，所以从未被隔离：

| 模型 | 寿命 | 重数语义 | 结果 |
|:--|:--|:--|:--|
| `simulations/zero_sum_closure_exit.py` | 有 $L{=}4$ | `set`（= 每词容量 1） | 熄灭 |
| `L2_period.py:52-92` | **无** | `Counter`（无容量） | 逐周期指数发散 |

=> "哪一个因素造成了差别"**在语料里答不出来**。本脚本把它拆成可独立开关的三项：

    CAP     每词容量 K ∈ {1, 2, 4, 8, 16, ∞}
    LIFETIME  到龄即入终端：L ∈ {8, 12, 16, 20} 或 无
    MEMORY    历史层保留的层数 k ∈ {2, 3, ∞}（D222 的有限记忆）

**字级精确、无粗粒化**（这是与 `L2_period.py` 的关键区别，后者把状态压在
(平衡, 奇偶) 桶上，指数增长有相当部分是粗粒化的产物）：

    长度 a 的状态 = 长度 2^a 的整数计数向量，下标即 ±1 词的位掩码。
    一步：new = repeat(counts, 2)      —— 全分支，精确整数重数，纯 numpy 向量化
    闭合：new 中 popcount == (a+1)/2 的分量（此时词和恰为 0）-> 写历史，退出活动层
    寿命：a+1 >= L 的分量 -> 终端 D
    容量：new = minimum(new, K)
    周期末：清空 -> 从历史重播种（每条闭合历史 w 给 w+ 与 w- 两个种子）

**判据**：跨周期看活动层总量的比值 r_n = |E_{n+1}|/|E_n| 是否 →1（稳态），
是 → 熄灭（0），还是 →const>1（发散）。
"""
from __future__ import annotations

import json
import os
import sys
from math import comb

import numpy as np

try:
    import pyopencl as cl
    from zcl import Engine
    HAVE_CL = True
except Exception:
    HAVE_CL = False

INITIAL_WORDS = ("+", "-", "++", "--", "+++", "---", "+-+", "-+-")


# --------------------------------------------------------------- 词 -> 位掩码
def encode(text: str) -> int:
    v = 0
    for i, c in enumerate(text):
        if c == "+":
            v |= 1 << i
    return v


# ------------------------------------------------------- 字级一步（纯 numpy）
def build_pc_tables(L: int):
    """每个长度的 popcount 表（闭合判定用）。"""
    tabs = {}
    for a in range(1, L + 1):
        if a % 2 == 0:
            tabs[a] = np.bitwise_count(np.arange(1 << a, dtype=np.uint32)) == (a // 2)
    return tabs


def run(L: int | None, cap: int | None, memory: int | None, cycles: int,
        seed_cap: int = 1 << 24):
    """
    L=None 表示无寿命（词长无上界；此时以 max_len 截断，靠 cap 或 memory 收敛）。
    返回每周期的统计。
    """
    max_len = L if L is not None else 24
    pc = build_pc_tables(max_len)

    # 初始：8 个站点各持一个初始词
    counts = [None] * (max_len + 1)
    for t in INITIAL_WORDS:
        a = len(t)
        if counts[a] is None:
            counts[a] = np.zeros(1 << a, dtype=np.int64)
        counts[a][encode(t)] += 1
    # 活动层状态另存一份（周期末要清空重播种）
    active = [None if c is None else c.copy() for c in counts]
    history = [None] * (max_len + 1)          # 每个长度上的闭合词计数向量
    stats = []

    for cyc in range(1, cycles + 1):
        closed_this_cycle = 0
        peak = 0
        dead_total = 0
        for step in range(max_len):
            nxt = [None] * (max_len + 1)
            for a in range(max_len + 1):
                c = active[a]
                if c is None or not c.any():
                    continue
                b = a + 1
                if b > max_len:
                    dead_total += int(c.sum())
                    continue
                rep = np.repeat(c, 2)                     # 全分支
                m = pc.get(b)
                if m is not None:                         # 闭合判定
                    cl_n = int(rep[m].sum())
                    if cl_n:
                        closed_this_cycle += cl_n
                        if history[b] is None:
                            history[b] = np.zeros_like(rep)
                        history[b] += rep * m
                        rep = rep.copy()
                        rep[m] = 0
                if cap is not None:
                    np.minimum(rep, cap, out=rep)
                if b >= max_len:                          # 到龄 -> 终端
                    dead_total += int(rep.sum())
                else:
                    if nxt[b] is None:
                        nxt[b] = rep
                    else:
                        nxt[b] += rep
            active = nxt
            tot = sum(int(c.sum()) for c in active if c is not None)
            peak = max(peak, tot)
        destroyed = sum(int(c.sum()) for c in active if c is not None)

        # ---- 有限记忆：只保留最高的 memory 层 ----
        if memory is not None:
            keep = sorted(i for i in range(max_len + 1) if history[i] is not None
                          and history[i].any())[-memory:]
            history = [history[i] if i in keep else None for i in range(max_len + 1)]

        h = sum(int(c.sum()) for c in history if c is not None)

        # ---- 重播种：每条闭合历史 w 给 w+ 与 w- ----
        active = [None] * (max_len + 1)
        for a in range(max_len + 1):
            c = history[a]
            if c is None or not c.any():
                continue
            b = a + 1
            if b > max_len:
                continue
            nz = np.nonzero(c)[0]
            vals = np.minimum(c[nz], seed_cap)
            seeds = np.repeat(nz, 2) * 2 + np.tile(np.array([0, 1]), len(nz))
            sv = np.repeat(vals, 2)
            if active[b] is None:
                active[b] = np.zeros(1 << b, dtype=np.int64)
            np.add.at(active[b], seeds, sv)

        stats.append(dict(cycle=cyc, closed=closed_this_cycle, peak=peak,
                          destroyed=destroyed, history=h,
                          active_after_reseed=int(sum(int(c.sum()) for c in active
                                                      if c is not None)),
                          dead=dead_total))
    return stats


def classify(stats):
    """按末段活动层峰值比值判稳态/熄灭/发散。"""
    peaks = [s["peak"] for s in stats]
    tail = peaks[-max(3, len(peaks) // 4):]
    if all(p == 0 for p in tail):
        return "熄灭"
    ratios = [tail[i + 1] / tail[i] for i in range(len(tail) - 1) if tail[i] > 0]
    if not ratios:
        return "?"
    r = float(np.median(ratios))
    if r < 1.02:
        return "稳态"
    return f"发散(r≈{r:.2f})"


def main():
    cycles = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print("=" * 100)
    print(f"L2 容量/寿命/记忆 三维扫描（字级精确，无粗粒化，{cycles} 个周期）")
    print("=" * 100)
    rows = []
    for L in (8, 12, 16, 20, None):
        for cap in (1, 2, 4, 16, 256, None):
            for memory in (2, 3, None):
                st = run(L, cap, memory, cycles)
                cls = classify(st)
                peaks = [s["peak"] for s in st]
                rows.append(dict(L=L, cap=cap, memory=memory, verdict=cls,
                                 peaks=peaks, closed=[s["closed"] for s in st],
                                 history=st[-1]["history"]))
                print(f"L={str(L):>4} cap={str(cap):>5} k={str(memory):>4} -> "
                      f"{cls:>14}  peak: {peaks[0]:>7} ... {peaks[-1]:>10}  "
                      f"history={st[-1]['history']:>9} closed={st[-1]['closed']:>9}", flush=True)
        print("-" * 100, flush=True)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "l2_capacity.json"), "w") as f:
        json.dump({"cycles": cycles, "rows": rows}, f, indent=2)
    print("\nsaved results/l2_capacity.json")

    # 汇总表
    print("\n" + "=" * 100)
    print("汇总：行 = (寿命 L, 容量 K)，列 = 记忆 k")
    for L in (8, 12, 16, 20, None):
        print(f"\n  寿命 L = {L}")
        print(f"    {'容量 K':>8} | {'k=2':>16} | {'k=3':>16} | {'k=∞':>16}")
        for cap in (1, 2, 4, 16, 256, None):
            cells = []
            for m in (2, 3, None):
                r = next(x for x in rows if x["L"] == L and x["cap"] == cap
                         and x["memory"] == m)
                cells.append(r["verdict"])
            print(f"    {str(cap):>8} | {cells[0]:>16} | {cells[1]:>16} | {cells[2]:>16}")


if __name__ == "__main__":
    main()
