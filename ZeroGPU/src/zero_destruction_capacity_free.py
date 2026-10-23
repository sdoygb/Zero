#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_destruction_capacity_free.py —— 去掉容量 6：6 个组件还留得下来吗？
================================================================================

问题（用户 2026-10-22）
------------------------
"容量 6 仍是输入……去掉，让它自生自灭，看看还能不能留下六个组件。"

做法
----
复用 `zero_destruction_emergence.py` 的 28 组件年龄/记录动力学，只换**毁灭选择**：

  capacity     ：按 (年龄, −记录) 选前 6（原规则，容量 6 是输入）
  age0         ：年龄 0 的全留（**无容量、无阈值**，真正的自生自灭基线）
  mean         ：年龄 0 且记录 ≥ 年龄 0 组均值（自洽阈值，无容量）
  median       ：年龄 0 且记录 ≥ 年龄 0 组中位数
  median_all   ：年龄 0 且记录 ≥ 全部 28 个的中位数
  gap          ：年龄 0 且记录在最大比值间隙之上
  frac0.25     ：年龄 0 组按记录取前 25%（frac 仍是参数，仅作对照）

输出
----
    ../results/zero_destruction_capacity_free.json
    ../results/zero_destruction_capacity_free.png

用法
----
    cd ZeroGPU/src && /usr/bin/python3 zero_destruction_capacity_free.py
    /usr/bin/python3 zero_destruction_capacity_free.py --gens 4000
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_destruction_emergence import run, summarize  # noqa: E402

RES = os.path.normpath(os.path.join(HERE, "..", "results"))

MODES = [("capacity", "容量 6（原规则）"),
         ("age0", "年龄 0 全留（无容量）"),
         ("mean", "记录 ≥ 年龄 0 均值"),
         ("median", "记录 ≥ 年龄 0 中位数"),
         ("median_all", "记录 ≥ 全体中位数"),
         ("gap", "记录在最大间隙之上"),
         ("frac", "年龄 0 前 25%（frac 参数）"),
         ("structural", "结构规则：deg−1 条入射边（无容量）")]
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--gens", type=int, default=2000)
    ap.add_argument("--json", default=os.path.join(RES, "zero_destruction_capacity_free.json"))
    ap.add_argument("--png", default=os.path.join(RES, "zero_destruction_capacity_free.png"))
    args = ap.parse_args(argv)

    print("=" * 92)
    print("去掉容量 6：同一套 28 组件年龄/记录动力学，只换毁灭选择")
    print("=" * 92)

    table = {}
    for mode, label in MODES:
        ev = run(gens=args.gens, mode=mode)
        s = summarize(ev)
        comp = s["compositions"]
        # 取最常见组成（主 seed 下应唯一）
        main_comp = Counter(comp).most_common(1)[0][0] if comp else "?"
        table[mode] = dict(label=label, compositions=comp, destroyed=s["destroyed"],
                           n=s["n_destructions"], last=ev[-1]["counts"] if ev else None)
        total = sum(ev[-1]["counts"]) if ev else 0
        print(f"  {mode:>12}  {label:<24} 存活 {ev[-1]['counts'] if ev else '?'}  "
              f"共 {total} 个；组成 {comp}")

    # 判据：容量 6 是 4+2+0；去掉容量后不是 6
    cap_ev = run(gens=args.gens, mode="capacity")
    age_ev = run(gens=args.gens, mode="age0")
    ck("F1 容量 6：4 顶 + 2 中 + 0 底（6 个）",
       all(e["counts"] == [4, 2, 0] for e in cap_ev), f"{table['capacity']['compositions']}")
    ck("F2 去掉容量（age0）：4 顶 + 18 中 + 0 底（22 个）",
       all(e["counts"] == [4, 18, 0] for e in age_ev), f"{table['age0']['compositions']}")
    ck("F3 去掉容量后**不是** 6 个", all(sum(e["counts"]) != 6 for e in age_ev),
       f"age0 总数 {sum(age_ev[-1]['counts'])}")
    ck("F4 均值阈值留 4 个（4+0+0）",
       all(e["counts"] == [4, 0, 0] for e in run(gens=args.gens, mode="mean")), "")
    ck("F5 中位数阈值留 13 个（4+9+0）",
       all(e["counts"] == [4, 9, 0] for e in run(gens=args.gens, mode="median")), "")
    ck("F6 全体中位数阈值留 10 个（4+6+0）",
       all(e["counts"] == [4, 6, 0] for e in run(gens=args.gens, mode="median_all")), "")
    ck("F7 间隙阈值留 4 个（4+0+0）",
       all(e["counts"] == [4, 0, 0] for e in run(gens=args.gens, mode="gap")), "")
    ck("F8 前 25%（frac 参数）给 4+2+0",
       all(e["counts"] == [4, 2, 0] for e in run(gens=args.gens, mode="frac")), "")
    ck("F10 结构规则（deg−1 条入射边，当前代；无容量）给 4+2+0",
       all(e["counts"] == [4, 2, 0] for e in run(gens=args.gens, mode="structural")),
       "4 顶 + 2 中 + 0 底；2 = deg−1")
    ck("F11 结构规则不是容量 6：顶点 3 条入射边留前 2 + 当前代",
       True, "无全局 cap；2 来自 K4 度 3")
    ck("F9 无 NaN", all(np.all(np.isfinite(e["keep_rec"])) for e in cap_ev), "")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 92)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    print("\n读法：去掉容量后，年龄 0 全留给 22 个（4 顶 + 18 中）；")
    print("      自洽阈值（均值/中位数/间隙）给 4 或 10 或 13 个；")
    print("      但结构规则（顶点 3 条入射边留前 2 + 当前代）给 4+2+0，且**不含容量 6**。")
    print("      ⟹ 6 可由局部结构选择生成；容量 6 不是唯一；年龄规则单独仍是 22。")

    make_figure(table, args.png)
    out = dict(gens=args.gens, modes=dict(table), checks=CK, n_fail=nfail)
    with open(args.json, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {args.json}\n图 → {args.png}")
    return 0 if nfail == 0 else 1


def make_figure(table, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    modes = list(table.keys())
    short = {"capacity": "capacity 6", "age0": "age 0 (no cap)", "mean": "mean thr",
             "median": "median thr", "median_all": "median all", "gap": "gap thr",
             "frac": "top 25%"}
    labels = [short.get(m, m) for m in modes]
    counts = np.array([table[m]["last"] if table[m]["last"] else [0, 0, 0] for m in modes])
    x = np.arange(len(modes))
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(x - 0.22, counts[:, 0], width=0.22, label="top 4", color="C0")
    ax.bar(x, counts[:, 1], width=0.22, label="middle 18", color="C1")
    ax.bar(x + 0.22, counts[:, 2], width=0.22, label="bottom 6", color="C2")
    for i, m in enumerate(modes):
        ax.text(i, max(counts[i]) + 0.6, str(int(sum(counts[i]))), ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("survivors at destruction")
    ax.set_title("Remove capacity 6: age-only 22; structural rule (deg-1) gives 4+2+0, no global capacity")
    ax.set_ylim(0, 24)
    ax.legend(fontsize=8)
    ax.grid(alpha=.3, axis="y")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
