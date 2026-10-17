#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
boundary_scan.py —— 扫出根目录 Z 篇的全部「边界条款」并成表
=============================================================
用途：为「逐篇修订过时边界条款」提供清单（判定仍由人做，机器只负责不漏项）。
判据词（用户指定）：不主张／不声称／不得主张／否证／墙／已证不可导出／未证／开放
输出：_staging/边界条款清单.md（按篇号排序；每篇一节）

运行：cd ZeroGPU/src && /usr/bin/python3 boundary_scan.py
"""
import os
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
OUT = os.path.join(ROOT, "_staging", "边界条款清单.md")

# 判据词（分三档：硬边界 / 状态词 / 否证词）
HARD = ["不主张", "不声称", "不得主张", "不得作为", "不能声称"]
STATE = ["未证", "仍开放", "开放项", "已排除", "不可导出", "不能导出", "挂账", "条件证成"]
NEG = ["否证", "墙", "no-go", "不可能", "排除"]

PATS = [(w, "硬") for w in HARD] + [(w, "态") for w in STATE] + [(w, "否") for w in NEG]

# 与内容无关的行（目录行、公式标签、参考文献、指针行）跳过
SKIP = re.compile(r"^\s*[-*]?\s*(§|##|表|图|\||\$\$|\\tag)")


def scan(path):
    hits = []
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            s = line.rstrip("\n")
            if not s.strip():
                continue
            for w, tag in PATS:
                if w in s:
                    hits.append((i, tag, w, s.strip()))
                    break
    return hits


def main():
    files = sorted([f for f in os.listdir(ROOT)
                    if re.match(r"Z\d+_.*\.md$", f)],
                   key=lambda x: int(re.match(r"Z(\d+)_", x).group(1)))
    lines = ["# 边界条款清单（机检，判定由人做）", "",
             f"扫描对象：根目录 Z 篇共 **{len(files)}** 篇。",
             "判据词：硬＝不主张/不声称/不得主张；态＝未证/仍开放/开放项/已排除/不可导出/挂账；否＝否证/墙/no-go/不可能/排除。",
             ""]
    tot = 0
    per = {}
    for f in files:
        h = scan(os.path.join(ROOT, f))
        per[f] = h
        tot += len(h)
    lines.append(f"**合计命中 {tot} 行**。")
    lines.append("")
    lines.append("## 逐篇计数")
    lines.append("")
    lines.append("| 篇 | 命中 | 硬 | 态 | 否 |")
    lines.append("|:--|--:|--:|--:|--:|")
    for f, h in per.items():
        z = int(re.match(r"Z(\d+)_", f).group(1))
        a = sum(1 for x in h if x[1] == "硬")
        b = sum(1 for x in h if x[1] == "态")
        c = sum(1 for x in h if x[1] == "否")
        lines.append(f"| Z{z:03d} | {len(h)} | {a} | {b} | {c} |")
    lines.append("")
    lines.append("## 逐篇明细")
    for f, h in per.items():
        if not h:
            continue
        z = int(re.match(r"Z(\d+)_", f).group(1))
        lines.append("")
        lines.append(f"### Z{z:03d}　`{f}`")
        lines.append("")
        for i, tag, w, s in h:
            s = s if len(s) <= 150 else s[:147] + "…"
            lines.append(f"- L{i} 〔{tag}·{w}〕 {s}")
    with open(OUT, "w", encoding="utf-8") as g:
        g.write("\n".join(lines) + "\n")
    print(f"扫描 {len(files)} 篇，命中 {tot} 行 → {os.path.relpath(OUT, os.getcwd())}")


if __name__ == "__main__":
    main()
