#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""oldcorpus_check.py —— Z138「引用纪律：正文不引旧库」的机械检查器。

纪律：旧库的推导与结果**可以照用**，但正文（根目录 Z*.md）里不得出现旧库的
      编号、篇名、章节，也不得用「旧库」这个词当来源。

例外：Z138 自身（要引述旧库编号才能定义它们）。

用法：
  python3 _staging/oldcorpus_check.py              # 全库检查（退出码 0=干净，1=有命中）
  python3 _staging/oldcorpus_check.py Z043_*.md    # 只查指定篇
  python3 _staging/oldcorpus_check.py --summary    # 只打印每篇命中数
"""
import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXEMPT = {"Z138_引用纪律：正文不引旧库.md"}

# ── 旧库编号（Z138 §三）：注意保留例外 ──
#   R1–R9 是形式化守则号、E1–E9 是核查标签 ⟹ 都要求两位以上数字
PATTERNS = [
    (r"\bGQR\d+\b", "旧库编号 GQR*"),
    (r"\bQMR\d+\b", "旧库编号 QMR*"),
    (r"\bG\d+\b", "旧库编号 G*"),
    (r"\bR\d{2,}\b", "旧库编号 R*（两位以上）"),
    (r"\bSMR\d+\b", "旧库编号 SMR*"),
    (r"\bSM_\w+", "旧库文件名 SM_*"),
    (r"\bU\d+\b", "旧库编号 U*（注意与 U(1) 区分）"),
    (r"\bD_arc\b|\bD\d{2,3}\b", "旧库/讨论稿编号 D*"),
    (r"\bE\d{2,}\b", "讨论稿编号 E*（两位以上）"),
]
WORDS = ["旧库", "旧语料", "旧文章", "旧理论", "第一代台账"]


def scan(path: pathlib.Path):
    hits = []
    for i, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
        found = []
        for pat, tag in PATTERNS:
            for m in re.finditer(pat, line):
                found.append((m.group(0), tag))
        for w in WORDS:
            if w in line:
                found.append((w, "来源指示词"))
        if found:
            toks = "、".join(sorted({t for t, _ in found}))
            hits.append((i, toks, line.strip()[:110]))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*", help="默认：根目录全部 Z*.md")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()

    files = [pathlib.Path(f) for f in a.files] or sorted(ROOT.glob("Z[0-9][0-9][0-9]_*.md"))
    total, dirty = 0, 0
    for p in files:
        if p.name in EXEMPT:
            if not a.summary:
                print(f"（例外，跳过）{p.name}")
            continue
        hits = scan(p)
        if not hits:
            continue
        dirty += 1
        total += len(hits)
        if a.summary:
            print(f"{len(hits):4d}  {p.name}")
        else:
            print(f"\n━━ {p.name}  （{len(hits)} 处）")
            for i, toks, frag in hits:
                print(f"  L{i:<5d} [{toks}]  {frag}")
    print(f"\n合计：{dirty} 篇 / {total} 处命中（检查 {len(files)} 篇；例外 {len(EXEMPT)} 篇）")
    print("判据：命中后人工判一次「是引述还是引用」（Z138 §三）")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
