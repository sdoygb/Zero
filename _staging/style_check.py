#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""style_check.py —— 工作纪律《文档纪律：正文不留工作痕迹》的机械检查器。

用途：扫根目录 Z*.md，命中"工作痕迹"词即报（行号 + 词 + 该行片段）。
纪律：正文只留两类内容 —— ①已证的事实与证明；②未证的开放项。
      违禁的是**版本比较标签**与**过程叙述**；**认识论标签**（已证／条件／具名输入／开放…）保留。

例外：该纪律文档自身（要引述违禁词才能定义它们）。

用法：
  python3 _staging/style_check.py              # 全库检查（退出码 0=干净，1=有命中）
  python3 _staging/style_check.py Z002_*.md    # 只查指定篇
  python3 _staging/style_check.py --summary    # 只打印每篇命中数
"""
import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXEMPT = {"文档纪律_正文不留工作痕迹.md"}

# ── 违禁词表（见该纪律文档 §五）：版本比较标签 ＋ 过程叙述 ＋ 内部台账 ──
FORBIDDEN = [
    # 版本比较标签
    "已收窄", "已被收窄", "已被取代", "升档", "降档", "不再承重", "仍成立",
    # 过程叙述
    "本次修订", "状态更新", "后期篇", "本轮", "第一版", "早先版本",
    "口径提醒", "据实", "不美化", "诚实边界", "不主张",
    # 内部台账
    "落点", "归属与闭合判据", "审计报告",
    # 注意：`档位` 是认识论词（结论的地位），不在违禁表；被禁的是它的版本比较用法（升档／降档／档位由…→）

    # 更正（对旧版本的）
    "一处更正", "更正",
]
# ── 保留词（认识论标签）：出现在同一行不改变判定，仅用于人工复核时对照 ──
KEEP = ["已证", "导出", "条件", "具名输入", "开放", "不可导出", "排除", "数值律", "数学定理", "约定", "档位"]


def scan(path: pathlib.Path):
    hits = []
    for i, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
        for w in FORBIDDEN:
            if w in line:
                hits.append((i, w, line.strip()[:110]))
                break
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
            for i, w, frag in hits:
                print(f"  L{i:<5d} [{w}]  {frag}")
    print(f"\n合计：{dirty} 篇 / {total} 处命中（检查 {len(files)} 篇；例外 {len(EXEMPT)} 篇）")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
