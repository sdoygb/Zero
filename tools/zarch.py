#!/usr/bin/env python3
"""zarch —— 跨旧树检索（五个副目录 ＋ 本库）。

为什么有这个工具：五棵旧理论树里有大量**没有迁进本库**的推导与计算，
`zfind` 只索引本库，只查本库就下结论会漏（2026-11-07 用户当面指出）。
开工第四步、以及任何"有没有／以前算过吗"的问题，先跑本脚本。

用法：
    python3 tools/zarch.py "极端 天气"           # 全部旧树 ＋ 本库，任一命中
    python3 tools/zarch.py "ENSO" --tree 0401    # 只查一棵树
    python3 tools/zarch.py "降尺度" --all        # 全部词都出现才算命中
    python3 tools/zarch.py "oni" --names         # 只列命中文件（不列片段）
    python3 tools/zarch.py "气候" --count        # 只要计数
    python3 tools/zarch.py "climate" --no-vault  # 不扫本库

口径：按文件内容做**词法**匹配（大小写不敏感），不是语义检索——
命中 0 只说明"这套词没出现"，不说明"没有内容"；换词再查一次。
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
VAULT = HERE.parent

# 五棵旧树（与 AGENTS.md §0 的四步、§6 的表同一口径）
ROOTS = {
    "0401": "/Volumes/download/from_system_20261106/A_archives/Zero260401",
    "0901": "/Volumes/download/from_system_20261106/A_archives/Zero260901",
    "1007": "/Volumes/download/from_system_20261106/A_archives/Zero261007",
    "1010": "/Volumes/download/from_system_20261106/A_archives/Zero261010",
    "1030": "/Volumes/download/from_system_20261106/A_archives/Zero261030",
}

EXTS = {".md", ".py", ".ipynb", ".tex", ".txt", ".lean", ".json", ".csv", ".yaml", ".yml"}

# 依赖、构建产物、导出副本：跳过多半是别人写的代码，不是推导
SKIP_DIRS = {
    ".git", ".venv", "venv", "site-packages", "node_modules", "__pycache__",
    "build", "dist", ".mypy_cache", ".pytest_cache", "pdf_export",
    ".pnpm-store", ".cache", "EFI回滚备份", "Users", ".idea", ".vscode",
    "检索台账",  # 台账里是查询词本身，扫进来会自我命中
}

MAX_BYTES = 4 * 1024 * 1024
SNIP = 110


def iter_files(root: Path, exts):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for name in filenames:
            if name in ("zarch.py",):
                continue
            p = Path(dirpath) / name
            if p.suffix.lower() not in exts:
                continue
            try:
                if p.stat().st_size > MAX_BYTES:
                    continue
            except OSError:
                continue
            yield p


def _rx(term: str) -> str:
    """ASCII 词加词边界（否则 enso 会命中 tensor）；含中文的词按子串。"""
    if term.isascii():
        return r"(?<![A-Za-z0-9_])" + re.escape(term) + r"(?![A-Za-z0-9_])"
    return re.escape(term)


def run_tree(key: str, root: Path, patterns, need_all, top, names_only, count_only):
    files = hits = 0
    snippet_bag = []
    name_hits = []
    for p in iter_files(root, EXTS):
        files += 1
        try:
            text = p.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        matched_terms = set()
        local = []
        for i, line in enumerate(text.splitlines(), 1):
            for term, rx in patterns:
                if rx.search(line):
                    matched_terms.add(term)
                    local.append((i, line.strip()))
                    break
        if not local:
            continue
        if need_all and len(matched_terms) < len(patterns):
            continue
        hits += len(local)
        name_hits.append((p, len(local)))
        if len(snippet_bag) < top:
            snippet_bag.append((p, local[:2]))
    return files, hits, name_hits, snippet_bag


def main() -> int:
    ap = argparse.ArgumentParser(description="跨旧树检索（五棵副目录 ＋ 本库）")
    ap.add_argument("query", help="空格分隔的关键词（默认任一命中）")
    ap.add_argument("--tree", action="append", choices=sorted(ROOTS) + ["vault"],
                    help="只查指定树，可重复；默认五棵旧树 ＋ 本库")
    ap.add_argument("--no-vault", action="store_true", help="不扫本库")
    ap.add_argument("--all", action="store_true", help="全部关键词都出现才算命中")
    ap.add_argument("--top", type=int, default=4, help="每棵树列多少条片段（默认 4）")
    ap.add_argument("--names", action="store_true", help="只列命中文件")
    ap.add_argument("--count", action="store_true", help="只报计数")
    args = ap.parse_args()

    terms = [t for t in args.query.split() if t]
    if not terms:
        ap.error("query 不能为空")
    patterns = [(t, re.compile(_rx(t), re.IGNORECASE)) for t in terms]

    targets = []
    for key, path in ROOTS.items():
        if args.tree and key not in args.tree:
            continue
        targets.append((key, Path(path)))
    if not args.no_vault and (not args.tree or "vault" in args.tree):
        targets.append(("本库", VAULT))

    print(f"zarch · 词={terms} · 模式={'AND' if args.all else 'OR'} · 树={len(targets)}")
    total_hits = 0
    zero_trees = []
    for key, root in targets:
        if not root.exists():
            print(f"\n[{key}] 不存在：{root}")
            continue
        files, hits, name_hits, snips = run_tree(
            key, root, patterns, args.all, args.top, args.names, args.count)
        total_hits += hits
        flag = "0 命中" if hits == 0 else f"{hits} 行 / {len(name_hits)} 文件"
        print(f"\n[{key}] {root}\n  扫 {files} 文件 → {flag}")
        if hits == 0:
            zero_trees.append(key)
            continue
        if args.count:
            continue
        if args.names:
            for p, _n in name_hits[: args.top * 3]:
                print(f"  · {p.relative_to(root)}")
            if len(name_hits) > args.top * 3:
                print(f"  · … 另 {len(name_hits) - args.top * 3} 个文件")
            continue
        for p, loc in snips:
            print(f"  · {p.relative_to(root)}")
            for i, text in loc:
                print(f"      L{i}: {text[:SNIP]}")

    print(f"\n合计 {total_hits} 行命中。")
    if zero_trees:
        print(f"0 命中的树：{', '.join(zero_trees)} —— **0 命中 ≠ 没有**，换同义词再跑一次。")
    return 0 if total_hits else 1


if __name__ == "__main__":
    sys.exit(main())
