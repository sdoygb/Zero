#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
latex_fix.py -- 修复 Markdown 中的 LaTeX 定界符，使 KaTeX(rehype-katex) 能渲染。

修复项：
  F1  \\(...\\)  ->  $...$      （行内；跨行时转为 $$ 块）
  F2  \\[...\\]  ->  $$...$$    （块级）
  F3  表格行内数学里的 |  -> \\vert （避免切断表格列）
  F4  \\[ ... \\] 之类残留检查

只处理围栏代码块之外的内容。
默认 dry-run；加 --apply 才写盘。
"""

import argparse
import glob
import os
import re
import sys

FENCE = re.compile(r"(```.*?```|~~~.*?~~~)", re.S)


def split_fences(text):
    """返回 [(is_code, chunk)]。"""
    out, last = [], 0
    for m in FENCE.finditer(text):
        out.append((False, text[last:m.start()]))
        out.append((True, m.group(0)))
        last = m.end()
    out.append((False, text[last:]))
    return out


def fix_inline(text):
    """\\(...\\) -> $...$；跨行时 -> $$ 块。"""
    def repl(m):
        body = m.group(1)
        if "\n" in body:
            return "$$\n" + body.strip() + "\n$$"
        return "$" + body + "$"
    return re.sub(r"\\\((.*?)\\\)", repl, text, flags=re.S)


def fix_bracket(text):
    """\\[...\\] -> $$ 块。"""
    return re.sub(r"\\\[(.*?)\\\]",
                  lambda m: "$$\n" + m.group(1).strip() + "\n$$",
                  text, flags=re.S)


def fix_table_pipes(text):
    """表格行内的行内数学含 | 时，把 | 换成 \\vert。"""
    out = []
    for line in text.split("\n"):
        if line.strip().startswith("|") and "$" in line:
            def repl(m):
                body = m.group(1)
                if "|" in body and "\\vert" not in body and "\\mid" not in body:
                    body = body.replace("|", "\\vert ")
                return "$" + body + "$"
            line = re.sub(r"\$([^$\n]+)\$", repl, line)
        out.append(line)
    return "\n".join(out)


def process(path, apply=False):
    raw = open(path, encoding="utf-8").read()
    n1 = len(re.findall(r"\\\(", raw))
    n2 = len(re.findall(r"\\\[", raw))
    n3 = len(re.findall(r"^\s*\|.*\$[^$\n]*\|[^$\n]*\$", raw, flags=re.M))

    chunks = split_fences(raw)
    new = []
    for is_code, chunk in chunks:
        if is_code:
            new.append(chunk)
        else:
            c = fix_inline(chunk)
            c = fix_bracket(c)
            c = fix_table_pipes(c)
            new.append(c)
    fixed = "".join(new)

    changed = fixed != raw
    if apply and changed:
        open(path, "w", encoding="utf-8").write(fixed)
    return changed, (n1, n2, n3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--root", default=os.path.expanduser("~/Downloads/lh"))
    ap.add_argument("--glob", default="*.md")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.root, args.glob)))
    tot = [0, 0, 0]
    nch = 0
    for f in files:
        changed, (a, b, c) = process(f, apply=args.apply)
        tot[0] += a
        tot[1] += b
        tot[2] += c
        if changed:
            nch += 1
    print("%s：扫描 %d 个文件" % ("已写盘" if args.apply else "dry-run", len(files)))
    print("  含 \\( 的文件行数计 %d 处；含 \\[ 的 %d 处；表格内含 | 的数学 %d 处"
          % (tot[0], tot[1], tot[2]))
    print("  将被修改（或已修改）的文件：%d" % nch)
    return 0


if __name__ == "__main__":
    sys.exit(main())
