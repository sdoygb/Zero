#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""audit_latex.py —— 对任意 markdown 目录做 GitHub 渲染**穷举体检**

用法：
    python3 audit_latex.py <目录> [文件名通配]

检查 9 类问题（与 Zero-CSS 的 verify/check_latex.py 同一套规则）：
  E1  $$ 块边界（开块前 / 闭块后需空行）
  E2  表格行内数学含裸 |
  E3  已知被 GitHub 拒绝的宏（\\boxed \\operatorname \\llbracket ...）
  E4  控制字符（TAB / 退格）
  E5  表格列数不一致
  E6  行内公式含裸 _（会被 Markdown 当强调吃掉）
  E7  行内公式含 [[ ]] \\[ \\]
  E8  行内公式用了最小安全集之外的宏
  E9  孤立数学行（含 LaTeX 但不被 $...$ 或 $$ 包围）
"""
import os
import re
import sys
import glob

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from safe_macros import SAFE, BANNED  # noqa: E402


def e1(lines):
    out, inb = [], False
    for i, l in enumerate(lines):
        if l.strip() != "$$":
            continue
        prev = lines[i - 1].strip() if i > 0 else ""
        nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
        if not inb:
            if prev and not prev.startswith("|") and not prev.startswith("#"):
                out.append((i + 1, "开块 $$ 前缺空行: " + prev[:45]))
            inb = True
        else:
            if nxt and not nxt.startswith("|"):
                out.append((i + 1, "闭块 $$ 后缺空行: " + nxt[:45]))
            inb = False
    if inb:
        out.append((len(lines), "存在未闭合的 $$ 块"))
    return out


def _inline_segs(l):
    """返回行内数学片段（排除 $$ 覆盖区间）。"""
    dd = [(m.start(), m.end()) for m in re.finditer(r"\$\$.+?\$\$", l)]
    res = []
    for m in re.finditer(r"\$`([^`\n]+)`\$", l):
        res.append((m.start(), m.group(1)))
    for m in re.finditer(r"(?<!\$)\$([^$\n]+)\$(?!\$)", l):
        if any(a <= m.start() < b for a, b in dd):
            continue
        res.append((m.start(), m.group(1)))
    return res


def e2(lines):
    out = []
    for i, l in enumerate(lines, 1):
        if not l.lstrip().startswith("|") or l.lstrip().startswith("|--"):
            continue
        for _, seg in _inline_segs(l):
            if "|" in seg:
                out.append((i, seg.strip()[:60]))
    return out


def e3(lines):
    out = []
    for i, l in enumerate(lines, 1):
        for mac, why in BANNED.items():
            if re.search(r"\\" + re.escape(mac) + r"(?![a-zA-Z])", l):
                out.append((i, "\\%s (%s)" % (mac, why)))
    return out


def e4(text):
    out = []
    for m in re.finditer(r"[\t\b]", text):
        ln = text[: m.start()].count("\n") + 1
        out.append((ln, "TAB" if m.group() == "\t" else "BACKSPACE"))
    return out


def e5(lines):
    out, tbl = [], []

    def flush():
        if len(tbl) > 1:
            from collections import Counter
            mode = Counter(c for _, c in tbl).most_common(1)[0][0]
            for ln, c in tbl:
                if c != mode:
                    out.append((ln, "| 数 %d（同表众数 %d）" % (c, mode)))

    for i, l in enumerate(lines, 1):
        if l.lstrip().startswith("|"):
            tbl.append((i, l.count("|")))
        else:
            flush()
            tbl = []
    flush()
    return out


def e6(lines):
    out = []
    for i, l in enumerate(lines, 1):
        for _, seg in _inline_segs(l):
            if re.search(r"(?<!\\)_", seg):
                out.append((i, seg.strip()[:60]))
    return out


def e7(lines):
    out = []
    for i, l in enumerate(lines, 1):
        for _, seg in _inline_segs(l):
            for pat in ("[[", "]]", r"\[", r"\]"):
                if pat in seg:
                    out.append((i, "%s 出现在行内公式: %s" % (pat, seg.strip()[:50])))
                    break
    return out


def e8(lines):
    out = []
    for i, l in enumerate(lines, 1):
        for _, seg in _inline_segs(l):
            for mm in re.finditer(r"\\([a-zA-Z]+)", seg):
                if mm.group(1) not in SAFE:
                    out.append((i, "\\%s 不在最小安全集: %s" % (mm.group(1), seg[:45])))
                    break
    return out


def e9(lines):
    out, inb, inc = [], False, False
    for i, l in enumerate(lines, 1):
        if l.strip().startswith("```"):
            inc = not inc
            continue
        if inc:
            continue
        if l.strip() == "$$":
            inb = not inb
            continue
        if inb or not l.strip():
            continue
        if l.lstrip().startswith(("|", "#", ">")):
            continue
        t = re.sub(r"\$`[^`\n]+`\$", "", l)
        t = re.sub(r"(?<!\$)\$[^$\n]+\$(?!\$)", "", t)
        t = re.sub(r"`[^`]*`", "", t)
        if re.search(r"\\[a-zA-Z]{2,}", t):
            out.append((i, l.strip()[:65]))
    return out


NAMES = {1: "$$ 块边界", 2: "表格内裸 |", 3: "被拒宏", 4: "控制字符",
         5: "表格列数", 6: "行内裸 _", 7: "行内方括号", 8: "不安全宏", 9: "孤立数学行"}


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    pat = sys.argv[2] if len(sys.argv) > 2 else "**/*.md"
    files = sorted(glob.glob(os.path.join(root, pat), recursive=True))
    files = [f for f in files if ".git/" not in f and "node_modules" not in f]
    print("=" * 84)
    print("GitHub 渲染穷举体检 —— 目录 %s" % root)
    print("文件数 %d" % len(files))
    print("=" * 84)
    grand = {k: 0 for k in NAMES}
    dirty = []
    for f in files:
        text = open(f, encoding="utf-8", errors="replace").read()
        lines = text.split("\n")
        res = {1: e1(lines), 2: e2(lines), 3: e3(lines), 4: e4(text),
               5: e5(lines), 6: e6(lines), 7: e7(lines), 8: e8(lines), 9: e9(lines)}
        n = sum(len(v) for v in res.values())
        if n:
            dirty.append((f, res, n))
            for k in NAMES:
                grand[k] += len(res[k])
    for f, res, n in sorted(dirty, key=lambda x: -x[2]):
        print("\n%s  (%d 处)" % (os.path.relpath(f, root), n))
        for k in sorted(NAMES):
            for ln, txt in res[k][:4]:
                print("   E%d L%-5d %s | %s" % (k, ln, NAMES[k], txt))
            if len(res[k]) > 4:
                print("   E%d ... 另有 %d 处" % (k, len(res[k]) - 4))
    print("\n" + "=" * 84)
    print("按类别统计：")
    for k in sorted(NAMES):
        print("   E%d %-14s %5d 处" % (k, NAMES[k], grand[k]))
    print("   合计 %d 处，涉及 %d / %d 个文件" % (sum(grand.values()), len(dirty), len(files)))
    return 1 if sum(grand.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
