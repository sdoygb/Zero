#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""fix_latex.py —— 对 markdown 目录做 GitHub 渲染修复（分阶段，可逐步验证）

阶段：
  1  控制字符：TAB -> \t 的反斜杠，退格 -> \b 的反斜杠（修复被吃掉的 \t \b \tau \beta）
  2  \boxed{...} -> 去壳（GitHub 不加载 bbox 扩展）
  3  \operatorname(*) -> \text{}；\mathrm -> \text{}
  4  \xrightarrow{X} -> \overset{X}{\longrightarrow}
  5  行内公式：裸 _ -> \_ ；[[ ]] -> [ ] ；\[ \] 去掉转义
  6  $$ 块边界补空行
  7  孤立数学行 -> 补 $$ 定界符（仅纯数学行；含 markdown 粗体的打印出来人工处理）

用法：python3 fix_latex.py <目录> <阶段号> [--apply]
不带 --apply 时只报告，不写文件。
"""
import os
import re
import sys
import glob

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def files_of(root, pat="**/*.md"):
    fs = sorted(glob.glob(os.path.join(root, pat), recursive=True))
    return [f for f in fs if ".git/" not in f and "node_modules" not in f]


def s1_control(text):
    """TAB -> 反斜杠 + 't' 的位置无法自动判定，故只做**已知模式**的安全替换。

    实际损坏模式是 `$\tau_i$` 被写成 `$<TAB>au_i$`，即 \\t 变 TAB。
    同理 \\b -> 退格（\\beta -> <BS>eta）。
    修法：把 TAB 还原为 `\\t`，退格还原为 `\\b`。这是**唯一合理**的还原，
    因为 markdown 正文里不应出现裸 TAB（缩进用空格）/退格。
    """
    n = text.count("\t") + text.count("\b")
    text = text.replace("\t", "\\t").replace("\b", "\\b")
    return text, n


def s2_boxed(text):
    """去掉 \\boxed{...} 外壳（配对花括号）。"""
    out, i, n = [], 0, 0
    while True:
        k = text.find("\\boxed", i)
        if k < 0:
            out.append(text[i:])
            break
        out.append(text[i:k])
        j = k + len("\\boxed")
        while j < len(text) and text[j] in " \t":
            j += 1
        if j >= len(text) or text[j] != "{":
            out.append(text[k:j])
            i = j
            continue
        depth, m = 0, j
        while m < len(text):
            if text[m] == "{":
                depth += 1
            elif text[m] == "}":
                depth -= 1
                if depth == 0:
                    break
            m += 1
        out.append(text[j + 1 : m])
        n += 1
        i = m + 1
    return "".join(out), n


def s3_operatorname(text):
    """\\operatorname(*) / \\mathrm -> \\text（保留正体且不吃相邻空格）。"""
    n = 0
    for pat in (r"\\operatorname\*?\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}",
                r"\\mathrm\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}"):
        text, k = re.subn(pat, r"\\text{\1}", text)
        n += k
    return text, n


def s4_xrightarrow(text):
    text, n = re.subn(r"\\xrightarrow\s*\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}",
                      r"\\overset{\1}{\\longrightarrow}", text)
    return text, n


def _protect_blocks(line):
    blocks = []

    def stash(m):
        blocks.append(m.group(0))
        return "\x00%d\x00" % (len(blocks) - 1)

    return re.sub(r"\$\$.+?\$\$", stash, line), blocks


def _restore(line, blocks):
    return re.sub(r"\x00(\d+)\x00", lambda m: blocks[int(m.group(1))], line)


def s5_inline(text):
    """行内公式：裸 _ -> \\_ ；[[ ]] -> [ ] ；\\[ \\] -> [ ]。"""
    lines = text.split("\n")
    out = []
    n = 0
    for line in lines:
        work, blocks = _protect_blocks(line)

        def conv(m):
            nonlocal n
            seg = m.group(1)
            new = re.sub(r"(?<!\\)_", r"\\_", seg)
            new = new.replace("[[", "[").replace("]]", "]")
            new = new.replace("\\[", "[").replace("\\]", "]")
            if new != seg:
                n += 1
            return "$" + new + "$"

        work = re.sub(r"(?<!\$)\$([^$\n]+)\$(?!\$)", conv, work)
        out.append(_restore(work, blocks))
    return "\n".join(out), n


def s6_blank_lines(text):
    """$$ 块前后补空行。"""
    lines = text.split("\n")
    out, n, inb = [], 0, False
    for i, l in enumerate(lines):
        if l.strip() == "$$":
            if not inb:
                if out and out[-1].strip() != "":
                    out.append("")
                    n += 1
                out.append(l)
                inb = True
            else:
                out.append(l)
                inb = False
                nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
                if nxt:
                    out.append("")
                    n += 1
            continue
        out.append(l)
    return "\n".join(out), n


def s7_orphans(text):
    """孤立数学行 -> 补 $$。只处理**纯数学**行；
    含 markdown 标记（** 或 | 或 #）的返回待人工列表。"""
    lines = text.split("\n")
    wrapped, manual = [], []
    inb = inc = False
    for i, l in enumerate(lines):
        if l.strip().startswith("```"):
            inc = not inc
            wrapped.append(l)
            continue
        if inc:
            wrapped.append(l)
            continue
        if l.strip() == "$$":
            inb = not inb
            wrapped.append(l)
            continue
        if inb or not l.strip() or l.lstrip().startswith(("|", "#", ">")):
            wrapped.append(l)
            continue
        t = re.sub(r"\$`[^`\n]+`\$", "", l)
        t = re.sub(r"(?<!\$)\$[^$\n]+\$(?!\$)", "", t)
        t = re.sub(r"`[^`]*`", "", t)
        if re.search(r"\\[a-zA-Z]{2,}", t):
            if "**" in l or "|" in l:
                manual.append(i + 1)
                wrapped.append(l)
            else:
                wrapped.append("$$")
                wrapped.append(l)
                wrapped.append("$$")
    return "\n".join(wrapped), manual


STAGES = {1: ("控制字符", s1_control, True), 2: ("\\boxed", s2_boxed, True),
          3: ("\\operatorname/\\mathrm", s3_operatorname, True),
          4: ("\\xrightarrow", s4_xrightarrow, True),
          5: ("行内 _ [[ ]]", s5_inline, True),
          6: ("$$ 块空行", s6_blank_lines, True),
          7: ("孤立数学行", s7_orphans, False)}


def main():
    root = sys.argv[1]
    stage = int(sys.argv[2])
    apply_ = "--apply" in sys.argv
    name, fn, returns_int = STAGES[stage]
    print("阶段 %d：%s   目录 %s   %s" % (stage, name, root, "APPLY" if apply_ else "DRY-RUN"))
    total = 0
    manual_all = []
    for f in files_of(root):
        txt = open(f, encoding="utf-8", errors="replace").read()
        res = fn(txt)
        if returns_int:
            new, n = res
        else:
            new, manual = res
            n = len(manual)
            if manual:
                manual_all.append((f, manual))
        if n and apply_ and new != txt:
            open(f, "w", encoding="utf-8").write(new)
        if n:
            print("   %-58s %4d" % (os.path.relpath(f, root), n))
            total += n
    print("合计 %d 处" % total)
    if manual_all:
        print("\n需人工处理的孤立块（含 markdown 标记）：")
        for f, ln in manual_all:
            print("   %s  行 %s" % (os.path.relpath(f, root), ln))


if __name__ == "__main__":
    main()
