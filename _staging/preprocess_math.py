#!/usr/bin/env python3
"""Preprocess Chinese-science Markdown before pandoc->xelatex.

Fixes the pervasive shorthand  `x_\\mathcal{Z}`  (letter-arg math
commands used as bare subscripts without an enclosing brace group),
which breaks under unicode-math (pandoc's xelatex path loads it
unconditionally).  The fix wraps the argument in braces:

    x_\\mathcal{Z}      ->  x_{\\mathcal{Z}}
    x_\\mathcal Z       ->  x_{\\mathcal{Z}}
    x_{\\mathcal{Z}}    ->  unchanged (already safe)

Only math segments ($...$, $$...$$, \\(...\\), \\[...\\]) are touched,
so body text is never altered.  Idempotent.
"""
import re
import sys

ARG_CMDS = (
    "mathcal", "mathbb", "mathscr", "mathfrak", "mathds", "mathbf",
    "mathsf", "mathtt", "mathit", "mathrm", "bm", "boldsymbol",
)
_CMD_ALT = "|".join(ARG_CMDS)
# Match `_`/`^` then `\cmd` where next non-space char is `{` or a letter.
_PAT = re.compile(r"[_^](\\(" + _CMD_ALT + r"))(?=\s*[{A-Za-z])")


def _skip_ws(s, i):
    while i < len(s) and s[i] in " \t":
        i += 1
    return i


def fix_math(seg: str) -> str:
    out = []
    pos = 0
    while True:
        m = _PAT.search(seg, pos)
        if not m:
            out.append(seg[pos:])
            break
        pre = m.group(0)[0]            # '_' or '^'
        cmd = m.group(1)               # backslash + name
        # emit everything before, then "_" or "^"
        out.append(seg[pos:m.start()])
        i = m.end()
        i = _skip_ws(seg, i)
        if i >= len(seg):
            out.append(seg[m.start():])
            break
        if seg[i] == "{":
            depth = 0
            k = i
            while k < len(seg):
                if seg[k] == "{":
                    depth += 1
                elif seg[k] == "}":
                    depth -= 1
                    if depth == 0:
                        break
                k += 1
            arg = seg[i:k + 1]
            out.append(pre + "{" + cmd + arg + "}")
            pos = k + 1
        elif seg[i].isalpha():
            out.append(pre + "{" + cmd + "{" + seg[i] + "}}")
            pos = i + 1
        else:
            out.append(seg[m.start():m.end()])
            pos = m.end()
    return "".join(out)

# known command-glue typos in the corpus (missing space after command)
_GLUE_FIXES = [
    (r"\\Lambda k", r"\\Lambda k"),
    (r"\\Lambdak", r"\\Lambda k"),
    (r"\\Thetak", r"\\Theta k"),
    (r"\\xiA", r"\\xi A"),
    (r"\\xiB", r"\\xi B"),
    (r"\\xiC", r"\\xi C"),
    (r"\\etaA", r"\\eta A"),
    (r"\\etaB", r"\\eta B"),
    (r"\\etaC", r"\\eta C"),
    (r"\\etaF", r"\\eta F"),
    (r"\\lambdaA", r"\\lambda A"),
    (r"\\lambdaB", r"\\lambda B"),
    (r"\\lambdaC", r"\\lambda C"),
    (r"\\muA", r"\\mu A"),
    (r"\\muB", r"\\mu B"),
    (r"\\muC", r"\\mu C"),
    (r"\\alphaA", r"\\alpha A"),
    (r"\\betaA", r"\\beta A"),
    (r"\\gammaA", r"\\gamma A"),
    (r"\\gammaB", r"\\gamma B"),
    (r"\\gammaC", r"\\gamma C"),
    (r"\\thetaA", r"\\theta A"),
    (r"\\thetaB", r"\\theta B"),
    (r"\\thetaC", r"\\theta C"),
    (r"\\phiA", r"\\phi A"),
    (r"\\psiA", r"\\psi A"),
    (r"\\deltaA", r"\\delta A"),
    (r"\\deltaB", r"\\delta B"),
    (r"\\per\^", r"\\perp^"),
    (r"\\per\{", r"\\perp{"),
    (r"\\per[^p]", r"\\perp "),
    # 数学段内的 \* 是 LaTeX 的「可选断字」命令，出现在上标位（如 w_2^\*）会报
    # Missing { inserted ... \discretionary。作者本意是上标星号 ⇒ 一律换成 \ast。
    (r"\\\*", r"\\ast"),
]

def _esc_text_block(m):
    """转义 \\text{...} / \\textbf{...} 等文本命令内的裸 _ 与 %。

    旧实现用 (\\text\\{)([^{}]*?)_([^{}]*?)(\\}) 一次替换，惰性组会一路吞到
    右花括号，导致块内**只有第一个 _ 被转义**，第二个起漏掉
    （例：\\text{—— \\`build_169\\` 的 \\`Z_Λ\\` 律} ⇒ Z_Λ 未转义 ⇒ Missing $ inserted）。
    另外 \\textbf{} / \\emph{} 等也是文本模式，裸 _ 同样非法
    （例：\\textbf{\\`l_P\\`（等价 \\`G\\`）} ⇒ Missing $ inserted）。
    这里统一整块处理：块内所有裸 _ / % 全部转义。
    """
    cmd, inner = m.group(1), m.group(2)
    inner = re.sub(r"(?<!\\)_", r"\\_", inner)
    inner = re.sub(r"(?<!\\)%", r"\\%", inner)
    return "\\" + cmd + "{" + inner + "}"


# 文本模式的命令（裸 _ / % 必须转义）——涵盖 pandoc 从 markdown 加粗/斜体生成的命令
_TEXT_CMDS = r"text|textbf|textit|textrm|textsf|texttt|textsc|textup|emph|underline"


def fix_glue(seg: str) -> str:
    for pat, rep in _GLUE_FIXES:
        seg = re.sub(pat, rep, seg)
    # 转义文本类命令 {} 内所有裸 _ / %（否则 xelatex 报 Missing $ inserted）
    seg = re.sub(r"\\(%s)\{([^{}]*)\}" % _TEXT_CMDS, _esc_text_block, seg)
    return seg


def escape_underscores_in_body(text: str) -> str:
    """Escape bare underscores in non-math body text ('T_c' in prose).

    Math segments ($...$, $$...$$, \\(...\\), \\[...\\]) are masked first so
    underscores inside them are left for the math fixer to handle.
    """
    import re as _re
    segs = []
    def stash(m):
        segs.append(m.group(0))
        return "\x00SEG%d\x00" % (len(segs) - 1)
    masked = _re.sub(
        r"\$\$.*?\$\$|\$[^$\n]{1,500}?\$|\\\(.*?\\\)|\\\[.*?\\\]",
        stash, text, flags=_re.S)
    # escape '_' in remaining body text (not already escaped)
    masked = _re.sub(r"(?<!\\)_", r"\\_", masked)
    for i, seg in enumerate(segs):
        masked = masked.replace("\x00SEG%d\x00" % i, seg)
    return masked

def process(text: str) -> str:
    # ---- pass 1: stash every math segment
    segs = []
    def stash_math(m):
        segs.append(m.group(0))
        return "\x00M%d\x00" % (len(segs) - 1)

    masked = re.sub(
        r"\$\$.*?\$\$|\$[^$\n]{1,500}?\$|\\\(.*?\\\)|\\\[.*?\\\]",
        stash_math, text, flags=re.S)

    # ---- pass 2: escape remaining bare underscores in body text only
    masked = re.sub(r"(?<!\\)_", r"\\_", masked)

    # ---- pass 3: restore each math segment, with math-level fixes applied
    for i, seg in enumerate(segs):
        fixed = fix_glue(fix_math(seg))
        masked = masked.replace("\x00M%d\x00" % i, fixed)
    return masked

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: preprocess_math.py <in.md> <out.md>")
        sys.exit(1)
    with open(sys.argv[1], encoding="utf-8") as f:
        txt = f.read()
    with open(sys.argv[2], "w", encoding="utf-8") as f:
        f.write(process(txt))
    print(f"processed {sys.argv[1]} -> {sys.argv[2]}")
