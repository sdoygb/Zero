#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""latex_check.py —— 根目录 Z*.md 的 LaTeX 语法体检。

背景（2026-10-17）：用户报"很多 LaTeX 公式多了很多大于号"。定位到两处：
  ① **引用块里的行间公式**：`> $$ … > $$` 共 679 行／49 篇——`> ` 是 Markdown 引用标记，
     在把公式排进独立段落时会被当成字面字符（渲染里就多出 `>`）。**已全库改为公式移出引用块**。
  ② **`$` 定界符失衡**：Z118 一处 `…=2$$` 多一个 `$`（把后面的 `$…$` 吃掉，并伪造出 `$$`）。

本检查器把它们固化为回归护栏，并加几条常见 LaTeX 错型。
2026-10-18 加 ⑥：**双反斜杠**（`\\text` 这类"命令被复制成两份"）——渲染会变成字面反斜杠，
而 `$`／花括号检查都查不出来。矩阵类环境（pmatrix 等）里的 `\\` 是行分隔，按环境区间白名单跳过。
2026-10-18 加 ⑦：**控制字符**（`\a`／`\b`／`\v` 被转义吃掉后残留的 BEL／BS／VT）——实测在 4 个文件里抓到过 9 处。

用法：
  python3 _staging/latex_check.py              # 全库
  python3 _staging/latex_check.py Z008_*.md    # 指定篇
  python3 _staging/latex_check.py --summary
"""
import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent


def strip_code(t: str) -> str:
    """去掉行内 code 与围栏 code：它们里面的 $ / { / > 不参与判定。"""
    t = re.sub(r"```.*?```", "", t, flags=re.S)
    return re.sub(r"`[^`\n]*`", "", t)


def scan(path: pathlib.Path):
    raw = path.read_text(encoding="utf-8")
    t = strip_code(raw)
    hits = []

    # ① 行间公式被引用块包住（`> $$`）
    for i, line in enumerate(raw.split("\n"), 1):
        if re.match(r"^> *\$\$", line):
            hits.append((i, "行间公式在引用块内（`> $$`：渲染会多出 `>`）", line.strip()[:100]))

    # ② `$$` 必须成对（奇数 = 有未闭合的行间公式）
    n2 = t.count("$$")
    if n2 % 2:
        hits.append((0, f"`$$` 数为奇数（{n2}）——有未闭合的行间公式", ""))

    # ③ 行内 `$` 必须成对（排除 `$$`）
    n1 = len(re.findall(r"(?<!\$)\$(?!\$)", t))
    if n1 % 2:
        hits.append((0, f"行内 `$` 数为奇数（{n1}）——有未闭合的行内公式", ""))

    # ④ 行间公式块内花括号不平衡（逐块查）
    for m in re.finditer(r"\$\$(.*?)\$\$", t, flags=re.S):
        body = m.group(1)
        if body.count("{") != body.count("}"):
            ln = t[: m.start()].count("\n") + 1
            hits.append((ln, f"行间公式花括号不平衡（{{ {body.count('{')} / }} {body.count('}')}）",
                         body.strip().replace("\n", " ")[:100]))

    # ⑤ 环境配对：\begin{X} 必须在同一篇里有 \end{X}
    #   （原第 ⑤ 条"命令缺花括号"已删：`\mathbb Z`／`\mathcal R`／`\operatorname{span}` 这类
    #    **单 token 参数是合法 LaTeX**，加它会产出大量假阳性——实测 89 处命中里绝大多数是它。）
    for m in re.finditer(r"\\begin\{([a-z*]+)\}", t):
        env = m.group(1)
        if f"\\end{{{env}}}" not in t:
            ln = t[: m.start()].count("\n") + 1
            hits.append((ln, f"\\begin{{{env}}} 无配对 \\end{{{env}}}", ""))

    # ⑦ 控制字符（`\a`／`\b`／`\v` 被转义吃掉后的残留：BEL／BS／VT）
    ctrl = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
    for m in ctrl.finditer(t):
        ln = t[: m.start()].count("\n") + 1
        code = ord(m.group(0))
        hits.append((ln, f"控制字符 0x{code:02x}（`\\a`／`\\b`／`\\v` 被转义吃掉）",
                     t[max(0, m.start() - 24): m.start() + 16].replace("\n", " ")))

    # ⑥ 双反斜杠：`\\`+字母 基本是"命令被复制成两份"（渲染成字面的 `\text` 等）。
    #    **矩阵类环境里的 `\\` 是行分隔**（如 `1&V\\V&1`），按环境区间白名单跳过。
    MATRIX_ENVS = ("matrix", "pmatrix", "bmatrix", "Bmatrix", "vmatrix", "Vmatrix",
                   "smallmatrix", "array", "aligned", "alignedat", "gathered", "split",
                   "cases", "subarray")
    protected = []
    for env in MATRIX_ENVS:
        for m in re.finditer(rf"\\begin\{{{env}\}}.*?\\end\{{{env}\}}", t, flags=re.S):
            protected.append((m.start(), m.end()))

    def in_protected(pos):
        return any(a <= pos < b for a, b in protected)

    for m in re.finditer(r"\\\\([A-Za-z])", t):
        if in_protected(m.start()):
            continue
        ln = t[: m.start()].count("\n") + 1
        hits.append((ln, "双反斜杠（命令被复制成两份：`\\\\text` 等，渲染成字面反斜杠）",
                     t[m.start(): m.start() + 24].replace("\n", " ")))

    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*", help="默认：根目录全部 Z*.md")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()

    files = [pathlib.Path(f) for f in a.files] or sorted(ROOT.glob("Z[0-9][0-9][0-9]_*.md"))
    total, dirty = 0, 0
    for p in files:
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
                where = f"L{i:<5d}" if i else "      "
                print(f"  {where} [{w}]  {frag}")
    print(f"\n合计：{dirty} 篇 / {total} 处命中（检查 {len(files)} 篇）")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
