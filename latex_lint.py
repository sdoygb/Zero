#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
latex_lint.py -- 扫描项目内所有 Markdown 文章的 LaTeX 语法问题。

分类：
  C1  $$ 未闭合（行间分隔符成对失败）
  C2  行内 $ 计数为奇（未配对）
  C3  数学内 { } 不平衡
  C4  \\left 与 \\right 数量不等
  C5  \\begin{...} 与 \\end{...} 不匹配
  C6  表格行内的数学含 | （会切断表格列）
  C7  使用了 \\(...\\) 或 \\[...\\]（与项目 $ 风格不一致，渲染器可能不支持）
  C8  $$ 块内含空行（部分渲染器会提前结束）
  C9  $$ 与正文同行（部分渲染器只认独立行的 $$）
  C10 \\text{}/\\mbox{} 内含中日韩字符（KaTeX 常显示为缺字）
  C11 可疑命令（黑名单）
  C12 嵌套顺序错误（} 关闭了未结束的环境）
  C13 数学内含 markdown 结构（链接/反引号/强调）——渲染器会先解析成 HTML，
      KaTeX 随即失败并【回退显示原始 LaTeX】

用法：
  python3 latex_lint.py               # 扫描默认根，打印摘要
  python3 latex_lint.py --detail      # 追加打印每条问题的位置
"""

import os
import re
import sys
from collections import Counter, defaultdict

CJK = re.compile(r"[\u3000-\u303f\u4e00-\u9fff\uff00-\uffef]")

SUSPECT_CMDS = [
    r"\begin{align}", r"\begin{equation}", r"\begin{gather}",
    r"\begin{eqnarray}", r"\begin{multline}", r"\begin{alignat}",
    r"\eqref", r"\label{", r"\tag{", r"\hbox", r"\vbox",
    r"\def", r"\newcommand", r"\renewcommand", r"\DeclareMathOperator",
    r"\intertext", r"\subequations",
]


def strip_code(text):
    """去掉围栏代码块与行内代码，但保留换行与字符偏移（行号才不会错位）。"""
    def blank(m):
        return re.sub(r"[^\n]", " ", m.group(0))
    text = re.sub(r"```.*?```", blank, text, flags=re.S)
    text = re.sub(r"~~~.*?~~~", blank, text, flags=re.S)
    text = re.sub(r"`[^`\n]*`", blank, text)
    return text


def line_of(text, pos):
    return text.count("\n", 0, pos) + 1


def brace_balance(s):
    depth = 0
    i = 0
    while i < len(s):
        c = s[i]
        if c == "\\":
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth < 0:
                return False
        i += 1
    return depth == 0


def nest_ok(s):
    """栈式检查：} 不能关闭尚未 \end 的环境（花括号平衡≠顺序正确）。"""
    st = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == "\\":
            m = re.match(r"\\begin\{([^}]*)\}", s[i:])
            if m:
                st.append(("env", m.group(1)))
                i += m.end()
                continue
            m = re.match(r"\\end\{([^}]*)\}", s[i:])
            if m:
                nm = m.group(1)
                if not st or st[-1][0] != "env" or st[-1][1] != nm:
                    return False
                st.pop()
                i += m.end()
                continue
            i += 2
            continue
        if c == "{":
            st.append(("brace", None))
        elif c == "}":
            if not st or st[-1][0] != "brace":
                return False
            st.pop()
        i += 1
    return all(k == "brace" for k, _ in st)


def begin_end_ok(s):
    starts = re.findall(r"\\begin\{([^}]*)\}", s)
    ends = re.findall(r"\\end\{([^}]*)\}", s)
    if len(starts) != len(ends):
        return False
    for a, b in zip(starts, ends):
        if a != b:
            return False
    return True


def collect_math(text):
    """返回 [(kind, content, pos)]，kind ∈ {display, inline, paren, bracket}。"""
    out = []
    # 行间 $$...$$
    for m in re.finditer(r"\$\$(.*?)\$\$", text, flags=re.S):
        out.append(("display", m.group(1), m.start()))
    # 去掉 $$ 区间后再找行内
    no_disp = re.sub(r"\$\$.*?\$\$", lambda m: " " * len(m.group(0)), text, flags=re.S)
    for m in re.finditer(r"\$([^$\n]+)\$", no_disp):
        out.append(("inline", m.group(1), m.start()))
    for m in re.finditer(r"\\\((.*?)\\\)", text, flags=re.S):
        out.append(("paren", m.group(1), m.start()))
    for m in re.finditer(r"\\\[(.*?)\\\]", text, flags=re.S):
        out.append(("bracket", m.group(1), m.start()))
    return out


def lint_file(path):
    raw = open(path, encoding="utf-8", errors="replace").read()
    return lint_text(raw)


def lint_text(raw):
    """对【文本】做体检（与 lint_file 同一套规则，供正控制复用）。"""
    text = strip_code(raw)
    issues = defaultdict(list)

    # C1  $$ 闭合
    if text.count("$$") % 2 != 0:
        issues["C1"].append("$$ 计数为 %d（未闭合）" % text.count("$$"))

    # C2  行内 $
    no_disp = re.sub(r"\$\$.*?\$\$", lambda m: " " * len(m.group(0)), text, flags=re.S)
    if no_disp.count("$") % 2 != 0:
        issues["C2"].append("行内 $ 计数为 %d（未配对）" % no_disp.count("$"))

    spans = collect_math(text)

    for kind, content, pos in spans:
        ln = line_of(text, pos)
        # C3 花括号
        if not brace_balance(content):
            issues["C3"].append("L%d %s：花括号不平衡" % (ln, kind))
        # C4 left/right
        nl = len(re.findall(r"\\left(?![a-zA-Z])", content))
        nr = len(re.findall(r"\\right(?![a-zA-Z])", content))
        if nl != nr:
            issues["C4"].append("L%d %s：\\left=%d \\right=%d" % (ln, kind, nl, nr))
        # C5 begin/end
        if not begin_end_ok(content):
            issues["C5"].append("L%d %s：begin/end 不匹配" % (ln, kind))
        if not nest_ok(content):
            issues["C12"].append("L%d %s：嵌套顺序错误（} 关闭了未结束的环境）" % (ln, kind))
        # C10 CJK in \text
        for m in re.finditer(r"\\(?:text|mbox|mathrm)\{([^{}]*)\}", content):
            if CJK.search(m.group(1)):
                issues["C10"].append("L%d %s：\\text 内含中文" % (ln, kind))
                break
        # C13 a) 数学内含 markdown 结构（链接/反引号/强调）
        m13 = re.search(r"\[[^\]]*\]\([^)]*\)|`|\*\*", content)
        if m13:
            issues["C13"].append("L%d %s：数学内含 markdown（%s）"
                                 % (ln, kind, m13.group(0)[:24]))
        # C13 b) \text{}/\mbox{}/\mathrm{} 内含 _ ^ 或反引号（KaTeX 硬报错）
        for m in re.finditer(r"\\(?:text|mbox|mathrm)\{([^{}]*)\}", content):
            bad = re.search(r"(?<!\\)[_^`]", m.group(1))   # 只抓【未转义】的（\_ 是合法的）
            if bad:
                issues["C13"].append("L%d %s：\\text 内含 '%s'（KaTeX 硬报错）"
                                     % (ln, kind, bad.group(0)))

        # C11 可疑命令
        for cmd in SUSPECT_CMDS:
            if cmd in content:
                issues["C11"].append("L%d %s：%s" % (ln, kind, cmd))
                break

        # C6 表格行内数学含 |
        if kind in ("inline", "paren") and "|" in content:
            line = raw.split("\n")[ln - 1] if ln - 1 < len(raw.split("\n")) else ""
            if line.strip().startswith("|") and line.count("$") == 2:
                issues["C6"].append("L%d：表格行内数学含 |" % ln)

    # C7 \(...\) \[...\]
    np_ = len(re.findall(r"(?<!\\)\\\(", text))
    nb_ = len(re.findall(r"(?<!\\)\\\[", text))
    if np_ or nb_:
        issues["C7"].append("\\( 出现 %d 次，\\[ 出现 %d 次" % (np_, nb_))

    # C8  $$ 块内空行
    for m in re.finditer(r"\$\$(.*?)\$\$", text, flags=re.S):
        if "\n\n" in m.group(1):
            issues["C8"].append("L%d：$$ 块内含空行" % line_of(text, m.start()))

    # C9  $$ 与正文同行
    for m in re.finditer(r"\$\$(.*?)\$\$", text, flags=re.S):
        ln = line_of(text, m.start())
        line = raw.split("\n")[ln - 1] if ln - 1 < len(raw.split("\n")) else ""
        if line.strip() not in ("$$",) and not line.strip().startswith("$$"):
            issues["C9"].append("L%d：$$ 与正文同行" % ln)

    return issues


ROOTS = [
    ("lh", os.path.expanduser("~/Downloads/lh")),
    ("modular-equilibrium", os.path.expanduser("~/Downloads/modular-equilibrium")),
]


def main():
    detail = "--detail" in sys.argv
    files = []
    for tag, root in ROOTS:
        if not os.path.isdir(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in ("__pycache__", ".git")]
            for fn in filenames:
                if fn.endswith(".md"):
                    files.append(os.path.join(dirpath, fn))
    files.sort()

    total = Counter()
    per_file = {}
    for f in files:
        iss = lint_file(f)
        if iss:
            per_file[f] = iss
            for k, v in iss.items():
                total[k] += len(v)

    print("扫描 Markdown 文件：%d" % len(files))
    print("有问题的文件：%d" % len(per_file))
    print()
    names = {
        "C1": "$$ 未闭合", "C2": "行内 $ 未配对", "C3": "花括号不平衡",
        "C4": "\\left/\\right 不等", "C5": "begin/end 不匹配",
        "C6": "表格内数学含 |", "C7": "用了 \\( \\) / \\[ \\]",
        "C8": "$$ 内含空行", "C9": "$$ 与正文同行",
        "C10": "\\text 内含中文", "C11": "可疑命令",
        "C12": "嵌套顺序错误", "C13": "数学内含 markdown/硬错字符",
    }
    for k in sorted(names):
        print("  %-4s %-22s %5d 处，涉及 %3d 个文件"
              % (k, names[k], total.get(k, 0),
                 sum(1 for v in per_file.values() if k in v)))

    if detail:
        print("\n" + "=" * 72)
        for f, iss in sorted(per_file.items()):
            print("\n%s" % f)
            for k in sorted(iss):
                for item in iss[k][:4]:
                    print("    [%s] %s" % (k, item))
    else:
        print("\n按文件排序（问题数最多的前 20）：")
        ranked = sorted(per_file.items(),
                        key=lambda kv: -sum(len(v) for v in kv[1].values()))
        for f, iss in ranked[:20]:
            n = sum(len(v) for v in iss.values())
            cats = ",".join(sorted(iss))
            print("  %3d  %-70s %s" % (n, os.path.relpath(f, os.path.expanduser("~")), cats))

    return 0


if __name__ == "__main__":
    sys.exit(main())
