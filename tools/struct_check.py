#!/usr/bin/env python3
"""结构检查：编号与引用的一致性（全库回归）。

七项：
  C1 条目重号     同篇内 (定理|命题|推论|定义|引理|界定|设定|记账规则) X.Y 重复
  C2 条目在节内   条目号的前缀必须等于它所在的 `## §N` 的 N（位置化编号）
  C3 条目顺序     同一节内条目号递增（撇号算 +0.5，故 4.2′ 排在 4.2 之后）
  C4 公式号重复   \\tag{X.Y} 同篇重复
  C5 公式号顺序   同一节内 \\tag 的编号递增（本条之前常被漏检）
  C6 悬空条目引用 正文引用了本篇不存在的条目号
  C7 悬空公式引用 正文引用了本篇不存在的公式号
  C8 跨篇悬空     ZbNNN 引用的篇不存在

用法：
    python3 tools/struct_check.py            # 全库
    python3 tools/struct_check.py Zb047      # 单篇
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ITEM_RE = re.compile(r"\*\*(定理|命题|推论|定义|引理|界定|设定|记账规则) ?(\d+)\.(\d+)(′?)")
TAG_RE = re.compile(r"\\tag\{(\d+)\.(\d+)([a-z]?)\}")
SEC_RE = re.compile(r"^## §(\d+(?:\.\d+)?)")


def notes(pattern=None):
    fs = sorted(ROOT.glob("Zb[0-9][0-9][0-9]_*.md"))
    if pattern:
        fs = [f for f in fs if f.name.startswith(pattern)]
    return fs


def check(path, have):
    s = path.read_text(encoding="utf-8")
    lines = s.split("\n")
    issues = []

    # ---- 收集条目与公式，带所在节
    items = []   # (行号, 类型, 大, 小, 撇, 节)
    tags = []    # (行号, 大, 小, 后缀, 节)
    sec = None
    for i, ln in enumerate(lines):
        m = SEC_RE.match(ln)
        if m:
            sec = m.group(1)
            continue
        m = ITEM_RE.search(ln)
        if m and not ln.lstrip().startswith("|"):
            items.append((i, m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), sec))
        m = TAG_RE.search(ln)
        if m:
            tags.append((i, int(m.group(1)), int(m.group(2)), m.group(3), sec))

    names = {(t, a, b, c) for _, t, a, b, c, _ in items}
    ids = {f"{a}.{b}{c}" for _, _, a, b, c, _ in items}
    tagids = {f"{a}.{b}{c}" for _, a, b, c, _ in tags}

    # ---- C1 条目重号
    seen = {}
    for i, t, a, b, c, _ in items:
        k = (t, a, b, c)
        if k in seen:
            issues.append((i + 1, "C1", f"{t} {a}.{b}{c} 重号（首见 L{seen[k]+1}）"))
        else:
            seen[k] = i

    # ---- C2 条目在节内
    for i, t, a, b, c, sc in items:
        if sc is None:
            continue
        head = sc.split(".")[0]
        if a != int(head):
            issues.append((i + 1, "C2", f"{t} {a}.{b}{c} 出现在 §{sc}"))

    # ---- C3 条目顺序（节内递增，撇号 +0.5）
    prev = {}
    for i, t, a, b, c, sc in items:
        if sc is None:
            continue
        key = int(b) * 2 + (1 if c else 0)
        if sc in prev and key <= prev[sc][0]:
            issues.append((i + 1, "C3", f"{t} {a}.{b}{c} 排在 {prev[sc][1]} 之后（乱序）"))
        prev[sc] = (key, f"{a}.{b}{c}")

    # ---- C4 公式号重复
    seen = {}
    for i, a, b, c, _ in tags:
        k = f"{a}.{b}{c}"
        if k in seen:
            issues.append((i + 1, "C4", f"\\tag{{{k}}} 重号（首见 L{seen[k]+1}）"))
        else:
            seen[k] = i

    # ---- C5 公式号顺序（节内递增）——本轮补的检查
    prev = {}
    for i, a, b, c, sc in tags:
        if sc is None:
            continue
        key = (a, b)
        if sc in prev and key < prev[sc][0]:
            issues.append((i + 1, "C5", f"\\tag{{{a}.{b}{c}}} 排在 ({prev[sc][1]}) 之后（乱序）"))
        prev[sc] = (key, f"{a}.{b}{c}" if not c else f"{a}.{b}{c}")

    # ---- C6 悬空条目引用
    for i, ln in enumerate(lines):
        for m in re.finditer(r"(定理|命题|推论|定义|引理)\s*(\d+)\.(\d+)(′?)", ln):
            pre = ln[max(0, m.start() - 40):m.start()]
            post = ln[m.end():m.end() + 12]
            if re.search(r"Zb\d{3}[^。]{0,20}$", pre) or re.match(r"\s*[（(]\s*Zb\d{3}", post):
                continue          # 跨篇引用（前者「Zb015 命题 4.2」；后者「推论 2.2（Zb044）」）由 C8 管
            if (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)) not in names:
                issues.append((i + 1, "C6", f"引用了不存在的 {m.group(0)}"))

    # ---- C7 悬空公式引用
    for i, ln in enumerate(lines):
        for m in re.finditer(r"\((\d+)\.(\d+)([a-z]?)\)", ln):
            k = f"{m.group(1)}.{m.group(2)}{m.group(3)}"
            if k in tagids:
                continue
            if f"{m.group(1)}.{m.group(2)}" in ids or m.group(2).isdigit() and len(m.group(2)) > 2:
                continue          # 条目号或数值（如 1836.153），不是公式号
            issues.append((i + 1, "C7", f"引用了不存在的公式 ({k})"))

    # ---- C8 跨篇悬空
    for i, ln in enumerate(lines):
        for m in re.finditer(r"Zb(\d{3})", ln):
            if m.group(0) not in have and path.name[:5] != m.group(0):
                issues.append((i + 1, "C8", f"引用了不存在的篇 {m.group(0)}"))

    return issues


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    have = {f.name[:5] for f in notes()}
    files = notes(arg)
    total = 0
    for f in files:
        iss = check(f, have)
        if iss:
            print(f"\n{f.name}")
            for ln, code, msg in iss:
                print(f"  L{ln:<4} [{code}] {msg}")
            total += len(iss)
    print(f"\n━━ 结构检查｜{len(files)} 篇｜{total} 处")
    return 0 if total == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
