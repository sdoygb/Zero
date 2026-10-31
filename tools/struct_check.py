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
  C9 控制字符     控制字符/零宽字符混进正文（如 \\x08 冒充 \\b，终端里看不出来）
  C10 自检表覆盖  正文的每条定理/命题/推论/引理/界定 都应在 `## 自检` 表里有一行
  C11 跨篇条目     ZbNNN 命题 X.Y 的条目必须在那一篇里真的存在（C8 只查篇号）
  C12 行内公式     行内 `$...$` 不许跨行断开（`\\rho` 被吃成换行时就会这样）

c10 的口径：自检表覆盖是**存量债**（2026-10-31 首测 146 处），默认只统计不判死；
基线写死在本文件里，只许降不许升——新写的内容漏了自检表会被算成硬失败。

用法：
    python3 tools/struct_check.py              # 全库（C1–C9 硬失败；C10 只统计）
    python3 tools/struct_check.py --strict     # C10 也计入退出码
    python3 tools/struct_check.py Zb047        # 单篇
    python3 tools/struct_check.py Zb047 --strict
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ITEM_RE = re.compile(r"\*\*(定理|命题|推论|定义|引理|界定|设定|记账规则) ?(\d+)\.(\d+)(′?)")
TAG_RE = re.compile(r"\\tag\{(\d+)\.(\d+)([a-z]?)\}")
SEC_RE = re.compile(r"^## §(\d+(?:\.\d+)?)")
C10_BASELINE = 98   # 自检表覆盖的存量欠账；只许降不许升
# 零宽与双向控制字符：肉眼不可见，会把"看起来对的"文本弄坏
ZERO_WIDTH = {"\u200b", "\u200c", "\u200d", "\u200e", "\u200f",
              "\u202a", "\u202b", "\u202c", "\u202d", "\u202e", "\ufeff"}


def notes(pattern=None):
    fs = sorted(ROOT.glob("Zb[0-9][0-9][0-9]_*.md"))
    if pattern:
        fs = [f for f in fs if f.name.startswith(pattern)]
    return fs


def build_item_index():
    """每篇有哪些条目：{「ZbNNN」: {(类型, 编号), ...}}。给 C11 用。"""
    idx = {}
    for f in notes():
        s = f.read_text(encoding="utf-8")
        idx[f.name[:5]] = {(m.group(1), f"{m.group(2)}.{m.group(3)}{m.group(4)}")
                           for m in ITEM_RE.finditer(s)}
    return idx


def check(path, have, index=None):
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

    # ---- C9 控制字符 / 零宽字符
    # 白名单只留 \n（换行）与 \t（缩进）。其余控制字符一律报——`\x08` 冒充 `\b`
    # 这类错在终端里与正确写法长得一样，只能靠字节层查。
    for i, ln in enumerate(lines):
        for j, ch in enumerate(ln):
            o = ord(ch)
            if o >= 32 and o != 127:
                continue
            issues.append((i + 1, "C9", f"第 {j + 1} 列有控制字符 U+{o:04X}（{ln[max(0,j-14):j]!r}▸）"))
            break
        else:
            for j, ch in enumerate(ln):
                if ch in ZERO_WIDTH:
                    issues.append((i + 1, "C9", f"第 {j + 1} 列有零宽字符 U+{ord(ch):04X}"))
                    break

    # ---- C10 自检表覆盖
    # 查「正文有、表里没有」。反过来（表里有、正文没有）由 C6 管。
    hit = None
    for i, ln in enumerate(lines):
        if re.match(r"^## 自检", ln):
            hit = i
            break
    if hit is None:
        issues.append((1, "C10", "本篇没有 `## 自检` 表"))
    else:
        tail = "\n".join(lines[hit:])
        # 表里允许紧凑写法：「命题 3.2／3.3／3.4」「定理 5.1／5.2／推论 5.3」
        # 故一个类型词之后要把整串数字都收下，不能只认第一个。
        covered = set()
        for m in re.finditer(
                r"(定理|命题|推论|定义|引理|界定|设定|记账规则)\s*"
                r"(\d+\.\d+)((?:[／/、]\s*\d+\.\d+)*)", tail):
            t, first, rest = m.group(1), m.group(2), m.group(3)
            covered.add((t, first))
            for n in re.findall(r"\d+\.\d+", rest):
                covered.add((t, n))
        need = {"定理", "命题", "推论", "引理", "界定"}      # 定义/设定/记账规则不强制进表
        for i, t, a, b, c, _ in items:
            if t not in need:
                continue
            if (t, f"{a}.{b}{c}") not in covered:
                issues.append((i + 1, "C10", f"{t} {a}.{b}{c} 不在自检表里"))

    # ---- C11 跨篇条目引用
    # C8 只查篇号存在；「Zb015 命题 4.2」里 Zb015 在、但它的命题 4.2 不在（那是推论 4.2）。
    if index:
        for i, ln in enumerate(lines):
            for m in re.finditer(r"(Zb\d{3})[^。\n|]{0,4}?(定理|命题|推论|引理|界定)\s*(\d+\.\d+)(′?)", ln):
                tgt, t, num = m.group(1), m.group(2), m.group(3) + m.group(4)
                if tgt == path.name[:5] or tgt not in index:
                    continue
                if (t, num) not in index[tgt]:
                    issues.append((i + 1, "C11", f"{tgt} 没有 {t} {num}（引用了不存在的对象）"))

    # ---- C12 行内公式跨行
    # 行内 `$...$` 必须在一行内闭合。`\rho` 被 Python 的 `\r` 吃掉再规范化成换行时，
    # 公式会从中间断开——这种损伤 C9 查不到（没有控制字符残留），只能查配对。
    in_display = False
    for i, ln in enumerate(lines):
        parts = ln.split("$$")
        # 每段 `$$` 交替翻转显示态
        for k, seg in enumerate(parts):
            if k % 2 == 1:
                in_display = not in_display
                continue
            if in_display:
                continue
            body = re.sub(r"\\\$", "", seg)          # 去掉转义的 \$
            if body.count("$") % 2:
                issues.append((i + 1, "C12", "行内 $ 未配对（公式可能被断行）"))
                break

    return issues


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    strict = "--strict" in sys.argv
    arg = args[0] if args else None
    have = {f.name[:5] for f in notes()}
    files = notes(arg)

    index = build_item_index()
    hard = 0                 # C1–C9：硬失败
    soft = []                 # C10：存量债，默认不判死
    for f in files:
        iss = check(f, have, index)
        h = [(ln, c, m) for ln, c, m in iss if c != "C10"]
        s10 = [(ln, c, m) for ln, c, m in iss if c == "C10"]
        if h:
            print(f"\n{f.name}")
            for ln, code, msg in h:
                print(f"  L{ln:<4} [{code}] {msg}")
        hard += len(h)
        if s10:
            soft.append((f.name, s10))

    if soft:
        tot = sum(len(v) for _, v in soft)
        print(f"\n── 自检表覆盖（C10，存量债）：{len(soft)} 篇 / {tot} 处")
        if strict:
            for name, v in soft:
                print(f"  {name[:5]}: {len(v)} 处")
        else:
            print("   （默认只统计；加 --strict 才计入退出码）")
        if strict:
            hard += tot                      # 严格模式：一处都不能少
        elif tot > C10_BASELINE:
            print(f"   ⚠ 比基线 {C10_BASELINE} 多 {tot - C10_BASELINE} 处——新写的内容漏了自检表")
            hard += tot - C10_BASELINE
        elif tot < C10_BASELINE:
            print(f"   ✓ 比基线 {C10_BASELINE} 少 {C10_BASELINE - tot} 处（可下调 C10_BASELINE）")

    print(f"\n━━ 结构检查｜{len(files)} 篇｜硬失败 {hard} 处"
          + (f"｜自检表待补 {sum(len(v) for _, v in soft)} 处" if soft else ""))
    return 0 if hard == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
