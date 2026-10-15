#!/usr/bin/env python3
"""生成 Zero 形式化**覆盖矩阵**（报告层的机器部分，可复跑）。

判据（三分归属，2026-10-14 定）：
  每条承重结论必须落在且只落在一个盒子里 ——
    ① 已证        ：kernel 里有对应定理（附依赖清单，由 check.py 给出）
    ② 具名输入    ：建模选择，在账上点名，其推论可进 kernel
    ③ 不可形式化  ：叙述/哲学/文献，写明理由
  **不存在第四种状态**（"看起来证了其实是重言式"那种）。

输入（只读）：
  ZeroLean/Zero/*.lean                 —— 实际进 kernel 的定理（按篇号前缀分组）
  _staging/FORMALIZATION_COVERAGE.md   —— 逐篇覆盖账（F/A/N 分类）
输出：
  _staging/FORMALIZATION_MATRIX.md     —— 每篇一行 + 汇总 + 缺口清单

用法：python3 ZeroLean/scripts/coverage_matrix.py
"""

import re
import sys
import pathlib
import collections

ROOT = pathlib.Path(__file__).resolve().parents[2]
LEAN = ROOT / "ZeroLean" / "Zero"
COVERAGE = ROOT / "_staging" / "FORMALIZATION_COVERAGE.md"
BACKLOG = ROOT / "_staging" / "LEAN_BACKLOG.md"
OUT = ROOT / "_staging" / "FORMALIZATION_MATRIX.md"


def theorems_by_note():
    """篇号 -> [定理名]（从 .lean 文本提取；只认 `theorem zNNN_*` / `lemma zNNN_*`）。"""
    pat = re.compile(r"^\s*(?:theorem|lemma)\s+(z(\d{3})_[A-Za-z0-9_']*)", re.M)
    out = collections.defaultdict(list)
    if not LEAN.is_dir():
        return out
    for f in sorted(LEAN.glob("*.lean")):
        txt = f.read_text(encoding="utf-8", errors="ignore")
        for m in pat.finditer(txt):
            out[int(m.group(2))].append(m.group(1))
    return out


def coverage_by_note():
    """篇号 -> Counter({'F':n,'A':n,'N':n})（从覆盖账提取；行式 `Znnn | 结论 | 类别 | ...`）。"""
    cov = collections.defaultdict(collections.Counter)
    if not COVERAGE.exists():
        return cov
    for line in COVERAGE.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = re.match(r"^\s*\*{0,2}Z(\d{3})\*{0,2}\s*\|(.+)$", line)
        if not m:
            continue
        n = int(m.group(1))
        parts = [p.strip().strip("*` ") for p in m.group(2).split("|")]
        if len(parts) >= 2 and parts[1][:1] in ("F", "A", "N"):
            cov[n][parts[1][:1]] += 1
    return cov


def all_notes():
    """库内全部 Z 篇号（从根目录文件名取）。"""
    notes = {}
    for f in sorted(ROOT.glob("Z*.md")):
        m = re.match(r"^Z(\d{3})_", f.name)
        if m:
            notes[int(m.group(1))] = f.name
    return notes


def status(thm_count, cov):
    if thm_count > 0:
        return "① 已证（部分）" if cov.get("F", 0) > thm_count else "① 已证"
    if cov.get("F", 0) > 0:
        return "⏳ 待形式化（F 类未做）"
    if cov.get("A", 0) > 0:
        return "② 仅具名输入"
    if cov.get("N", 0) > 0:
        return "③ 不可形式化"
    return "❔ 未记账"


def main():
    thms = theorems_by_note()
    cov = coverage_by_note()
    notes = all_notes()
    backlog_txt = BACKLOG.read_text(encoding="utf-8", errors="ignore") if BACKLOG.exists() else ""

    lines = []
    lines.append("# Zero 形式化覆盖矩阵（机器生成，勿手改）\n")
    lines.append("> 由 `ZeroLean/scripts/coverage_matrix.py` 从 kernel 实际状态 + 覆盖账生成，可复跑。")
    lines.append("> 判据：三分归属（① 已证 ／ ② 具名输入 ／ ③ 不可形式化），**不存在第四种状态**。")
    lines.append("> 依赖清单（每条定理用了哪些公理）见 `zero_lean` 的输出或 `ZeroLean/.lean-report.json`。\n")

    tally = collections.Counter()
    lines.append("| 篇 | kernel 定理 | 覆盖账 F/A/N | 状态 |")
    lines.append("|---|---|---|---|")
    for n in sorted(notes):
        t = thms.get(n, [])
        c = cov.get(n, collections.Counter())
        st = status(len(t), c)
        tally[st] += 1
        shown = "、".join(t[:3]) + (f" …(+{len(t)-3})" if len(t) > 3 else "") if t else "—"
        lines.append(f"| Z{n:03d} | {shown} | {c.get('F',0)}/{c.get('A',0)}/{c.get('N',0)} | {st} |")

    lines.append("\n## 汇总\n")
    lines.append(f"- 库内篇数：**{len(notes)}**")
    lines.append(f"- 有 kernel 定理的篇：**{sum(1 for n in notes if thms.get(n))}**")
    lines.append(f"- 手写定理总数：**{sum(len(v) for v in thms.values())}**（口径：`.lean` 里的 `theorem/lemma zNNN_*`）")
    for k, v in tally.most_common():
        lines.append(f"  - {k}：{v} 篇")

    # 缺口清单：F 类未做 / 未记账
    gaps = [n for n in sorted(notes) if not thms.get(n) and cov.get(n, {}).get("F", 0) > 0]
    lines.append(f"\n## 缺口（F 类有内容但 kernel 里没有）—— {len(gaps)} 篇\n")
    lines.append("、".join(f"Z{n:03d}" for n in gaps) if gaps else "（无）")

    lines.append(f"\n## 参考\n- 覆盖账：`_staging/FORMALIZATION_COVERAGE.md`")
    lines.append(f"- 待办账：`_staging/LEAN_BACKLOG.md`（{len(backlog_txt.splitlines())} 行）")
    lines.append("- 审计报告：`_staging/FORMALIZATION_AUDIT.md`")

    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"已写出 {OUT.relative_to(ROOT)}")
    print(f"篇数 {len(notes)}｜有定理 {sum(1 for n in notes if thms.get(n))}｜手写定理 {sum(len(v) for v in thms.values())}｜缺口 {len(gaps)} 篇")
    for k, v in tally.most_common():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    sys.exit(main())
