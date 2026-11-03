#!/usr/bin/env python3
"""生成 Zb 新库的形式化覆盖矩阵（报告层，可复跑）。

输入：
  _staging/FORMALIZATION_COVERAGE_Zb.md   —— 65 篇逐篇覆盖账（状态 + 说明）
  Zero/ZbNNN_*.lean                        —— 实际进 kernel 的 `zbNNN_*` 定理
输出：
  _staging/FORMALIZATION_MATRIX_Zb.md      —— 每篇一行 + 汇总 + 缺口清单

口径：
  · ① 已证（部分）：该篇模块里至少有一条 `theorem/lemma zbNNN_*`；
  · ② 具名输入/诠释；③ 数值/拟合；⏳ 待形式化：由覆盖账给出，脚本只搬运。
"""
import re, pathlib, collections

ROOT = pathlib.Path(__file__).resolve().parent.parent
VAULT = ROOT.parent
COVERAGE = VAULT / "_staging" / "FORMALIZATION_COVERAGE_Zb.md"
OUT = VAULT / "_staging" / "FORMALIZATION_MATRIX_Zb.md"

def load_coverage():
    rows = {}
    if not COVERAGE.exists():
        return rows
    for line in COVERAGE.read_text(encoding="utf-8").splitlines():
        m = re.match(r"^\|\s*Zb(\d{3})\s*\|\s*(.*?)\s*\|\s*(.*?)\s*\|\s*(.*?)\s*\|$", line)
        if m:
            rows[int(m.group(1))] = {"title": m.group(2), "status": m.group(3), "note": m.group(4)}
    return rows

def load_theorems():
    thms = collections.defaultdict(list)
    mods = {}
    for p in sorted((ROOT / "Zero").glob("Zb[0-9][0-9][0-9]_*.lean")):
        m = re.match(r"Zb(\d{3})_", p.name)
        if not m:
            continue
        n = int(m.group(1))
        mods[n] = p.name
        txt = p.read_text(encoding="utf-8", errors="ignore")
        for mm in re.finditer(r"^\s*(?:theorem|lemma)\s+(zb\d{3}_[A-Za-z0-9_']*)", txt, re.M):
            thms[n].append(mm.group(1))
    return thms, mods

def main():
    cov = load_coverage()
    thms, mods = load_theorems()
    notes = sorted(set(cov) | set(thms))
    lines = []
    lines.append("# Zero 形式化覆盖矩阵（Zb 新库，机器生成）\n")
    lines.append("> 由 `ZeroLean/scripts/coverage_matrix.py` 从 kernel 实际状态 + 覆盖账生成，可复跑。")
    lines.append("> 判据：① 已证（kernel 定理）／② 具名输入或诠释／③ 数值或拟合／⏳ 待形式化。")
    lines.append("> **注意**：⏳ 不是「已证」，只是「有数学核、本轮未做」；① 也只覆盖该篇的**部分**命题。\n")

    tally = collections.Counter()
    lines.append("| 篇 | 模块 | kernel 定理 | 状态 | 说明 |")
    lines.append("|---|---|---:|---|---|")
    for n in notes:
        row = cov.get(n, {"title": "(未登记)", "status": "❔", "note": ""})
        t = thms.get(n, [])
        st = row["status"]
        if t and not st.startswith("①"):
            st = "① 已证（部分）"
        tally[st] += 1
        shown = "、".join(t[:2]) + (f" …(+{len(t)-2})" if len(t) > 2 else "") if t else "—"
        lines.append(f"| Zb{n:03d} | {mods.get(n, '—')} | {len(t)} | {st} | {row['note']} |")

    lines.append("\n## 汇总\n")
    lines.append(f"- 篇数：**{len(notes)}**")
    lines.append(f"- 有 kernel 定理的篇：**{sum(1 for n in notes if thms.get(n))}**")
    lines.append(f"- `zbNNN_*` 定理数：**{sum(len(v) for v in thms.values())}**")
    for k, v in tally.most_common():
        lines.append(f"  - {k}：{v} 篇")

    gaps = [n for n in notes if cov.get(n, {}).get("status", "").startswith("⏳")]
    lines.append(f"\n## 待形式化清单（⏳）—— {len(gaps)} 篇\n")
    lines.append("、".join(f"Zb{n:03d}" for n in gaps) if gaps else "（无）")
    lines.append("\n## 参考\n- 覆盖账：`_staging/FORMALIZATION_COVERAGE_Zb.md`")
    lines.append("- 机器检查：`python3 ZeroLean/scripts/check.py --coverage`")
    lines.append("- 构建：`cd ZeroLean && export PATH=\"$HOME/.local/lean4/bin:$PATH\" && lake build`")

    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"已写出 {OUT}")
    print(f"篇数 {len(notes)}｜有定理 {sum(1 for n in notes if thms.get(n))}｜定理 {sum(len(v) for v in thms.values())}")
    for k, v in tally.most_common():
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
