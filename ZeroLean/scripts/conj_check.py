#!/usr/bin/env python3
"""共扼覆盖率检查（原型）。

把 `Zero/Conjugation.lean` 的 conjRegistry 当作「已登记共扼检验」的清单，
把 `Zero/*.lean` 的所有 theorem/lemma 当作「待检验」的清单。

- 已登记：有 ConjEntry（proof 已在 Lean 里通过）；
- 未登记：普通 Lean 可能过，但**没有共扼检验** → 这是「缺的定理」候选清单；
- 若某条 ConjEntry 的 proof 失败，Lean build 会直接报错 → 缺的共扼引理暴露。

用法：python3 ZeroLean/scripts/conj_check.py
"""
from __future__ import annotations
import glob
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZERO = os.path.join(ROOT, "Zero")

# 1) 所有 theorem/lemma
theorems: list[tuple[str, str]] = []
for path in sorted(glob.glob(os.path.join(ZERO, "*.lean"))):
    if os.path.basename(path).startswith("Conjugation"):
        continue  # 框架/示例文件自身不计入库定理
    text = open(path, encoding="utf-8").read()
    for m in re.finditer(r"^\s*(?:theorem|lemma)\s+([A-Za-z0-9_\.\']+)", text, re.M):
        theorems.append((m.group(1), os.path.basename(path)))

# 2) conjRegistry 里的登记名（字符串字面量）
registered: list[str] = []
for rp in sorted(glob.glob(os.path.join(ZERO, "Conjugation*.lean"))):
    reg_text = open(rp, encoding="utf-8").read()
    # 只取 conjRegistry / conjRegistryReal 块里的名字；候选（ConjCandidate）不算登记
    for m in re.finditer(
        r"(?:conjRegistry\w*)\s*:\s*List ConjEntry\s*:=\s*\[(.*?)\n\]",
        reg_text,
        re.S,
    ):
        registered += re.findall(r'⟨"([^"]+)"', m.group(1))

# 3) 报告
print(f"theorem/lemma 总数: {len(theorems)}")
print(f"conjRegistry 登记数: {len(registered)}  {registered}")
unregistered = [t for t in theorems if t[0] not in registered]
print(f"未登记共扼检验: {len(unregistered)}")
print("前 20 条未登记（下一批登记候选）：")
for name, f in unregistered[:20]:
    print(f"  - {name}  ({f})")

# 4) JSON 落盘
out = {
    "theorems": len(theorems),
    "registered": registered,
    "unregistered": [{"name": n, "file": f} for n, f in unregistered],
}
with open(os.path.join(ROOT, "scripts", "conj_report.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("→ scripts/conj_report.json")
