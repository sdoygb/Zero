#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z0_axiom_hygiene_check.py -- 全库"公理家谱"卫生检查（防 A0–A5 污染回流）

背景（用户裁决 2026-10-03）：Zero 的公理**只有一条** = `Z0`（零不断乱动）。
`A0–A5` 是旧理论的公理命名，已降为定理；全库行文必须：
  · 当前时态表述里**不得**再把 A0–A5 说成公理；
  · "公理"二字只许出现在 Z0 的语境里；
  · A 标号只允许作为**历史命名**出现（`A0–A5（历史命名）`、`（历史标号 A5）` 这类降级形态）。

本脚本是**回归闸门**：任何文件重新写出"A 是公理"的表述就报错。

  F1  禁止短语（当前时态）：独立公理 A2/A3/A5、"公理 A3"、A0–A5 公理、公理表/公理集（指 A0–A5 时）……
  F2  "公理"出现的上下文逐条分类：只允许 Z0 语境或明确的"历史/已撤回"标注
  F3  家谱映射表在两处在位（G0 §0.1 与 Z0）
  F4  旧文件名 G39_is_A3_* 不得再出现
"""
from __future__ import annotations

import io
import os
import re
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   " + d) if d else ""))


FILES = sorted(f for f in os.listdir(HERE) if f.endswith(".md"))
PYF = sorted(f for f in os.listdir(HERE) if f.endswith(".py") and f != os.path.basename(__file__))
SELF = os.path.basename(__file__)


def txt(f):
    return io.open(os.path.join(HERE, f), encoding="utf-8", errors="replace").read()


# ----------------------------------------------------------------------
head("F1  禁止短语：当前时态不得再把 A0–A5 说成公理")

FORBIDDEN = [
    r"独立公理\s*\**\s*A2",
    r"三条独立公理",
    r"3 条独立公理",
    r"公理\s*A[0-5]\b",
    r"A[0-5]\s*公理",
    r"A0[–-]A5\s*公理",
    r"公理只有\s*(?!\**\s*\[?`?Z0)",
    r"六条公理",
    # ↓ 2026-10-03 补：任何"多于一条"的公理计数（旧句式"三条独立公理"的泛化）
    #   注意："三条定理/三条性质/三条尺子/三条路" 等**不是**公理计数，不在此列
    r"[一二三四五六七八九十兩两0-9]+\s*(?:条|个)\s*公理",
    r"公理\s*(?:有[一二三四五六七八九十兩两]|[一二三四五六七八九十兩两]\s*条)",
    r"(?:Z1|Z2|Z3|Z4|Z5)[^。；\n]{0,12}(?:条|个)\s*公理",
]
HIST = re.compile(r"旧|历史|原|曾|撤回|当时|此前|migrate|A2Z")
# 合法提及：唯一公理就是 Z0 的各种说法
LEGIT = re.compile(
    r"公理只有一条|公理只有 \[?`?Z0|唯一公理.*Z0|公理只有 Z0"
    # 历史重构：描述"当年有几条"（如"从 5 条公理压到 1 条"）
    r"|从.{0,12}\d+\s*条公理"
    r"|压到\s*(?:一|1|两|二)\s*条(?:\s*公理)?"      # 仅"压到 1 条[公理]"，不放行"压到 … 三条独立公理"
    r"|\d+\s*条公理压到"
    # 反事实/否定：如"Z0 之外的第二条公理会违反纪律"
    r"|之外的第[二两三四五六七八九十]条公理|第二条公理\)|第[二两]条公理\）"
)
hits = []
for f in FILES + PYF:
    if f in (SELF, "A2Z_MIGRATION_RECORD.md", "A2Z_MIGRATION_SPEC.md", "migrate_a_to_z.py"):
        continue
    body = txt(f)
    for pat in FORBIDDEN:
        for m in re.finditer(pat, body):
            ln = body[:m.start()].count("\n") + 1
            line = body.split("\n")[ln - 1]
            # 已明确标注为历史/旧记述的行允许保留旧措辞（如『旧记"3 条独立公理"』）
            if HIST.search(line) or LEGIT.search(line):
                continue
            hits.append((f, ln, line.strip()[:110]))
check("全库无『独立公理 A2／A3／A5』『公理 A< n>』『A0–A5 公理』『六条公理』",
      not hits, ("%d 处，例如 %s" % (len(hits), hits[:3])) if hits else "")

# ----------------------------------------------------------------------
head("F1b  牙齿自测：合成违规样本必须被拦、合法写法必须放行")

TEETH = [
    ("三条公理全过", True),
    ("Zero 有三条公理", True),
    ("Z1 定理 1、Z2、Z3 三条公理", True),
    ("三条独立公理（A2、A3、A5）", True),
    ("六条公理", True),
    ("两条公理", True),
    ("公理 A3", True),
    ("A0–A5 公理", True),
    ("公理只有 Z0 一条", False),
    ("公理只有一条：零不断乱动", False),
    ("从 5 条公理压到 1 条", False),
    ("把它们升为公理（即 Z0 之外的第二条公理）会违反纪律", False),
    ("三条尺子夹住 pi", False),
    ("三条定理（Z1 定理 1、Z2、Z3）", False),
    ("三条性质全过（非自反、对称、对序稳定）", False),
]
_teeth_bad = []
for _t, _should in TEETH:
    _flag = any(re.search(_pp, _t) for _pp in FORBIDDEN) and not (HIST.search(_t) or LEGIT.search(_t))
    if _flag != _should:
        _teeth_bad.append((_t, _should, _flag))
check("牙齿自测 %d 例（合成违规必拦 / 合法必放）" % len(TEETH), not _teeth_bad,
      ("误判 %s" % _teeth_bad) if _teeth_bad else "全部符合预期")
if _teeth_bad:
    print("      ^ 注意：测试串里若含 CJK，请在模式里用 [ \\t]* 而不是 \\s*（Python3 的 \\s 不匹配零宽）")

# ----------------------------------------------------------------------
head("F2  『公理』上下文分类：只许 Z0 语境或历史标注")

ALLOW = [
    r"Z0",                     # Z0 自身
    r"历史", r"已撤回", r"旧", r"原公理", r"曾", r"当时",
    r"扩充条款", r"具名输入", r"不设概率", r"识别",
    r"唯一公理", r"公理只有", r"一条公理", r"1 条公理",
    r"公理层之外", r"公理集完备性",
    # 否定式（"不是公理""非公理"）与"不增公理"政策语
    r"不是[^。]{0,8}公理", r"非公理", r"不[是再]\s*公理", r"不算公理",
    r"不增公理", r"不新增公理", r"不引入[^。]{0,6}公理", r"没有[^。]{0,6}公理",
    r"不是公理", r"公理层", r"公理化的", r"公理体系", r"公理化",
]
bad_ctx = []
for f in FILES:
    if f in (SELF, "A2Z_MIGRATION_SPEC.md"):
        continue
    body = txt(f)
    for m in re.finditer(r"公理", body):
        s = max(0, m.start() - 60)
        e = m.end() + 60
        ctx = body[s:e].replace("\n", " ")
        if not any(re.search(a, ctx) for a in ALLOW):
            ln = body[:m.start()].count("\n") + 1
            bad_ctx.append((f, ln, ctx.strip()[:130]))
check("【报告项，不阻断】每处『公理』的上下文分类", True,
      "%d 处落入「Z0 语境／历史标注／否定式」之外，需人工看；"
      "注意 D 系列另有自己的『公理表』（该文登记表），与本条无关" % len(bad_ctx))
if bad_ctx:
    print("      待复核样例：")
    for f, ln, ctx in bad_ctx[:8]:
        print("        %s:%d  %s" % (f, ln, ctx[:100]))

# ----------------------------------------------------------------------
head("F3  家谱映射表在位（G0 §0.1 与 Z0）")

g0 = txt("G0_bottom_layer_and_derivation_route.md")
z0 = txt("Z0_zero_never_rests_single_axiom.md")
check("G0 §0.1 家谱映射表在位", ("§0.1" in g0) and ("历史命名" in g0) and ("Z1 定理 1" in g0))
check("G0 声明公理只有 Z0", ("公理只有" in g0) and ("Z0" in g0))
check("G0 仍保留历史锚点（脚本依赖）", "唯一的基数量" in g0 and "不设概率" in g0)
check("Z0 仍是唯一公理文档", "唯一公理" in z0 or "公理只有一条" in z0)
check("Z0 给出 Z1–Z5 的导出", all(("**Z%d**" % k) in z0 for k in range(1, 6)))

# ----------------------------------------------------------------------
head("F4  旧文件名不得回流")

check("G39 旧文件名 G39_is_A3_* 不存在",
      not os.path.exists(os.path.join(HERE, "G39_is_A3_designed_for_the_theorems.md")))
check("G39 新文件名在位", os.path.exists(os.path.join(HERE, "G39_is_Z2_designed_for_the_theorems.md")))
EXEMPT = {SELF, "A2Z_MIGRATION_SPEC.md", "A2Z_MIGRATION_RECORD.md"}
old_refs = [f for f in FILES + PYF
            if f not in EXEMPT and "G39_is_A3_designed_for_the_theorems" in txt(f)]
check("全库无对旧文件名的引用", not old_refs, ("%s" % old_refs[:4]) if old_refs else "")

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
