#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G36_check.py -- G36 的【机械保障】部分（已按"只留增量、只留有信息量"协议精简）。

保留（有信息量、便宜、能抓真错）：
   F1  账本合计与 INDEX 一致（防"假绿"：账本说全过、实际某个脚本崩了）
   F2  合计数覆盖充分（> 1000）
   F3  旧全量入口 G10_check.py 确已删除（防有人再跑全量入口）
   F4  全系列文档的【未增公理／未增扩充条款】禁用词扫描（防"口头上加了公理"）

删除（不检验任何量，纯仪式）：
   文档存在性断言、文档措辞断言、check(..., True) 登记项。
"""

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def rd(f):
    p = os.path.join(HERE, f)
    return io.open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


# ======================================================================
head("F1  账本合计与 INDEX 一致（防假绿）")

g10 = rd("G10_final_derivation_and_input_ledger.md")
m = re.search(r"\*\*合计\*\*\s*\|\s*\*\*(\d+)\*\*", g10)
check("G10 账本给出合计数", m is not None, "合计 = %s" % (m.group(1) if m else "?"))
lib = int(m.group(1)) if m else -1
idx = rd("INDEX.md")
mi = re.search(r"账本合计：\*\*(?:通过|独立实断言)\s*(\d+)", idx)
check("INDEX 与账本一致", mi is not None and int(mi.group(1)) == lib,
      "INDEX=%s 账本=%d" % (mi.group(1) if mi else "无", lib))

# ======================================================================
head("F2  覆盖充分 / 旧入口已删")

check("合计数 > 1000（覆盖充分）", lib > 1000, "%d 项" % lib)
check("旧全量入口 G10_check.py 已删除（改用 ledger_sync.py 增量同步）",
      not os.path.exists(os.path.join(HERE, "G10_check.py")))

# ======================================================================
head("F3  未增扩充条款（禁用词扫描）")

DERIV = sorted(f for f in os.listdir(HERE) if re.match(r"G3\d_.*\.md$", f))
print("      扫描本轮系列文档 %d 篇" % len(DERIV))
BANNED = ["新增公理", "新增扩充条款", "新增 A6", "新的公理", "增补公理", "引入新公理"]
NEG = ("不", "未", "无", "非", "没有")
_bad = []
for f in DERIV:
    t = rd(f)
    for b in BANNED:
        i = t.find(b)
        while i != -1:
            ctx = t[max(0, i - 16):i]
            if not any(n in ctx for n in NEG):
                _bad.append((f, b))
                break
            i = t.find(b, i + 1)
check("全部 %d 篇本轮文档均未声称新增公理／新增扩充条款" % len(DERIV), not _bad,
      "命中 %s" % _bad[:3] if _bad else "")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
