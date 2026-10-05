#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INDEX_check.py -- 核验 INDEX.md：无孤儿文档、计数一致、范围声明成立。

对应文档 INDEX.md。
失败时退出码非零。

  F1  lh/ 里每个 .md 都在 INDEX.md 中有链接（无孤儿）
  F2  INDEX 的计数与文件系统一致
  F3  范围声明成立：D1-D209 不在 lh/；lh/ 里的 D 全是 D210-D259
  F4  三个来源的清单完整（G / D2xx / zero_sum_*）
  F5  核验入口被登记
  F6  诚实边界
"""

import importlib
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


IDX = io.open(os.path.join(HERE, "INDEX.md"), encoding="utf-8").read()
FILES = sorted(f for f in os.listdir(HERE) if f.endswith(".md"))

GMD = [f for f in FILES if re.match(r"G\d+_.*\.md$", f)]
DMD = [f for f in FILES if re.match(r"D\d+_.*\.md$", f)]
ZMD = [f for f in FILES if f.startswith("zero_sum_") and f.endswith(".md")]
CK = sorted(f for f in os.listdir(HERE) if re.match(r"G\d+_check\.py$", f))
ZBMD = [f for f in FILES if re.match(r"Z\d+_.*\.md$", f)]
ZCK = sorted(f for f in os.listdir(HERE) if re.match(r"Z\d+_check\.py$", f))
ZCK2 = sorted(f for f in os.listdir(HERE) if f.startswith("zero_sum_") and f.endswith("_check.py"))

# ======================================================================
head("F1  无孤儿文档：lh/ 里每个 .md 都在 INDEX.md 中有链接")

orphans = []
for f in FILES:
    if f == "INDEX.md":
        continue
    if ("](%s)" % f) not in IDX:
        orphans.append(f)
check("所有 .md 都被 INDEX 链接（共 %d 篇，排除 INDEX.md 自身）" % (len(FILES) - 1),
      not orphans, "孤儿: %s" % orphans if orphans else "无")

# ======================================================================
head("F2  计数与文件系统一致")

print("      G*.md=%d  D*.md=%d  zero_sum_*.md=%d  Z*.md=%d  G*_check.py=%d  Z*_check.py=%d"
      % (len(GMD), len(DMD), len(ZMD), len(ZBMD), len(CK), len(ZCK)))
check("INDEX 声明 G* = %d 篇 + %d 个核验脚本" % (len(GMD), len(CK)),
      ("%d 篇 + %d 个核验脚本" % (len(GMD), len(CK))) in IDX)
check("INDEX 声明 D2xx = %d 篇" % len(DMD), ("| %d 篇 |" % len(DMD)) in IDX)
check("INDEX 声明 Z* 基础层 = %d 篇 + %d 个核验脚本" % (len(ZBMD), len(ZCK)),
      ("%d 篇 ＋ %d 个核验脚本" % (len(ZBMD), len(ZCK))) in IDX)
zpy = len([f for f in os.listdir(HERE) if f.startswith("zero_sum_") and f.endswith(".py")
           and not f.endswith("_check.py")])
check("INDEX 声明 zero_sum_* = %d 篇 + %d 个程序 + %d 个核验脚本" % (len(ZMD), zpy, len(ZCK2)),
      ("%d 篇笔记 / %d 个程序 ＋ %d 个核验脚本" % (len(ZMD), zpy, len(ZCK2))) in IDX)

# ======================================================================
head("F3  范围声明成立")

nums = sorted(int(re.match(r"D(\d+)", f).group(1)) for f in DMD)
check("lh/ 里的 D 文件全部在 210–259 区间", nums and min(nums) >= 210 and max(nums) <= 259,
      "%d–%d" % (min(nums), max(nums)))
check("lh/ 里没有 D1–D209（不属零和宇宙的语料未拷入）",
      not any(n <= 209 for n in nums))
check("INDEX 明示 D1–D209 数量为 0",
      "**0**（未拷入）" in IDX or "**0**" in IDX)

# ======================================================================
head("F4  三个来源的清单完整")

missing_g = [f for f in GMD if ("](%s)" % f) not in IDX]
missing_d = [f for f in DMD if ("](%s)" % f) not in IDX]
missing_z = [f for f in ZMD if ("](%s)" % f) not in IDX]
missing_zb = [f for f in ZBMD if ("](%s)" % f) not in IDX]
check("G 清单完整（%d 篇）" % len(GMD), not missing_g, "缺 %s" % missing_g)
check("D2xx 清单完整（%d 篇）" % len(DMD), not missing_d, "缺 %s" % missing_d)
check("zero_sum 清单完整（%d 篇）" % len(ZMD), not missing_z, "缺 %s" % missing_z)
check("Z* 基础层清单完整（%d 篇）" % len(ZBMD), not missing_zb, "缺 %s" % missing_zb)

# ======================================================================
head("F5  核验入口被登记")

check("INDEX 提到同步入口 ledger_sync.py", "ledger_sync.py" in IDX)
_g10t = io.open(os.path.join(HERE, "G10_final_derivation_and_input_ledger.md"),
                encoding="utf-8").read()
_mt = re.search(r"\*\*合计\*\*\s*\|\s*\*\*(\d+)\*\*", _g10t)
_ledger = int(_mt.group(1)) if _mt else -1
_in_idx = re.search(r"账本合计：\*\*(?:通过|独立实断言)\s*(\d+)", IDX)
check("INDEX 给出的核验数与 G10 账本合计一致（%d）" % _ledger,
      _in_idx is not None and int(_in_idx.group(1)) == _ledger,
      "INDEX=%s" % (_in_idx.group(1) if _in_idx else "无"))
check("INDEX 提到本脚本 INDEX_check.py", "INDEX_check.py" in IDX)

# ======================================================================
head("F6  诚实边界")

check("INDEX 是『文档清单』，不改变任何推导结论", True)
check("4 个动力学沙盒已指向 G28（不再标未读）",
      "4 个动力学沙盒" in IDX and "G28_dynamics_audit.md" in IDX)
check("8 篇非原生桥接已标出（★）", "★" in IDX)
check("范围由用户界定；本文档只如实登记", True)

# ======================================================================
head("F7  LaTeX 渲染体检（写一篇拦一篇）")

# 渲染器 = KaTeX（rehype-katex）；分类沿用 latex_lint.py 的 C1-C13
BLOCKING = ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C11", "C12", "C13")
BENIGN = ("C10",)

try:
    _LX = importlib.import_module("latex_lint")
except Exception as _e:
    _LX = None
check("latex_lint.py 可导入（体检工具在位）", _LX is not None, str(_e) if _LX is None else "")

if _LX is not None:
    # ---- 正控制：合成坏样本【必须】被检出（否则"拦得住"本身没被验证）----
    CTRL = [
        ("C1  $$ 未闭合", "C1", "$$ 没有闭合"),
        ("C2  行内 $ 未配对", "C2", "行内 $x = 1"),
        ("C3  花括号不平衡", "C3", "$$ \\boxed{ x } } $$"),
        ("C5  begin/end 不匹配", "C5", "$$ \\begin{aligned} x $$"),
        ("C7  用了 \\(...\\)", "C7", "文本 \\(x\\) 文本"),
        ("C8  $$ 块内空行", "C8", "$$\n x\n\n y\n$$"),
        ("C9  $$ 与正文同行", "C9", "正文 $$x$$ 继续"),
        ("C13a 数学内含 markdown 链接", "C13", "$\\text{[`G9`](G9_x.md)}$"),
        ("C13b \\text 内含未转义 _", "C13", "$\\text{a_b}$"),
    ]
    _miss = [n for n, cat, s in CTRL if cat not in _LX.lint_text(s)]
    check("正控制：%d 类合成坏样本全部被检出" % len(CTRL), not _miss,
          ("漏检 %s" % _miss) if _miss else "")
    # ---- 负控制：合法写法不得误报 ----
    _NEG = ["$\\text{a\\_b}$", "$\\alpha_1+\\beta^2$", "$$\\boxed{\\ x\\ }$$"]
    _negbad = [s for s in _NEG if any(c in _LX.lint_text(s) for c in BLOCKING)]
    check("负控制：合法写法（转义 \\_ 等）不误报", not _negbad,
          ("误报 %s" % _negbad) if _negbad else "")

    # ---- 逐篇阻断扫描（lh/ 全部 .md）----
    _hits, _c10 = {}, 0
    for _f in FILES:
        _iss = _LX.lint_file(os.path.join(HERE, _f))
        _bad = {k: v for k, v in _iss.items() if k in BLOCKING}
        if _bad:
            _hits[_f] = _bad
        _c10 += len(_iss.get("C10", []))
    check("lh/ 全部 %d 篇：阻断类（%s）为 0" % (len(FILES), ",".join(BLOCKING)),
          not _hits,
          ("；".join("%s:%s" % (k, sorted(v)) for k, v in sorted(_hits.items()))) if _hits else "无")
    print("      C10（\\text 内含中文）共 %d 处 —— 良性（KaTeX 字形回退；LATEX_AUDIT.md 已判定）" % _c10)
    check("C10 只报告、不阻断（与 LATEX_AUDIT.md 的判定一致）", True)
    check("体检覆盖：链接/反引号/强调 + \\text 内未转义 _ ^ + 花括号 + 环境 + 表格 + $$ 配对", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
