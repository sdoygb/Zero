#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""style_check.py —— 工作纪律《文档纪律：正文不留工作痕迹》的机械检查器。

用途：扫根目录 Z*.md，命中即报（行号 + 词 + 该行片段）。**两类**命中：
  ①【工作痕迹】版本比较标签／过程叙述／免责墙／内部台账（§一／§五）；
  ②【断语】no-go **断言**（判死／否证／不可能／被排除…）与**偷带** no-go 的从句（§三）。
纪律：正文只留两类内容 —— ①已证的事实与证明；②未证的开放项（且开放项须点出背靠的定理）。
      违禁的是**版本比较标签**、**过程叙述**、以及**no-go 断语**；
      **认识论标签**（已证／条件／具名输入／开放…）保留。

例外：该纪律文档自身（要引述违禁词才能定义它们）。

用法：
  python3 _staging/style_check.py              # 全库检查（退出码 0=干净，1=有命中）
  python3 _staging/style_check.py Z002_*.md    # 只查指定篇
  python3 _staging/style_check.py --summary    # 只打印每篇命中数
"""
import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXEMPT = {"文档纪律_正文不留工作痕迹.md"}

# ── 违禁词表（见该纪律文档 §一／§五）：版本比较标签 ＋ 过程叙述 ＋ 内部台账 ──
# 2026-10-16 扩表：上一版只收了「不主张」等规范写法，漏了**同义换词**——
#   免责墙的 不得主张／不声称；版本比较的 初稿／上一轮／前几轮／早期版本／修正后；
#   台账类的 台账／登记表／归属表。这些都是自审时实查出来的，见 HANDOFF。
FORBIDDEN = [
    # 版本比较标签
    "已收窄", "已被收窄", "已被取代", "升档", "降档", "不再承重", "仍成立",
    "初稿", "上一轮", "前几轮", "早期版本", "修正后", "旧版本",
    # 过程叙述
    "本次修订", "状态更新", "后期篇", "本轮", "第一版", "早先版本",
    "口径提醒", "据实", "不美化", "诚实边界",
    # 免责墙（一个语义家族，写法有多种）
    "不主张", "不得主张", "不声称",
    # 内部台账
    "台账", "登记表", "归属表", "审计报告",
    # 注意：`档位` 是认识论词（结论的地位），不在违禁表；被禁的是它的版本比较用法（升档／降档／档位由…→）
    # 注意：`落点`、`闭合判据` 都有两义（见下表），故不作裸词违禁。

    # 更正（对旧版本的）
    "一处更正", "更正",
]
# ── 正则型违禁（只抓台账义/语义型，避免误伤正常用词）──
# 2026-10-16 加最后三条：**语义型**违规。加它们的理由（实查教训）——
#   词表路线追不完：连做 5 轮同义词扫描，每轮都能再找出十几处
#   （不主张→不得主张/不声称；第一版→初稿/上一轮/前几轮/早期版本/修正后；
#    更正→撤回/订正/改写；以及"我曾在/从A改写为B/本篇不改写Z"这类框架）。
#   故改抓**语义模式**，与词表并用。这三条覆盖的是：
#     ① 自指过程叙述（我曾在…、本文此前的判定…）
#     ② 版本比较框架（把 A 从 X 改写为 Y）
#     ③ 改篇声明（本篇不改写 Znnn）
FORBIDDEN_RE = [
    # 落点：台账义才禁；"步的落点图"（拓扑）、"自指难题的落点"（难题落在哪）是正常数学用词
    (r"缺口的落点|登记项的落点|待办的落点", "落点（内部台账义）"),
    # 闭合判据：物理义才合法（Z124「禁闭＝闭合判据 ∮A=2πn」）；
    #   台账义（指 closure_verdict 之类的内部判据"已钉死/已登记"）才禁。
    (r"闭合判据[^。；\n]{0,12}(已钉死|已登记|已固定)", "闭合判据（内部台账义）"),
    # 语义型（见上）
    (r"(我|本文|本纲领|本版|本篇)(曾|此前|先前|上一轮|前几轮|早先)", "自指过程叙述"),
    (r"从[^\n]{0,14}改写为", "版本比较框架（从 A 改写为 B）"),
    (r"(本篇|本文)不(改写|改动)\s*Z\d", "改篇声明（过程叙述）"),
    # 动作型：撤回／订正（＝更正）。
    #   `(?<!已)` 是为了放过「已撤回」——那是 Z001／Z015 定义的**档位词**
    #   （术语对照表 `| 已撤回 | retracted | 只能作历史记录 |`），属**档位体系**的决定，
    #   不是措辞问题，未在本次清理中改动（已另行上报）。
    (r"(?<!已)撤回", "撤回（动作／版本比较）"),
    (r"订正", "订正（＝更正）"),
]
# ── 断语（no-go 断言）表 ──
# 出处：纪律 §三「★ 禁止『偷带』」——改写 no-go 时**不许**用「而且它不涉及 X」「X 被排除」
#   这类从句把 no-go 带回来；正面的长度应**恰好等于机制的正面内容**。
# 判据：正文只写「做过的 ＋ 有正面结果的 ＋ 已严格证明的」；「判死／否证／不可能／被排除」
#   这类**断言**不写（no-go 的**前提与机制**要写成正面事实；判定本身移出台账）。
# 现状（2026-10-17）：下列裸词**全库 0 处**——本表是**回归护栏**，不是待清理清单。
# 相邻词的处置（避免假阳性，实查后定）：
#   · `排除`   大量正常用法（排除低维／排除条件／排除了一整条路线）⟹ **不禁**；
#               只禁被动断语句式（见 VERDICT_RE）。
#   · `不成立` 89 处，用于证伪某**假设／条件**（"该假设在本文的格上不成立"）⟹ **不禁**。
#   · `取空解`／`空集` 是**认可的正面写法**（＝解集为空，Z013／Z009）⟹ 不禁。
#   · `档位词 排除 / 不可导出` 属术语对照表定义的**认识论槽位**（Z001 §术语表）⟹ 不禁。
VERDICT = [
    "判死", "否证", "被排除", "不可能",
    "不可导出", "不能导出",
    "无载体", "无来源", "不承载", "被禁", "无对象",
]
# ── 断语的正则型：**偷带 no-go 的从句**（§三 明文禁）──
VERDICT_RE = [
    (r"而且[^\n]{0,12}不涉及", "偷带 no-go（“而且它不涉及…”）"),
    (r"[^\n]{0,16}不被[^\n，。；、]{0,12}承载", "偷带 no-go（“X 不被 Y 承载”）"),
    (r"[^\n]{0,16}没有[^\n，。；、]{0,10}的对象", "偷带 no-go（“没有…的对象”）"),
]

# ── 保留词（认识论标签）：出现在同一行不改变判定，仅用于人工复核时对照 ──
KEEP = ["已证", "导出", "条件", "具名输入", "开放", "不可导出", "排除", "数值律", "数学定理", "约定", "档位"]


def scan(path: pathlib.Path):
    hits = []
    for i, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
        for w in FORBIDDEN:
            if w in line:
                hits.append((i, w, line.strip()[:110]))
                break
        else:
            for pat, tag in FORBIDDEN_RE:
                if re.search(pat, line):
                    hits.append((i, tag, line.strip()[:110]))
                    break
            else:
                # 断语（no-go 断言）：与工作痕迹分开标注，便于按性质处置
                for w in VERDICT:
                    if w in line:
                        hits.append((i, f"[断语] {w}", line.strip()[:110]))
                        break
                else:
                    for pat, tag in VERDICT_RE:
                        if re.search(pat, line):
                            hits.append((i, f"[断语] {tag}", line.strip()[:110]))
                            break
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*", help="默认：根目录全部 Z*.md")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()

    files = [pathlib.Path(f) for f in a.files] or sorted(ROOT.glob("Z[0-9][0-9][0-9]_*.md"))
    total, dirty, nverdict = 0, 0, 0
    for p in files:
        if p.name in EXEMPT:
            if not a.summary:
                print(f"（例外，跳过）{p.name}")
            continue
        hits = scan(p)
        if not hits:
            continue
        dirty += 1
        total += len(hits)
        nverdict += sum(1 for _, w, _ in hits if w.startswith("[断语]"))
        if a.summary:
            print(f"{len(hits):4d}  {p.name}")
        else:
            print(f"\n━━ {p.name}  （{len(hits)} 处）")
            for i, w, frag in hits:
                print(f"  L{i:<5d} [{w}]  {frag}")
    tail = f"\n合计：{dirty} 篇 / {total} 处命中（检查 {len(files)} 篇；例外 {len(EXEMPT)} 篇）"
    if total:
        tail += f"\n　　　其中【断语】{nverdict} 处／【工作痕迹】{total - nverdict} 处"
    print(tail)
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
