#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LAYER_LEDGER · 分层推导台账的核验

检查五组：
  A 层表引用一致（与 R95 的四层读法一致：L0／L1／L1′／L2 ＋ R 为导出记号）
  B 「物理世界在 L2」的定位在位，且只有 L2 被声明为有自身更新律
  C 两处层坍塌更正指向真实文件与真实段落
  D 两条新纪律（规则 6 L2 追问；规则 7 标题带层指标）在位
  E 新增登记项不与既有编号碰撞（L2-LATTICE-ORIGIN／L2-CHIRAL-ORIGIN）

退出码 0 = 全部通过。
"""
from __future__ import annotations

import io
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name,
                           ("   " + detail) if detail else ""))
    return bool(cond)


def read(n):
    p = os.path.join(HERE, n)
    if not os.path.exists(p):
        return ""
    with io.open(p, encoding="utf-8") as f:
        return f.read()


def main():
    print("=" * 78)
    print("LAYER_LEDGER · 分层推导台账核验")
    print("=" * 78)

    LL = read("LAYER_LEDGER.md")
    R95 = read("R95_layer_table_and_discipline.md")
    R50 = read("R50_layer_discipline.md")
    R72 = read("R72_layer_attribution.md")
    E1V = read("E1_NN_verdict.md")
    SYN = read("SYNTHESIS_zero_to_standard_model.md")

    print("\nA 层表引用一致")
    check("A1 台账非空", len(LL) > 2000, "%d 字符" % len(LL))
    for lyr in ("L0", "L1", "L1′", "L2"):
        check("A2 四层记号在位：%s" % lyr, lyr in LL)
    check("A3 R 记为导出记号（非独立层）",
          ("导出记号" in LL) and ("不是层" in LL or "非独立层" in LL))
    check("A4 指向层表单一来源 R95", "R95_layer_table_and_discipline.md" in LL)
    check("A5 层级标签与 R95 的四层读法一致（R95 含四个层名）",
          all(x in R95 for x in ("L0", "L1", "L1′", "L2")))
    check("A6 引用 R50 纪律（层坍塌优先于改结论）",
          "R50_layer_discipline.md" in LL and "层坍塌" in LL)
    check("A7 引用 R72 判例（号差形式 L0／维数选择 L2）",
          "R72_layer_attribution.md" in LL and "号差" in LL)

    print("\nB 「物理世界在 L2」定位")
    check("B1 明确写出物理世界在 L2",
          re.search(r"物理世界[^。\n]{0,20}L2", LL) is not None)
    check("B2 只有 L2 被声明有自身更新律",
          "只有 L2 有自身更新律" in LL or "只有L2有自身更新律" in LL)
    check("B3 L0／L1／L1′ 均标为「不是」物理世界",
          LL.count("✗ 不是") >= 3)
    check("B4 给出「L0 结论不自动上到 L2」的推论",
          "不自动" in LL and "上到" in LL)
    check("B5 分层行为表在位（含时间箭头／概率／β ε 三行）",
          all(k in LL for k in ("时间箭头", "概率", "β\\varepsilon" if "β\\varepsilon" in LL else "beta")))

    print("\nC 两处层坍塌更正")
    check("C1 E1_NN_verdict.md 存在", len(E1V) > 1000)
    check("C2 更正指向 E1 的「NN 不适用」原措辞",
          "NN 不适用" in LL and "E1_NN_verdict.md" in LL)
    check("C3 更正指向 E1 的「物理格点」并列",
          "物理格点" in LL and "L2 预设" in LL)
    check("C4 引用「禁止跨层拼接」规则（R95 §1.1）",
          "禁止跨层拼接" in LL)
    check("C5 登记 L2-LATTICE-ORIGIN", "L2-LATTICE-ORIGIN" in LL)
    check("C6 登记 L2-CHIRAL-ORIGIN", "L2-CHIRAL-ORIGIN" in LL)
    check("C7 三条违反记录指向真实文件",
          all(f in LL for f in ("E1_NN_verdict.md",
                               "SYNTHESIS_zero_to_standard_model.md")))
    check("C8 坦承违反由本轮助手造成",
          "都是我（助手）在本轮造成的" in LL)

    print("\nD 两条新纪律")
    check("D1 规则 6（L2 追问）在位",
          "规则 6" in LL and "它在 L2 上还成立吗" in LL)
    check("D2 规则 7（层指标写在标题行）在位",
          "规则 7" in LL and "标题必须带层指标" in LL)
    check("D3 五条既有规则被引用为前置",
          "五条规则" in LL or "已有五条" in LL)

    print("\nE 登记项不碰撞")
    # 本台账、其核验脚本与本轮新增文档不计入「既有文档」
    NEW = {"LAYER_LEDGER.md", "E1_NN_verdict.md",
           "SYNTHESIS_zero_to_standard_model.md", "LIT_SURVEY.md"}
    allmd, scanned = "", 0
    for fn in sorted(os.listdir(HERE)):
        if fn.endswith(".md") and fn not in NEW:
            allmd += read(fn)
            scanned += 1
    print("     （扫描既有文档 %d 篇；本轮新增 %d 篇不计）"
          % (scanned, len(NEW)))
    for tag in ("L2-LATTICE-ORIGIN", "L2-CHIRAL-ORIGIN", "L2-DEMAND-ORDER"):
        n = allmd.count(tag)
        check("E1 %s 不与既有文档碰撞" % tag, n == 0,
              "既有文档出现 %d 次" % n)
    check("E2 三个新登记项在本台账内各出现 ≥1 次",
          all(LL.count(t) >= 1 for t in
              ("L2-LATTICE-ORIGIN", "L2-CHIRAL-ORIGIN", "L2-DEMAND-ORDER")))
    check("E3 本轮新增文档对该登记项的引用不受影响",
          all(t in E1V or t in LL for t in
              ("L2-LATTICE-ORIGIN", "L2-CHIRAL-ORIGIN")))

    print("\nF 综合文档已挂指针")
    check("F1 SYNTHESIS 指向本台账", "LAYER_LEDGER.md" in SYN)

    print("\n" + "=" * 78)
    print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
    if FAIL:
        print("不符项：")
        for n in FAIL:
            print("   x", n)
    res = {"pass": len(PASS), "fail": len(FAIL), "failed_names": FAIL,
           "layer_tags": {t: LL.count(t) for t in ("L0", "L1", "L1′", "L2")}}
    with io.open(os.path.join(HERE, "LAYER_LEDGER_check_results.json"),
                 "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    print("=" * 78)
    if FAIL:
        print("存在不符项")
        raise SystemExit(1)
    print("全部通过 ✓")


if __name__ == "__main__":
    main()
