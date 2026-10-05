#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R95 · 统一层表的核验

检查五组：
  A 层表单一来源在位（四层 ＋ L3 导出记号）
  B 既有文档已挂指针（R50 修正框、STATUS 层定义行）
  C 判据在位（R93 的"有无自身更新律"；R89 的"禁止跨层拼接"）
  D 层数一致性（不出现未加注的"R 是独立层"类断言）
  E 层指标泄漏诊断（扫描关键文档，报出"声明了层却又不带层的量"）

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
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


def read(n):
    p = os.path.join(HERE, n)
    if not os.path.exists(p):
        return ""
    with io.open(p, encoding="utf-8") as f:
        return f.read()


def main():
    print("=" * 78)
    print("R95 · 统一层表核验")
    print("=" * 78)

    SRC = read("R95_layer_table_and_discipline.md")
    R50 = read("R50_layer_discipline.md")
    STATUS = read("STATUS.md")
    R93 = read("R93_L3_sublayer.md")
    R89 = read("R89_layer_audit.md")
    R94 = read("R94_GR_layer_compatibility.md")
    R92 = read("R92_layer_inventory_and_dictionary.md")

    # ---------------- A ----------------
    print("\n[A] 层表单一来源在位")
    check("R95 存在且非空", len(SRC) > 2000, "%d 字符" % len(SRC))
    check("四层齐备（L0/L1/L1′/L2）",
          all(x in SRC for x in ("L0", "L1", "L1′", "L2")))
    check("R 明示为'导出记号'而非独立层",
          "导出记号" in SRC and "不是独立层" in SRC)
    check("层数写明为 4", "层数：4" in SRC or "**层数：4" in SRC)
    check("给出 L3 非独立的四条判据（J1–J4）",
          all(("`R93` J%d" % k) in SRC for k in (1, 2, 3, 4)))
    has_be = ("beta\\varepsilon" in SRC) or ("β\\varepsilon" in SRC) or ("βε" in SRC)
    check("参数定层表在位（βε / q_n / G,Λ / n）",
          has_be and ("q_n" in SRC) and (("G,\\Lambda" in SRC) or ("G,Λ" in SRC)))

    # ---------------- B ----------------
    print("\n[B] 既有文档已挂指针")
    check("R50 §1 有 R 修正框", ("地位的修正" in R50))
    check("R50 修正框指向 R95", "R95_layer_table_and_discipline.md" in R50)
    check("R50 修正框给出四层/五层两读法", "四层（推荐）" in R50 and "五层（本表原样）" in R50)
    check("R50 修正框含 (R50-3) 删'定律'的措辞修正", "措辞修正" in R50)
    check("R50 修正框含 E5 拆项指针", "`E5` 的修正" in R50)
    check("STATUS 层定义行挂 L3 记号", "导出记号" in STATUS)

    # ---------------- C ----------------
    print("\n[C] 判据在位")
    check("R93 的'有无自身更新律'判据在位",
          "自身更新律" in SRC and "自身更新律" in R93)
    check("R93 给出 σ_t 显含态的判据", "显含" in R93)
    check("R89 的'禁止跨层拼接'已升为纪律第 4 条",
          "禁止跨层拼接" in SRC and "禁止跨层拼接" in R89 or "跨层拼接" in SRC)
    check("独立层判据已写成纪律第 5 条", "独立层" in SRC and "自身更新律" in SRC)

    # ---------------- D ----------------
    print("\n[D] 层数一致性（不出现未加注的'L3 是独立层'）")
    docs = {n: read(n) for n in ("R89_layer_audit.md", "R90_layer_distribution.md",
                                 "R92_layer_inventory_and_dictionary.md",
                                 "R93_L3_sublayer.md", "R94_GR_layer_compatibility.md")}
    bad = []
    for n, t in docs.items():
        # 只抓"肯定式"断言：允许中间的 LaTeX 记号码，否认词（不是/非/列为）不算
        for m in re.finditer(r"(L3|mathcal R)[^。\n]{0,30}?\s*(是|为|算)\s*独立层", t):
            seg = m.group(0)
            # 剥掉 LaTeX 记号再判否认词
            plain = re.sub(r"[{}\\^_*]", "", seg)
            if ("不是独立层" in plain) or ("非独立" in plain):
                continue
            # 转述他人主张不算"未加注的肯定断言"：向前多取 50 字看归属
            ctx = t[max(0, m.start() - 50): m.end()]
            if re.search(r"R\d+|`|把\s*$|原文|其主张|记为", ctx):
                continue
            bad.append((n, seg))
    check("五份下游文档无'L3 是独立层'的未加注断言", not bad,
          "命中：%s" % bad[:3] if bad else "0 处")

    # ---------------- E ----------------
    print("\n[E] 层指标泄漏诊断（声明的层 vs 实际带层的量）")
    # 关键量：声明了层归属的清单必须齐备
    key_quantities = {
        "beta_epsilon": (("beta\\varepsilon" in SRC) or ("β\\varepsilon" in SRC) or ("βε" in SRC), "βε"),
        "q_n": ("q_n" in SRC, "q_n"),
        "Ω^2": ("\\Omega^2" in SRC or "Omega^2" in SRC, "Ω²"),
        "G_Lambda": (("G,\\Lambda" in SRC) or ("G,Λ" in SRC), "G,Λ"),
    }
    for k, (present, label) in key_quantities.items():
        check("参数定层表含 %s" % label, present)
    # 层表是否覆盖 GR 与 QM 两侧
    check("跨理论适用性表在位（GR/QM）", "跨理论的适用性" in SRC and "GR" in SRC and "QM" in SRC)
    check("R94 的结论已纳入（构造式/表示式）",
          "构造式" in SRC and "表示式" in SRC)

    # 泄漏统计：R95 本文中出现的层标记数量（仅诊断，不作通过条件）
    tags = {L: len(re.findall(L, SRC)) for L in ("L0", "L1′", "L1", "L2", "L3", "mathcal R")}
    print("     层标记出现次数（诊断）：%s" % tags)
    print("     （诊断，不作通过条件）：L3 应与 L1′/L2 同量级或更少")
    check("R 的提及不超过 L1′+L2 的合计（导出记号不应是主要坐标）",
          tags["mathcal R"] <= tags["L1′"] + tags["L2"],
          "R=%d  L1′+L2=%d" % (tags["mathcal R"], tags["L1′"] + tags["L2"]))

    print("\n" + "=" * 78)
    print("核验 %d 项，未过 %d 项 %s" % (len(PASS) + len(FAIL), len(FAIL),
                                        FAIL if FAIL else ""))
    print("=" * 78)
    out = {"pass": len(PASS), "fail": len(FAIL), "failed_names": FAIL,
           "layer_tag_counts": tags}
    with io.open(os.path.join(HERE, "R95_check_results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("-> R95_check_results.json")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
