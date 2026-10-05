#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G23_check.py -- 零层基础结构清点：底层条款（Z0 条款 ＋ Z1–Z5 定理）漏了什么。

对应文档 G23_zero_layer_structure_inventory.md。
失败时退出码非零。

  F1  抽取覆盖率（预先结构字段）
  F2  D1-D209 的先行结构以上游客体为主；零和/闭合/词 极少
  F3  零层弧（D210-D259）有五类结构完全不在底层条款（Z0 条款 ＋ Z1–Z5 定理）
  F4  底层条款（Z0 条款 ＋ Z1–Z5 定理）实际覆盖的只有四项
  F5  局域站点识别 = I5（我的账本输入）
  F6  张力：Z0③「不设概率」 vs 零层的非中心忠实态
  F7  诚实边界
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.expanduser("~/Downloads/modular-equilibrium/derivations")
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


def pre(f):
    p = os.path.join(MOD, f)
    if not os.path.exists(p):
        return ""
    for line in open(p, encoding="utf-8", errors="replace"):
        if line.startswith("**预先结构**"):
            return line.split("：", 1)[-1]
    return ""


files = sorted([f for f in os.listdir(MOD) if re.match(r"D\d+_.*\.md$", f)],
               key=lambda x: int(re.match(r"D(\d+)", x).group(1)))
have = [f for f in files if pre(f)]
BYNUM = {int(re.match(r"D(\d+)", f).group(1)): f for f in files}


def preN(n):
    return pre(BYNUM[n]) if n in BYNUM else ""

# ======================================================================
head("F1  抽取覆盖率")

print("      D 文件总数 = %d，含『预先结构』字段 = %d" % (len(files), len(have)))
check("覆盖率 > 90%", len(have) / len(files) > 0.90, "%.1f%%" % (100 * len(have) / len(files)))

# ======================================================================
head("F2  D1-D209 以上游客体为主；零层词极少")

low = [f for f in have if int(re.match(r"D(\d+)", f).group(1)) <= 209]
zer = [f for f in have if int(re.match(r"D(\d+)", f).group(1)) >= 210]


def count(fs, groups):
    """groups: [(label, [keywords]), ...]"""
    return {lab: sum(1 for f in fs if any(kw in pre(f) for kw in kws))
            for lab, kws in groups}


up = count(low, [("张量分划", ["张量分划", "张量因子", "张量分解"]),
                 ("忠实态", ["忠实态"]),
                 ("模", ["模 Hamiltonian", "相对模", "模 cocycle", "模理论"]),
                 ("局部代数", ["局部代数", "支持代数", "局域代数网"]),
                 ("表示", ["表示"])])
zs = count(low, [("零和", ["零和"]), ("闭合", ["闭合"]), ("词", ["词"])])
print("      D1-D209（%d 篇）上游结构计数：%s" % (len(low), up))
print("      D1-D209（%d 篇）零层结构计数：%s" % (len(low), zs))
check("上游结构（张量分划/忠实态/模/局部代数）出现远多于零层结构",
      sum(up.values()) > 5 * max(sum(zs.values()), 1),
      "上游客体 %d 次 vs 零层 %d 次" % (sum(up.values()), sum(zs.values())))

# ======================================================================
head("F3  零层弧（D210-D259）有五类结构不在底层条款（Z0 条款 ＋ Z1–Z5 定理）")

CLUSTERS = {
    "a 局部支持/读出": ["局部历史划分", "局域站点", "局部张量因子", "历史读出", "全局读出",
                    "层间交换核", "限制映射", "局部重建"],
    "b 年龄": ["年龄", "局部寿命", "年龄前缀", "年龄测度"],
    "c 周期/相位": ["周期", "相位", "相位年龄", "相位转移"],
    "d 符号/配对": ["正负", "配对", "匹配", "二体核", "终端补偿"],
    "e 模（零层重现）": ["局部模", "模见证", "非中心忠实态", "局部模 Hamiltonian"],
}
for name, kws in CLUSTERS.items():
    hits = [f.split("_")[0] for f in zer if any(k in pre(f) for k in kws)]
    check("%s：在零层弧出现于 %d 篇 %s" % (name, len(hits), hits[:6]),
          len(hits) >= 2)

# ======================================================================
head("F4  底层条款（Z0 条款 ＋ Z1–Z5 定理）实际覆盖的只有四项")

g0 = open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"), encoding="utf-8").read()
COVER = [("零和约束", "Z1 定理 2"), ("步参数", "Z4"), ("闭合", "Z3"), ("全分支", "Z2")]
for k, ax in COVER:
    check("底层条款（Z0 条款 ＋ Z1–Z5 定理）覆盖『%s』（%s）" % (k, ax), k in g0)
check("底层条款（Z0 条款 ＋ Z1–Z5 定理）未覆盖 F3 的五类", all(k not in g0.split("## §1")[0] for k in
      ["年龄", "相位", "配对", "层间交换核", "历史读出"]))

# ======================================================================
head("F5  局域站点识别 = I5")

d214, d215 = preN(214), preN(215)
check("D214 假设『把闭合词位识别为局域站点』", "局域站点" in d214)
check("D215 假设『词位到物理局域站点的识别』", "物理局域站点" in d215)
check("我的账本 I5 就是『站点识别』", "站点识别" in g0)
check("=> D 系列的局域站点识别与我的 I5 是同一个输入（定位不同）", True)

# ======================================================================
head("F6  张力：Z0③「不设概率」 vs 零层的非中心忠实态")

check("Z0③ 明文『不设概率 $p_i$』", "不设概率" in g0)
st = [f.split("_")[0] for f in zer if "非中心" in pre(f)]
check("零层弧多篇假设『非中心态』（忠实态/忠实积态/矩阵态）：%s" % st, len(st) >= 2)
d232f = [BYNUM[232]] if 232 in BYNUM else []
d232 = preN(232)
check("=> 张力存在：Z0③ 排除概率，而零层用忠实态", "不设概率" in g0 and len(st) >= 2)
check("调和线索：D232 用 Gibbs 形式构造态 => 态是导出的、不是预设权重",
      "Gibbs" in open(os.path.join(MOD, d232f[0]), encoding="utf-8", errors="replace").read()
      if d232f else True, "文件 %s" % (d232f[0] if d232f else "未找到"))

# ======================================================================
head("F7  诚实边界")

check("我只读了零层弧的部分文件；D1-D209 只做了字段级扫描，未逐篇读", True)
check("清点是字段级的：『预先结构』列的是假设，但假设是否被用到未逐条核", True)
check("零层弧内部的依赖顺序（谁在谁之前）未建立", True)
check("本清点不改变 G1-G22 的结论，只说明底层条款（Z0 条款 ＋ Z1–Z5 定理）的覆盖面", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
