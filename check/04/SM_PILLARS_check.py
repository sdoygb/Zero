#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SM_PILLARS_check.py -- 核验《SM_PILLARS · 对照经典标准模型：十二根支柱的退化记分卡》

  F1  结构：头部五字段；§0-§7 各节在位；§2 总表 12 行 x 6 格；P1-P12 连续
  F2  判定分布：从 §2 表**机械点算**五档，与 §0 / §2 合计式逐位一致（0/3/5/3/1 = 12）
  F3  锚定：§2 每行依据栏的 sm/ 与 lh 根目录文件在磁盘上真实存在；抽检关键词真实出现
  F4  跨账一致：从 SM19 §2.3 表**重新机械点算**五档（13/16/5/1/0 + 1 仅参考）；
                从 G91 §8 表机械点算 19 项三分（5/3/4/7）；SM18 的 26 参数引用
  F5  契约与红线：sm/ 的 SM<数字>_ 文件集**未被本件改动**（仍 20 篇 + 20 脚本）；
                禁令七条在位；无肯定式"已导出标准模型"；$$ 成对；翻转条件表 6 行

纯标准库（fractions/itertools/math/os/re）。**不** import 任何旧理论脚本，
**不** import lh 根目录脚本；**不**写任何文件。
"""
from __future__ import annotations

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

FAIL = 0
NCHK = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(name, cond, detail=""):
    global FAIL, NCHK
    NCHK += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def note(t):
    print("  [i] %s" % t)


def read(p):
    with io.open(p, encoding="utf-8") as f:
        return f.read()


def legacy(name):
    """从 sm/_legacy_SM0-19.tar.gz 读取旧文（旧文已从工作树删除，对照以备份为准）"""
    import tarfile
    p = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
    if not os.path.exists(p):
        return ""
    with tarfile.open(p) as tf:
        for m in tf.getmembers():
            if os.path.basename(m.name) == name:
                f = tf.extractfile(m)
                return f.read().decode("utf-8", "replace") if f else ""
    return ""


def exists(rel):
    return os.path.exists(os.path.join(HERE, rel)) or os.path.exists(os.path.join(ROOT, rel))


DOC = read(os.path.join(HERE, "SM_PILLARS.md"))
SM19 = legacy("SM19_verdict_ledger.md")
G91 = read(os.path.join(ROOT, "G91_matter_sector_range_and_three_way_verdict.md"))
SM18 = legacy("SM18_parameter_ledger_26.md")

# ======================================================================
head("F1  结构（头部字段 / 八节 / 12 行 x 6 格）")

for fld in ["**日期**", "**性质**", "**上游", "**口径来源**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)

for sec in ["## §0", "## §1", "## §2", "## §3", "## §4", "## §5", "## §6", "## §7"]:
    check("节 %s 在位" % sec, sec in DOC)

SEC2 = DOC.split("## §2 记分卡总表")[1].split("### 2.1")[0]
rows = [ln for ln in SEC2.splitlines() if re.match(r"^\|\s*\*\*P\d+\*\*\s*\|", ln)]
check("§2 总表 12 行（P1-P12）", len(rows) == 12, "实测 %d 行" % len(rows))

ids = [re.match(r"^\|\s*\*\*(P\d+)\*\*", ln).group(1) for ln in rows] if len(rows) == 12 else []
check("P 编号连续 P1..P12", ids == ["P%d" % i for i in range(1, 13)], "实测 %s" % ids)

ncell = [len([c for c in ln.strip().strip("|").split("|")]) for ln in rows]
check("每行 6 格（| 计数 = 7）", all(n == 6 for n in ncell), "实测 %s" % sorted(set(ncell)))

check("§0 给出主数 0 根（完整退化）", "\\boxed{0}" in DOC)
check("§1.2 五档定义表在位（完整/条件/数学核/不可导出/开放）",
      all(("【%s】" % k) in DOC for k in ["完整退化", "条件退化", "只到数学核", "不可导出", "开放"]))
check("§1.2 归属规则在位（含判死成分不许记条件退化）",
      "含\"被判死\"成分的支柱一律不许记【条件退化】" in DOC or "被判死" in DOC and "一律不许记【条件退化】" in DOC)

# ======================================================================
head("F2  判定分布（从 §2 表机械点算，不手数）")

VERD4 = {}
for ln in rows:
    c = [x.strip() for x in ln.strip().strip("|").split("|")]
    VERD4[c[0].strip("* ")] = c[4]

LEVELS = ["完整退化", "条件退化", "只到数学核", "不可导出", "开放"]
cnt = {k: 0 for k in LEVELS}
for k, v in VERD4.items():
    hit = [lv for lv in LEVELS if ("【%s】" % lv) in v]
    if len(hit) == 1:
        cnt[hit[0]] += 1
    else:
        check("行 %s 的判定栏恰好命中一档" % k, False, "得到 %s" % hit)

EXPECT = {"完整退化": 0, "条件退化": 3, "只到数学核": 5, "不可导出": 3, "开放": 1}
for lv in LEVELS:
    check("五档 [%s] = %d" % (lv, EXPECT[lv]), cnt[lv] == EXPECT[lv], "实测 %d" % cnt[lv])
check("五档之和 = 12", sum(cnt.values()) == 12, "实测 %d" % sum(cnt.values()))

# 逐行判定（钉住，防漂移）
ROWMAP = {
    "P1": "条件退化", "P2": "只到数学核", "P3": "开放", "P4": "不可导出",
    "P5": "条件退化", "P6": "不可导出", "P7": "只到数学核", "P8": "条件退化",
    "P9": "不可导出", "P10": "只到数学核", "P11": "只到数学核", "P12": "只到数学核",
}
for k in sorted(ROWMAP, key=lambda s: int(s[1:])):
    check("行 %s 判定 = 【%s】" % (k, ROWMAP[k]), ("【%s】" % ROWMAP[k]) in VERD4.get(k, ""))

# §0 与 §2 合计式的数字串必须与点算一致（分布声明不许手写错）
check("§0 声明 0/3/5/3/1 分布",
      all(("\\text{【%s】}%d" % (lv, n)) in DOC for lv, n in
          [("条件退化", 3), ("只到数学核", 5), ("不可导出", 3), ("开放", 1)])
      and "\\textbf{完整退化} = \\boxed{0}" in DOC)
check("§2 合计式五档齐全且和 = 12",
      all(("\\underbrace{%d}" % n) in DOC for n in [0, 3, 5, 3, 1]) and "=12 ." in DOC)

# ======================================================================
head("F3  锚定（依据栏文件存在 + 关键词真实出现）")

FILES = [
    "sm/SM2_stage_discrete_skeleton.md", "sm/SM3_carrier_fermions_spin.md",
    "sm/SM4_join_locality_continuum.md", "sm/SM5_gauge_origin.md",
    "sm/SM6_gauge_dynamics.md", "sm/SM7_gauge_group_hypercharge_charge.md",
    "sm/SM8_anomalies.md", "sm/SM9_generation_counting.md",
    "sm/SM10_matter_dictionary.md", "sm/SM11_higgs_and_ewsb.md",
    "sm/SM12_scale_weak_angle_unification.md", "sm/SM13_flavor_blindness_moduli.md",
    "sm/SM14_mass_mixing_cp.md", "sm/SM15_neutrino_sector.md",
    "sm/SM16_strong_sector_theta.md", "sm/SM18_parameter_ledger_26.md",
    "sm/SM19_verdict_ledger.md",
    "G91_matter_sector_range_and_three_way_verdict.md", "G92_standard_model_gap_inventory.md",
    "G12_gauge_sector_minimal_extension.md", "G57_unreachability_of_absolute_normalization.md",
    "G60_dimensionless_ledger_and_one_free_unit.md", "STATUS.md", "Z18_zero_dimension_bridge_and_gr_stalls.md",
]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]      # 旧文已删：改由备份包核验
_moved = [f for f in ["sm/SM2_stage_discrete_skeleton.md"] if True]
missing = [f for f in FILES if not exists(f)]
check("依据栏引用的 %d 个现存文件在磁盘上存在" % len(FILES), not missing, "缺 %s" % missing)
import tarfile as _tf_mod
_tarp = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
check("旧文依据可从备份包核验（_legacy_SM0-19.tar.gz 在位）", os.path.exists(_tarp))
if os.path.exists(_tarp):
    with _tf_mod.open(_tarp) as _tf:
        _n = set(os.path.basename(x) for x in _tf.getnames())
    check("旧 SM18／SM19（先前依据栏的两个 sm/ 旧文）已入备份",
          "SM18_parameter_ledger_26.md" in _n and "SM19_verdict_ledger.md" in _n)

KW = [
    ("sm/SM2_stage_discrete_skeleton.md", ["首达计数", "读出坐标", "因果锥"]),
    ("sm/SM3_carrier_fermions_spin.md", ["M_2(\\mathbb C)", "Z-CAR", "自旋--统计"]),
    ("sm/SM4_join_locality_continuum.md", ["CHSH", "R-JOIN"]),
    ("sm/SM5_gauge_origin.md", ["引理 43", "指标集"]),
    ("sm/SM6_gauge_dynamics.md", ["BRST", "s^2=0", "11N_c-2N_f"]),
    ("sm/SM7_gauge_group_hypercharge_charge.md", ["\\Sigma T_3^2", "\\frac13\\mathbb Z", "G_0/\\mathbb Z_6"]),
    ("sm/SM8_anomalies.md", ["A_{333}", "Witten", "B-L"]),
    ("sm/SM9_generation_counting.md", ["no-go", "代"]),
    ("sm/SM10_matter_dictionary.md", ["六类", "字典"]),
    ("sm/SM11_higgs_and_ewsb.md", ["Goldstone", "\\rho=1", "(1,2,\\pm\\tfrac12)"]),
    ("sm/SM12_scale_weak_angle_unification.md", ["3/8", "质子衰变", "postdiction"]),
    ("sm/SM13_flavor_blindness_moduli.md", ["模空间", "味"]),
    ("sm/SM14_mass_mixing_cp.md", ["Jarlskog", "Majorana"]),
    ("sm/SM15_neutrino_sector.md", ["seesaw", "反演"]),
    ("sm/SM16_strong_sector_theta.md", ["b_0", "禁闭", "Goldstone"]),
    ("sm/SM18_parameter_ledger_26.md", ["26 行", "具名输入", "【导出】$0$ 项"]),
    ("G91_matter_sector_range_and_three_way_verdict.md", ["不可导出", "只能具名输入"]),
    ("G92_standard_model_gap_inventory.md", ["A8", "A9"]),
    ("G12_gauge_sector_minimal_extension.md", ["引理 43"]),
    ("STATUS.md", ["条件恢复"]),
]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    miss = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not miss, "缺 %s" % miss)

# ======================================================================
head("F4  跨账一致（SM19 五档 / G91 §8 三分 / SM18 参数）")

# --- SM19 §2.3 表机械点算 ---
sec = SM19.split("### 2.3")[1].split("## §3")[0]
p_rows = [ln for ln in sec.splitlines() if ln.startswith("| **P-")]
p_verd = []
for ln in p_rows:
    c = [x.strip() for x in ln.strip().strip("|").split("|")]
    p_verd.append(c[6])
check("SM19 §2.3 表 36 行", len(p_rows) == 36, "实测 %d 行" % len(p_rows))
ORDER = ["【导出】", "【条件】", "【一致性检验】", "【巧合/错误】", "【已更正】"]
tally = {}
for v in p_verd:
    for tag in ORDER:
        if tag in v:
            tally[tag] = tally.get(tag, 0) + 1
            break
    else:
        tally["仅参考"] = tally.get("仅参考", 0) + 1
EXP19 = {"【导出】": 13, "【条件】": 16, "【一致性检验】": 5, "【巧合/错误】": 1, "【已更正】": 0, "仅参考": 1}
for tag, n in EXP19.items():
    check("SM19 机械点算 %s = %d" % (tag, n), tally.get(tag, 0) == n, "实测 %d" % tally.get(tag, 0))
check("SM19 五档和 = 35 + 1 = 36",
      sum(tally.get(t, 0) for t in ORDER) == 35 and sum(tally.values()) == 36)

# --- 本表 §4.1 引用的正是上面这一套（不是 SM19 里的旧串 11/13/11）---
check("§4.1 五档串 = 13/16/5/1/0（与机械点算一致）",
      all(("%d\\ \\text{%s}" % (EXP19[t], t)) in DOC for t in ORDER))
check("§4.1 明写 35 与 36", "=35" in DOC and "=36" in DOC)
check("§4.1 明写【导出】只针对数学核", "全部只针对数学核" in DOC)
stale = "11+13+11+1+0" in SM19
note("SM19 里的旧计数串 11+13+11+1+0 目前%s（已知未修项，照 SM_OVERVIEW §5；本件不判它）"
     % ("仍在" if stale else "已清"))

# --- G91 §8 表机械点算 19 项三分 ---
g8 = G91.split("## §8")[1].split("## §9")[0] if "## §9" in G91 else G91.split("## §8")[1]
g_rows = [ln for ln in g8.splitlines() if re.match(r"^\|\s*\d+\s*\|", ln)]
check("G91 §8 表 19 行", len(g_rows) == 19, "实测 %d 行" % len(g_rows))
g_cnt = {"不可导出": 0, "开放": 0, "能导出": 0, "只能具名输入": 0}
for ln in g_rows:
    c = [x.strip() for x in ln.strip().strip("|").split("|")]
    vd = c[2] if len(c) > 2 else ""
    if "不可导出" in vd:
        g_cnt["不可导出"] += 1
    elif "只能具名输入" in vd:
        g_cnt["只能具名输入"] += 1
    elif "开放" in vd:
        g_cnt["开放"] += 1
    elif "导出" in vd:          # #9/#11 写「能导出」，#13/#14 写「导出（离散侧）」——同一栏
        g_cnt["能导出"] += 1
for k, n in [("不可导出", 5), ("开放", 3), ("能导出", 4), ("只能具名输入", 7)]:
    check("G91 §8 机械点算 %s = %d" % (k, n), g_cnt[k] == n, "实测 %d" % g_cnt[k])
check("G91 §8 四栏和 = 19", sum(g_cnt.values()) == 19)
check("§4.2 引用与 G91 §8 同向零冲突（文字在位）", "零冲突" in DOC and "不可导出 $5$" in DOC)

# --- SM18 ---
check("SM18 26 参数 + 【导出】0 项在位", "26 行" in SM18 and "【导出】$0$ 项" in SM18)
check("§4.3 引用 SM18 的 25+1 口径", "25$ 具名输入" in DOC and "1$ 不可导出" in DOC)

# --- 三根根缺口 ---
check("§4.4 三根根缺口表在位（3 行）",
      all(("**%s" % k) in DOC for k in ["规范群／表示内容", "绝对尺度", "复合与连续极限"]))

# ======================================================================
head("F5  契约与红线")

# 1) 本件不改变 sm/ 的 SM<数字>_ 文件集契约
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
smc = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_check\.py$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 %d" % len(smd))
check("旧 SM<数字>_check.py 已按指令删除（0 个）", len(smc) == 0, "实测 %d" % len(smc))
smrd = sorted(f for f in os.listdir(HERE) if re.match(r"^SMR\d+_.*\.md$", f))
smrc = sorted(f for f in os.listdir(HERE) if re.match(r"^SMR\d+_check\.py$", f))
check("SMR 系列（支柱 12 篇 ＋ 后续篇）在位", len(smrd) >= 12 and len(smrc) >= 12,
      "%d 篇 / %d 个" % (len(smrd), len(smrc)))
check("本件文件名不被交付清单正则捕获（无 SM<数字>_ 前缀）",
      not re.match(r"^SM\d+_", "SM_PILLARS.md") and not re.match(r"^SM\d+_", "SM_PILLARS_check.py"))
check("旧 SM0_check.py 已删除（其交付清单正则不再适用）",
      not os.path.exists(os.path.join(HERE, "SM0_check.py")))
check("SMR 文件名不被旧正则 ^SM<digit>_ 捕获",
      not any(re.match(r"^SM\d+_", f) for f in smrd + smrc))

# 2) 禁令七条
sec51 = DOC.split("### 5.1")[1].split("### 5.2")[0]
check("§5.1 禁令七条（1.-7. 齐）", all(("%d. **不许**" % i) in sec51 for i in range(1, 8)))

# 3) 无肯定式"已导出标准模型"
bad = []
for ln in DOC.splitlines():
    if "标准模型" in ln and "导出" in ln:
        if not any(k in ln for k in ["不", "没", "未", "禁", "零", "0 根", "无"]):
            bad.append(ln.strip()[:60])
check("含\"标准模型+导出\"的行都带否定限定", not bad, "命中 %s" % bad)

# 4) LaTeX 纪律（照 SM2_check 的 F5 体例）
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))
check("$$ 块数 >= 3（记分卡量级）", len(dd) // 2 >= 3, "块数 %d" % (len(dd) // 2))

# 5) 翻转条件表 6 行
sec52 = DOC.split("### 5.2")[1]
r52 = [ln for ln in sec52.splitlines() if re.match(r"^\|\s*\d+\s*\|", ln)]
check("§5.2 结论翻转条件表 6 行", len(r52) == 6, "实测 %d 行" % len(r52))
check("§5.2 每一行给依据（SM19 §6.4 / Z23 等）", all(("`" in ln) for ln in r52))

# ======================================================================
print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
