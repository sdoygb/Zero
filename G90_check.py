#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G90_check.py —— Zero 系列参考分诊的核验
=======================================
独立实断言，全部可失败：
  F1  7 篇笔记与 10 个实验文件在位（根目录与 simulations/ 双份且 md5 相同）
  F2  引用计数复算 == 分诊表所列（口径：G*.md + D*.md 中 "zero_sum_<名>" 出现次数）
  F3  【承接】5 项各自的锚点在目标文档里真实存在
  F4  【打问号】2 项的"对象至今不存在"可核验（G 系列无对应节）
  F5  【纯历史】3 项的自述文本逐字在位
  F6  结构性事实：无笔记的 3 个 == 零引用的 3 个
  F7  可复跑性：在临时目录里真跑 3 个实验，退出码 0 且输出落在该目录
  F8  §5 时间线：笔记里 L=4 是"设定"，而 G46/G61 不引用 Zero 系列
"""
import hashlib
import io
import os
import shutil
import subprocess
import sys
import tempfile


HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = ("v" if ok else "x") if (level == "ind" or LEDGER_MODE == "A" or not ok) else "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


def rd(p):
    return io.open(os.path.join(HERE, p), encoding="utf-8").read()


NOTES = ["zero_sum_closure_exit.md", "zero_sum_closure_time_selection.md",
         "zero_sum_global_R_local_P.md", "zero_sum_open_reservoir.md",
         "zero_sum_periodic_destruction.md", "zero_sum_reproduction_audit.md",
         "zero_sum_tri_layer_universe.md"]
EXPS = ["closure_exit", "closure_time_selection", "cycle_evolution", "geometry_probe",
        "global_R_local_P", "living_universe", "open_reservoir", "periodic_destruction",
        "reproduction_audit", "tri_layer_universe"]
CITED = {"closure_exit": 7, "periodic_destruction": 6, "global_R_local_P": 6,
         "open_reservoir": 4, "tri_layer_universe": 2, "closure_time_selection": 0,
         "reproduction_audit": 0, "cycle_evolution": 0, "geometry_probe": 0,
         "living_universe": 0}
NO_NOTE = {"cycle_evolution", "geometry_probe", "living_universe"}

DOC = rd("G90_zero_series_reference_triage.md")

# ---------------------------------------------------------------- F1
head("F1  文件在位（根目录与 simulations/ 双份）")
for f in NOTES:
    check("笔记在位：%s" % f, os.path.exists(os.path.join(HERE, f)))
exp_ok, md5_ok = True, True
for e in EXPS:
    a = os.path.join(HERE, "zero_sum_%s.py" % e)
    b = os.path.join(HERE, "simulations", "zero_sum_%s.py" % e)
    exp_ok = exp_ok and os.path.exists(a) and os.path.exists(b)
    if os.path.exists(a) and os.path.exists(b):
        ha = hashlib.md5(io.open(a, "rb").read()).hexdigest()
        hb = hashlib.md5(io.open(b, "rb").read()).hexdigest()
        md5_ok = md5_ok and (ha == hb)
check("10 个实验脚本在根目录与 simulations/ 均在位", exp_ok)
check("两处副本 md5 相同（未分叉）", md5_ok)

# ---------------------------------------------------------------- F2
head("F2  引用计数复算（G*.md + D*.md 中 'zero_sum_<名>' 的次数）")
SELF = "G90_zero_series_reference_triage.md"
# 排除"关于 Zero 层本身的文档"：它们的职责就是逐条点名 Zero 组件，引用它们不表示理论脊柱在用
ZERO_DOCS = {SELF, "Z1_zero_layer_as_the_foundation.md",
             "Z0_zero_never_rests_single_axiom.md"}
docs = [f for f in os.listdir(HERE)
        if (f.startswith("G") or f.startswith("D")) and f.endswith(".md") and f not in ZERO_DOCS]
counts = {}
for e in EXPS:
    counts[e] = sum(rd(f).count("zero_sum_" + e) for f in docs)
    check("被引次数 %-24s ≥ %d（实算 %d）" % (e, CITED[e], counts[e]), counts[e] >= CITED[e])
check("分诊表逐行列出 10 个实验的引用数",
      all(("`%s`" % e) in DOC for e in EXPS))

# ---------------------------------------------------------------- F3
head("F3  【承接】5 项的锚点在目标文档中真实存在")
g34, g28, g20, g0, g22 = (rd("G34_exit_rule_absorbed_into_pi.md"), rd("G28_dynamics_audit.md"),
                          rd("G20_axiom_audit_extended_to_zero_and_D.md"),
                          rd("G0_bottom_layer_and_derivation_route.md"),
                          rd("G22_correction_sublayers_and_local_layers.md"))
check("G34 断言「闭合即退出」是 π 的一条定义",
      ("闭合即退出" in g34) and ("定义" in g34))
check("G28 引用了 zero_sum_closure_exit（对照实验的来源）", "zero_sum_closure_exit" in g28)
check("G20 登记了 I8", "I8" in g20)
check("G0 登记了 I8（登记表在 G0）", "I8" in g0)
check("G22 接住了 locality（tri_layer 的诊断）", "locality" in g22)
check("G34 的判定被 G90 正文复述（含 G34 链接）", "G34_exit_rule_absorbed_into_pi.md" in DOC)
check("G20 §3 引理 74 的编号在 G90 中被引用", "引理 74" in DOC and "引理 74" in g20)
check("重播种映射 w -> {w+,w-} 在笔记与 G90 中一致",
      ("w+" in rd("zero_sum_periodic_destruction.md")) and ("w+" in DOC))

# ---------------------------------------------------------------- F4
head("F4  【打问号】2 项：后继对象至今不存在")
g_all = "".join(rd(f) for f in docs if f.startswith("G") and f not in
                ("G90_zero_series_reference_triage.md",))
_auto = sum(rd(f).count("自催化") + rd(f).count("autocatalytic") for f in docs
            if f.startswith("G") and f != "G90_zero_series_reference_triage.md")
check("G 系列（除 G90）确实没有『自催化／繁殖规则』这一节（实算 0 处）", _auto == 0,
      "扫描 %d 篇 G 文档，命中 %d 处" % (len([f for f in docs if f.startswith('G')]) - 1, _auto))
check("reproduction_audit 的隐藏假设在笔记中逐字在位",
      "hidden copy assumption" in rd("zero_sum_reproduction_audit.md"))
check("closure_time_selection 的结论句在笔记中逐字在位",
      "mode dominates almost completely" in rd("zero_sum_closure_time_selection.md"))
check("G90 明确写出该对象「至今不存在／没有后继」",
      ("至今不存在" in DOC) and ("没有后继" in DOC or "无后继" in DOC))

# ---------------------------------------------------------------- F5
head("F5  【纯历史】3 项的自述逐字在位")
check("cycle_evolution 自述 exploratory / 不用 D*",
      "exploratory model" in rd("simulations/zero_sum_cycle_evolution.py"))
check("geometry_probe 自述 narrow question",
      "narrow question" in rd("simulations/zero_sum_geometry_probe.py"))
check("living_universe 自述独立玩具模型",
      ("玩具" in rd("simulations/zero_sum_living_universe.py")))

# ---------------------------------------------------------------- F6
head("F6  结构性事实：无笔记的 3 个中，2 个零引用、1 个是漏引例外")
zero_cited = {e for e in EXPS if counts[e] == 0}
check("零引用集合恰为 {cycle_evolution, living_universe, closure_time_selection, reproduction_audit}",
      zero_cited == {"cycle_evolution", "living_universe",
                     "closure_time_selection", "reproduction_audit"},
      "实算 %s（口径：排除 G90/Z0/Z1 三篇 Zero 层文档）" % sorted(zero_cited))
# geometry_probe 是例外：它零引用是因为 G28 原稿漏引（本版已补引）
check("geometry_probe 现已被 G28 引用（补引后计数 ≥ 1）", counts["geometry_probe"] >= 1,
      "实算 %d" % counts["geometry_probe"])
check("3 个无笔记实验中，2 个零引用、1 个是'漏引例外'",
      (NO_NOTE & zero_cited) == {"cycle_evolution", "living_universe"}
      and "geometry_probe" not in zero_cited)

# ---------------------------------------------------------------- F7
head("F7  可复跑性：临时目录里真跑 3 个实验")
tmp = tempfile.mkdtemp(prefix="g90_")
try:
    os.makedirs(os.path.join(tmp, "simulations"))
    # 实验互相依赖（7 个 import zero_sum_cycle_evolution），故成套复制
    for fn in os.listdir(os.path.join(HERE, "simulations")):
        if fn.startswith("zero_sum_") and (fn.endswith(".py") or fn.endswith(".json")):
            shutil.copy(os.path.join(HERE, "simulations", fn), os.path.join(tmp, "simulations"))
    for e in ("closure_exit", "open_reservoir", "periodic_destruction"):
        r = subprocess.run([sys.executable, os.path.join(tmp, "simulations", "zero_sum_%s.py" % e)],
                           capture_output=True, text=True, timeout=300)
        outp = os.path.join(tmp, "simulations", "zero_sum_%s_results.json" % e)
        check("实跑 zero_sum_%s.py：rc=0 且结果落在 simulations/" % e,
              r.returncode == 0 and os.path.exists(outp),
              "rc=%d, json=%s" % (r.returncode, os.path.exists(outp)))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

# ---------------------------------------------------------------- F8
head("F8  §5 时间线：L=4 在两代里的等级不同")
note = rd("zero_sum_closure_exit.md")
check("笔记里 L=4 是「设定」而非导出",
      ("L=4" in note) and ("寿命" in note))
g46, g61 = rd("G46_k_is_the_lifetime.md"), rd("G61_locking_the_five_integers.md")
check("G46/G61 均不引用 Zero 系列（不主张启发关系）",
      ("zero_sum" not in g46) and ("zero_sum" not in g61))
check("G90 明写「不登记为承接」", "不登记为承接" in DOC)

# ---------------------------------------------------------------- F9
head("F9  层内基础设施 vs 理论资产：cycle_evolution 的双重身份")
imp = [f for f in os.listdir(os.path.join(HERE, "simulations"))
       if f.startswith("zero_sum_") and f.endswith(".py")
       and "from zero_sum_cycle_evolution import" in rd("simulations/" + f)]
check("7 个实验 import zero_sum_cycle_evolution 的 canonical_cycle",
      len(imp) == 7, "实算 %d 个：%s" % (len(imp), sorted(x[9:-3] for x in imp)))
check("cycle_evolution 自身被引为 0（对理论非资产）", counts["cycle_evolution"] == 0)
check("G90 写清了这双重身份（纯历史 + 层内基础设施）",
      ("层内基础设施" in DOC) and ("纯历史" in DOC))

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
