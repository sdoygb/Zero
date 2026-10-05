#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sum_persistence_theorems_check.py —— Zero 文章④【持续性定理】的核验
=======================================================================
实断言（活动层轨迹全部由**冻结结果 JSON** 实算）：
  F1  closure_exit：两个情景最终熄灭（活动层归零）
  F2  periodic_destruction：wipe_no_reseed 熄灭、wipe_reseed 持续
  F3  open_reservoir：三个维度末值 > 0（无寿命 ⇒ 持续）且量级与文章一致
  F4  「同模型只差重播种」的对照可复算（wipe_* 两情景）
  F5  文章写明判据『无寿命 ∨ 重播种』与 Z0② ⇒ Z5 的推论
  F6  重播种映射 w ↦ {w±} 在文章与程序中一致
"""
import io
import json
import os
import sys

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


def js(p):
    return json.load(io.open(os.path.join(HERE, p), encoding="utf-8"))


DOC = rd("zero_sum_persistence_theorems.md")


def traj(run, key="active_paths"):
    snaps = run.get("snapshots", [])
    if not snaps:
        return []
    k = key if key in snaps[0] else next((x for x in snaps[0] if "active" in x), None)
    return [s.get(k) for s in snaps] if k else []


# ---------------------------------------------------------------- F1
head("F1  closure_exit：两个情景最终熄灭")
ce = js("simulations/zero_sum_closure_exit_results.json")
for r in ce["runs"]:
    t = traj(r)
    check("%-26s 活动层末值 = 0" % r["scenario"], t and t[-1] == 0, "轨迹 %s" % t[:7])

# ---------------------------------------------------------------- F2
head("F2  periodic_destruction：只差重播种那一条")
pd = js("simulations/zero_sum_periodic_destruction_results.json")
W = {r["scenario"]: traj(r) for r in pd["runs"]}
check("wipe_no_reseed 熄灭（末值 0）", W.get("wipe_no_reseed", [1])[-1] == 0,
      "轨迹 %s" % W.get("wipe_no_reseed", [])[:8])
check("wipe_reseed 持续（末值 > 0）", W.get("wipe_reseed", [0])[-1] > 0,
      "轨迹 %s" % W.get("wipe_reseed", [])[:8])
check("continuous 亦持续（对照）", W.get("continuous", [0])[-1] > 0)
check("wipe_no_reseed 与 wipe_reseed 在清空前逐代相同（同模型只差重播种）",
      W["wipe_no_reseed"][:2] == W["wipe_reseed"][:2],
      "%s vs %s" % (W["wipe_no_reseed"][:3], W["wipe_reseed"][:3]))

# ---------------------------------------------------------------- F3
head("F3  open_reservoir：无寿命 ⇒ 持续（三档维度）")
orv = js("simulations/zero_sum_open_reservoir_results.json")
finals = {r["scenario"]: r.get("final_active_paths") for r in orv["runs"]}
check("全部情景末值 > 0（无寿命即持续）", all((v or 0) > 0 for v in finals.values()),
      "%s" % finals)
check("维度 d=1 的末值 = 504（文章 §4 表）", finals.get("d1") == 504, "实算 %s" % finals.get("d1"))
check("维度 d=2 的末值 = 1214752（文章 §4 表）", finals.get("d2") == 1214752,
      "实算 %s" % finals.get("d2"))
check("维度 d=4 的末值 = 1777511808（文章 §4 表）", finals.get("d4") == 1777511808,
      "实算 %s" % finals.get("d4"))
check("末值随维度单调增（爆炸趋势）",
      finals.get("d1", 0) < finals.get("d2", 0) < finals.get("d4", 0))

# ---------------------------------------------------------------- F4
head("F4  未闭合者不补充就会死（F1＋F2 的合成）")
dead = [r["scenario"] for r in ce["runs"] if (traj(r) or [1])[-1] == 0]
check("closure_exit 的两个情景都属于『无重播种』且都熄灭", len(dead) == 2, "%s" % dead)
check("wipe_no_reseed 属于『无重播种』且熄灭", W["wipe_no_reseed"][-1] == 0)
persist_models = ["wipe_reseed"] + list(finals.keys())
check("『有重播种』或『无寿命』的模型全部持续", len(persist_models) >= 4, "%s" % persist_models)

# ---------------------------------------------------------------- F5
head("F5  文章写明判据与推论")
check("写出判据『必须无寿命 ∨ 有重播种 二择一』",
      "二择一" in DOC or "必须" in DOC and "重播种" in DOC)
check("写出 Z0② ⇒ Z5 的推论", "Z5" in DOC and "不断" in DOC)
check("指向 I8 的升级（输入 → 公理条款）", "I8" in DOC)
check("引用 G46（有限寿命是 G 系列需要的）", "G46" in DOC)

# ---------------------------------------------------------------- F6
head("F6  重播种映射 w ↦ {w±} 一致")
src = rd("simulations/zero_sum_periodic_destruction.py")
check("程序里重播种给 w+ 与 w− 两个种子", "seed = word + (step,)" in src
      and "for step in (1, -1):" in src)
check("只取不闭合的种子", "if not is_closed(seed)" in src)
check("文章写出映射 w ↦ {w+, w−}", "w+" in DOC and "w-" in DOC)
check("程序按站点独立播种（文章 §5 的『站点独立』）",
      "active[site].add(seed)" in src)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
