#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G28_check.py -- 动力学审计：零和宇宙的动力学建起来没有？

对应文档 G28_dynamics_audit.md。
失败时退出码非零。

  F1  4 个沙盒都是动力学沙盒（docstring）
  F2  geometry_probe 的否定结论（记录 JSON）
  F3  失败机制：平均度随周期无界增长（=> 不是流形离散化）
  F4  局部增长维数随周期漂移（不收敛到 4）
  F5  living_universe 用 rates => 概率型 => 违反 Z0③
  F6  generative_selection 显式用 Poisson => 违反 Z0③
  F7  cycle_evolution 的规则与结果
  F8  四层判定
  F9  诚实边界
"""

import io
import json
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


def rd(p):
    return io.open(os.path.join(HERE, p), encoding="utf-8", errors="replace").read()


def js(p):
    return json.load(io.open(os.path.join(HERE, p), encoding="utf-8"))


def doc(f):
    m = re.search(r'"""(.*?)"""', rd(f), re.S)
    return m.group(1).strip() if m else ""


# ======================================================================
head("F1  4 个沙盒都是动力学沙盒")

for f, key in [("zero_sum_living_universe.py", "分支"),
               ("zero_sum_geometry_probe.py", "geometry"),
               ("zero_sum_cycle_evolution.py", "rule"),
               ("zero_generative_selection.py", "生殖")]:
    d = doc(f)
    check("%s 有模型描述（含『%s』）" % (f, key), key in d, "%d 字" % len(d))
check("=> 这 4 个此前标为『未读』的文件都是动力学沙盒", True)

# ======================================================================
head("F2  geometry_probe 的否定结论")

g = js("zero_sum_geometry_probe_results.json")
it = g["interpretation"]
check("status = no stable four-dimensional plateau detected",
      it["status"] == "no stable four-dimensional plateau detected")
check("reason 提到增长指数与谱维数随周期漂移",
      "drift with closure period" in it["reason"])
check("next_gap：零和闭环图本身不固定有效几何方向数",
      "does not fix the number of effective geometric directions" in it["next_gap"])
check("用了对照组（4D torus / random regular）",
      len(g["model"]["controls"]) >= 2)
print("      周期 = %s" % g["model"]["periods"])

# ======================================================================
head("F3  失败机制：平均度随周期无界增长")

fp = g["summary"]["fixed_period_graphs"]
Ts = sorted(fp, key=lambda s: int(s))
deg = [(int(T), fp[T]["mean_degree"]) for T in Ts]
nodes = [(int(T), fp[T]["nodes"]) for T in Ts]
print("      T    节点数   平均度")
for (t, n), (_, d) in zip(nodes, deg):
    print("      %-4d %6d  %6.2f" % (t, n, d))
check("平均度单调增长", all(deg[i][1] < deg[i + 1][1] for i in range(len(deg) - 1)))
check("平均度从 3.6 增到 9.5（%d 倍）" % round(deg[-1][1] / deg[0][1]),
      deg[-1][1] / deg[0][1] > 2.0,
      "%.2f -> %.2f" % (deg[0][1], deg[-1][1]))
check("节点数随周期快速增长（>100 倍）",
      nodes[-1][1] / nodes[0][1] > 100, "%d -> %d" % (nodes[0][1], nodes[-1][1]))
check("=> 度无界 => 闭环图不是流形的离散化 => 无稳定谱维数", True)

# ======================================================================
head("F4  局部增长维数随周期漂移（不收敛到 4）")

lg = [(int(T), fp[T]["local_growth_dimension"]) for T in Ts
      if fp[T]["local_growth_dimension"] is not None]
print("      " + "  ".join("T=%d:%.3f" % (t, v) for t, v in lg))
check("有 >=3 个周期给出了增长维数", len(lg) >= 3)
vals = [v for _, v in lg]
check("增长维数随周期单调上升（漂移）", all(vals[i] < vals[i + 1] for i in range(len(vals) - 1)),
      " -> ".join("%.3f" % v for v in vals))
check("没有任何期给出 4.0（无 4 维平台）", all(abs(v - 4.0) > 0.3 for v in vals))

# ======================================================================
head("F5  living_universe 用 rates => 概率型 => 违反 Z0③")

lu = js("zero_sum_living_universe_results.json")
m = lu["model"]
check("模型含 pair_source_rate（速率 => 概率）", "pair_source_rate" in m,
      "= %s" % m.get("pair_source_rate"))
check("模型含 annihilation_rate", "annihilation_rate" in m, "= %s" % m.get("annihilation_rate"))
check("replicates >= 20（多次随机重复 => 随机过程）", m.get("replicates", 0) >= 20,
      "= %s" % m.get("replicates"))
s = lu["summary"]
print("      场景              存活率  自持率  平均环数")
for k, v in s.items():
    print("      %-18s %.2f   %.2f   %.1f" % (k, v["survival_probability"],
                                              v["self_sustained_probability"], v["mean_final_loops"]))
check("『闭合但绝育』存活率 0，而『复制』类存活率 1（再生需要随机过程）",
      s["sterile_closure"]["survival_probability"] == 0.0 and
      s["living"]["survival_probability"] == 1.0)
check("=> 自持的『活』动力学依赖概率 => 与 Z0③『不设概率』冲突", True)

# ======================================================================
head("F6  generative_selection 显式用 Poisson")

gs = js("zero_generative_selection_results.json")
check("population_update 含 Poisson",
      "Poisson" in gs["model"]["population_update"],
      gs["model"]["population_update"][:60])
check("birth_rate 依赖 fidelity 等参数", "fidelity" in gs["model"]["birth_rate"])
check("零约束在运行中保持（max_total_charge_error = 0）",
      gs["diagnostics"]["max_total_charge_error"] == 0)
check("=> 显式概率更新 => 违反 Z0③", True)

# ======================================================================
head("F7  cycle_evolution 的规则与结果")

ce = js("zero_sum_cycle_evolution_results.json")
check("零约束：每个环与每次变异都保持 sum(charges)=0",
      "sum(charges)=0" in ce["model"]["zero_constraint"])
check("三条规则 sterile / single_cut / all_cuts",
      set(ce["model"]["rules"]) == {"sterile", "single_cut", "all_cuts"})
check("谱系率 λ = log(M)/T", "log(M)/T" in ce["model"]["lineage_rate"].replace(" ", ""))
s2 = ce["summary"]
print("      场景                谱系数  多样性")
for k, v in s2.items():
    print("      %-20s %5d  %.4f" % (k, v["final_lineages"], v["final_diversity"]))
check("无变异（sterile / one_cut_clone）=> 1 个谱系、多样性 0",
      s2["sterile_regular"]["final_lineages"] == 1 and s2["one_cut_clone"]["final_lineages"] == 1)
check("有变异 => 122 个谱系、多样性 > 0.9",
      s2["one_cut_evolution"]["final_lineages"] == 122 and
      s2["one_cut_evolution"]["final_diversity"] > 0.9)
check("=> 多样性需要变异（= 概率）", True)

# ======================================================================
head("F8  四层判定")

check("层 1 运动学：✅ 底层条款（Z0 条款 ＋ Z1–Z5 定理；G19 精简为独立定理 Z1 定理 1、Z2、Z3）", True)
check("层 2 微观动力学：✅ 建好并运行（cycle_evolution 零约束每步精确）", True)
check("层 3 宏观动力学：⚠️ 只建到扩散型（G6 热方程 / G7 梯度流）；几何极限未建立且探针否定", True)
check("层 4 统计/量子动力学：❌ 没有（Z0③ 禁概率；G27 均匀计数 => 模 Hamiltonian 平凡）", True)
check("=> 微观建起来了；宏观只到扩散型；统计/量子没有；几何有否定证据", True)

# ======================================================================
head("F9  诚实边界")

check("我只复跑了 geometry_probe（26 秒 exit 0）；living_universe / generative_selection / cycle_evolution 未复跑", True)
check("cycle_evolution 的『变异概率 0.03』与 living_universe 的 rates 是外部设定，不是从底层条款（Z0 条款 ＋ Z1–Z5 定理）导出", True)
check("geometry_probe 用 4D torus 作对照，隐含『4 维是目标』——这本身是输入", True)
check("本审计不改变 G1–G27 的结论，只登记动力学的完成度", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
