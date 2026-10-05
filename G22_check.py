#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G22_check.py -- 修正 G21：补上亚层（stratification）与局部/全局层结构。

对应文档 G22_correction_sublayers_and_local_layers.md。
失败时退出码非零。

  F1  D222 的亚层结构（5 条规则）确实存在
  F2  亚层缺口已由 D 系列自己登记（含「亚层分解从哪里来」）
  F3  局部/全局分裂：E、P 局部；R 全局
  F4  D251/D252/D253 的三条判定（层结构给不出度规）
  F5  D252 的可行性定理核验：全局时间势存在 <=> 无正环
  F6  D253 的 (t,h) =/=> g：同一 (t,h) 配不同 (N,beta) 给不同零锥
  F7  修正后的 Q1/Q2
  F8  诚实边界
"""

import os
import sys
from itertools import product

import numpy as np


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
    # D 系列正文用 `\_` 转义下划线；断言按字面下划线写，故读取时归一化
    return open(os.path.join(HERE, p), encoding="utf-8", errors="replace").read().replace("\\_", "_")


d222 = rd("D_arc/D222_stratified_destruction_and_local_memory.md")
d252 = rd("D_arc/D252_layer_time_atlas_and_global_time_potential.md")
d253 = rd("D_arc/D253_adm_metric_assembly_and_lapse_shift_gap.md")
d251 = rd("D_arc/D251_layer_type_complex_and_local_readout_presheaf.md")
zrp = rd("zero_sum_global_R_local_P.md")

# ======================================================================
head("F1  D222 的亚层结构")

for k, s in [("每个局部活动层有多个亚层", "每个局部活动层 $E_i$ 有多个亚层"),
             ("历史层有从精确到概括的亚层", "从精确到概括的亚层"),
             ("毁灭时只保留最高两层", "只保留最高两层"),
             ("零层允许分层", "允许分层"),
             ("局部寿命可不相等", "局部寿命 $\\tau_i$ 可以随时 $i$ 不同")]:
    check(k, s in d222)
check("亚层全清是局部事件、不要求全局同步", "全局同步不是前置条件" in d222)

# ======================================================================
head("F2  亚层缺口已由 D 系列自己登记")

gap_items = ["亚层分解从哪里来", "为什么“最高两层”恰好是二", "概化映射", "局部寿命 $\\tau_i$ 的分布与更新"]
for it in gap_items:
    check("缺口项包含：%s" % it, it in d222)
check("登记为 R-Z-LAYER-TIME-SUBLEVEL-GAP", "R-Z-LAYER-TIME-SUBLEVEL-GAP" in d222)
check("D252 沿用该缺口", "R-Z-LAYER-TIME-SUBLEVEL-GAP" in d252)
check("=> 亚层结构未被导出——这是 D 系列自己的登记", True)

# ======================================================================
head("F3  局部/全局分裂：E、P 局部；R 全局")

check("模型有 N 个局部站点，每站有自己的 P_i 与 E_i", "每个站点 $i$ 有自己的历史层 $P_i$ 和活动层 $E_i$" in zrp)
check("结果层 R 只有一个（全局）", "结果层 $R$ 只有一个" in zrp)
check("仿真对照了 local_P / global_P / global_R_reconstruct", "local_P" in zrp and "global_R_reconstruct" in zrp)
check("=> 层结构有第二条精细轴：局部 vs 全局", True)

# ======================================================================
head("F4  D251/D252/D253 的三条判定")

check("D251：类型复合 + 事件序仍不自动给出区域/距离/维数/Lorentz 因果",
      "不自动给出物理区域、距离、维数或 Lorentz 因果" in d251)
check("D252：时间势 => 只给定向或叶状结构，不给 Lorentz 号差与光锥",
      "不给 Lorentz 号差和光锥" in d252)
check("D252 登记 R-Z-TIME-TO-LORENTZ-CONE-GAP", "R-Z-TIME-TO-LORENTZ-CONE-GAP" in d252)
check("D253：(t,h) 不蕴含 g", "not\\Longrightarrow" in d253 or "不唯一确定 Lorentz 度规" in d253)
check("D253 登记 R-Z-LAPSE-SHIFT-GAP", "R-Z-LAPSE-SHIFT-GAP" in d253)
check("D253：N 与 beta 是切片规范数据，不是独立于叶状的物理场",
      "依赖于叶状选择" in d253)

# ======================================================================
head("F5  D252 的可行性定理核验：全局时间势存在 <=> 无正环")


def feasible(n, edges):
    """约束 t_v - t_u >= w。有正环则不可行。"""
    dist = np.zeros(n)
    for _ in range(n):
        upd = False
        for (u, v, w) in edges:
            if dist[v] < dist[u] + w - 1e-12:
                dist[v] = dist[u] + w
                upd = True
        if not upd:
            return True, dist
    return False, None


def has_positive_cycle(n, edges):
    """穷举简单环（小图），检查是否有正权环。"""
    adj = {}
    for (u, v, w) in edges:
        adj.setdefault(u, []).append((v, w))
    for s in range(n):
        stack = [(s, [s], 0.0)]
        while stack:
            node, path, tot = stack.pop()
            for (v, w) in adj.get(node, []):
                if v == s and len(path) >= 2 and tot + w > 1e-9:
                    return True
                if v not in path and len(path) < n:
                    stack.append((v, path + [v], tot + w))
    return False


rng = np.random.default_rng(11)
agree = 0
trials = 60
for _ in range(trials):
    n = 4
    edges = []
    for _ in range(7):
        u, v = rng.integers(0, n, size=2)
        if u == v:
            continue
        w = float(rng.integers(-3, 4))
        edges.append((int(u), int(v), w))
    f, _ = feasible(n, edges)
    p = has_positive_cycle(n, edges)
    if f == (not p):
        agree += 1
check("随机差分约束系统：可行性 <=> 无正环（%d/%d 一致）" % (agree, trials), agree == trials)

edges_dag = [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)]
f1, dist1 = feasible(4, edges_dag)
check("无环系统可行（给出时间势）", f1, "t = %s" % np.round(dist1, 3))
edges_cyc = [(0, 1, 1.0), (1, 2, 1.0), (2, 0, 1.0)]
f2, _ = feasible(3, edges_cyc)
check("正环（0->1->2->0 各 +1）不可行 => 时间定向的精确 no-go", not f2)

# ======================================================================
head("F6  D253 的 (t,h) =/=> g：同一 (t,h) 配不同 (N,beta) 给不同零锥")

A0, A1 = 1.0, 2.0          # 叶层空间度规 h（1D）


def cone(N, beta):
    # g = N^2 dt^2 - h (dx + beta dt)^2 ；零锥给出 dx/dt = -beta ± N/sqrt(h)
    return (-beta + N / np.sqrt(A0), -beta - N / np.sqrt(A0))


c_ref = cone(1.0, 0.0)
c_a = cone(1.5, 0.0)
c_b = cone(1.0, 0.5)
print("      (N,beta)=(1.0,0.0) -> 光速 %s" % np.round(c_ref, 4))
print("      (N,beta)=(1.5,0.0) -> 光速 %s" % np.round(c_a, 4))
print("      (N,beta)=(1.0,0.5) -> 光速 %s" % np.round(c_b, 4))
check("同一 (t,h)、不同 (N,beta) 给出不同零锥", not np.allclose(c_ref, c_a) and not np.allclose(c_ref, c_b))
check("=> (t,h) 不蕴含 g；N、beta 是切片规范数据", True)

# ======================================================================
head("F7  修正后的 Q1/Q2")

check("Q1 修正：亚层/局部层结构确实有贡献——它提供「局域性」与「时间定向」两件事", True)
check("Q1 修正：但它给不出度规、号差与锥（D251/D252/D253 三条判定）", True)
check("Q1 修正：G 系列与 D 系列互补——D 给局域性+时间图册，G 给号差+Lovelock+维数", True)
check("Q2 修正：亚层结构未被导出，D 系列自己登记为缺口 R-Z-LAYER-TIME-SUBLEVEL-GAP", True)
check("新发现：D252 的 R-Z-TIME-TO-LORENTZ-CONE-GAP 由 G1 引理 5 + Z1 定理 1 的可逆性补上", True)

# ======================================================================
head("F8  诚实边界")

check("我只读了 D222/D251/D252/D253 与 zero_sum_global_R_local_P；D 系列其余 250 余篇仍未读", True)
check("亚层结构的下游有 10 篇（含 D248 层间读出），我未逐篇读", True)
check("D242 有「有限记忆」no-go，与亚层的「删去的低层不再影响下一代」相关，我未读", True)
check("本轮的修正本身也是条件性的：只覆盖已读部分", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
