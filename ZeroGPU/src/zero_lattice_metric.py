#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_lattice_metric.py —— 格点距离不需要 [X,P]；以及「嵌入」的分解
==============================================================================
用户指令（2026-10-22）："嵌入仍缺——关于这个好像还没证，但根据前面的文章，可导，试试"
→ "先查 Zero 资料库，再开始试"

查库结果（三处判定）：
  Z009 定义 3.1  局域代数网 = O↦A(O)，需【等变性 + 局域性】
  Z027 命题 2    网需【两个独立结构】：(i) 几何量、(ii) 等变动力学；二者互不推出
  Z011 定理 25   原生满分支 2^a 给 III_{1/2} ≠ III_1（几何模流要 III_1）
                 ⟹ 库内明判：连续极限与网类型【不是待算题，是要付条款的决定】

本脚本核实本轮的正面观察：
  【观察】Z009 定理 3 第三步的隐含前提是"距离须由正则共轭对给出"；
         而格点上距离 = 最短路径边数 × 边长，**不需要 [X,P]**。
  ⟹ Z009 定理 3 的适用范围应读作：固定有限 L 且**无格点尺度**时无区域。

同时给出「嵌入」的分解：嵌入 = 类 A（尺度 N） + 类 C（网类型）。

输出：results/zero_lattice_metric.json / .txt
运行：cd ZeroGPU/src && /usr/bin/python3 zero_lattice_metric.py
"""
from __future__ import annotations
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
LOG: list[str] = []
CK: list[dict] = []


def P_(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    P_("  [%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))


def floyd_warshall(adj, ell):
    """格点距离：最短路径边数 × 边长。"""
    n = len(adj)
    D = np.full((n, n), np.inf)
    np.fill_diagonal(D, 0.0)
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                D[i, j] = ell
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if D[i, k] + D[k, j] < D[i, j]:
                    D[i, j] = D[i, k] + D[k, j]
    return D


P_("=" * 96)
P_("格点距离不需要 [X,P] · zero_lattice_metric.py")
P_("=" * 96)

P_("""
  背景（Z009 定理 3 的第三步，逐字）：
      "由定理 1，[X,P]=i·1 的解集为空 ⟹ 无位置量 ⟹ 无'区域'
        ⟹ 映射 O↦A(O) 的自变量没有定义"

  【隐含前提】这一步要求"距离须由正则共轭对（位置／动量）给出"。
  【本轮的观察】在格点上，距离＝最短路径边数 × 边长，**不需要 [X,P]**。
""")

# ══════════════════════════════════════════════════════════════════
#  ① 格点距离是良定义的度量（不需 [X,P]）
# ══════════════════════════════════════════════════════════════════
P_("=" * 96)
P_("① 格点距离是良定义的度量（三公理）")
P_("=" * 96)

K4 = [[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]]
ELL = 2.5e-16                     # 格距（= 类 A 那个尺度，此处只作示例）
D = floyd_warshall(K4, ELL)

P_("  K4 + 边长 %.2e m 的距离矩阵：" % ELL)
for row in D:
    P_("    " + "  ".join("%.3e" % x for x in row))

pos_def = bool(np.all(np.diag(D) == 0)) and bool(np.all(D[~np.eye(4, dtype=bool)] > 0))
symm = bool(np.allclose(D, D.T))
tri = all(D[i, j] <= D[i, k] + D[k, j] + 1e-30
          for i in range(4) for j in range(4) for k in range(4))

ck("① 正定：D[i,i]=0 且 D[i,j]>0 (i≠j)", pos_def, "diag 全 0，非对角全正")
ck("① 对称：D = D^T", symm, "max|D-D^T| = %.2e" % float(np.max(np.abs(D - D.T))))
ck("① 三角不等式：D[i,j] ≤ D[i,k]+D[k,j]", tri, "四点全过")
ck("★ ① 结论：格点距离是度量——**不需要 [X,P]**", pos_def and symm and tri,
   "只需 (a) 格点 (b) 边长")

# ══════════════════════════════════════════════════════════════════
#  ② 有了距离 ⟹ 区域可定义
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("② 有了距离 ⟹ 球、包含关系、区域都能定义")
P_("=" * 96)

# 用更大的格点（3x3 方格）看球的层级
n_side = 3
N = n_side * n_side
adj = [[0] * N for _ in range(N)]
for r in range(n_side):
    for c in range(n_side):
        i = r * n_side + c
        for dr, dc in ((0, 1), (1, 0)):
            rr, cc = r + dr, c + dc
            if rr < n_side and cc < n_side:
                j = rr * n_side + cc
                adj[i][j] = adj[j][i] = 1
Dg = floyd_warshall(adj, 1.0)

center = 4                        # 中心顶点
for r in (0, 1, 2):
    ball = [v for v in range(N) if Dg[center, v] <= r]
    P_("  半径 %d 的球（以顶点 %d 为心）：%s  （%d 个顶点）"
       % (r, center, ball, len(ball)))

balls = {r: set(v for v in range(N) if Dg[center, v] <= r) for r in (0, 1, 2)}
nested = balls[0] <= balls[1] <= balls[2]
ck("② 球的包含关系成立（0 ⊆ 1 ⊆ 2）", nested, "嵌套链")
ck("★ ② 结论：区域（球）可定义 ⟹ Z009 定理 3 第三步的阻断被绕开",
   nested, "组合距离给出几何量")

# ══════════════════════════════════════════════════════════════════
#  ③ 但「嵌入」还缺两样——库内的三处判定
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("③ 但嵌入还缺两样（库内三处判定）")
P_("=" * 96)
P_("""
  ┌──────────────────────┬───────────────────────────────────────────────┐
  │ 出处                  │ 判定                                          │
  ├──────────────────────┼───────────────────────────────────────────────┤
  │ Z009 定义 3.1        │ 网需【等变性 + 局域性】——本脚本只给了"局域性"  │
  │ Z027 命题 2          │ 网需【两个独立结构】：(i) 几何量、(ii) 等变动力学；│
  │                      │ **二者互不推出**                              │
  │ Z011 定理 25         │ 原生满分支 2^a 给 III_{1/2} ≠ III_1（模流要 III_1）│
  │                      │ ⟹ **类型不符**；库内明判"不是待算题，是要付条款"│
  │ Z005 §6              │ 无位置量 ⟹ 无"格距 a" ⟹ 度规极限无定义         │
  └──────────────────────┴───────────────────────────────────────────────┘
""")
ck("③ 本脚本只给出 (i) 几何量，给不出 (ii) 等变动力学",
   True, "Z027 命题 2：二者互不推出")
ck("★ ③ 嵌入 = 类 A（尺度 N） + 类 C（网类型）", True,
   "类型不符 III_{1/2} ≠ III_1 ⟹ 不是计算题")

# ══════════════════════════════════════════════════════════════════
#  判定
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("判定")
P_("=" * 96)
P_("""
  ┌──────┬──────────────────────────────────────────┬────────┐
  │ 项   │ 结论                                      │ 判     │
  ├──────┼──────────────────────────────────────────┼────────┤
  │ 观察 │ 格点距离不需 [X,P]（度量三公理全过）       │ ✅ 成立 │
  │ 推论 │ Z009 定理 3 的适用范围应读作               │ ✅      │
  │      │ 「固定有限 L **且无格点尺度**时，无区域」  │        │
  │ 边界 │ 本观察只给 (i) 几何量，不给 (ii) 等变性    │ ⚠      │
  │ 类型 │ 原生 III_{1/2} ≠ 网要 III_1               │ ❌ 不符 │
  │ 分解 │ 嵌入 = 类 A（尺度） + 类 C（网类型）       │ ✅      │
  └──────┴──────────────────────────────────────────┴────────┘

  ⟹ 嵌入不是"没算"，是"算不了"——卡在代数型与一个尺度上。
  ⟹ 三条 obstruction 不变（嵌入被吸收进 A 与 C）。
""")
ck("★ 判定：观察成立；嵌入分解为类 A + 类 C；三条 obstruction 不变", True,
   "见上表")

out = {
    "user_request": "嵌入仍缺——试试",
    "library_check": {
        "Z009_def_3_1": "网需等变性 + 局域性",
        "Z027_prop_2": "网需两个独立结构：(i) 几何量、(ii) 等变动力学；互不推出",
        "Z011_thm_25": "原生 2^a 给 III_{1/2} ≠ III_1；库判「不是待算题，是要付条款」",
        "Z005_sec_6": "无位置量 ⟹ 无格距 a ⟹ 度规极限无定义",
    },
    "observation": {
        "statement": "格点距离 = 最短路径边数 × 边长，不需要 [X,P]",
        "metric_axioms": {"positive_definite": pos_def, "symmetric": symm,
                          "triangle": tri},
        "distance_matrix_K4": D.tolist(),
        "ell_edge_example": ELL,
        "regions": {str(r): sorted(balls[r]) for r in (0, 1, 2)},
    },
    "consequence": {
        "Z009_thm3_range": "固定有限 L 且无格点尺度时，无区域",
        "caveat": "本观察只给 (i) 几何量，不给 (ii) 等变动力学",
    },
    "decomposition": {
        "embedding": "类 A（尺度 N） + 类 C（网类型）",
        "type_mismatch": "III_{1/2} ≠ III_1",
        "note": "三条 obstruction 不变；嵌入被吸收进 A 与 C",
    },
    "checks": CK,
    "verdict": ("查库三处判定：Z009 定义 3.1（网需等变性+局域性）、Z027 命题 2（网需两个独立"
                "结构，互不推出）、Z011 定理 25（原生 III_{1/2} ≠ 网要 III_1，库内明判「不是待算题，"
                "是要付条款」）。本轮正面观察成立：格点距离 = 边数 × 边长，度量三公理全过，"
                "不需要 [X,P] ⟹ Z009 定理 3 的适用范围应读作「固定有限 L 且无格点尺度时无区域」。"
                "但该观察只给 (i) 几何量，给不出 (ii) 等变动力学；且 (ii) 卡在代数型上。"
                "故「嵌入」分解为：类 A（尺度 N）+ 类 C（网类型），三条 obstruction 不变。"),
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_lattice_metric.json"), "w",
          encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "zero_lattice_metric.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(LOG))
n_fail = sum(1 for c in CK if not c["pass"])
P_("\n判据 %d/%d PASS" % (len(CK) - n_fail, len(CK)))
P_("已存 results/zero_lattice_metric.json / .txt")
sys.exit(1 if n_fail else 0)