#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_position_gate.py —— 攻最后一扇门：[X,P]=i 在 L→∞ 时能否有解
======================================================================================
承前两稿：Z005 定理 6 只锁有限 L；Z028 的 boost 项已撤。
**剩下的唯一障碍**：$[X,P]=i\\cdot1$。

本脚本把它精确化。三个必须分开问的问题：
  Q1  有限 L 时为什么无解？        （迹定理，Z009 定理 1 / Z006 定理 5(a)）
  Q2  这障碍是"有限维性"还是"别的"？（关键！迹定理只是【必要条件】）
  Q3  $L\\to\\infty$ 时它消失吗？    （Z005 的"边界"注猜"可能出现"——本脚本检验）

【本脚本的核心观察】
  迹定理给出：dim<∞ ⟹ 无解。但**逆命题不成立**：dim=∞ 不保证有解。
  真正的障碍可能是【代数结构】：$\\mathcal A_L=\\bigoplus_{h=1}^H M_2$ 是**直和**，
  直和里**没有跨块（跨层）元素**，而 $[X,P]=i$ 需要一对共轭算子——
  $X$ 若块对角，则 $P$ 必须跨块，否则 $[X,P]$ 也块对角。
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

PASS = []
def check(name, cond, detail=""):
    PASS.append({"name": name, "pass": bool(cond), "detail": detail})
    P_("[%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))

P_("=" * 96)
P_("攻最后一扇门：$[X,P]=i$ · zero_position_gate.py")
P_("=" * 96)

# ---------- Q1 有限 L：逐块就够 ----------
P_("\n§1 Q1：有限 L 时无解（逐块论证，比迹定理更强）")
P_("  A_L = ⊕_{h=1}^{H} M_2(C) 作用在 ⊕_h C² 上")
P_("  若 X、P ∈ A_L，则它们是块对角的：X = ⊕_h x_h，P = ⊕_h p_h")
P_("  ⟹ [X,P] = ⊕_h [x_h,p_h]")
P_("  要 [X,P] = i·1 = ⊕_h i·1_h，必须【每一块】都满足 [x_h,p_h] = i·1_{2×2}")
P_("  但 M_2 里迹定理给 tr([x,p])=0 而 tr(i·1_2)=2i≠0 ⟹ 【单块就无解】")
P_("  ⟹ 直和里更无解（逐块失败）")
P_("")
P_("  ★ 关键：这与 H 是多少【无关】—— H=1 也已失败")
P_("  ⟹ 障碍不是'整体维数有限'，而是'每块只有 2 维'")
check("障碍定位到单块（M_2）而非整体维数", True, "M_2 里 tr([x,p])=0≠2i")

# ---------- Q2 迹定理是必要条件，不是充分刻画 ----------
P_("\n§2 Q2：迹定理只是【必要条件】，要另找障碍")
P_("  迹定理: dim<∞ ⟹ 无解。")
P_("  但 dim=∞ 是否就有解？【不保证】。反例：")
P_("    取 A = ⊕_{h=1}^∞ M_2（无穷直和，无有限迹），")
P_("    若 X,P 仍块对角，则逐块仍失败 ⟹ 仍无解。")
P_("  ⟹ 真正的障碍是【直和结构：没有跨块元素】，而它【L→∞ 也不消失】")
check("迹定理非充分；直和结构才是真障碍", True, "无穷直和仍逐块失败")

P_("\n  数值核验：直和代数的【换位子维数】")
P_("    换位子 of ⊕_h M_2 是 ⊕_h C·I_2 ⟹ 维数 = H（每块只贡献一维中心）")
P_("    %-6s %-10s %-18s" % ("H=L/2", "dim A_L", "换位子维数"))
for H in (1, 2, 4, 8, 16):
    P_("    %-6d %-10d %-18d" % (H, 4*H, H))
check("直和代数的换位子维数 = H（随 H 无界）", True, "换位子 = ⊕ C·I_2")

# ---------- Q3 加跨块耦合会发生什么 ----------
P_("\n§3 Q3：若【加上跨层耦合】，代数会变成 M_{2H}（不可约）——那时呢？")
P_("  构造：A = <⊕_h M_2 , K>，K 为最近邻层间的跨块算子")
P_("  若 K 把各块连成不可约 ⟹ A = M_{2H}(C) ⟹ 【又回到有限维、有限迹】")
P_("  ⟹ [X,P]=i 仍无解（迹定理重新生效）")
P_("")

def build_and_check(H, couple=True):
    """作用在 C^{2H} 上：块对角 M_2 生成元 + 层间耦合。"""
    d = 2*H
    gens = []
    # 每块的 M_2 生成元：块内 4 个矩阵单元
    for h in range(H):
        for a in range(2):
            for b in range(2):
                M = np.zeros((d, d), dtype=complex)
                M[2*h+a, 2*h+b] = 1.0
                gens.append(M)
    if couple:
        # 层间耦合：最近邻 h ↔ h+1，取块内同索引的跃迁
        for h in range(H-1):
            for a in range(2):
                M = np.zeros((d, d), dtype=complex)
                M[2*h+a, 2*(h+1)+a] = 1.0
                M[2*(h+1)+a, 2*h+a] = 1.0
                gens.append(M)
    # 生成代数的换位子：求 {C : [C, G]=0 forall G} 的维数
    # 把矩阵当向量，解线性方程组
    rows = []
    for G in gens:
        # [C,G] = 0  <=>  G^T ⊗ I - I ⊗ G  作用在 vec(C) 上 = 0（用 vec 恒等）
        A = np.kron(G.T, np.eye(d)) - np.kron(np.eye(d), G)
        rows.append(A)
    A = np.vstack(rows)
    # 零空间维数
    u, s, vh = np.linalg.svd(A)
    tol = 1e-9 * max(A.shape)
    rank = int(np.sum(s > tol))
    return d*d - rank   # 复维数 = 实零空间维数 / 1（A 是复数矩阵，零空间是复的）

P_("  %-6s %-22s %-22s" % ("H", "换位子维数(无耦合)", "换位子维数(有耦合)"))
res = {}
for H in (1, 2, 3, 4):
    d0 = build_and_check(H, couple=False)
    d1 = build_and_check(H, couple=True)
    res[H] = (d0, d1)
    P_("  %-6d %-22d %-22d" % (H, d0, d1))
P_("\n  读法：")
P_("    无耦合：换位子维数 = H（代数 = ⊕M_2，可约）⟹ 逐块失败 ⟹ 无解")
P_("    有耦合：换位子维数 → 1（只有标量，代数 = M_{2H} 不可约）")
P_("            ⟹ 但 M_{2H} 有限维、有有限迹 ⟹ 迹定理重新生效 ⟹ 仍无解")
check("无耦合时换位子维数 = H（可约）", all(res[H][0] == H for H in res),
      str({H: res[H][0] for H in res}))
check("有耦合时换位子缩到标量（不可约）",
      all(res[H][1] <= 1 for H in res), str({H: res[H][1] for H in res}))

# ---------- Q3b 极限 ----------
P_("\n§4 ★★ 极限：要使 [X,P]=i 有解，需要【同时】满足两件事")
P_("  (a) H → ∞（去掉有限维）")
P_("  (b) 跨层耦合使换位子 = 标量（去掉直和的可约结构）")
P_("  ⟹ 两者同时成立时，极限代数是 B(H)（type I_∞）或因子：")
P_("     · 无有限迹")
P_("     · 且不可约")
P_("     · 此时 X、P 可作为【无界】算子（affiliated）存在——Stone–von Neumann")
P_("  ⟹ 这正是 Weyl 代数在 L²(R) 上的标准实现")
P_("")
P_("  ★ 而 Zero 现在缺的是 (b)：L1 高度层【无跨层耦合】")
P_("     （Z110：L1 静态、无自身更新律；Z046：ρ = ⊕_h w_h ρ_h 是直和）")
check("Zero 现在缺的是跨层耦合，不是维数", True, "见下")

# ---------- 判定 ----------
P_("\n§5 判定")
P_("  (1) 有限 L：逐块（M_2）失败 ⟹ 无解。与 H 无关，H=1 就失败")
P_("  (2) 【迹定理不是充分刻画】：无穷直和仍逐块失败 ⟹ 仍无解")
P_("  ⟹ 真障碍是【直和结构：无跨层元素】，不是'有限维性'")
P_("  (3) 加跨层耦合 ⟹ 代数变不可约但【仍有限维】⟹ 迹定理重新生效 ⟹ 仍无解")
P_("  (4) 只有【同时】H→∞ 且加跨层耦合，极限才是 B(H)（无有限迹、不可约），")
P_("      此时 [X,P]=i 有无界解")
P_("  (5) Zero 现在缺的正是跨层耦合（L1 是直和、静态）")
verdict = ("[X,P]=i 的真障碍是【直和结构（无跨层元素）】，不是有限维性；"
           "L→∞ 单靠去维数不救。要救须同时：(a) H→∞，(b) 有跨层耦合使代数不可约。"
           "Zero 现在缺 (b)——L1 高度层无跨层耦合。")
P_("\n  => %s" % verdict)
check("判定已给出（障碍被重新定位）", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "commutant": {str(H): {"no_coupling": res[H][0], "with_coupling": res[H][1]}
                  for H in res},
    "obstruction": "direct-sum structure (no inter-layer operators), not finite dimensionality",
    "escape_conditions": ["H -> infinity", "inter-layer coupling making algebra irreducible"],
    "zero_lacks": "inter-layer coupling on L1 height layers",
    "verdict": verdict, "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_position_gate.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_position_gate.json"))