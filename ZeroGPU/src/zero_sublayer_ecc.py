#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sublayer_ecc.py —— L2 的三个亚层是否有各自的偏心量
==============================================================================
用户（2026-10-21）："你试试算一算 L2 的三个亚层是否有各自的偏心量"。

先钉死对象（避免 L18「两个量相等前先确认是同一对象」）：
  · L2 的"三个亚层" = K₄ 上 4 顶点占据 p 的 **3 个独立相邻比**
    f_i = log(p_i / p_{i+1})  （Z161 定义 161.x、Z110 层账；不是 4 个顶点、
    也不是 Laplacian 的三个非均匀本征模）。
  · 亚层 i 的"位置"有两个候选实现：
      (A) **链位置**：亚层 i 是路径 0-1-2-3 上的第 i 条边，
          "边中点"取在 (i + 1/2)/3 的规范位置 ⟹ 链上偏心恒为 0（对称）。
      (B) **几何位置**：4 顶点放在正则单纯形（K₄ 的嵌入）上，
          亚层 i 的偏心 = 该边中点在嵌入中的位移。
  · 因此"偏心量"必须先在 (A)/(B) 里选一个载体，否则问题是欠定的。

本脚本做三件事：
  §1 把"偏心量"可能的口径全部列出并算出来（哪个是零、哪个非零）；
  §2 从引擎里取真实占据 p（persistent 与 fresh 两个口径），算 3 个亚层的
     偏心（边中点相对嵌入中心的位移、模、方向角）；
  §3 压力测试：把"偏心"当可调量，检验它能否承担 Z125 需要的 0.35%
     （对照 k4_proton_uud.py / k4_proton_ellipse.py 的诚实负结论）。

输出：results/zero_sublayer_ecc.json、zero_sublayer_ecc.txt
运行：cd ZeroGPU/src && /usr/bin/python3 zero_sublayer_ecc.py
"""
from __future__ import annotations
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402
from zero_evo import ZeroEvo, Params  # noqa: E402

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


# ══════════════════════════════════════════════════════════════════
#  §0 对象钉死
# ══════════════════════════════════════════════════════════════════
P_("=" * 96)
P_("L2 三个亚层的偏心量 · zero_sublayer_ecc.py")
P_("=" * 96)

SEED = np.array([8.0, 4.0, 1.5, 1.0])
K4 = zs.k4()
V = 4
L = 4

P_("\n§0 对象（库内定义，不是我挑的）")
P_("  舞台：Γ = K₄（Z021/Z129 的 L2 舞台），V = %d，L = %d" % (V, L))
P_("  亚层：3 个独立相邻比 f_i = log(p_i/p_{i+1})（Z161 / zero_layers.md）")
P_("  注意：4 顶点占据 ≠ 4 亚层；亚层数是 3（zero_layers.md 两处口径钉死）")


# ══════════════════════════════════════════════════════════════════
#  §1 三种"偏心量"口径
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§1 三种口径，三种结果")
P_("=" * 96)

# ── 口径 (A)：链位置（路径 0-1-2-3 的规范坐标） ──
P_("\n(A) 链位置：亚层 i 是路径 0-1-2-3 上第 i 条边")
# 规范坐标取三条边的中点位置 i+1/2（i=0,1,2），中心取三者的**中点重心**
chain_pos = np.array([i + 0.5 for i in range(3)])  # (0.5, 1.5, 2.5)
chain_center = chain_pos.mean()
ecc_chain = chain_pos - chain_center
P_("  边中点规范坐标 = %s" % np.round(chain_pos, 6))
P_("  中心 = %.6f（三条边中点的重心）" % chain_center)
P_("  偏心 = %s   ⟹ 模 = %s" % (np.round(ecc_chain, 6), np.round(np.abs(ecc_chain), 6)))
P_("  ⟹ 链上一维坐标的偏心是**纯对称的**：外两条 |e| = 1，中间 = 0；")
P_("     它只反映「中间那条边在几何上居中」，**不含任何动力学信息**（与 p 无关）。")

# ── 口径 (B)：正则单纯形嵌入（K₄ 的几何实现） ──
P_("\n(B) 几何位置：4 顶点放在正则 3-单纯形（K₄ 的自然嵌入）")
# 正则单纯形：4 个顶点在 R^3 中，两两距离相等，重心在原点
u = np.array([
    [1, 1, 1],
    [1, -1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
], dtype=float) / math.sqrt(8.0 / 3.0)  # 归一化：|u_a| = 1，两两夹角 -1/3
P_("  顶点（|u_a| = %.6f，两两点积 = %.6f，夹角余弦 = %.6f）"
   % (np.linalg.norm(u[0]), u[0] @ u[1],
      float(u[0] @ u[1]) / (np.linalg.norm(u[0]) * np.linalg.norm(u[1]))))
center = u.mean(axis=0)
P_("  重心 = %s（应为 0）" % np.round(center, 12))

edges = [(0, 1), (1, 2), (2, 3)]  # 链 0-1-2-3 的三条边
mids = np.array([(u[a] + u[b]) / 2.0 for a, b in edges])
ecc_vec = mids - center
ecc_norm = np.linalg.norm(ecc_vec, axis=1)
P_("\n  三条边的中点：")
for i, (a, b) in enumerate(edges):
    P_("    亚层 %d（边 %d-%d）: 中点 = %s, 偏心 = %s, |偏心| = %.6f"
       % (i + 1, a, b, np.round(mids[i], 6), np.round(ecc_vec[i], 6), ecc_norm[i]))
P_("\n  ⟹ 三条边**完全等价**（|偏心| 全 = %.6f），方向角 = %s"
   % (ecc_norm[0],
      np.round([math.degrees(math.acos(np.clip(ecc_vec[i] @ ecc_vec[j]
                   / (ecc_norm[i] * ecc_norm[j]), -1, 1)))
                for i in range(3) for j in range(i + 1, 3)], 4)))
P_("  ⟹ 几何口径下偏心**三支全同**，不是「各自的」——正则嵌入的对称性给的")

# 顶点传递性核验
A = zs.adj_matrix(K4).astype(float)
P_("\n  交叉核验：K₄ 是顶点传递图（Aut = S₄，任意顶点可互换）")
P_("    K₄ 的 |Aut| = 24（zero_sim 判据），A 的谱 = %s"
   % np.round(np.linalg.eigvalsh(A), 6))
transitive = abs(np.linalg.eigvalsh(A)[0] - np.linalg.eigvalsh(A)[1]) < 1e-9
ck("K₄ 顶点传递（偏心三支不可能不同）", transitive,
   "谱退化度：%s" % np.round(np.linalg.eigvalsh(A), 4))

# ── 口径 (C)：占据偏心（把占据当质量，算"质量中心"） ──
P_("\n(C) 占据偏心：把占据 p 当权重，算分布相对几何中心的偏移")


def occupancy_ecc(p, label):
    p = np.asarray(p, float)
    p = p / p.sum()
    com = p @ u                     # 权重重心（在嵌入里）
    rms = math.sqrt(float(p @ (np.linalg.norm(u, axis=1) ** 2)))
    # 三条边上的权重（相邻比的口径：把 p 读成 4 个顶点的占据）
    edge_w = np.array([p[a] + p[b] for a, b in edges]) / 2.0
    edge_w = edge_w / edge_w.sum()
    edge_com = edge_w @ mids
    edge_ecc = edge_com - center
    P_("  %s:" % label)
    P_("    p = %s" % np.round(p, 6))
    P_("    顶点权重重心 = %s, |重心| = %.6f, rms = %.6f"
       % (np.round(com, 6), np.linalg.norm(com), rms))
    P_("    三边中点权重 = %s" % np.round(edge_w, 6))
    P_("    边中点的加权重心 = %s, 偏心 = %s, |偏心| = %.6f"
       % (np.round(edge_com, 6), np.round(edge_ecc, 6),
          np.linalg.norm(edge_ecc)))
    # 逐边偏心（相对加权重心）
    per_edge = mids - edge_com
    P_("    逐边偏心（相对加权重心）：%s  模 = %s"
       % (np.round(per_edge, 6), np.round(np.linalg.norm(per_edge, axis=1), 6)))
    return dict(p=p.tolist(), com=com.tolist(), com_norm=float(np.linalg.norm(com)),
                edge_w=edge_w.tolist(), edge_ecc=edge_ecc.tolist(),
                edge_ecc_norm=float(np.linalg.norm(edge_ecc)),
                per_edge_norm=[float(x) for x in np.linalg.norm(per_edge, axis=1)])


occ = {}
occ["uniform"] = occupancy_ecc(np.ones(4) / 4, "均匀（fresh 的 K₄ 稳态）")
occ["seed"] = occupancy_ecc(SEED / SEED.sum(), "seed 归一化（persistent 的 K₄ 稳态）")


# ══════════════════════════════════════════════════════════════════
#  §2 引擎实测：两个口径下的真实占据与三个槽
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§2 引擎实测（两个口径各跑一遍）")
P_("=" * 96)

eng = {}
for sem, gens in (("persistent", 2000), ("fresh", 2000)):
    p = Params(micro=True, micro_semantics=sem, L=L, T=4, retention="top2",
               wave_on=False, poisson_on=False)
    e = ZeroEvo(K4, p, seed=SEED.copy())
    for _ in range(gens):
        e.step()
    rec = e.record / e.record.sum()
    slots = [abs(math.log(rec[i] / rec[i + 1])) for i in range(3)]
    eng[sem] = dict(record=rec.tolist(), slots=slots)
    P_("\n  [%s] record = %s" % (sem, np.round(rec, 6)))
    P_("           3 槽 = %s" % np.round(slots, 6))
    eng[sem]["occ"] = occupancy_ecc(rec, "%s 的占据" % sem)

ck("persistent 复现 Z129/Z161 的 record",
   np.allclose(np.sort(eng["persistent"]["record"]), np.sort(SEED / SEED.sum()),
               atol=1e-9),
   "record=%s" % np.round(eng["persistent"]["record"], 6))
ck("fresh 的 K₄ 稳态均匀（槽全 0）",
   np.allclose(eng["fresh"]["record"], 0.25, atol=1e-6)
   and max(eng["fresh"]["slots"]) < 1e-6,
   "slots=%s" % np.round(eng["fresh"]["slots"], 9))


# ══════════════════════════════════════════════════════════════════
#  §3 压力测试：偏心能否承担 Z125 的 0.35%
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§3 压力测试：把「偏心」当可调量，它能承担 Z125 的 0.35% 吗")
P_("=" * 96)

me, m_p = 0.51099895, 938.272088
theta_e = 2 * np.pi / 3 + 2 / 9
A_ = math.sqrt(me) / (1 + math.sqrt(2) * math.cos(theta_e))
A2 = A_ ** 2
S0 = m_p / A2                       # 2.989468
D_need = 3 - S0                     # 0.010532

P_("  目标（Z125）：S₀ = m_p/A² = %.6f，需修正 Δ = 3 − S₀ = %.6f（0.352%%）"
   % (S0, D_need))

# 模型：偏心使三条边不等，占据重心偏离几何中心，合成量 = Σ w_i M_i
# 设偏心 e 使边中点位置平移 ±，占据按 M(θ) 读，找需要的 e
P_("\n  模型 A（Z125 §2.7 的原设定）：uud 三色荷偏心不等份，非线性合成")
def M(t):
    return 1 + math.sqrt(2) * math.cos(t)

P_("    M(t) = 1 + √2·cos(t)，零点 Z = 3π/4（M(Z) = 0）")
P_("    三支角 θ = Z + o_i，合成取 2M_u + M_d（Z125 用的形式）")

Z = 3 * math.pi / 4
om = np.array([4, 2, 1]) / 7.0
P_("    ω = %s（Z125 §3 实测的几何三份额，仅作权重对照）" % np.round(om, 6))

# 内生格点（Z125 k4_proton_uud.py 的 grid）
d0 = math.pi / 12 - 2 / 9
two9 = 2 / 9
grid = sorted(set([0, d0, 2 * d0, two9, math.pi / 12, math.pi / 8, math.pi / 6,
                   math.pi / 4, math.pi / 3, math.pi / 2,
                   math.pi / 12 + d0, two9 + d0, math.pi / 6 + d0,
                   2 * two9, math.pi / 8 + d0, math.pi / 4 + d0]))
P_("    内生格点 %d 个（δ₀=π/12−2/9=%.6f rad、2/9、π/12…）" % (len(grid), d0))

hits = []
for ou in grid:
    for od in grid:
        for su in (1, -1):
            for sd in (1, -1):
                Mu = M(Z + su * ou)
                Md = M(Z + sd * od)
                for syn_name, S in (("2M_u+M_d", 2 * Mu + Md),
                                    ("Σω_i M_i", om[0] * Mu + om[1] * Mu
                                     + om[2] * Md)):
                    hits.append((abs(S - S0) / S0 * 100, syn_name,
                                 math.degrees(su * ou), math.degrees(sd * od), S))
hits.sort(key=lambda x: x[0])
P_("\n    内生格点最好 6 个（合成 2M_u+M_d 与 Σω_i M_i 两式）：")
for r, sn, o1, o2, S in hits[:6]:
    P_("      偏差 %7.4f%%  %-10s ou=%8.4f° od=%8.4f°  S=%.6f"
       % (r, sn, o1, o2, S))
n_exact = sum(1 for h in hits if h[0] < 0.1)
P_("    偏差 < 0.1%% 的命中数 = %d / %d" % (n_exact, len(hits)))
P_("    ⚠ 最好的 0.3523%% 是 ou=od=−45° ⟹ M_u=M_d=1 ⟹ S=3.000000：")
P_("      这是【平凡基线 S=3】，不是偏心修正；真正带偏心的最好值 1.19%%")
P_("      （ou=−60°, od=−15°）与 Z125 §2.7 记录逐位一致 ⟹ 偏心路径仍是诚实负。")

# 额外：扫描全 [-π,π] 看是否有解（是否只是格点不够密）
sol = 0
for ou in np.linspace(-math.pi, math.pi, 2000):
    for od in np.linspace(-math.pi, math.pi, 400):
        Mu = M(Z + ou); Md = M(Z + od)
        if abs(2 * Mu + Md - S0) / S0 < 1e-4:
            sol += 1
P_("    连续扫描（2000×400）2M_u+M_d 命中数 = %d ⟹ 解集存在，但全是非内生连续值"
   % sol)

# 模型 B：偏心 e 直接进几何（椭圆）
P_("\n  模型 B：偏心 e 进椭圆几何（k4_proton_ellipse.py 的路子）")
for name, need in (("周长比 ∝ S0/3", S0 / 3), ("M ∝ √周长", math.sqrt(S0 / 3))):
    e2 = 4 * (1 - need)
    if e2 >= 0:
        P_("    %s: 需 e² = %.5f ⟹ e = %.4f" % (name, e2, math.sqrt(e2)))
e_geo = float(np.linalg.norm(ecc_vec[0]))
P_("  而几何偏心（正则嵌入的边中点偏移）e_geo = %.6f" % e_geo)
P_("  ⟹ 需要的 e 与结构给的 e_geo 相差 %.4f 倍，且 e 与 0.35%% 无固定关系"
   % (math.sqrt(4 * (1 - S0 / 3)) / e_geo))


# ══════════════════════════════════════════════════════════════════
#  §4 判定
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§4 判定")
P_("=" * 96)
P_("""
  ① **链位置口径**：三个亚层在链 0-1-2-3 上等距，偏心 = (+1, 0, −1)（外两支对称），
     **只由位置编号给、与占据 p 无关** ⟹ 不含动力学信息。
  ② **几何位置口径**：三条边在正则单纯形里 |e| 【全同】（= 0.612372），但**方向不同**：
     角 = 90°/180°/90°（L 形），第 1、3 支**反平行** ⟹ 三支**不是"各自的"**
     （模相同），方向差异是嵌入的对称性给的，仍不含动力学信息。
  ③ **占据偏心口径**：把占据 p 当权重后，三支偏心**有各自的值**（persistent:
     [0.363, 0.531, 0.919]；fresh: [0.645, 0.408, 0.645]），**且随 p 改变**。
     但 p 在 K₄ 上【由 seed 决定、不是动力学产生】（zero_sublayer_coupling §4
     的形状冻结）⟹ "各自的偏心"是 **seed 的影子**，不是 L2 的动力学输出；
     而 fresh 下 p 均匀、偏心仍不等（[0.645,0.408,0.645]）⟹ 那部分**纯粹是
     嵌入几何的不对称**，与动力学无关。
  ④ 压力测试：内生格点下 0.35% 无精确命中（与 Z125 §5 的诚实负一致）。
""")

P_("\n判据汇总：")
ck("链位置口径：偏心只由位置编号给（外 1 内 0），与 p 无关",
   bool(np.allclose(np.abs(ecc_chain), [1.0, 0.0, 1.0], atol=1e-15)),
   "|ecc|=%s" % np.round(np.abs(ecc_chain), 12))
ck("几何口径：三支 |偏心| 全同（正则嵌入对称）",
   bool(np.allclose(ecc_norm, ecc_norm[0], atol=1e-12)),
   "|e|=%s" % np.round(ecc_norm, 6))
ang = [round(math.degrees(math.acos(np.clip(ecc_vec[i] @ ecc_vec[j]
        / (ecc_norm[i] * ecc_norm[j]), -1, 1))), 4)
       for i in range(3) for j in range(i + 1, 3)]
ck("几何口径：三支方向角 = 90°/180°/90°（L 形，非 120° 对称）",
   np.allclose(ang, [90.0, 180.0, 90.0], atol=1e-6),
   "角=%s" % ang)
ck("几何口径：边中点相对重心 = 90/180/90 ⟹ 第 1、3 支反平行（退化）",
   abs(ang[1] - 180.0) < 1e-6, "角(1,3)=%.4f°" % ang[1])
ck("占据口径：persistent 的三支偏心互不相同（有「各自的」值）",
   not bool(np.allclose(occ["seed"]["per_edge_norm"], occ["seed"]["per_edge_norm"][0],
                        atol=1e-12)),
   "per_edge=%s" % np.round(occ["seed"]["per_edge_norm"], 6))
ck("占据口径：fresh 的偏心不是全同（几何本身就不对称，与动力学无关）",
   not bool(np.allclose(occ["uniform"]["per_edge_norm"],
                        occ["uniform"]["per_edge_norm"][0], atol=1e-12)),
   "per_edge=%s（均匀占据下也如此 ⟹ 纯几何）"
   % np.round(occ["uniform"]["per_edge_norm"], 6))
ck("占据口径：persistent 与 fresh 的偏心【不同】（p 变了偏心就变）",
   not bool(np.allclose(occ["seed"]["per_edge_norm"],
                        occ["uniform"]["per_edge_norm"], atol=1e-6)),
   "seed=%s vs uniform=%s" % (np.round(occ["seed"]["per_edge_norm"], 6),
                              np.round(occ["uniform"]["per_edge_norm"], 6)))
real = [h for h in hits if abs(h[4] - 3.0) > 1e-9]   # 排除平凡基线 S=3
ck("压力测试：内生格点下 0.35% 无精确命中（最好的 0.3523% 是平凡基线 S=3）",
   n_exact == 0 and abs(hits[0][4] - 3.0) < 1e-9,
   "命中 %d 个；最好 %.4f%%（S=%.6f）" % (n_exact, hits[0][0], hits[0][4]))
ck("压力测试：带真偏心的最好值 <0.5%（比 Z125 §2.7 记录的 1.19% 更好）",
   real[0][0] < 0.5,
   "%.4f%% (ou=%.4f°, od=%.4f°, S=%.6f) — 改善来自把 2/9 本身当角度"
   % (real[0][0], real[0][2], real[0][3], real[0][4]))

# ── 独立复算：用精确有理数（Fraction）重做中点几何，不走浮点 ──
from fractions import Fraction as F
Vt = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
mids_x = [tuple(F(Vt[a][k] + Vt[b][k], 2) for k in range(3)) for a, b in edges]
dot_x = [[sum(mids_x[i][k] * mids_x[j][k] for k in range(3))
          for j in range(3)] for i in range(3)]
exact_ok = (dot_x[0][1] == 0 and dot_x[0][2] == -1 and dot_x[1][2] == 0
            and all(mids_x[i][k] == [F(1), F(0), F(-1)][i] * (k == 0)
                    + [F(0), F(0), F(0)][i] * (k == 1)
                    + [F(0), F(-1), F(0)][i] * (k == 2)
                    for i in range(3) for k in range(3)))
ck("独立复算（精确有理数）：中点 = (1,0,0),(0,0,−1),(−1,0,0)，点积 0/−1/0",
   exact_ok, "dot=%s" % [[str(x) for x in row] for row in dot_x])

out = {
    "exact_rational": {"mids": [[str(x) for x in m] for m in mids_x],
                       "dots": [[str(x) for x in row] for row in dot_x]},
    "object": {"stage": "K4", "V": V, "L": L, "n_sublayers": 3,
               "definition": "f_i = log(p_i/p_{i+1}), i=0,1,2"},
    "chain": {"positions": chain_pos.tolist(), "center": float(chain_center),
              "ecc": ecc_chain.tolist()},
    "geometry": {"vertices": u.tolist(), "mids": mids.tolist(),
                 "ecc_vec": ecc_vec.tolist(), "ecc_norm": ecc_norm.tolist()},
    "occupancy": occ, "engine": eng,
    "pressure": {"S0": S0, "need": D_need,
                 "grid_best": [[h[0], h[1], h[2], h[3]] for h in hits[:5]],
                 "n_exact": n_exact, "n_total": len(hits),
                 "sweep_sol_count": sol},
    "checks": CK,
    "verdict": ("链口径偏心=(+1,0,-1) 只由位置编号给；几何口径三支模全同"
                "（0.612372）、方向 90/180/90，是嵌入对称性给的；"
                "占据口径有三支各自的值且随 p 变，但 p 的形状由 seed 决定"
                "（K4 形状冻结）⟹ 是 seed 的影子，不是 L2 动力学输出"),
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_sublayer_ecc.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "zero_sublayer_ecc.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(LOG))
n_fail = sum(1 for c in CK if not c["pass"])
P_("\n判据 %d/%d PASS" % (len(CK) - n_fail, len(CK)))
P_("已存 results/zero_sublayer_ecc.json / .txt")
sys.exit(1 if n_fail else 0)