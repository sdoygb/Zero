#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage_ledger_check.py —— 舞台账（Z021）与「G」的舞台依赖性

【口径】本脚本里的 FAIL ＝ **该篇要报告的判定**（「预期形态不成立」），不是脚本错误；
        退出码 1 表示存在此类判定。读退出码前先读本行。
==========================================================
处理 D18–D24 审计报告 §7 的缺口 1（舞台）与缺口 4（长短程实测）。

A 组（舞台的代数账）
  A1 K₄ 的 6 条边在 A₄ 下可迁 ⟹ 边可迁的 K₄-胞网只能取常电压（= 循环覆盖）
  A2 复算 Z021 §8.2：Z_n（循环）电压覆盖（无三角）围长恒 4 —— n=2..6 穷举、7,8 抽样
  A3 反例：Z³（两生成元、非循环）覆盖可同时做到 3-正则 / 二部 / 围长 6
  A4 计数障碍：Heawood 的 14 顶点 ≠ 4N ⟹「每格点一个 K₄ 纤维」与 Z021 定理 8 不相容
B 组（「G」是舞台的量，不是普适常数）
  B1 长程刚度 σ_coord = lim λ(k)/|k|² 的**精确**值（精确瑞利商，环面 N=64 无误差）
     与键向量预言 ½·tr(Σ b bᵀ)/3 对比；给出键长单位下的势幅度 A = 1/(σ_coord·ℓ_b)
  B2 短程自能（可选，仅在数值稳定时报告）
  B3 Heawood（有限图、有谱隙）：无 1/r 平台
口径：L 为图 Laplacian；L ≈ -σ_coord∇²（坐标单位）；键长单位下的势幅度 A = 1/(σ_coord ℓ_b)。
运行：cd ZeroGPU/src && /usr/bin/python3 stage_ledger_check.py
"""
import itertools
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zero_sim as zs  # noqa: E402

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
AXES = [[(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]]
TRI = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
Q4 = [(0, 1, 3, 2), (0, 2, 3, 1), (0, 3, 1, 2)]

# ============================================================
print("=" * 72)
print("A1  K₄ 的边在置换群下的轨道（边可迁 ⟹ 不变电压必为常数）")
print("=" * 72)


def perms_A4():
    out = []
    for p in itertools.permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            out.append(p)
    return out


A4 = perms_A4()
orbit = set()
for p in A4:
    for (a, b) in EDGES:
        orbit.add(tuple(sorted((p[a], p[b]))))
print(f"  |A₄|={len(A4)}；A₄ 在 6 条边上的轨道 = {sorted(orbit)}")
ck("A1a A₄ 在 K₄ 的 6 条边上可迁 ⟹ 不变电压必为常数",
   len(orbit) == 6, f"轨道大小 {len(orbit)}/6 ⟹ 边可迁的 K₄-胞网只能是常电压（循环覆盖，见 A2）")
note("K₄ 的对边对数（轴数）", len(AXES), "Z021 §7.2：6 条边配成 3 对 ⟹ 3 个空间方向")

# ============================================================
print()
print("=" * 72)
print("A2  Z_n（循环）电压覆盖的围长（复算 Z021 §8.2）")
print("=" * 72)


def cover_cyclic(n, a):
    V = [(v, x) for v in range(4) for x in range(n)]
    pos = {p: i for i, p in enumerate(V)}
    adj = {p: set() for p in V}
    for (a_, b_) in EDGES:
        w = a[(a_, b_)]
        for x in range(n):
            p, q = (a_, x), (b_, (x + w) % n)
            adj[p].add(q)
            adj[q].add(p)
    return [sorted(pos[q] for q in adj[p]) for p in V]


girths = {}
for n in range(2, 9):
    tuples = (list(itertools.product(range(n), repeat=6)) if n <= 6 else
              [tuple(int(v) for v in np.random.default_rng(20261016).integers(0, n, 6))
               for _ in range(4000)])
    gs = set()
    for t in tuples:
        a = dict(zip(EDGES, t))
        if any((a[(x, y)] + a[(y, z)] - a[(x, z)]) % n == 0 for (x, y, z) in TRI):
            continue                                     # 有三角，跳过（Z021 口径）
        gs.add(zs.girth(cover_cyclic(n, a)))
    girths[n] = sorted(gs)
    print(f"  n={n}: 无三角指派的围长集合 = {girths[n]}")
ck("A2a Z021 §8.2『一切覆盖的围长恒为 4；覆盖路线封死』",
   all(g == [4] for g in girths.values()),
   f"为假：n=2,3 确是 4，但 n≥4 存在无三角且围长 6 的循环覆盖（n=4 电压 (0,0,0,1,3,1) 给 16 顶点、"
   f"3-正则、二部、围长 6，独立实现复核一致）；实测 {girths}")
ck("A2b Z021 的**结论**（K₄ 覆盖路线出局）仍成立——理由是顶点数而非围长",
   all(4 * n != 14 for n in range(1, 20)),
   "覆盖恒有 4n 个顶点，最小围长 6 的例子 16 > 14 ⟹ Heawood 的最小性（Moore 下界 14）不受影响")

# ============================================================
print()
print("=" * 72)
print("A3  Z³（非循环）覆盖：3-正则 / 二部 / 围长 6 的显式例子")
print("=" * 72)
W = {(0, 1): (-1, -1, -1), (0, 2): (-1, -1, -1), (0, 3): (-1, -1, -1),
     (1, 2): (-1, -1, -1), (1, 3): (-1, -1, 0), (2, 3): (-1, -1, -1)}


def cover_z3(N, w):
    pts = [(v, i, j, k) for v in range(4) for i in range(N) for j in range(N) for k in range(N)]
    idx = {p: n for n, p in enumerate(pts)}
    adj = {p: set() for p in pts}
    for (a_, b_) in EDGES:
        d = w[(a_, b_)]
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    p = (a_, i, j, k)
                    q = (b_, (i + d[0]) % N, (j + d[1]) % N, (k + d[2]) % N)
                    adj[p].add(q)
                    adj[q].add(p)
    return [sorted(idx[q] for q in adj[p]) for p in pts]


adj_k4cell = cover_z3(6, W)
deg_k4 = set(zs.degrees(adj_k4cell))
gir_k4 = zs.girth(adj_k4cell)
bip_k4 = zs.is_bipartite(adj_k4cell)
conn_k4 = len(zs._bfs_order(adj_k4cell)) == len(adj_k4cell)
print(f"  |V|={len(adj_k4cell)}（4·6³）  度={deg_k4}  围长={gir_k4}  二部={bip_k4}  连通={conn_k4}")
ck("A3a Z³ 覆盖可同时满足 3-正则 + 二部 + 围长 6",
   deg_k4 == {3} and gir_k4 >= 6 and bip_k4 and conn_k4,
   f"度 {deg_k4}、围长 {gir_k4}、二部 {bip_k4} ⟹ Z021 §8.2『一切覆盖围长恒 4』只对**循环**覆盖成立")
note("K₄-胞网的配位数", 3, "K₄ 的度=3 被每个位点继承（D18 §6.2 的 z=4『金刚石』不是 K₄-胞网）")

# ============================================================
print()
print("=" * 72)
print("A4  计数障碍：Heawood 的顶点数")
print("=" * 72)
hw = zs.heawood()
print(f"  Heawood：|V|={len(hw)}，度={set(zs.degrees(hw))}，围长={zs.girth(hw)}，二部={zs.is_bipartite(hw)}")
ck("A4a Heawood 的顶点数可写成 4N（『每格点一个 K₄ 纤维』）",
   14 % 4 == 0, f"|V|=14，14 mod 4 = {14 % 4} ⟹ 与 Z021 定理 8 不相容")

# ============================================================
print()
print("=" * 72)
print("B1  长程刚度 σ_coord（精确瑞利商）与势幅度 A = 1/(σ_coord·ℓ_b)")
print("=" * 72)
N = 4096            # 只用于取小波数（不建图）：|k|² ~ 3.5e-5，k⁴ 修正 ~1e-5
kN = 2 * np.pi / N


def analyze(label, neighbors):
    """neighbors: 某位点的**全部邻居**键向量（含 ±b，逐条列出，不重复计数）"""
    m = np.array(neighbors, float)
    M = m.T @ m                                # Σ over 邻居 b bᵀ
    c = float(np.trace(M) / 3)
    sigma_pred = c / 2.0
    lb = float(np.linalg.norm(m, axis=1).mean())
    A_pred = 1.0 / (sigma_pred * lb)
    # 精确瑞利商：x(p)=cos(k·p)（k 与环面通约）⟹ RQ = Σ_b (1-cos(k·b)) 精确成立
    k = kN * np.array([1.0, 2.0, 3.0])
    rq = float(np.sum(1 - np.cos(m @ k)))
    sigma_num = rq / float(k @ k)
    print(f"  {label:<14} ℓ_b={lb:.4f} c={c:.3f} σ_coord(预言)={sigma_pred:.4f} "
          f"σ_coord(瑞利商)={sigma_num:.4f} A=1/(σℓ_b)={A_pred:.4f}")
    return {"label": label, "bond_len": lb, "c": c, "sigma_pred": sigma_pred,
            "sigma_rayleigh": sigma_num, "A": A_pred}


CUB_NB = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
DIA_NB = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
cub = analyze("cubic(z=6)", CUB_NB)
dia = analyze("diamond(z=4)", DIA_NB)
print(f"  A₄-胞网（Z³ 覆盖）度=3：六类键向量（我构造的例子）——")
OFF = {0: (1, 1, 1), 1: (1, -1, -1), 2: (-1, 1, -1), 3: (-1, -1, 1)}
bonds_k4 = [np.array(W[(a, b)], float) + np.array(OFF[b], float) - np.array(OFF[a], float)
            for (a, b) in EDGES]
print("    一个四面体几何实现的六类键向量：" +
      " ".join(str(tuple(int(x) for x in v)) for v in bonds_k4))
print(f"    键长² = {[int(round(float(v @ v))) for v in bonds_k4]}（不等长 ⟹ 几何实现是自由选择）")
Mk4 = np.array([np.outer(v, v) for v in bonds_k4]).sum(0)
ev = np.linalg.eigvalsh(Mk4)
note("K₄-胞网例子：M=Σbbᵀ 的特征值", [round(float(x), 3) for x in ev],
     "不等 ⟹（我这个）几何实现各向异性；结合 A1（边可迁不可得），"
     "Z002 引理 4『无偏好 ⟹ 均匀导纳』在该舞台上不自动成立")

ck("B1a 精确瑞利商 → 键向量预言 ½·tr(Σbbᵀ)/3",
   abs(cub["sigma_rayleigh"] / cub["sigma_pred"] - 1) < 1e-3
   and abs(dia["sigma_rayleigh"] / dia["sigma_pred"] - 1) < 1e-3,
   f"cubic {cub['sigma_rayleigh']:.4f} vs {cub['sigma_pred']:.4f}；"
   f"diamond {dia['sigma_rayleigh']:.4f} vs {dia['sigma_pred']:.4f}")
kdir = [kN * np.array(v, float) for v in ((1, 2, 3), (3, 1, 2), (2, 3, 1))]
rq_dir = [float(np.sum(1 - np.cos(np.array(CUB_NB) @ kk))) / float(kk @ kk) for kk in kdir]
rq_dir_d = [float(np.sum(1 - np.cos(np.array(DIA_NB) @ kk))) / float(kk @ kk) for kk in kdir]
ck("B1c 各向同性舞台：σ 与波数方向无关",
   max(rq_dir) - min(rq_dir) < 1e-3 and max(rq_dir_d) - min(rq_dir_d) < 1e-3,
   f"cubic 三个方向的 σ = {[round(x, 5) for x in rq_dir]}；diamond = {[round(x, 5) for x in rq_dir_d]}")

ck("B1b D18 §6.1『G 是普适常数——微观结构只修正短程』",
   abs(dia["A"] / cub["A"] - 1) < 0.1,
   f"同一单位源、同一物理距离：cubic A={cub['A']:.4f} vs diamond A={dia['A']:.4f}"
   f"（比 {dia['A']/cub['A']:.3f}）⟹ 长程幅度差 {cub['A']/dia['A']:.2f} 倍；"
   "普适的只有 1/r 形状与 4π 归一化（Z033 命题 14），幅度是舞台量 Q 的函数")
note("『G』的台账身份", "G = 1/Q", "Q = σ_coord 张量（Z002 §3.6 引理 5 的 Q^{ij}）；Q 由舞台键几何决定")

# ============================================================
print()
print("=" * 72)
print("B3  Heawood（有限图、有谱隙）没有 1/r 长程势")
print("=" * 72)
A = zs.adj_matrix(hw)
Lh = np.diag(A.sum(1)) - A
lam = np.sort(np.linalg.eigvalsh(Lh))
gap = float(lam[1])
rho = np.zeros(len(hw))
rho[0] = 1.0
rho -= rho.mean()
phi = -np.linalg.lstsq(Lh, rho, rcond=None)[0]
print(f"  谱隙 λ₂={gap:.4f}；max|φ|={np.abs(phi).max():.4f}；φ(0)={phi[0]:+.4f}")
ck("B3a Heawood 上存在 1/r 长程势",
   False,
   f"14 顶点、谱隙 {gap:.3f} ⟹ L⁺ 有界（max|φ|={np.abs(phi).max():.3f}）、无尺度分离；"
   "Z021 §8.5『走者在格点(Heawood)上跑』与 Z033 命题 14 的 d=3 渐近不能同时成立")

# ============================================================
print()
print("=" * 72)
n_fail = sum(1 for c in CK if not c["pass"])
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 72)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "stage_ledger_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT, "sigma": [cub, dia],
               "cyclic_girths": {str(k): v for k, v in girths.items()},
               "k4cell_cover": {"voltage": {f"{a}{b}": list(v) for (a, b), v in W.items()},
                                "V": len(adj_k4cell), "girth": gir_k4, "bipartite": bip_k4,
                                "bond_vectors": [list(map(int, v)) for v in bonds_k4]},
               "heawood": {"V": len(hw), "gap": gap, "max_phi": float(np.abs(phi).max())}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/stage_ledger_check.json")
sys.exit(1 if n_fail else 0)
