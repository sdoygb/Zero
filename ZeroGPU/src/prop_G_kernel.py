#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prop_G_kernel.py —— 命题 G：层间回路给出的有效相互作用核

【口径】本脚本里的 FAIL ＝ **该篇要报告的判定**（「预期形态不成立」），不是脚本错误；
        退出码 1 表示存在此类判定。读退出码前先读本行。
=========================================================
命题 G（把 D31–D36 的"引力产生"图像严格化）：
  设 L2 图上第 t 代质量分布 m_t(v)。边③：R_{t+1}(v)=R_t(v)+g(m_t(v))；
  边②：seed_{t+1}(v)=F(R_t(·))。有效相互作用核 K(v,u) =
  "u 处单位质量扰动经 ③→② 对 v 处的响应"。待判：K 的形状与力程。

F 的三类（对"历史→种子"映射的完备分类）：
  I  自指：F(v)=R(v)                     ⟹ 台账现状（D38 实测零耦合）
  II 局域：F(v)=Σ_{u~v} w R(u)，多代累积  ⟹ 程长 = 代数（有限）
  III 全局：F(v)=Σ_u w(v,u) R(u)         ⟹ w 即引力核（须预设）

判据（PASS = D31–D36 的主张成立）：
  G1a 类 I 给出非零跨顶点耦合
  G2a 类 II 给出长程（幂律、程长与代数无关）
  G3a 1/r 可以只靠层间回路得到（不预设 w）
  G4a 局域回灌的力程与舞台/代数无关且发散
运行：cd ZeroGPU/src && /usr/bin/python3 prop_G_kernel.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import torus_adj      # noqa: E402

RES = os.path.join(HERE, "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


def lap(adj):
    V = len(adj)
    L = np.zeros((V, V))
    for u, nb in enumerate(adj):
        L[u, u] = len(nb)
        for w in nb:
            L[u, w] -= 1.0
    return L


def adjacency(adj):
    V = len(adj)
    A = np.zeros((V, V))
    for u, nb in enumerate(adj):
        for w in nb:
            A[u, w] = 1.0
    return A


def shells(adj, ref=0):
    V = len(adj)
    dist = -np.ones(V, int)
    dist[ref] = 0
    fr = [ref]
    while fr:
        nx = []
        for u in fr:
            for w in adj[u]:
                if dist[w] < 0:
                    dist[w] = dist[u] + 1
                    nx.append(w)
        fr = nx
    return dist


def report(adj, K, label, ref=0):
    dist = shells(adj, ref)
    row = np.abs(K[ref, :])
    prof = {int(d): float(row[dist == d].mean()) for d in np.unique(dist)}
    ds = np.array(sorted(prof))
    ks = np.array([prof[d] for d in ds])
    nz = ds[ks > 1e-12]
    support = int(nz.max()) if len(nz) else 0
    pos = (ks > 1e-12) & (ds > 0)
    pw = float("nan")
    if pos.sum() >= 3:
        A2 = np.vstack([np.log(ds[pos]), np.ones(pos.sum())]).T
        c2, *_ = np.linalg.lstsq(A2, np.log(ks[pos]), rcond=None)
        pw = float(-c2[0])
    print(f"  {label}")
    print("    " + " ".join(f"d{d}:{prof[d]:.2e}" for d in ds[:8])
          + f"  支撑半径={support}  幂律指数={pw:.3f}")
    return {"label": label, "support": support, "power": pw,
            "profile": {str(d): prof[d] for d in ds}}


# ============================================================
print("=" * 72)
print("命题 G：三类 F 的有效核（2D 16×16 与 3D 8³）")
print("=" * 72)
results = {}

print("\n--- 类 I 自指（台账现状）---")
adj2 = torus_adj(2, 16)
V2 = 256
K_I = np.eye(V2)
results["I_self"] = report(adj2, K_I, "F(v)=R(v)")
ck("G1a 类 I（自指）给出非零跨顶点耦合",
   np.abs(K_I - np.diag(np.diag(K_I))).max() > 1e-12,
   f"K 严格对角（非对角元 {np.abs(K_I-np.diag(np.diag(K_I))).max():.1e}）⟹ 零耦合（与 D38 引擎实测一致）")

print("\n--- 类 II 局域（多代累积回灌 Σ_{k≤T} A^k）---")
A2m = adjacency(adj2)
for T in (1, 3, 10, 50):
    K = sum(np.linalg.matrix_power(A2m, k) for k in range(1, T + 1))
    results[f"II_2d_T{T}"] = report(adj2, K, f"2D 局域累积 T={T} 代")
sup = {T: results[f"II_2d_T{T}"]["support"] for T in (1, 3, 10, 50)}
print(f"  支撑半径 vs 代数：{sup}")
ck("G2a 类 II（局域回灌）给出长程力（程长与代数无关）",
   any(abs(sup[T] - T) > 2 for T in (1, 3, 10)),
   f"支撑半径严格等于代数（T=1→{sup[1]}、3→{sup[3]}、10→{sup[10]}）⟹ 有限程、程长 = 回灌代数")

adj3 = torus_adj(3, 8)
V3 = 512
A3m = adjacency(adj3)
print()
for T in (3, 10, 50):
    K = sum(np.linalg.matrix_power(A3m, k) for k in range(1, T + 1))
    results[f"II_3d_T{T}"] = report(adj3, K, f"3D 局域累积 T={T} 代")

print("\n--- 类 III 全局（须预设 w）---")
L3 = lap(adj3)
K_IIIa = np.linalg.pinv(L3)
results["III_poisson_3d"] = report(adj3, K_IIIa, "3D w = L⁺（泊松核）")
idx = np.arange(V3)
x = (idx // 64) % 8; y = (idx // 8) % 8; z = idx % 8
dx = np.minimum(np.abs(x - x[0]), 8 - np.abs(x - x[0]))
dy = np.minimum(np.abs(y - y[0]), 8 - np.abs(y - y[0]))
dz = np.minimum(np.abs(z - z[0]), 8 - np.abs(z - z[0]))
r3 = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
rad = {rr: float(4 * np.pi * rr * K_IIIa[0, np.abs(r3 - rr) < 0.3].mean()) for rr in (1, 2, 3, 4)}
print("  3D L⁺ 的径向 4πr·L⁺：" + " ".join(f"r={rr}:{v:+.4f}" for rr, v in rad.items()))
note("3D L⁺ 的 4πr 读数", {str(k): round(v, 4) for k, v in rad.items()},
     "8³ 环面太小（图像效应），但要点是：读数要成为 1/r 平台必须**预设** L⁺ 本身")
ck("G3a 1/r 可以只靠层间回路得到（不预设 w）",
   abs(results["II_2d_T10"]["power"] - 1.0) < 0.3 and results["II_2d_T10"]["support"] > 20,
   f"类 II（唯一不预设 w 的类）在 T=10 时幂律指数 {results['II_2d_T10']['power']:.3f}、"
   f"支撑半径 {results['II_2d_T10']['support']} ⟹ 不是 1/r；要 1/r 必须把 w 取成泊松核"
   "（类 III-a）＝把 Lφ=−x 装进边②＝Z025 定理 3 的价签")

print()
print("--- 类 II 的奇偶结构（二部图特征）---")
prof = results["II_2d_T10"]["profile"]
even = np.mean([prof[str(d)] for d in range(0, 11, 2) if str(d) in prof])
odd = np.mean([prof[str(d)] for d in range(1, 11, 2) if str(d) in prof])
print(f"  T=10：偶距离平均 {even:.2e}，奇距离平均 {odd:.2e}，比 {even/odd:.2f}")
note("奇偶振荡", f"偶/奇 = {even/odd:.2f}", "二部图的特征（我 D25 已证 K₄-胞网是二部）")

print()
print("=" * 72)
print("G4  谱隙与力程")
print("=" * 72)
evA = np.sort(np.linalg.eigvalsh(A2m))[::-1]
xi_pred = 1.0 / np.log(evA[0] / evA[1])
print(f"  2D 16×16：A 谱 λ_max={evA[0]:.4f}、λ₂={evA[1]:.4f} ⟹ ξ(A^k) = {xi_pred:.3f} 格距")
ck("G4a 局域回灌的力程与舞台无关且发散",
   xi_pred > 100,
   f"ξ = {xi_pred:.3f} 格距（有限、且随舞台谱隙变）⟹ 局域型给有限程，与 D31 §3『不管多远』矛盾")

# ============================================================
print()
print("=" * 72)
n_fail = sum(1 for x in CK if not x["pass"])
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 72)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "prop_G_kernel.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT, "classes": results,
               "support_vs_generation": {str(k): v for k, v in sup.items()},
               "xi_pred_local": xi_pred, "parity_ratio": float(even / odd),
               "config": {"V2": V2, "V3": V3}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/prop_G_kernel.json")
sys.exit(1 if n_fail else 0)