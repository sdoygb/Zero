"""
z0_bound_states.py --- 粒子机制：Γ 上的缺陷束缚态。

**为什么算这个**

前几轮我说"最小系统没有粒子"，依据是**均匀格点**上的谱密度连续、gapless。
但"Γ 是均匀的"是我**默认**的——这又是一个预设。Z1 只说"存在有限连通图 Γ"，没说是均匀的。

判据（标准散射/微扰论）：
    有限秩微扰**不改变本质谱** => 宿主连续谱之外的任何本征值都是**真束缚态**。
    宿主取长度 N 的链，其邻接算子连续谱 = [-2, 2]（N -> 无穷）。
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import networkx as nx


def attach_cycle(G: nx.Graph, c: int, k: int, nid: int):
    """在顶点 c 上挂一个 k-环（顶点编号 nid..nid+k-1）。"""
    for j in range(k):
        G.add_edge(c if j == 0 else nid + j - 1, nid + j if j < k - 1 else c)
    return G


def bound_states(G: nx.Graph, edge: float = 2.0, tol: float = 1e-9):
    """返回连续谱 [-edge, edge] 之外的本征值及其逆参与比。"""
    A = nx.to_numpy_array(G)
    w, V = np.linalg.eigh(A)
    ipr = (V ** 4).sum(axis=0)
    out = (w < -edge - tol) | (w > edge + tol)
    idx = np.nonzero(out)[0]
    return [(float(w[i]), float(ipr[i])) for i in idx]


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    res = {}

    print("=" * 92)
    print(f"Γ 上的束缚态（宿主：长度 {N} 的链；连续谱 = [-2, 2]；有限秩微扰不改变本质谱）")
    print("=" * 92)
    base = nx.path_graph(N)
    bs = bound_states(base)
    print(f"  {'干净链（无缺陷）':<22} 带外能级 {len(bs)} 个")
    res["clean"] = {"count": len(bs)}

    print()
    table = []
    for k in (3, 4, 5, 6, 8, 12):
        G = nx.path_graph(N)
        attach_cycle(G, N // 2, k, N)
        b = bound_states(G)
        table.append({"k": k, "count": len(b),
                      "energies": [round(e, 6) for e, _ in sorted(b)],
                      "loc_sites": [round(1.0 / p, 2) for _, p in sorted(b)]})
        loc = ", ".join(f"{1/p:.1f}格点@λ={e:+.4f}" for e, p in sorted(b))
        print(f"  {'链 + ' + str(k) + '-环':<22} 带外能级 {len(b)} 个   {loc}")
    res["single_defect"] = table

    print()
    print("=" * 92)
    print("两个缺陷：束缚态能量会不会移动？（= 有没有相互作用）")
    print("=" * 92)
    G1 = nx.path_graph(N)
    attach_cycle(G1, N // 2, 4, N)
    b1 = sorted(e for e, _ in bound_states(G1))
    print(f"  单个四环                       能量 = {[round(x,6) for x in b1]}")
    inter = []
    for d in (4, 5, 8, 12, 20, 60, 150):
        G2 = nx.path_graph(N)
        attach_cycle(G2, N // 2, 4, N)
        attach_cycle(G2, N // 2 + d, 4, N + 10)
        b2 = bound_states(G2)
        pos = sorted(e for e, _ in b2 if e > 0)
        # 与单缺陷正能量的最大偏差（配对比较：两个缺陷 -> 两个正能级）
        ref = [x for x in b1 if x > 0][0]
        dev = max(abs(x - ref) for x in pos) if pos else float("nan")
        inter.append({"d": d, "positive_energies": [round(x, 8) for x in pos],
                      "max_shift_vs_single": dev})
        print(f"  两个四环 d={d:<4}                正能量 = "
              f"{[round(x,6) for x in pos]}   最大位移 = {dev:.3e}")
    res["two_defects"] = inter
    print("""
  判读（看完数据才写的）：
    ① d 小（4,5,8）时能级**分裂**（2.383 -> 2.402/2.360）=> 两个局域模**有重叠**
    ② d >= 20 时能级**与单缺陷完全简并**（位移 < 1e-6）=> **无长程相互作用**
    ③ 力程 ~ 局域化长度（~5 格点）。这是**短程有效相互作用**，
       机制是局域波函数的指数小重叠——正是束缚态物理的标准图像。""")

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "z0_bound_states.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("\nsaved results/z0_bound_states.json")


if __name__ == "__main__":
    main()
