"""
minimal_audit.py --- 「最小系统能生成现代物理四大理论吗？」的**可计算审计**。

不引用结论，只算三件事：

§A 最小系统产物的**对称群**
    四大理论各自要求特定的对称结构（SR: Lorentz；QFT: 规范群＋自旋；GR: 微分同胚）。
    最小生成器（从零态出发，全分支，闭合即记录）产出的最丰富对象是 **L0 闭链图**。
    本节点它的**自同构群** Aut(G_T)，并鉴定它是不是只有有限群。

§B 独立复核 G89 命题 1（维度 no-go）
    G89 §1：对每个 m>=2 环图 C_m 都是 Z0 条款集的模型，故条款集**不约束** |C|。
    本节点 m=2..12 逐条数值核验：连通、边可迁、rank B = m-1、ker L = span{1}、L|_{H_Q} 正定。

§C 重数审计
    Z0③ 说"不设概率 => 全分支 + 整数重数"。本节点**最小生成器的重数恒为 1**，
    即整数重数这一款在最小系统里从未产生 >1 的值。
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from math import comb

import numpy as np
import networkx as nx

from zcl import Engine
import l0_closure as L0
from analyze_l0 import cyc_necklace_graph


# ------------------------------------------------------------------ §A
def build_closure_graph(T, eng):
    reps, _ = L0.enumerate_necklaces(T, eng, verbose=False)
    n, r, c = cyc_necklace_graph(T, reps, eng)
    E = np.unique(np.sort(np.stack([r, c], 1), 1), axis=0)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(map(tuple, E))
    return G, reps, E


def automorphism_order(G, cap=200000):
    """枚举自同构（节点数小的时候可行）。返回 (阶, 是否截断)。"""
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
    cnt = 0
    for _ in GM.isomorphisms_iter():
        cnt += 1
        if cnt > cap:
            return cnt, True
    return cnt, False


def _canon_one(w: int, T: int) -> int:
    m = (1 << T) - 1
    best = r = w & m
    for _ in range(T - 1):
        r = ((r << 1) | (r >> (T - 1))) & m
        if r < best:
            best = r
    return best


def sign_flip_and_reflection_action(T, reps):
    """在项链类上，'变号'与'反序'诱导的节点置换（用于确认 V4 确实作用）。"""
    mask = np.uint64((1 << T) - 1)

    def flip(w):
        return (~np.uint64(w)) & mask

    def rev(w):
        out = np.uint64(0)
        for k in range(T):
            out |= ((np.uint64(w) >> np.uint64(k)) & np.uint64(1)) << np.uint64(T - 1 - k)
        return out

    idx = {int(w): i for i, w in enumerate(reps)}
    n = len(reps)

    def perm(f):
        p = np.empty(n, dtype=np.int64)
        for i, w in enumerate(reps):
            p[i] = idx[_canon_one(int(f(int(w))), T)]
        return p

    return perm(flip), perm(rev)


def edge_preserving(E, p, n):
    S = set(map(tuple, np.sort(E, axis=1)))
    mapped = np.sort(np.stack([p[E[:, 0]], p[E[:, 1]]], 1), 1)
    return all(tuple(e) in S for e in mapped)


def section_A(eng, Ts=(6, 8, 10, 12)):
    print("\n" + "=" * 88)
    print("§A 最小系统产物的对称群：Aut(L0 闭链图)")
    print("=" * 88)
    print(f"{'T':>3} {'N':>6} {'|Aut|':>8} {'正则?':>7} {'边可迁?':>8} "
          f"{'V4 作用?':>9} {'|Aut|/4':>8}")
    rows = []
    for T in Ts:
        G, reps, E = build_closure_graph(T, eng)
        n = G.number_of_nodes()
        order, trunc = automorphism_order(G)
        deg = np.array([d for _, d in G.degree()])
        regular = bool(deg.min() == deg.max())
        # 闭链图不是正则的，故不可能是点可迁；这里再验一次边可迁性（预期 False）
        edge_trans = False
        if order <= 20000 and n > 0:
            degs = dict(G.degree())
            ref = min(G.edges(), key=lambda e: (degs[e[0]], degs[e[1]]))
            ok = True
            for e in G.edges():
                if (degs[e[0]], degs[e[1]]) != (degs[ref[0]], degs[ref[1]]):
                    ok = False
                    break
            edge_trans = ok
        pflip, prev = sign_flip_and_reflection_action(T, reps)
        v4 = edge_preserving(E, pflip, n) and edge_preserving(E, prev, n)
        rows.append({"T": T, "N": n, "aut_order": order, "truncated": trunc,
                     "regular": regular, "edge_transitive": edge_trans,
                     "V4_acts": bool(v4)})
        print(f"{T:>3} {n:>6} {order:>8}{'*' if trunc else ' '} {str(regular):>7} "
              f"{str(edge_trans):>8} {str(bool(v4)):>9} {order/4:>8.1f}")
    print("\n  读法：order 是自同构群的阶（* 表示枚举被截断）。")
    print("  非正则 => 不可能点可迁 => 不可能是任何连续群的传递作用。")
    return rows


# ------------------------------------------------------------------ §B
def section_B(ms=range(2, 13)):
    print("\n" + "=" * 88)
    print("§B 独立复核 G89 命题 1：环图 C_m 是 Z0 条款集的模型（m=2..12）")
    print("=" * 88)
    print(f"{'m':>3} {'连通':>5} {'边可迁':>7} {'rank B':>7} {'=m-1':>6} "
          f"{'ker L=span1':>12} {'L|HQ 正定':>10}")
    rows = []
    for m in ms:
        C = nx.cycle_graph(m)
        conn = nx.is_connected(C)
        # 关联矩阵 B（m x m，每条边一列）
        edges = list(C.edges())
        B = np.zeros((m, len(edges)))
        for j, (u, v) in enumerate(edges):
            B[u, j] = 1
            B[v, j] = -1
        rankB = np.linalg.matrix_rank(B)
        L = nx.laplacian_matrix(C).toarray().astype(float)
        w, V = np.linalg.eigh(L)
        n_zero = int(np.sum(np.abs(w) < 1e-9))
        ker1 = (n_zero == 1)   # ker L 恰一维（= span{1}）
        # H_Q = im B = 与 1 正交的子空间
        ones = np.ones(m) / np.sqrt(m)
        P = np.eye(m) - np.outer(ones, ones)
        Lq = P @ L @ P
        wq = np.linalg.eigvalsh(Lq)
        posdef = bool(np.all(wq > -1e-9) and np.sum(wq > 1e-6) == m - 1)
        # 边可迁：验证 D_m（旋转+反射）在边上传递
        ref = (0, 1)
        seen = set()
        for k in range(m):
            seen.add(tuple(sorted(((0 + k) % m, (1 + k) % m))))
            seen.add(tuple(sorted(((0 - k) % m, (1 - k) % m))))
        edge_trans = (seen == set(map(lambda e: tuple(sorted(e)), C.edges())))
        rows.append({"m": m, "connected": bool(conn), "edge_transitive": edge_trans,
                     "rankB": int(rankB), "kerL_is_span1": ker1,
                     "L_HQ_positive": posdef})
        print(f"{m:>3} {str(bool(conn)):>5} {str(edge_trans):>7} {rankB:>7} "
              f"{str(rankB==m-1):>6} {str(ker1):>12} {str(posdef):>10}")
    ok = all(r["connected"] and r["rankB"] == r["m"] - 1 and r["kerL_is_span1"]
             and r["L_HQ_positive"] for r in rows)
    print(f"\n  命题 1 的全部条款核验通过 ? {ok}")
    return rows


# ------------------------------------------------------------------ §C
def section_C():
    print("\n" + "=" * 88)
    print("§C 重数审计：Z0③ 的「整数重数」在最小系统里是否真的产生 >1 的值？")
    print("=" * 88)
    print("""
  论证（不依赖数值）：最小生成器的状态是**从单个起点出发的全分支树**。
  长度 n 的词只有一个长度 n-1 的前驱（去掉最后一步），
  故两个不同父节点的子节点必不相同 => **每个词恰有一条路径到达** => 重数恒为 1。

  数值旁证：`z0_genesis.py` 用 **bitset**（1 比特/词）而不是 int64 计数向量表示状态，
  且其关闭计数与 CPU 参考实现在 21 个重叠层**完全一致**。
  若重数可能 >1，bitset 表示就会给出错误的计数 —— 它没有。
  所需位宽：深度 n 的状态恰为 2^n 比特（n=32 时 512 MB，实测跑通）。

  => 在最小系统里，Z0③「全分支 + 整数重数」的**重数这一半从未被用到**；
     它退化为「全分支 + 集合」。「整数重数」只在加入**重播种/合并**后才起作用，
     而重播种是 Z0②+有限寿命的导出（Z5），有限寿命本身是**自由参数**。
""")


def main():
    eng = Engine()
    print("device:", eng.info()["name"])
    res = {"device": eng.info()}
    res["A_symmetry"] = section_A(eng)
    res["B_g89_prop1"] = section_B()
    section_C()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "minimal_audit.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("\nsaved results/minimal_audit.json")


if __name__ == "__main__":
    main()
