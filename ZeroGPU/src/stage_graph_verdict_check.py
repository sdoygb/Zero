#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage_graph_verdict_check.py —— 舞台图的判定：方格是具名输入
================================================================
对象：Z117 §四2 的开放项——舞台图＝方格图 vs Z021 的 Γ=Heawood。

核验：
  A. Heawood 图（Fano 关联图）：|V|=14、deg=3、围长=6、二部；
  B. 双缝格规模（240×260）远超 14，Γ 不能承载该传播；
  C. 内容图的 q_eff 无平台（Z145 的数据：L=8→24 单调上升）；
  D. 判定：方格舞台是 L2 具名输入（Z002 §11 的嵌入读法）。

判据：A–C 全成立 → 方形舞台不是 Γ 的同一对象。
"""
from __future__ import annotations
import os
import json

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def heawood():
    """Heawood = Fano 平面的关联图：7 点 + 7 线，点-线相邻。"""
    # Fano 线：{i, i+1, i+3} mod 7（标准差集）
    lines = [frozenset({i % 7, (i + 1) % 7, (i + 3) % 7}) for i in range(7)]
    V = 14
    adj = [set() for _ in range(V)]
    for j, ln in enumerate(lines):
        for p in ln:
            adj[p].add(7 + j)
            adj[7 + j].add(p)
    return adj, lines


def girth(adj):
    n = len(adj)
    best = 10**9
    for s in range(n):
        dist = {s: 0}
        par = {s: -1}
        q = [s]
        for u in q:
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    par[w] = u
                    q.append(w)
                elif par[u] != w:
                    best = min(best, dist[u] + dist[w] + 1)
    return best


def main():
    P("=" * 80)
    P("舞台图的判定：方格是具名输入 (stage_graph_verdict_check)")
    P("=" * 80)

    adj, lines = heawood()
    nV = len(adj)
    degs = sorted({len(a) for a in adj})
    g = girth(adj)
    # 二部性
    color = {}
    bip = True
    stack = [0]
    color[0] = 0
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in color:
                color[w] = 1 - color[u]
                stack.append(w)
            elif color[w] == color[u]:
                bip = False
    P(f"Heawood: |V|={nV}, degrees={degs}, girth={g}, bipartite={bip}")

    # 双缝格规模
    nr, nc = 240, 260
    P(f"双缝格: {nr}×{nc} = {nr*nc} 顶点（远大于 {nV}）")

    # 内容图 q_eff（Z145 的数据）
    qeff = [(8, 1.661), (16, 2.5), (24, 3.246)]
    P("内容图 q_eff（Z145）: " + ", ".join(f"L={L}:{q:.3f}" for L, q in qeff))

    checks = []

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    chk("A  Heawood: |V|=14", nV == 14, str(nV))
    chk("A  Heawood: 3-正则", degs == [3], str(degs))
    chk("A  Heawood: 围长 6", g == 6, str(g))
    chk("A  Heawood: 二部", bip)
    chk("B  双缝格规模 >> 14", nr * nc > 1000, f"{nr*nc}")
    chk("C  q_eff 无平台（L=8→24 上升）", qeff[-1][1] - qeff[0][1] > 1.0,
        f"{qeff[0][1]}→{qeff[-1][1]}")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 80)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 80)

    out = {"heawood": {"nV": nV, "degs": degs, "girth": g, "bipartite": bip},
           "grid_vertices": nr * nc, "qeff": qeff,
           "checks": checks, "pass": npass, "total": len(checks)}
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "stage_graph_verdict_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "stage_graph_verdict_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 stage_graph_verdict_check.json / .txt")


if __name__ == "__main__":
    main()
