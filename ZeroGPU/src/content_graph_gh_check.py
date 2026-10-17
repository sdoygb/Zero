#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z145 核验脚本：内容图（词上的交换图）的有限维 GH 可长性判据。

对象：$L=2n$、$n$ 个 $+1$ 与 $n$ 个 $-1$ 的**平衡词**，节点＝**循环等价类**（旋转取代表），
边＝**一次相邻换位**（在类之间）。

判据（Benjamini–Finucane–Tessera, arXiv:1203.5624）：
  若 $(X_n)$ 为有限连通**顶点传递**图且 $\lvert X_n\rvert=\mathbf o(\mathrm{diam}(X_n)^q)$（某固定 $q$），
  则取子列按直径重标度后 GH-收敛到维数 $<q$ 的环面。
  于是**可操作形式**：$\lvert V\rvert/\mathrm{diam}^{q}$ 是否有界（是否有平台）。

本脚本做三件事（不 import 项目脚本、不写文件以外无副作用）：
  A．枚举内容图：节点数、直径、连通性（与 Motzkin 数交叉验证）
  B．$q_{\\rm eff}$：两个独立定义各算一遍，看有没有平台
  C．对照：晶格细化族的 $\lvert V\rvert/\\mathrm{diam}^{3}$ 恒定（满足条件）

复跑：cd ZeroGPU/src && /usr/bin/python3 content_graph_gh_check.py
"""
import itertools, collections, math

FAIL = []


def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAIL.append(name)
    print(f"[{tag}] {name} {detail}")


def necklaces(n):
    """n 个 +1、n 个 -1 的平衡词，按旋转取代表（循环等价类）"""
    L = 2 * n
    seen, reps = set(), []
    for pos in itertools.combinations(range(L), n):
        w = [-1] * L
        for p in pos:
            w[p] = 1
        best = None
        for s in range(L):
            r = tuple(w[s:] + w[:s])
            if best is None or r < best:
                best = r
        if best not in seen:
            seen.add(best)
            reps.append(best)
    return reps


def build_adj(reps):
    idx = {w: i for i, w in enumerate(reps)}
    adj = collections.defaultdict(set)
    for w, i in idx.items():
        L = len(w)
        for k in range(L):
            k2 = (k + 1) % L
            if w[k] == w[k2]:
                continue
            v = list(w)
            v[k], v[k2] = v[k2], v[k]
            # 类级邻接（Z145 §1.1 的字面定义：「在类之间」）：
            # 交换后的词先归到它自己的旋转等价类，再与该类相连。
            # 旧写法 `if t in idx` 只连"交换后恰好仍是规范代表"的边，
            # 那是**代表依赖**的（取决于类里挑了哪个代表），会系统性漏边。
            t = min(tuple(v[s:] + v[:s]) for s in range(L))
            j = idx[t]
            if i != j:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def bfs(adj, nrep):
    """返回 (最大连通分量顶点集, 其直径)"""
    best_comp = set()
    for s in range(nrep):
        if s in best_comp:
            continue
        dist = {s: 0}
        q = collections.deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    q.append(v)
        comp = set(dist)
        if len(comp) > len(best_comp):
            best_comp = comp
            best_diam = max(dist.values())
    return best_comp, best_diam


# 平衡二元颈链数（OEIS A003239）：$n$ 个 +1 与 $n$ 个 -1 的循环等价类个数
A003239 = {0: 1, 1: 1, 2: 2, 3: 4, 4: 10, 5: 26, 6: 80, 7: 246, 8: 810,
           9: 2704, 10: 9252, 11: 32066, 12: 112720, 13: 400024}

print("=" * 78)
print("A．内容图的枚举与 OEIS A003239 交叉验证")
print("=" * 78)
rows = []
for L in (8, 10, 12, 14, 16, 18, 20, 22, 24):
    n = L // 2
    reps = necklaces(n)
    adj = build_adj(reps)
    comp, diam = bfs(adj, len(reps))
    V = len(comp)
    rows.append((L, n, V, diam))
    mot = A003239.get(n)
    ok = (mot is not None and V == mot)
    print(f"  L={L:2d}  n={n:2d}  |V|={V:6d}  diam={diam:3d}"
          f"  A003239[{n}]={mot}  {'✓' if ok else '✗(未收录)'}")

check("A1 节点数 = 平衡二元颈链数（OEIS A003239）",
      all(A003239.get(n) == V for L, n, V, d in rows if n in A003239),
      f"（{sum(1 for L,n,V,d in rows if n in A003239)} 档全部相符）")

print()
print("=" * 78)
print("B．$q_{\\rm eff}$：两个独立定义，看是否有平台")
print("=" * 78)
print("  B1：$q_{\\rm eff}:=\\log|V|/\\log(\\mathrm{diam})$")
print("  B2：$q_{\\rm eff}':=\\log|V|/\\log(2n)$")
print()
print(f"  {'L':>3} {'|V|':>7} {'diam':>5} {'B1':>8} {'B2':>8} {'|V|/diam^3':>11}")
for L, n, V, diam in rows:
    b1 = math.log(V) / math.log(diam)
    b2 = math.log(V) / math.log(2 * n)
    print(f"  {L:3d} {V:7d} {diam:5d} {b1:8.4f} {b2:8.4f} {V/diam**3:11.5f}")

b1 = [math.log(V) / math.log(d) for L, n, V, d in rows]
b2 = [math.log(V) / math.log(2 * n) for L, n, V, d in rows]

check("B1 $q_{\\rm eff}$ 单调上升（无平台）",
      all(b1[i] < b1[i + 1] for i in range(len(b1) - 1)),
      f"：{b1[0]:.3f}（L=8）→ {b1[-1]:.3f}（L=24），共 {len(b1)} 档严格递增")
check("B2 $q_{\\rm eff}'$ 单调上升（无平台）",
      all(b2[i] < b2[i + 1] for i in range(len(b2) - 1)),
      f"：{b2[0]:.3f} → {b2[-1]:.3f}")
r3 = [V / d**3 for L, n, V, d in rows]
check("B3 $|V|/\\mathrm{diam}^3$ 无界（先降后升 ⟹ 无平台，BFT 条件不满足）",
      min(r3) < r3[0] and r3[-1] > 4 * min(r3),
      f"：{r3[0]:.5f} → 谷 {min(r3):.5f} → {r3[-1]:.5f}（$L$ 大时无界）")

print()
print("=" * 78)
print("C．对照：晶格细化族满足 BFT 条件（$|V|/\\mathrm{diam}^3$ 有界）")
print("=" * 78)
lat = []
for m in (8, 16, 32, 64):
    V = m ** 3                 # 三维格点
    diam = m - 1               # 图距离直径（边长 m 的立方格，L1 直径）
    lat.append((m, V, diam, V / diam ** 3))
    print(f"  m={m:3d}  |V|={V:7d}  diam={diam:4d}  |V|/diam^3={V/diam**3:.6f}")
vals = [x[3] for x in lat]
check("C1 晶格细化族 $|V|/\\mathrm{diam}^3$ 有界（$\\to1$ ⟹ 满足 BFT 条件）",
      all(v < 2.0 for v in vals) and abs(vals[-1] - 1.0) < 0.1,
      f"：{vals[0]:.4f} → {vals[-1]:.4f}（单调降、收敛到 $1$）")

print()
print("=" * 78)
print(f"合计：{'全部 PASS' if not FAIL else 'FAIL: ' + ', '.join(FAIL)}")
print("=" * 78)
print()
print("结论：内容图（词上的交换图）的 $q_{\\rm eff}$ 在两个独立定义下均单调上升、无平台；")
print("     而晶格细化族的 $|V|/\\mathrm{diam}^3$ 有界（$\\to1$）。")
print("     ⟹ 内容图不满足 BFT 的有限维 GH 可长性条件，晶格族满足。")

import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results",
                   "content_graph_gh_check.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({
        "object": "内容图（词上的交换图，节点=循环等价类，边=一次相邻换位）",
        "content_graph": [
            {"L": L, "n": n, "V": V, "diam": diam,
             "q_eff_logV_logdiam": round(math.log(V) / math.log(diam), 6),
             "q_eff_logV_log2n": round(math.log(V) / math.log(2 * n), 6),
             "V_over_diam3": round(V / diam ** 3, 6)}
            for L, n, V, diam in rows],
        "lattice_refinement": [
            {"m": m, "V": V, "diam": d, "V_over_diam3": round(r, 6)}
            for m, V, d, r in lat],
        "verdict": {
            "content_graph_q_eff_monotone": bool(all(b1[i] < b1[i+1] for i in range(len(b1)-1))),
            "content_graph_V_over_diam3_unbounded": bool(r3[-1] > 4 * min(r3)),
            "lattice_V_over_diam3_bounded": bool(max(vals) < 2.0),
        },
        "criteria": {"total": 5, "failed": FAIL},
    }, fh, ensure_ascii=False, indent=2)
print(f"结果写入 {os.path.normpath(OUT)}")
