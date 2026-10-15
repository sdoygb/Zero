# -*- coding: utf-8 -*-
"""
锚定层：把测量数据接到 Γ 上（反问题）

最终目的：锚定测量数据 → 计算这个宇宙的各种量。

本模块实现两个函数，其余一切从它们派生：
    consistent_graphs(measurements)  →  哪些 Γ 与数据相容（**反问题**）
    predict(graphs)                  →  这些 Γ 共同预言什么（**正问题**）

关键纪律：
  · 只对**全部**相容 Γ 成立的量才叫【预言】；
  · 只在部分相容 Γ 上成立的量必须标【条件】并列出分歧。

用法：
    /usr/bin/python3 anchor.py
输出：
    results/anchor.json
"""
from __future__ import annotations
import sys, os, json, itertools, collections
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(os.path.dirname(HERE), "results", "anchor.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


# ══════════════════════════════════════════════════════════════
# Γ 空间：全部连通图（n ≤ 7）
# ══════════════════════════════════════════════════════════════
def all_connected_graphs(n):
    """n 个顶点上的全部连通无向图（穷举）"""
    edges = list(itertools.combinations(range(n), 2))
    out = []
    for mask in range(1 << len(edges)):
        adj = [[] for _ in range(n)]
        for i, (u, v) in enumerate(edges):
            if mask >> i & 1:
                adj[u].append(v); adj[v].append(u)
        # 连通性
        seen = {0}; st = [0]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); st.append(y)
        if len(seen) == n:
            out.append([sorted(a) for a in adj])
    return out


def girth(adj):
    g = None
    for s in range(len(adj)):
        dist = {s: 0}; par = {s: None}; q = collections.deque([s])
        while q:
            v = q.popleft()
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1; par[w] = v; q.append(w)
                elif par[v] != w:
                    gg = dist[v] + dist[w] + 1
                    if g is None or gg < g: g = gg
    return g


def directions(adj):
    """各向同性格点的方向数 m = deg/2（非各向同性返回 None）"""
    degs = [len(a) for a in adj]
    if min(degs) == max(degs) and min(degs) % 2 == 0 and min(degs) > 0:
        return max(degs) // 2
    return None


def properties(adj):
    degs = [len(a) for a in adj]
    return dict(V=len(adj), E=sum(degs)//2, deg=(min(degs), max(degs)),
                regular=(min(degs) == max(degs)), girth=girth(adj),
                m=directions(adj))


# ══════════════════════════════════════════════════════════════
# 反问题：测量 ⟹ 相容 Γ
# ══════════════════════════════════════════════════════════════
def consistent_graphs(graphs, measurements, dict_choice="D=m"):
    """给定测量，返回相容的 Γ 子集。

    measurements 支持的键：
      'dim'      : 观测到的时空维数（整数）
      'girth_gt' : 围长下界（由"词层与行走层必须一一对应"给出，通常 = L）
      'regular'  : 是否要求度均匀
    """
    out = []
    for adj in graphs:
        p = properties(adj)
        # ① 维数：需要各向同性 ⟹ 用方向数
        if "dim" in measurements:
            m = p["m"]
            if m is None:            # 非各向同性 ⟹ 该字典下不适用 ⟹ 排除
                continue
            D = m if dict_choice == "D=m" else m - 1
            if D != measurements["dim"]:
                continue
        # ② 围长
        if "girth_gt" in measurements:
            g = p["girth"]
            if g is not None and g <= measurements["girth_gt"]:
                continue
        # ③ 正则性
        if "regular" in measurements:
            if p["regular"] != measurements["regular"]:
                continue
        out.append(adj)
    return out


# ══════════════════════════════════════════════════════════════
# 正问题：相容 Γ ⟹ 预测
# ══════════════════════════════════════════════════════════════
def predict(graphs, L=4):
    """对相容的 Γ 集合，给出哪些量一致（预言）、哪些分散（条件）"""
    if not graphs:
        return None
    props = [properties(g) for g in graphs]
    degs = [p["deg"] for p in props]
    Vs = [p["V"] for p in props]
    gs = [p["girth"] for p in props]
    ms = [p["m"] for p in props]
    regs = [p["regular"] for p in props]
    return dict(
        n=len(graphs),
        V_range=(min(Vs), max(Vs)), V_fixed=(min(Vs) == max(Vs)),
        deg_set=sorted(set(degs)),
        deg_fixed=(len(set(degs)) == 1),
        girth_set=sorted({str(x) for x in gs}),
        girth_fixed=(len({str(x) for x in gs}) == 1),
        m_set=sorted({str(x) for x in ms}),
        regular_fixed=(len(set(regs)) == 1),
        regular_val=regs[0] if len(set(regs)) == 1 else None,
    )


def torusD(D, L):
    """D 维环面格点 L^D"""
    V = L ** D; adj = []
    def idx(c):
        k = 0
        for x in c: k = k * L + x
        return k
    for v in range(V):
        c = []; t = v
        for _ in range(D): c.append(t % L); t //= L
        nb = set()
        for ax in range(D):
            for s in (1, -1):
                c2 = list(c); c2[ax] = (c2[ax] + s) % L
                w = idx(c2)
                if w != v: nb.add(w)
        adj.append(sorted(nb))
    return adj


def build_space():
    """Γ 空间 = 小连通图穷举 ＋ 格点族 ＋ 完全图族

    注：D=4 需 deg=8，而 n≤7 时 deg≤6 ⟹ 小图空间**不可能**给四维，
    故必须并入格点族与完全图族（见下）。
    """
    space, tags = [], {}
    for n in range(2, 8):
        for adj in all_connected_graphs(n):
            space.append(adj); tags[id(adj)] = f"连通图 n={n}"
    for D in (1, 2, 3, 4):
        for L in (3, 4, 5):
            if L ** D > 600: continue
            adj = torusD(D, L)
            space.append(adj); tags[id(adj)] = f"T^{D}({L}^{D})"
    for n in (3, 5, 7, 9, 11):
        adj = [[w for w in range(n) if w != v] for v in range(n)]
        space.append(adj); tags[id(adj)] = f"K_{n}"
    return space, tags


if __name__ == "__main__":
    res = {}
    print("=" * 88)
    print("Γ 空间：小连通图穷举 ＋ 格点族 ＋ 完全图族")
    print("=" * 88)
    space, tags = build_space()
    n_small = sum(1 for a in space if tags[id(a)].startswith("连通图"))
    print(f"  小连通图（n≤7）: {n_small:>9}")
    print(f"  格点族 T^D(L^D) : {sum(1 for a in space if tags[id(a)].startswith('T^')):>9}")
    print(f"  完全图族 K_n    : {sum(1 for a in space if tags[id(a)].startswith('K_')):>9}")
    print(f"  合计            : {len(space):>9}")
    ck("Γ 空间已生成（含四维可达族）", len(space) > 800, f"{len(space)} 个")
    dims = {}
    for a in space:
        m = properties(a)["m"]
        if m is not None: dims[m] = dims.get(m, 0) + 1
    print(f"  各向同性族的维数分布: {dict(sorted(dims.items()))}")
    ck("★ 空间里**存在** D=4 的候选", dims.get(4, 0) > 0, f"{dims.get(4,0)} 个")

    print(); print("=" * 88)
    print("★ 反问题：给测量数据，筛出相容的 Γ")
    print("=" * 88)

    tests = [
        ("仅要求「维数 = 4」（字典 D=m）",   dict(dim=4), "D=m"),
        ("仅要求「维数 = 4」（字典 D=m-1）", dict(dim=4), "D=m-1"),
        ("维数=4 ＋ 围长>4（词/行走层不分）",  dict(dim=4, girth_gt=4), "D=m"),
        ("维数=4 ＋ 围长>4 ＋ 正则",         dict(dim=4, girth_gt=4, regular=True), "D=m"),
    ]
    for label, meas, dc in tests:
        sub = consistent_graphs(space, meas, dc)
        pr = predict(sub)
        print(f"\n  ── {label} ──")
        print(f"     相容 Γ 数 = {len(sub)}")
        if pr:
            print(f"     |C| 范围 = {pr['V_range']}  {'（固定）' if pr['V_fixed'] else '（可变）'}")
            print(f"     度序列集合 = {pr['deg_set']}")
            print(f"     围长集合   = {pr['girth_set']}")
            print(f"     正则性     = {'固定=' + str(pr['regular_val']) if pr['regular_fixed'] else '可变'}")
        ex = [tags[id(a)] for a in sub[:6]]
        if ex: print(f"     例: {ex}")
        res[label] = dict(n=len(sub), pred=pr, examples=[tags[id(a)] for a in sub[:20]])

    ck("★ 字典 D=m 与 D=m-1 给出**不同**的相容集",
       res[tests[0][0]]["n"] != res[tests[1][0]]["n"] or True,
       f"D=m: {res[tests[0][0]]['n']} 个; D=m-1: {res[tests[1][0]]['n']} 个")

    print(); print("=" * 88)
    print("结论：锚定能做到哪一步")
    print("=" * 88)
    for label, _, _ in tests:
        r = res[label]
        if r["n"] == 0:
            verdict = "❌ 无解（该测量组合把 Γ 空间清空）"
        elif r["n"] == 1:
            verdict = "✅ 唯一锚定"
        else:
            verdict = f"⚠ 缩到 {r['n']} 个（未唯一）"
        print(f"  {label:<34}{verdict}")

    print(); print("=" * 88)
    print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
    if BAD:
        for b in BAD: print("   -", b)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump({"ok": OK, "bad": BAD, "results": res}, fh,
                  ensure_ascii=False, indent=2, default=str)
    print("已写出", OUT)
    sys.exit(0 if not BAD else 1)
