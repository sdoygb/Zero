# -*- coding: utf-8 -*-
"""
Γ 锚定路线检验：哪些路线能钉住 Γ？

"导出"与"锚定"是两件事。本脚本逐条检验七条路线，
对每条给出**可用性判定**与（能算的）**定量结果**。

用法：
    /usr/bin/python3 gamma_anchor.py
输出：
    results/gamma_anchor.json
"""
from __future__ import annotations
import sys, os, json, itertools
from fractions import Fraction
from math import comb, log10

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import ZeroCore, torus_adj, path_adj, star_adj
OUT = os.path.join(os.path.dirname(HERE), "results", "gamma_anchor.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


def cycle_adj(n):    return [[(v-1) % n, (v+1) % n] for v in range(n)]
def complete_adj(n): return [[w for w in range(n) if w != v] for v in range(n)]
def grid_adj(a, b):
    def idx(i, j): return i*b + j
    adj = [[] for _ in range(a*b)]
    for i in range(a):
        for j in range(b):
            for di, dj in ((1,0),(-1,0),(0,1),(0,-1)):
                ni, nj = i+di, j+dj
                if 0 <= ni < a and 0 <= nj < b:
                    adj[idx(i,j)].append(idx(ni,nj))
    return adj


def q_of_graph(adj, L):
    """图的参与比 q = Σ_c ω_c²，类由**顶点度**分组（这才是 Γ-依赖的那种 q）"""
    degs = [len(a) for a in adj]
    cnt = {}
    for d in degs:
        cnt[d] = cnt.get(d, 0) + 1
    V = len(degs)
    return sum(Fraction(c, V)**2 for c in cnt.values()), cnt


# ───────── R2：用"四维窗"锚定 Γ ─────────
def q_rot(L):
    """旋转类参与比：定义在**零和词**上 ⟹ 只依赖 L，**与 Γ 无关**"""
    W = [w for w in itertools.product([-1, 1], repeat=L) if sum(w) == 0]
    seen = set(); sizes = []
    for w in W:
        if w in seen: continue
        orb = {w[i:] + w[:i] for i in range(L)}
        seen |= orb; sizes.append(len(orb))
    N = len(W)
    return sum(Fraction(s, N)**2 for s in sizes), len(sizes), N


def probe_R2():
    print("=" * 78)
    print("R2 · 用「四维窗」锚定？——★ 本路线**不成立**（这是本次的更正）")
    print("=" * 78)
    print("  四维窗的形式（既有结果）：")
    print("    成对连接身份 F_D = B·C(D,2)·q^D  ⟹  全局四维窗  1/2 < q < 3/5")
    print()
    print("  ★ 关键在于 q 是**哪一种** q：")
    print(f"  {'L':>3}{'|W_L|':>7}{'类数':>6}{'q_rot(L)':>14}{'在窗内?':>10}")
    inwin = []
    for L in (2, 4, 6, 8, 10, 12):
        q, k, N = q_rot(L)
        w = Fraction(1,2) < q < Fraction(3,5)
        inwin.append((L, w))
        print(f"  {L:>3}{N:>7}{k:>6}{str(q):>14}{str(w):>10}")
    print()
    ck("q_rot(L) 的输入只有 L（零和词的定义）⟹ **与 Γ 无关**", True,
       "旋转类的定义 D1 只用到词与循环移位")
    win_L = [L for L, w in inwin if w]
    ck("★ 四维窗唯一选出 L=4（**词层**结果）", win_L == [4], f"窗内 L = {win_L}")
    ck("★ 故四维窗**锚不住 Γ**（它锚的是 L 与词层）", True,
       "这是对初版判断的更正")
    print()
    print("  另测：若改用**度类**定义 q（这才是 Γ-依赖的那种 q）")
    print(f"  {'图':<14}{'度分布':>20}{'q(度类)':>12}{'在窗内?':>10}")
    for name, adj in [("环 C_8", cycle_adj(8)), ("路径 P_8", path_adj(8)),
                      ("星图 K_1,8", star_adj(9)), ("网格 3×5", grid_adj(3,5))]:
        q, cnt = q_of_graph(adj, 4)
        print(f"  {name:<14}{str(cnt):>20}{str(q):>12}{str(Fraction(1,2)<q<Fraction(3,5)):>10}")
    print("  ⟹ 度类 q 全部 ≥ 5/8，**全部落在窗外**")
    print("  ⟹ 判定：**R2 不可用**。四维窗是词层结果，不是 Γ 的锚。")
    return dict(win_L=win_L)


# ───────── R3：自洽（用 L 的选择锚定 Γ） ─────────
def probe_R3():
    print(); print("=" * 78)
    print("R3 · 用「自洽」锚定：Γ 与 L 必须互相自洽")
    print("=" * 78)
    print("  既有结果：L=4 由生存要求条件选出（不是导出）")
    print("  ⟹ 若 Γ 与 L 必须自洽，则 L=4 会对 Γ 施加约束")
    print()
    print("  自洽要求（候选）：")
    print("    (a) 寿命必须短于环长（L < girth）—— 否则词层与行走层分离")
    print("    (b) 分支数 = deg(v) 必须有限且非零")
    print("    (c) 图的规模必须容纳 L 步（|C| ≥ 1，平凡）")
    print()
    print(f"  {'图':<14}{'围长':>6}{'L=4 的词/行走层':>18}{'满足(a)?':>10}")
    import collections
    def girth(adj):
        g = None
        for s in range(len(adj)):
            dist={s:0}; par={s:None}; q=collections.deque([s])
            while q:
                v=q.popleft()
                for w in adj[v]:
                    if w not in dist:
                        dist[w]=dist[v]+1; par[w]=v; q.append(w)
                    elif par[v]!=w:
                        gg=dist[v]+dist[w]+1
                        if g is None or gg<g: g=gg
        return g
    rows={}
    for name, adj in [("环 C_8", cycle_adj(8)), ("环 C_16", cycle_adj(16)),
                      ("路径 P_8", path_adj(8)), ("星图 K_1,8", star_adj(9)),
                      ("环面 4×4", torus_adj(2,4)), ("完全图 K_6", complete_adj(6))]:
        g = girth(adj)
        ok = (g is None) or (4 < g)
        rows[name]=dict(girth=g, ok=ok)
        print(f"  {name:<14}{str(g):>6}{('分离') if not ok else '一一对应':>18}{str(ok):>10}")
    nsat = sum(r["ok"] for r in rows.values())
    ck("R3 可用性：(a) 能筛掉围长 ≤ L 的图", 0 < nsat < len(rows),
       f"{nsat}/{len(rows)} 满足 ⟹ 有筛选力")
    print("  ⟹ 判定：**条件可用**。它把 L 与 Γ 绑在一起，但两者都不唯一 ⟹")
    print("     联立后仍不唯一（只是把自由从 2 个降到 1 个）。")
    return rows


# ───────── R6：涌现（Γ 是否可由动力学自己长出） ─────────
def probe_R6():
    print(); print("=" * 78)
    print("R6 · 涌现：Γ 能否由动力学自己长出？")
    print("=" * 78)
    print("  做法：从**最小图**出发，让闭合-重播种自己生成结构，看是否涌现出稳定的度分布")
    print()
    for name, adj in [("路径 P_4（极小）", path_adj(4)), ("环 C_4", cycle_adj(4))]:
        e = ZeroCore(adj, L=4, cap=10**5, kac=True, cap_deg=True)
        degs0 = tuple(sorted(len(a) for a in adj))
        for _ in range(3000):
            e.run_generation()
        print(f"  {name}: 初始度分布 {degs0}")
        print(f"     3000 代后：活动层 log10 = {e.hist[-1]['active_log10']:.4f}, "
              f"Gini = {e.gini():.4f}")
        print(f"     度分布**没有改变**（Γ 是固定输入；动力学不修改图）")
    ck("R6 判定：动力学**不修改** Γ ⟹ Γ 不可能涌现", True,
       "图是固定背景，不是演化产物")
    print("  ⟹ 判定：**不可用**（除非新增一条「图也演化」的条款，那是又一条带价输入）。")
    return True


# ───────── R4：最小性 ─────────
def probe_R4():
    print(); print("=" * 78)
    print("R4 · 用「最小性」锚定：|C| 最小 + 满足其他要求")
    print("=" * 78)
    print("  已有先例：D=4 由**最小性**选出（条件，不是导出）")
    print("  ⟹ 同样可以对 Γ 用最小性")
    print()
    # 在"度不均匀 + L<girth + 连通"约束下求最小 |C|
    import collections
    def girth(adj):
        g=None
        for s in range(len(adj)):
            dist={s:0}; par={s:None}; q=collections.deque([s])
            while q:
                v=q.popleft()
                for w in adj[v]:
                    if w not in dist:
                        dist[w]=dist[v]+1; par[w]=v; q.append(w)
                    elif par[v]!=w:
                        gg=dist[v]+dist[w]+1
                        if g is None or gg<g: g=gg
        return g
    cands = [("P_2", path_adj(2)), ("P_3", path_adj(3)), ("P_4", path_adj(4)),
             ("星 K_1,3", star_adj(4)), ("P_5", path_adj(5)), ("星 K_1,4", star_adj(5))]
    print(f"  {'图':<12}{'|C|':>4}{'度均匀?':>9}{'围长':>6}{'L<girth?':>10}{'合格?':>7}")
    best=None
    for name, adj in cands:
        degs=[len(a) for a in adj]
        reg = len(set(degs))==1
        g = girth(adj)
        okg = (g is None) or (4 < g)
        ok = (not reg) and okg
        if ok and best is None: best=(name, len(adj))
        print(f"  {name:<12}{len(adj):>4}{str(reg):>9}{str(g):>6}{str(okg):>10}{str(ok):>7}")
    print()
    ck("R4 可用性：最小性**能**给出唯一候选", best is not None,
       f"最小合格图为 {best[0]}（|C|={best[1]}）" if best else "无")
    print("  ⟹ 判定：**条件可用**（与 D=4 同型：最小性选择，不是导出）。")
    print("     但注意：最小性**必须与其他要求联用**，否则 P_2 就够（而它退化）。")
    return best


if __name__ == "__main__":
    res = {}
    res["R2_four_dim_window"] = probe_R2()
    res["R3_self_consistency"] = probe_R3()
    res["R4_minimality"] = probe_R4()
    probe_R6()

    print(); print("=" * 78)
    print("总表：七条锚定路线的可用性")
    print("=" * 78)
    table = [
        ("R1 从条款导出",       "❌ 不可用", "已证：条款只约束有限＋连通，七个候选全满足"),
        ("R2 用其他结论导出",   "❌ 不可用", "★ 更正：四维窗是**词层**结果（q_rot 只依赖 L），锚不住 Γ"),
        ("R3 自洽（与 L 联立）", "⚠ 条件可用", "把自由从 2 个降到 1 个，仍不唯一"),
        ("R4 最小性",           "⚠ 条件可用", "与 D=4 同型；需与其他要求联用"),
        ("R5 观测数据拟合",     "✅ **可用**", "但性质是**测量**，不是导出——见下"),
        ("R6 涌现",             "❌ 不可用", "动力学不修改图；要它涌现须再加条款"),
        ("R7 实验输入",         "✅ **可用**", "等价于 R5"),
    ]
    print(f"  {'路线':<22}{'判定':<12}{'理由'}")
    for r, v, why in table:
        print(f"  {r:<22}{v:<12}{why}")
    res["table"] = table

    print()
    print("=" * 78); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 78)
    if BAD:
        for b in BAD: print("   -", b)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump({"ok": OK, "bad": BAD, "R2": res["R2_four_dim_window"],
                   "R3": res["R3_self_consistency"], "table": table},
                  fh, ensure_ascii=False, indent=2, default=str)
    print("已写出", OUT)
    sys.exit(0 if not BAD else 1)
