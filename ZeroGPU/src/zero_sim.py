#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sim —— Zero 模拟器**公共核心**（结构 ＋ 算法 ＋ 判据自检）

为什么有这个文件：`ZeroGPU/src/` 里已有 130+ 个判定脚本，但公共结构被反复复制
（`fano_heawood()` 在 algebraic_capacity / capacity_destruction / saturation_alignment /
distribution_engine 里各写了一份，`girth()` 在 anchor 里）。本模块把**结构、算法、几何**
收拢成一层，新脚本 `from zero_sim import ...` 即可 —— 别再复制粘贴。

三层结构（与 Zero 的分层对应）：
  · **L0 结构层**：`k4` / `heawood`（＝Fano 关联图）/ `cycle` / `path` / `grid` / `double_cover`
  · **L0 算法层**：`degrees` `girth` `spectrum` `aut_order` `closed_walks` `return_prob`
  · **L1 词层**：`layer_distribution`（DP，L 任意）/ `balanced_words`（穷举，小 L）/ `orbit_count`
  · **耦合几何**：`simplex_vectors` `gram` `rho`

判据（`sim_selftest`，写 `results/zero_sim.json`）覆盖，全部可复跑：
  Heawood：|V|=14、度 3、二部、围长 6、|E|=21、谱 {3, √2⁶, (−√2)⁶, −3}、零化多项式
           (A−3)(A+3)(A²−2)=0、|Aut|=336、顶点稳定子 24、其于 3 条入射边的作用 = 6（满 S₃）、
           闭道闭式 7a_{2m}−3^{2m}=6·2^m（a_k=Tr A^k/14）、质量隙 m=−log(√2/3)=0.7520386984、
           关联长度 ξ=1.3297188059（Z060/Z061）
  词层：  Σ_h b_h = C(L,L/2)、b_h 恒偶、b_H=2、Z032(2.1) 镜像闭式 b_h=N_≤(h)−N_≤(h−1)（Z032/Z047/Z049）
  耦合：  正则 D-单纯形（D=1..6）：零和、等长、ρ=−1/(m−1)、Gram 秩 = D（Z021）
  组态图层：补偿移动图 = Johnson 图 J(L,L/2)：正则度 (L/2)²、连通、无向（时间可逆）；
           Z073 结构定理：上行永远可走、下行会被阻（L=8/10/12 计数 4/32/166，与篇内表相符）
  结构补充：Heawood 的 **3-边着色**（每色 7 边完美匹配 ⟹ 3 个轴，Z021 §8）、
           K₄ 的 **−1 空间在 (123) 循环下分解 1⊕ω⊕ω²**（Z090）、
           K₄ 与 Heawood 的 **Laplacian 谱**（Z002/Z025/Z033/Z102）

用法：
    /usr/bin/python3 zero_sim.py            # 自检 ＋ 写 results/zero_sim.json
    /usr/bin/python3 zero_sim.py --quiet    # 只写文件
    from zero_sim import heawood, spectrum  # 当库用
"""
from __future__ import annotations

import json
import math
import pathlib
import sys
from collections import Counter, deque
from fractions import Fraction

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]          # ZeroGPU/
OUT = ROOT / "results" / "zero_sim.json"

OK: list = []
BAD: list = []
ROWS: list = []


def ck(name: str, cond: bool, detail: str = "") -> bool:
    """登记一条判据（与库内其它脚本同风格）。"""
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    if not QUIET:
        print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ─────────────────────────── L0 结构层 ───────────────────────────

def fano_lines() -> list:
    """Fano 平面的 7 条线：循环差集 {j, j+1, j+3} (mod 7)。"""
    return [tuple(sorted({j % 7, (j + 1) % 7, (j + 3) % 7})) for j in range(7)]


def heawood() -> list:
    """Heawood 图 = Fano 平面的关联图。顶点 0..6 = 点，7..13 = 线（线 j ↔ 点 7+j）。"""
    adj = [set() for _ in range(14)]
    for j, line in enumerate(fano_lines()):
        L = 7 + j
        for p in line:
            adj[p].add(L)
            adj[L].add(p)
    return [sorted(s) for s in adj]


def k4() -> list:
    """K₄（Zero 的通道图；Γ 的唯一刻画，见 doc/gamma_origin_verdict.md）。"""
    return [sorted(v for v in range(4) if v != u) for u in range(4)]


def cycle(n: int) -> list:
    return [sorted({(u - 1) % n, (u + 1) % n}) for u in range(n)]


def path(n: int) -> list:
    return [sorted(v for v in (u - 1, u + 1) if 0 <= v < n) for u in range(n)]


def grid(n: int) -> list:
    """n×n 方格图（Zero 里作"舞台输入"用的那张图，见 Z105/Z115）。"""
    adj = [set() for _ in range(n * n)]
    for y in range(n):
        for x in range(n):
            u = y * n + x
            for dx, dy in ((1, 0), (0, 1)):
                xx, yy = x + dx, y + dy
                if xx < n and yy < n:
                    v = yy * n + xx
                    adj[u].add(v)
                    adj[v].add(u)
    return [sorted(s) for s in adj]


def double_cover(adj: list) -> list:
    """双重覆盖 G □ K₂（Z095 里 K₄□K₂ 的那种构造）：顶点 (v, ε)。"""
    n = len(adj)
    out = [set() for _ in range(2 * n)]
    for v in range(n):
        out[v].add(n + v)
        out[n + v].add(v)
        for w in adj[v]:
            out[v].add(w)
            out[n + v].add(n + w)
    return [sorted(s) for s in out]


# ─────────────────────────── L0 算法层 ───────────────────────────

def degrees(adj: list) -> list:
    return [len(a) for a in adj]


def is_bipartite(adj: list) -> bool:
    n = len(adj)
    color = [-1] * n
    for s in range(n):
        if color[s] != -1:
            continue
        color[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if color[v] == -1:
                    color[v] = 1 - color[u]
                    q.append(v)
                elif color[v] == color[u]:
                    return False
    return True


def girth(adj: list) -> int:
    """围长：逐边去掉后求两端最短距（BFS），取最小环长。无环返回 inf。"""
    best = math.inf
    S = [set(a) for a in adj]
    for u in range(len(adj)):
        for v in S[u]:
            if u > v:
                continue
            dist = {u: 0}
            q = deque([u])
            while q:
                x = q.popleft()
                for y in S[x]:
                    if (x == u and y == v) or (x == v and y == u):
                        continue
                    if y not in dist:
                        dist[y] = dist[x] + 1
                        q.append(y)
            if v in dist:
                best = min(best, dist[v] + 1)
    return best


def adj_matrix(adj: list) -> np.ndarray:
    """邻接矩阵（注意 `adj` 是邻接表，行长度不一，不能直接 np.array）。"""
    n = len(adj)
    A = np.zeros((n, n))
    for i, a in enumerate(adj):
        for j in a:
            A[i, j] = 1.0
    return A


def spectrum(adj: list) -> np.ndarray:
    """邻接矩阵的谱（实对称 ⟹ eigvalsh，升序）。"""
    return np.linalg.eigvalsh(adj_matrix(adj))


def spectrum_rounded(adj: list, nd: int = 6) -> list:
    return sorted(round(float(x), nd) for x in spectrum(adj))


def _bfs_order(adj: list, start: int = 0) -> list:
    """BFS 顶点序（让自同构回溯的约束尽早生效）。"""
    n = len(adj)
    S = [set(a) for a in adj]
    order, seen = [], set()
    for s in [start] + [v for v in range(n) if v != start]:
        if s in seen:
            continue
        seen.add(s); order.append(s)
        q = deque([s])
        while q:
            u = q.popleft()
            for w in sorted(S[u]):
                if w not in seen:
                    seen.add(w); order.append(w); q.append(w)
    return order


def _aut_search(adj: list, fixed: dict | None = None, node_cap: int = 20_000_000):
    """生成所有保邻接的双射 {顶点: 像}（可指定若干强制像）。"""
    n = len(adj)
    S = [set(a) for a in adj]
    order = _bfs_order(adj)
    img, used, nodes = {}, [False] * n, 0
    for v, u in (fixed or {}).items():
        img[v], used[u] = u, True

    def bt(i: int):
        nonlocal nodes
        nodes += 1
        if nodes > node_cap:
            raise RuntimeError("node_cap")
        if i == n:
            yield dict(img)
            return
        v = order[i]
        if v in img:
            yield from bt(i + 1)
            return
        for u in range(n):
            if used[u] or len(S[u]) != len(S[v]):
                continue
            if all((w in S[v]) == (wi in S[u]) for w, wi in img.items()):
                img[v], used[u] = u, True
                yield from bt(i + 1)
                used[u] = False
                del img[v]

    yield from bt(0)


def aut_order(adj: list, fixed: dict | None = None) -> int:
    """自同构群的阶（`fixed` 给强制像，例如 {0: 0} 求顶点稳定子阶）。"""
    return sum(1 for _ in _aut_search(adj, fixed))


def stabilizer_edge_action(adj: list, v0: int = 0) -> int:
    """顶点 v0 的稳定子在 v0 的入射边上作用的**像的大小**（Heawood：6 = 满 S₃）。"""
    base = list(adj[v0])
    perms = set()
    for img in _aut_search(adj, {v0: v0}):
        perms.add(tuple(base.index(img[w]) for w in base))
    return len(perms)


def adj_power(adj: list, k: int) -> np.ndarray:
    return np.linalg.matrix_power(adj_matrix(adj), k)


def adj_power_int(adj: list, k: int) -> list:
    """A^k 的**精确整数**版本（float 在 k≥40 已失真：Tr A^60 ~ 8×10^28）。"""
    n = len(adj)
    M = [[1 if j in adj[i] else 0 for j in range(n)] for i in range(n)]
    R = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    for _ in range(k):
        R = [[sum(R[i][t] * M[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    return R


def closed_walks(adj: list, k: int) -> int:
    """闭道总数 Tr A^k（精确整数）。"""
    R = adj_power_int(adj, k)
    return sum(R[i][i] for i in range(len(R)))


def return_prob(adj: list, v: int, k: int) -> Fraction:
    """简单随机行走 k 步后回到 v 的概率 = (A^k)[v,v] / deg(v)^k（精确有理数）。"""
    d = len(adj[v])
    return Fraction(closed_walks_at(adj, v, k), d ** k)


def closed_walks_at(adj: list, v: int, k: int) -> int:
    """从 v 出发长度 k 的闭道数 (A^k)[v,v]（精确整数）。"""
    return adj_power_int(adj, k)[v][v]


def bipartition(adj: list) -> tuple:
    """二部划分 (A, B)（非二部返回 (None, None)）。"""
    n = len(adj)
    color = [-1] * n
    for s0 in range(n):
        if color[s0] != -1:
            continue
        color[s0] = 0
        q = deque([s0])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if color[v] == -1:
                    color[v] = 1 - color[u]; q.append(v)
                elif color[v] == color[u]:
                    return None, None
    return [v for v in range(n) if color[v] == 0], [v for v in range(n) if color[v] == 1]


def max_matching(part_a: list, adj: list, allowed: set | None = None) -> dict:
    """二部最大匹配（Kuhn 增广路）。`allowed` 是允许使用的边集 (u, v)（u ∈ part_a）。"""
    match_b: dict = {}

    def try_kuhn(u, visited):
        for v in sorted(adj[u]):
            if allowed is not None and (u, v) not in allowed:
                continue
            if v in visited:
                continue
            visited.add(v)
            if v not in match_b or try_kuhn(match_b[v], visited):
                match_b[v] = u
                return True
        return False

    for u in part_a:
        try_kuhn(u, set())
    return match_b


def edge_color_3(adj: list) -> list | None:
    """3-正则二部图的 **3-边着色**（反复取完美匹配；König 线着色定理的构造版）。

    返回 `[(u, v, color), ...]`；非二部或匹配失败返回 None。
    用途：Z021 §8 的"3-边着色给 3 个轴"。
    """
    A, B = bipartition(adj)
    if A is None:
        return None
    remaining = {(u, v) for u in A for v in adj[u]}
    out = []
    for c in range(3):
        m = max_matching(A, adj, remaining)
        if len(m) != len(A):                     # 不是完美匹配 ⟹ 构造失败
            return None
        for v, u in m.items():
            out.append((u, v, c))
            remaining.discard((u, v))
    if remaining:
        return None
    return sorted(out)


def laplacian(adj: list) -> np.ndarray:
    """图 Laplacian L = D − A（Z002/Z025/Z033/Z102 用的对象）。"""
    A = adj_matrix(adj)
    return np.diag(A.sum(axis=1)) - A


def k4_minus_one_decomposition() -> tuple:
    """K₄ 邻接矩阵的 −1 特征空间在 (1 2 3) 循环置换下分解为 1 ⊕ ω ⊕ ω²（Z090）。

    返回（−1 空间的维数, 该 3×3 限制算子的特征值）。
    """
    A = adj_matrix(k4())
    vals, vecs = np.linalg.eigh(A)
    E = vecs[:, np.isclose(vals, -1)]                     # −1 特征空间（3 维）
    P = np.zeros((4, 4), dtype=complex)                   # 置换成 (1 2 3)：1→2→3→1
    P[2, 1] = P[3, 2] = P[1, 3] = 1.0
    P[0, 0] = 1.0
    # 把 P 限制到 E 上（E 是 A 的不变子空间，P 与 A 对易 ⟹ P 保持 E）
    M = np.linalg.solve(E.T @ E, E.T @ (P @ E))
    return E.shape[1], np.linalg.eigvals(M)


# ─────────────────────────── L1 词层 ───────────────────────────

def balanced_words(L: int) -> list:
    """长度 L、和为 0 的词（±1 序列），穷举（小 L 用，L ≤ 20 约 18 万）。"""
    out = []
    for mask in range(1 << L):
        s = 0
        w = []
        for i in range(L):
            step = 1 if (mask >> i) & 1 else -1
            s += step
            w.append(step)
        if s == 0:
            out.append(tuple(w))
    return out


def layer_of(w: tuple) -> int:
    """层高 h(w) = max |前缀和|（与 Z001/Z079/Z047 同口径）。"""
    s = hi = lo = 0
    for step in w:
        s += step
        hi = max(hi, s)
        lo = min(lo, s)
    return max(hi, -lo)


def layer_distribution(L: int) -> Counter:
    """层高分布（DP，状态 =(当前和, 最大前缀, 最小前缀)），L 任意大。"""
    dp = {(0, 0, 0): 1}
    for _ in range(L):
        ndp: dict = {}
        for (s, hi, lo), cnt in dp.items():
            for step in (1, -1):
                s2 = s + step
                key = (s2, max(hi, s2), min(lo, s2))
                ndp[key] = ndp.get(key, 0) + cnt
        dp = ndp
    c: Counter = Counter()
    for (s, hi, lo), cnt in dp.items():
        if s == 0:
            c[max(hi, -lo)] += cnt
    return c


def word_runs(w: tuple) -> int:
    """**极大常值段数**（线性口径）—— Z075 的"段结构" `b_r`。

    口径已钉死（2026-10-15 数值复算）：L=8 上按本定义的分布
    `{2:2, 3:6, 4:18, 5:18, 6:18, 7:6, 8:2}`（Σ=70=C(8,4)）与 Z075 正文表**逐位相符**；
    环上口径（`runs_cyclic`）给 `{2:8, 4:36, 6:24, 8:2}` —— 不符，故段结构是**线性**段数。
    """
    return 1 + sum(1 for i in range(1, len(w)) if w[i] != w[i - 1])


def runs_cyclic(w: tuple) -> int:
    """环上极大常值段数（Z074 的"扇区"是**正段数**，另一个统计量，别与 `word_runs` 混）。"""
    r = word_runs(w)
    if len(w) > 1 and w[0] == w[-1]:
        r -= 1
    return r


def run_distribution(L: int) -> Counter:
    return Counter(word_runs(w) for w in balanced_words(L))


def word_area(w: tuple) -> int:
    """词的**面积** `a(w) = Σ_k 前缀和_k`。

    对平衡词，下面三种写法**完全等价**（分部求和）：`Σ_k prefix_k`、`Σ_{k<L} prefix_k`、
    `−Σ_i i·w_i`（0 起或 1 起索引都行，差一个常数 `Σw = 0`）——2026-10-14 实测四种写法
    在 L ≤ 20 上给出**同一**零面积计数（2/8/58/526/5448），故 Z048 的定义无口径歧义。
    """
    s = 0
    total = 0
    for x in w:
        s += x
        total += s
    return total


def zero_area_counts(Lmax: int = 20) -> dict:
    """各偶长度 L 上**零面积平衡词**的个数（Z048：仅在 L ≡ 0 (mod 4) 时非零）。"""
    out = {}
    for L in range(2, Lmax + 1, 2):
        out[L] = sum(1 for w in balanced_words(L) if word_area(w) == 0)
    return out


def mirror_N_leq(L: int, h: int) -> int:
    """`N_≤(L,h)`：在带 [−h, h] 内从 0 回到 0 的 ±1 游走数（方法镜像法）。

    上下两壁在 ±(h+1)，镜像电荷交替排布，第 j 对给 (−1)^j C(L, L/2 + j(h+1))；
    **j 取遍 ℤ**（负 j 也贡献，例如 L=4,h=1 的 j=−1 给 −C(4,0)）——只取 j≥0 会算错。
    """
    total, step = 0, h + 1
    j = -(L // step + 2)
    while True:
        idx = L // 2 + j * step
        if idx > L:
            break
        if 0 <= idx <= L:
            total += (1 if j % 2 == 0 else -1) * math.comb(L, idx)
        j += 1
    return total


def mirror_shell(L: int, h: int) -> int:
    """Z032(2.1)：壳计数 `b_h = N_≤(L,h) − N_≤(L,h−1)`（h 为层高 `max(M⁺, M⁻)`）。"""
    return mirror_N_leq(L, h) - mirror_N_leq(L, h - 1)


def orbit_count(dist: Counter) -> dict:
    """每层二元轨道数 = b_h/2（Z049：取负 w↦−w 无不动点）。"""
    return {h: b // 2 for h, b in dist.items() if b % 2 == 0}


# ─────────────────── L2 组态图层：补偿移动（Z024/Z072/Z073） ───────────────────

def comp_neighbors(w: tuple) -> list:
    """**补偿移动**（Z072 §2.2）：在任意一对 `(+1, −1)` 之间搬一个单位（保长、保零和）。

    等价说法：平衡词 ↔ "+1 位置集合" `S`（`|S| = L/2`），补偿移动就是"换一个元素"
    —— 故组态图恰是 **Johnson 图 J(L, L/2)**，度为 `|S|·|Sᶜ| = (L/2)²`。
    """
    pos = [i for i, x in enumerate(w) if x == 1]
    neg = [i for i, x in enumerate(w) if x == -1]
    out = []
    for i in pos:
        for j in neg:
            x = list(w)
            x[i], x[j] = -1, 1
            out.append(tuple(x))
    return out


def move_graph_stats(L: int) -> dict:
    """补偿移动图的：度数集合、是否正则、是否连通、均匀测度是否平稳（正则 ⟺ 平稳）。"""
    ws = balanced_words(L)
    degs = Counter(len(comp_neighbors(w)) for w in ws)
    # 连通性（BFS，顶点用元组）
    seen = {ws[0]}
    q = deque([ws[0]])
    while q:
        u = q.popleft()
        for v in comp_neighbors(u):
            if v not in seen:
                seen.add(v)
                q.append(v)
    return {"n": len(ws), "degrees": dict(degs), "regular": len(degs) == 1,
            "degree": next(iter(degs)) if len(degs) == 1 else None,
            "connected": len(seen) == len(ws), "target_degree": (L // 2) ** 2}


def up_move_report(L: int) -> dict:
    """Z073 结构定理：(i) `h < H` 时上行永远可走；(ii) 下行会被阻（计数）。"""
    H = L // 2
    no_up = [w for w in balanced_words(L)
             if layer_of(w) < H and not any(layer_of(x) == layer_of(w) + 1
                                            for x in comp_neighbors(w))]
    blocked_down = [w for w in balanced_words(L)
                    if layer_of(w) > 1 and not any(layer_of(x) == layer_of(w) - 1
                                                   for x in comp_neighbors(w))]
    return {"no_up": len(no_up), "blocked_down": len(blocked_down),
            "blocked_layers": sorted({layer_of(w) for w in blocked_down}),
            "reversible": all(w in comp_neighbors(x) for w in balanced_words(L)
                              for x in comp_neighbors(w))}


# ─────────────────────────── 耦合几何（Z021） ───────────────────────────

def simplex_vectors(D: int) -> np.ndarray:
    """正则 D-单纯形的规范构造：u_a = e_a − 1/m·1（m = D+1），形状 (m, m)。"""
    m = D + 1
    U = np.eye(m) - np.ones((m, m)) / m
    return U


def gram(D: int) -> np.ndarray:
    U = simplex_vectors(D)
    return U @ U.T


def rho(D: int) -> Fraction:
    """成对夹角余弦的精确值 −1/D（D = m−1）。"""
    return Fraction(-1, D)


def peak_argmax(q: Fraction, Dmax: int = 40) -> int:
    """F_D = C(D,2)·q^D 的极大点（精确有理数比较）。"""
    best, bestv = 0, Fraction(0)
    for D in range(Dmax + 1):
        v = Fraction(math.comb(D, 2)) * q ** D
        if v > bestv:
            best, bestv = D, v
    return best


# ─────────────────────────── 自检 ───────────────────────────

def sim_selftest() -> dict:
    print("== L0 结构：Heawood（Fano 关联图） ==")
    H = heawood()
    lines = fano_lines()
    pairs = Counter()
    for ln in lines:
        for a in range(3):
            for b in range(a + 1, 3):
                pairs[tuple(sorted((ln[a], ln[b])))] += 1
    ck("Fano：7 条线、每对点恰在一条线上", len(lines) == 7 and len(pairs) == 21
       and set(pairs.values()) == {1})
    ck("Heawood：|V|=14、度 3、|E|=21",
       len(H) == 14 and set(degrees(H)) == {3} and sum(degrees(H)) // 2 == 21)
    ck("Heawood：二部（点的部分 ⟷ 线的部分）", is_bipartite(H))
    ck("Heawood：围长 = 6", girth(H) == 6, f"girth={girth(H)}")

    print("== L0 谱：{3, √2⁶, (−√2)⁶, −3}、零化多项式、闭道闭式、质量隙 ==")
    sp = spectrum(H)
    srt = spectrum_rounded(H)
    target = sorted([3.0] + [round(math.sqrt(2), 6)] * 6 + [round(-math.sqrt(2), 6)] * 6 + [-3.0])
    ck("谱 = {3, (√2)×6, (−√2)×6, −3}", srt == target, f"{sorted(set(srt))}")
    A = adj_matrix(H)
    I = np.eye(14)
    Z = (A - 3 * I) @ (A + 3 * I) @ (A @ A - 2 * I)
    ck("零化多项式 (A−3)(A+3)(A²−2)=0", float(np.abs(Z).max()) < 1e-8,
       f"max|Z|={float(np.abs(Z).max()):.2e}")
    n_aut = aut_order(H)
    ck("|Aut(Heawood)| = 336", n_aut == 336, f"|Aut|={n_aut}")
    n_stab = aut_order(H, fixed={0: 0})
    ck("顶点稳定子阶 = 24", n_stab == 24, f"|Stab|={n_stab}")
    ck("顶点稳定子在 3 条入射边上的作用 = 6（满 S₃）",
       stabilizer_edge_action(H, 0) == 6, f"像的个数={stabilizer_edge_action(H, 0)}")
    exact = all(7 * Fraction(closed_walks(H, 2 * m), 14) - 3 ** (2 * m) == 6 * 2 ** m
                for m in range(0, 31))
    ck("闭道闭式 7a_{2m} − 3^{2m} = 6·2^m（a_k = Tr A^k/14，m ≤ 30）", exact)
    mgap = -math.log(math.sqrt(2) / 3)
    ck("质量隙 m = −log(√2/3) = 0.7520386984", abs(mgap - 0.7520386984) < 1e-9,
       f"m={mgap:.10f}，ξ={1/mgap:.10f}")
    ck("返回概率偶数步精确式 p(2k) = (9^k + 6·2^k)/(7·9^k)（k ≤ 12）",
       all(return_prob(H, 0, 2 * k) == Fraction(9 ** k + 6 * 2 ** k, 7 * 9 ** k)
           for k in range(0, 13)))
    ck("偏差衰减率 = (√2/3)^{2k} ⟹ 质量隙 m = −log(√2/3)（k ≤ 12）",
       all(return_prob(H, 0, 2 * k) - Fraction(1, 7)
           == Fraction(6, 7) * Fraction(2, 9) ** k for k in range(0, 13)))
    ck("奇数步返回概率 = 0（二部图无奇闭道，k ≤ 11）",
       all(return_prob(H, 0, 2 * k + 1) == 0 for k in range(0, 12)))

    print("== L1 词层：b_h 分布、镜像闭式、二元轨道（Z032/Z047/Z049） ==")
    for L in (4, 8, 12, 16):
        dist = layer_distribution(L)
        brute = Counter(layer_of(w) for w in balanced_words(L)) if L <= 12 else None
        ck(f"L={L}：Σ b_h = C(L,L/2)", sum(dist.values()) == math.comb(L, L // 2),
           f"Σ={sum(dist.values())}")
        ck(f"L={L}：b_h 恒偶", all(b % 2 == 0 for b in dist.values()))
        ck(f"L={L}：b_H = 2（顶层）", dist[max(dist)] == 2)
        ck(f"L={L}：Z032(2.1) 镜像闭式与 DP 逐层相符",
           all(mirror_shell(L, h) == b for h, b in dist.items()),
           f"b={dict(sorted(dist.items()))}")
        if brute is not None:
            ck(f"L={L}：DP 与穷举逐层相符", dict(brute) == dict(dist))
    ck("L=4：b = (层 1 有 4 词、层 2 有 2 词)，层权 2/3 : 1/3",
       dict(layer_distribution(4)) == {1: 4, 2: 2}, f"{dict(layer_distribution(4))}")
    ck("L=16：b_h 与 Z032 数值逐位相符（顶层 2）",
       dict(layer_distribution(16)) == {h: mirror_shell(16, h)
                                        for h in sorted(layer_distribution(16))})

    print("== Z074/Z075：段结构（线性极大常值段数） ==")
    d8 = run_distribution(8)
    ck("L=8：段分布 {2:2, 3:6, 4:18, 5:18, 6:18, 7:6, 8:2}（Z075 表）",
       dict(sorted(d8.items())) == {2: 2, 3: 6, 4: 18, 5: 18, 6: 18, 7: 6, 8: 2},
       f"{dict(sorted(d8.items()))}")
    ck("L=8：Σ_r b_r = C(8,4) = 70", sum(d8.values()) == 70)
    ck("L=4：段分布 {2:2, 3:2, 4:2}，Σ = C(4,2) = 6",
       dict(sorted(run_distribution(4).items())) == {2: 2, 3: 2, 4: 2})
    ck("口径判据：环上段数分布与 Z075 表**不符**（故段结构取线性口径）",
       dict(sorted(Counter(runs_cyclic(w) for w in balanced_words(8)).items()))
         != {2: 2, 3: 6, 4: 18, 5: 18, 6: 18, 7: 6, 8: 2})

    print("== Z024/Z072/Z073：补偿移动图（= Johnson 图 J(L,L/2)） ==")
    for L in (4, 6, 8, 10):
        st = move_graph_stats(L)
        ck(f"L={L}：补偿移动图正则、度 = (L/2)² = {st['target_degree']}",
           st["regular"] and st["degree"] == st["target_degree"], f"{st['degrees']}")
        ck(f"L={L}：连通（Z024 的不可约性）", st["connected"], f"顶点 {st['n']}")
    # Z073 表：L=8/10/12 的受阻下行计数 = 4/32/166（与篇内表逐位相符）
    # ★ 精化（本条自检抓出，2026-10-15）：受阻词**并非全在 h=2** —— L=12 起层 3 也出现
    #   （L=12：2 个；L=14：26 个），故篇内"下行被阻"的正确读法是"层 ≥ 2 且随 L 增长"。
    expect = {8: 4, 10: 32, 12: 166}
    for L in (8, 10, 12):
        rep = up_move_report(L)
        ck(f"L={L}：上行永远可走（h<H 时必有一邻升高 1）", rep["no_up"] == 0)
        ck(f"L={L}：下行会被阻，计数 = {expect[L]}（Z073 表）", rep["blocked_down"] == expect[L],
           f"层分布 {rep['blocked_layers']}")
        ck(f"L={L}：边关系对称（无向 ⟹ 时间可逆）", rep["reversible"])
    rep14 = up_move_report(14)
    ck("L=14：上行永远可走、下行受阻 714（层 2 有 688、层 3 有 26 —— 非全在 h=2）",
       rep14["no_up"] == 0 and rep14["blocked_down"] == 714, f"{rep14}")

    print("== Z048：零面积词（面积口径的四种写法等价） ==")
    za = zero_area_counts(20)
    ck("零面积计数 L=4/8/12/16/20 = 2/8/58/526/5448",
       [za[L] for L in (4, 8, 12, 16, 20)] == [2, 8, 58, 526, 5448], f"{za}")
    ck("L ≡ 2 (mod 4) ⟹ 零面积词数为 0（L ≤ 20）",
       all(za[L] == 0 for L in range(2, 21, 2) if L % 4 == 2))
    # 分部求和：对平衡词 Σ_{k<L} prefix_k = −Σ_i i·w_i（0 起索引），故三种写法同值
    ck("三种面积写法在平衡词上等价（L ≤ 12 逐词）",
       all(sum(i * w[i] for i in range(L)) == -word_area(w)
           and sum((i + 1) * w[i] for i in range(L)) == -word_area(w)
           for L in (4, 8, 12) for w in balanced_words(L)))

    print("== 耦合几何：正则 D-单纯形（Z021） ==")
    for D in range(1, 7):
        m = D + 1
        U = simplex_vectors(D)
        G = gram(D)
        off = G[0, 1]
        ck(f"D={D}：零和（Σu_a=0）", float(np.abs(U.sum(axis=0)).max()) < 1e-12)
        ck(f"D={D}：等长（‖u_a‖²=1−1/m）",
           float(np.abs(np.diag(G) - (1 - 1 / m)).max()) < 1e-12)
        ck(f"D={D}：ρ = ⟨u_a,u_b⟩ = −1/m = −1/(D+1)（归一化后 = −1/D）",
           abs(off + 1 / m) < 1e-12, f"ρ={off:.6f}")
        ck(f"D={D}：Gram 秩 = D", int(np.linalg.matrix_rank(G, tol=1e-9)) == D,
           f"rank={int(np.linalg.matrix_rank(G, tol=1e-9))}")
    ck("四维窗 (1/2, 3/5) 含 q=5/9", Fraction(1, 2) < Fraction(5, 9) < Fraction(3, 5))
    ck("q=5/9：F_D = C(D,2)q^D 的唯一峰 D=4", peak_argmax(Fraction(5, 9)) == 4)
    ck("q=5/9：不落在三维窗 (1/3,1/2)",
       not (Fraction(1, 3) < Fraction(5, 9) < Fraction(1, 2)))
    ck("D=3：夹角余弦 = −1/3（四面体角 arccos(−1/3)）", rho(3) == Fraction(-1, 3))

    print("== L0 结构：K₄ 与双重覆盖（Z089/Z095） ==")
    K = k4()
    ck("K₄：4 顶点、度 3、|E|=6", set(degrees(K)) == {3} and sum(degrees(K)) // 2 == 6)
    ck("K₄：围长 3（故层对齐需换图）", girth(K) == 3)
    DC = double_cover(K)
    spdc = sorted(round(float(x), 6) for x in spectrum(DC))
    ck("K₄□K₂：8 顶点、4-正则",
       len(DC) == 8 and set(degrees(DC)) == {4})
    ck("K₄□K₂：谱 = {4, 2, 0³, (−2)³}",
       spdc == sorted([4.0, 2.0, 0.0, 0.0, 0.0, -2.0, -2.0, -2.0]), f"{spdc}")

    print("== Z021 §8 / Z090 / Z002：3-边着色、−1 空间分解、Laplacian ==")
    ec = edge_color_3(H)
    ok_ec = ec is not None and len(ec) == 21
    if ok_ec:
        by_color = {c: [(u, v) for u, v, cc in ec if cc == c] for c in range(3)}
        ends = {c: [x for e in es for x in e] for c, es in by_color.items()}
        ok_ec = all(len(es) == 7 and len(set(ends[c])) == 14 for c, es in by_color.items())
    ck("Heawood 有 3-边着色（每色 7 边 = 完美匹配，无共端冲突）", ok_ec,
       f"边数={len(ec) if ec else None}")
    dimE, ev = k4_minus_one_decomposition()
    omega = np.exp(2j * np.pi / 3)
    ck("K₄：−1 空间维数 = 3 且在 (123) 循环下分解为 {1, ω, ω²}",
       dimE == 3 and np.allclose(sorted(ev, key=lambda z: (round(z.real, 6), round(z.imag, 6))),
                                 sorted([1 + 0j, omega, np.conj(omega)],
                                        key=lambda z: (round(z.real, 6), round(z.imag, 6))),
                                 atol=1e-9),
       f"特征值={np.round(sorted(ev, key=lambda z: z.real), 6)}")
    lapk = sorted(round(float(x), 6) for x in np.linalg.eigvalsh(laplacian(k4())))
    ck("K₄ Laplacian 谱 = {0, 4, 4, 4}", lapk == [0.0, 4.0, 4.0, 4.0], f"{lapk}")
    laph = np.linalg.eigvalsh(laplacian(H))
    ck("Heawood Laplacian 谱 = {0, 3−√2×6, 3+√2×6, 6}",
       np.allclose(sorted(laph), sorted([0.0] + [3 - math.sqrt(2)] * 6 + [3 + math.sqrt(2)] * 6 + [6.0]),
                   atol=1e-9))

    res = {
        "module": "zero_sim",
        "assertions_ok": len(OK),
        "assertions_bad": len(BAD),
        "failed": BAD,
        "rows": ROWS,
        "heawood": {"|V|": 14, "degree": 3, "girth": girth(H), "|E|": 21,
                    "spectrum": sorted(set(srt)), "|Aut|": n_aut,
                    "vertex_stabilizer": n_stab,
                    "mass_gap": mgap, "xi": 1 / mgap},
        "layers": {L: dict(sorted(layer_distribution(L).items()))
                   for L in (4, 8, 12, 16)},
        "zero_area_counts": za,
        "run_distribution": {str(L): dict(sorted(run_distribution(L).items())) for L in (4, 8)},
        "move_graph": {str(L): move_graph_stats(L) for L in (4, 6, 8, 10)},
        "up_down": {str(L): up_move_report(L) for L in (8, 10, 12)},
        "edge_color_3": ec,
        "k4_minus_one": {"dim": int(dimE), "eigenvalues": [[float(z.real), float(z.imag)] for z in ev]},
    }
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n判据：{len(OK)} 过 / {len(BAD)} 败  →  {OUT}")
    if BAD:
        print("失败项：")
        for b in BAD:
            print(f"  ✗ {b}")
    return res


QUIET = False

if __name__ == "__main__":
    QUIET = "--quiet" in sys.argv
    r = sim_selftest()
    sys.exit(1 if r["assertions_bad"] else 0)