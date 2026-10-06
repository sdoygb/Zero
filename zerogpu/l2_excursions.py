"""
l2_excursions.py --- 零乱动的**远足分解** → 多个 L2 区域（语料原生机制）

**机制（来自语料 `UNIVERSE.md` §2，零参数）**
    Z0①「零从不停」⇒ 走回零后**从零继续** ⇒ 世界线自动断成一串**远足 (excursion)**；
    每条远足 = 一个"从零重新开始"的区域 = 一个**区域性 L1/L2**。
    远足数分布的闭式 $E(u)=1-\sqrt{1-4u}$，$E(1/4)=1$（临界分支过程）。

**Γ 补上之后的新问题（本脚本）**
    没有 Γ 时，远足序列 i.i.d. ⇒ 微正则理想气体 ⇒「无关联、无散射、无相互作用」（语料已证）。
    有了 Γ，**远足会在顶点上相撞**：同一个顶点被多条远足踩到 ⟹ 整数重数合并 ⟹ **相互作用**。
    于是"多个 L2"不再是形式分解，而是**由碰撞图（空间重叠）决定的连通块**：

        远足 = 节点；两条远足相邻 ⟺ 它们的足迹在**非原点**顶点上相交
        连通块 = **一个区域（一个 L2）**

    —— 分区**不手工**、**无自由参数**（除 Γ 与总长度 $n$）。

**测量**
    · 区域数、最大区域占比（渗流式问题）、区域大小分布
    · 每个区域：远足条数、总原时、足迹顶点数、回转半径、诱导子图的有效维数 d_s
    · 与 D 的关系（D=1,2,3,4 同一套代码）

**两条读法（HANDOFF 坑 4）**：本脚本用**不吸收**读法（Z0① 字面）——回到原点后继续走。

等级：【数值证据】。D、L、n 是输入；Γ 的来源（L2-LATTICE-ORIGIN）仍开放。

用法：/usr/bin/python3 l2_excursions.py
输出：results/l2_excursions.json
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def simulate(D, L, n_steps, n_walk, seed=0):
    """
    在 ℤ^D 环面上做**不吸收**随机走动（每条可用边等概率）。
    返回 per-walker 的远足列表：[(length, visited_sites_set), ...]（不含原点）。
    """
    rng = np.random.default_rng(seed)
    origin = 0
    excur = []          # excur[w] = [(len, frozenset(sites != origin))]
    for w in range(n_walk):
        pos = np.zeros(D, np.int64)
        step = rng.integers(0, 2 * D, size=n_steps)
        axis = step // 2
        sign = np.where(step % 2 == 0, 1, -1)
        cur_len = 0
        sites = set()
        lst = []
        for t in range(n_steps):
            a = int(axis[t])
            pos[a] = (pos[a] + int(sign[t])) % L
            cur_len += 1
            key = 0
            for d in range(D):
                key = key * L + int(pos[d])
            if key == origin:
                lst.append((cur_len, frozenset(sites)))
                cur_len = 0
                sites = set()
            else:
                sites.add(key)
        if cur_len > 0:                      # 收尾的未闭合段
            lst.append((cur_len, frozenset(sites)))
        excur.append(lst)
    return excur


def region_stats(excur, L, D, min_size=2):
    """碰撞图（非原点顶点相交）的连通块 = 区域。"""
    flat = []
    for w, lst in enumerate(excur):
        for e, (ln, sites) in enumerate(lst):
            flat.append((w, e, ln, sites))
    n_ex = len(flat)
    parent = list(range(n_ex))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    site_map: dict[int, int] = {}
    for i, (w, e, ln, sites) in enumerate(flat):
        for s in sites:
            j = site_map.get(s)
            if j is None:
                site_map[s] = i
            else:
                union(j, i)
    comp: dict[int, list[int]] = {}
    for i in range(n_ex):
        comp.setdefault(find(i), []).append(i)

    rows = []
    for r, members in sorted(comp.items(), key=lambda kv: -len(kv[1])):
        ln_tot = sum(flat[i][2] for i in members)
        support = set()
        for i in members:
            support |= flat[i][3]
        rows.append({"n_excursions": len(members), "total_length": int(ln_tot),
                     "support_sites": len(support), "mean_excursion_len":
                         round(float(ln_tot) / max(1, len(members)), 3),
                     "_support": support})
    return {"n_excursions": n_ex, "regions": rows}


def collision_stats(excur):
    """
    碰撞强度：同一顶点被 k 条远足踩到 ⟹ C(k,2) 对相互作用。
    语料在**无 Γ** 时证过"远足序列 i.i.d. ⇒ 无相互作用"；这里量 Γ 补回多少。
    """
    site_count = {}
    for lst in excur:
        for _ln, sites in lst:
            for s in sites:
                site_count[s] = site_count.get(s, 0) + 1
    n_ex = sum(len(l) for l in excur)
    pairs = sum(k * (k - 1) // 2 for k in site_count.values())
    colliding_sites = sum(1 for k in site_count.values() if k >= 2)
    # 相邻远足（同一条世界线上前后两条）是否相撞 —— 语料说无 Γ 时它们独立
    adj_pairs = 0
    adj_collide = 0
    for lst in excur:
        for i in range(len(lst) - 1):
            adj_pairs += 1
            if lst[i][1] & lst[i + 1][1]:
                adj_collide += 1
    return {"excursions": int(n_ex),
            "visited_site_slots": int(sum(site_count.values())),
            "distinct_sites": int(len(site_count)),
            "colliding_pairs": int(pairs),
            "pairs_per_excursion": round(pairs / max(1, n_ex), 4),
            "colliding_site_frac": round(colliding_sites / max(1, len(site_count)), 4),
            "adjacent_pairs": int(adj_pairs),
            "adjacent_collision_frac": round(adj_collide / max(1, adj_pairs), 4)}


def support_dimension(support, L, D, n_sample=200000, r_max=14, seed=0):
    """区域足迹的**生长维数**：从区域质心附近出发，vol(r) ~ r^{d_g}。"""
    if len(support) < 20:
        return None
    rng = np.random.default_rng(seed)
    coords = np.array([[ (s // (L ** (D - 1 - d))) % L for d in range(D)] for s in support], np.int64)
    if len(coords) > n_sample:
        coords = coords[rng.choice(len(coords), n_sample, replace=False)]
    c0 = coords.mean(0)
    seen = {tuple(c) for c in coords}
    start = tuple(np.round(c0).astype(np.int64) % L)
    if start not in seen:
        start = tuple(coords[0])
    frontier = {start}
    visited = {start}
    vol = []
    for _ in range(r_max):
        nxt = set()
        for c in frontier:
            for d in range(D):
                for sg in (1, -1):
                    cc = list(c)
                    cc[d] = (cc[d] + sg) % L
                    cc = tuple(cc)
                    if cc in seen and cc not in visited:
                        visited.add(cc)
                        nxt.add(cc)
        frontier = nxt
        vol.append(len(visited))
        if not frontier:
            break
    if len(vol) < 4:
        return None
    r = np.arange(1, len(vol) + 1, dtype=float)
    v = np.array(vol, float)
    s0 = max(1, int(0.4 * len(r)))
    return float(np.polyfit(np.log(r[s0:]), np.log(v[s0:]), 1)[0])


def main(argv):
    cfgs = [(1, 4001, 400, 3000), (2, 201, 400, 3000), (3, 61, 400, 1500), (4, 27, 400, 800)]
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "reading": "不吸收（Z0① 字面）",
           "region_definition": "远足 = 节点；两条远足相邻 ⟺ 足迹在非原点顶点相交；连通块 = 区域",
           "runs": []}
    print("=" * 88)
    print(f"{'D':>2} {'L':>5} {'steps':>6} {'walkers':>8} {'远足数':>9} {'区域数':>7} "
          f"{'最大区域占比':>12} {'top5 区域(远足数)':>26} {'d_g(top)':>9}")
    for D, L, n_steps, n_walk in cfgs:
        t0 = time.time()
        ex = simulate(D, L, n_steps, n_walk)
        st = region_stats(ex, L, D)
        cst = collision_stats(ex)
        regs = st["regions"]
        tot = max(1, st["n_excursions"])
        top = [r["n_excursions"] for r in regs[:5]]
        dg = support_dimension(regs[0]["_support"], L, D) if regs else None
        row = {"D": D, "L": L, "n_steps": n_steps, "n_walkers": n_walk,
               "n_excursions": st["n_excursions"], "n_regions": len(regs),
               "collisions": cst,
               "largest_region_excursion_frac": round(top[0] / tot, 4) if top else None,
               "regions": [{k: v for k, v in r.items() if k != "_support"} for r in regs[:40]],
               "d_growth_largest_region": None if dg is None else round(dg, 4),
               "seconds": round(time.time() - t0, 1)}
        out["runs"].append(row)
        print(f"     碰撞: 相撞顶点占比={cst['colliding_site_frac']} "
              f"每远足对数={cst['pairs_per_excursion']} 相邻远足相撞率={cst['adjacent_collision_frac']}")
        print(f"{D:>2} {L:>5} {n_steps:>6} {n_walk:>8} {st['n_excursions']:>9} {len(regs):>7} "
              f"{row['largest_region_excursion_frac']:>12} {str(top):>26} "
              f"{('%.3f' % dg) if dg else '   -  ':>9}  [{row['seconds']}s]", flush=True)
    with open(os.path.join(OUT, "l2_excursions.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 88)
    print("写出 results/l2_excursions.json")
    print("\n读法：区域数 > 1 ⟺ 远足在空间上**没有**全部相撞（Γ 引起的相互作用是局部的）。")
    print("      最大区域占比 → 1 ⟺ 所有远足连成一片（一个 L2）。")


if __name__ == "__main__":
    main(sys.argv[1:])
