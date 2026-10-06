"""
z0_loop.py --- 零参数**自持环**：封闭 L0 → 闭合 → 记录 → 毁灭/重播种 → 回灌

**为什么需要它**：到本轮为止所有模拟都是**单向枚举**（L0→L1 只写不回灌），
枚举完就结束 —— 世界不会自己接着活。用户的图景要求闭环。

**环的定义（尽量零参数）**

```
        ┌──────────────── 回灌（只读保留的记录） ────────────────┐
        ↓                                                       │
   活动层 E（走者集合，各带 (起点 u, 当前位置 v)）                  │
        │ 每一步：沿 Γ 走一条边（Z2：所有可用边都实例化）           │
        ↓                                                       │
   v == u ？ ── 否 ──→ 继续走                                     │
        │ 是（**闭合事件**，Σw=0）                                │
        ↓                                                       │
   写记录（远足长度等）→ 毁灭（移除该走者）→ **重播种 r 条** ──────┘
```

**唯一的外加量 = 重播种强度 $r$**（每次闭合新播几条）。但 $r$ 的地位是**可判定的**，不是自由的：

$$
\text{每条闭合移除 }1\text{ 条、重播 }r\text{ 条}
\ \Longrightarrow\
\begin{cases}
r<1 & \text{活动量指数**熄灭**} \\
r=1 & \text{活动量**精确守恒**（Z0① 字面：闭合不停，继续走）} \\
r>1 & \text{活动量**指数爆炸**}（Z2 分支读法：}r=\\deg v>1\text{）}
\end{cases}
$$

所以本脚本给出**相图**，并指出：**Z0① 的字面读法恰好落在 $r=1$（临界、不调参自持）**，
而 **Z2 的全分支读法 $r=\\deg v>1$ 落在爆炸相**（这就是"体积 $2^n$"的来源）。
两者正是语料自己的 L2-a（保守窗）／L2-b（产生分支）两个亚层（`R95:32`）。

**稳态里测什么**（这才是本脚本的价值所在）
  ① 活动量 A(t)：核验 $r=1$ 守恒、$r\\ne1$ 熄灭/爆炸
  ② **稳态远足长度分布** $\rho(L)$：有没有极点（粒子）还是只有支割
  ③ **稳态下的 L2 区域数**：远足碰撞图的连通块（复用 `l2_excursions.region_stats`）
  ④ 记录层的**增长率**（每秒记几条）——有没有出现"时钟"

**等级**：全部【数值证据】；$r=1$ 的守恒是【精确】（计数论证）。
用法：/usr/bin/python3 z0_loop.py
输出：results/z0_loop.json
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

RECORDS = []


def rec(sec, obj, readout, x, y, meta=None, note=None, kind="curve", maxpts=64):
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) > maxpts:
        idx = np.unique(np.r_[np.linspace(0, len(x) - 1, maxpts - 1).astype(int), len(x) - 1])
        x, y = x[idx], y[idx]
    RECORDS.append({"id": f"{sec}:{obj}:{readout}", "sector": sec, "object": obj,
                    "readout": readout, "kind": kind,
                    "x": [float(v) for v in x], "y": [float(v) for v in y],
                    "meta": meta or {}, **({"note": note} if note else {})})


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def run_loop(A, n0, steps, r=1.0, seed=0, collect_support=False):
    """
    自持环。返回统计。
      A            : Γ 的邻接（CSR，无自环）
      n0           : 初始走者数
      r            : 每次闭合重播种几条（可为小数：按期望取整）
      collect_support : 是否记录远足的足迹（做碰撞图用；费内存）
    """
    N = A.shape[0]
    indptr, indices = A.indptr, A.indices
    deg = np.diff(indptr).astype(np.int64)
    rng = np.random.default_rng(seed)

    start = rng.integers(0, N, n0).astype(np.int64)
    pos = start.copy()
    age = np.zeros(n0, np.int64)              # 每条走者的当前远足长度
    act_hist, closed_hist = [], []
    exlen = []                                # 闭合的远足长度
    supports = [] if collect_support else None
    visited = [set([int(s)]) for s in start] if collect_support else None
    n_closed_total = 0

    for t in range(steps):
        # --- 走一步（保守：每条走者正好走一条边；Z2 的"全分支"见 r 的讨论）---
        m = len(pos)
        off = (rng.random(m) * deg[pos]).astype(np.int64)
        pos = indices[indptr[pos] + off]
        age += 1
        if collect_support:
            for i in range(m):
                visited[i].add(int(pos[i]))
        # --- 闭合：回到起点 ---
        cl = np.nonzero(pos == start)[0]
        if len(cl):
            n_closed_total += len(cl)
            exlen.extend(age[cl].tolist())
            if collect_support:
                for i in cl:
                    supports.append((int(age[i]), frozenset(visited[i])))
            if r <= 0:
                keep = np.setdiff1d(np.arange(m), cl, assume_unique=False)
                start, pos, age = start[keep], pos[keep], age[keep]
                if collect_support:
                    visited = [visited[i] for i in keep]
            else:
                # 重播种：把闭合者重置为"起点=当前位置、年龄 0"，并按 r 追加新走者
                age[cl] = 0
                if collect_support:
                    for i in cl:
                        visited[i] = {int(pos[i])}
                k = int(np.floor(r)) + (1 if rng.random() < (r - np.floor(r)) else 0)
                if k > 1:
                    extra = np.repeat(cl, k - 1)
                    start = np.concatenate([start, pos[extra]])
                    pos = np.concatenate([pos, pos[extra]])
                    age = np.concatenate([age, np.zeros(len(extra), np.int64)])
                    if collect_support:
                        visited += [{int(pos[i])} for i in extra]
                elif k == 0:
                    keep = np.setdiff1d(np.arange(len(pos)), cl, assume_unique=False)
                    start, pos, age = start[keep], pos[keep], age[keep]
                    if collect_support:
                        visited = [visited[i] for i in keep]
        act_hist.append(len(pos))
        closed_hist.append(n_closed_total)
        if len(pos) == 0:
            act_hist.extend([0] * (steps - t - 1))
            closed_hist.extend([n_closed_total] * (steps - t - 1))
            break
        if len(pos) > 4_000_000:
            act_hist.extend([len(pos)] * (steps - t - 1))
            closed_hist.extend([n_closed_total] * (steps - t - 1))
            break

    return {"r": r, "n0": int(n0), "steps": steps,
            "activity": act_hist, "closed_cum": closed_hist,
            "excursion_lengths": exlen, "supports": supports}


def _first_return_1d(m):
    """1D 无限直线首次返回在 2m 步的精确概率：C(2m,m)/((2m-1)4^m)。"""
    import math
    return math.comb(2 * m, m) / ((2 * m - 1) * 4.0 ** m)


def excursion_structure(exlen, D=None, Lmax=None):
    """稳态远足长度分布的结构：$\rho(L)$ 的标度（处理二部奇偶、限拟合范围）。"""
    if len(exlen) < 50:
        return None
    Lall = np.array(exlen, float)
    Lall = Lall[Lall > 0]
    mean_all = float(Lall.mean()) if len(Lall) else float("nan")
    L = Lall[Lall <= Lmax] if Lmax else Lall
    if len(L) < 50:
        return None
    vals, cnt = np.unique(L, return_counts=True)
    # 结构分类：幂律？指数？
    m = cnt > 0
    out = {"n_excursions": int(len(L)), "mean_len_in_fit_range": round(float(L.mean()), 4),
           "mean_len_all": round(mean_all, 4), "fit_Lmax": Lmax,
           "max_len": int(L.max()), "lengths": vals.astype(int).tolist(),
           "counts": cnt.astype(int).tolist()}
    if len(vals) >= 6:
        lv, lc = np.log(vals[m]), np.log(cnt[m])
        b1 = np.polyfit(lv, lc, 1)
        r2_pow = 1 - ((lc - np.polyval(b1, lv)) ** 2).sum() / max(((lc - lc.mean()) ** 2).sum(), 1e-30)
        b2 = np.polyfit(vals[m], lc, 1)
        r2_exp = 1 - ((lc - np.polyval(b2, vals[m])) ** 2).sum() / max(((lc - lc.mean()) ** 2).sum(), 1e-30)
        out["power_law"] = {"exponent": round(float(b1[0]), 4), "r2": round(float(r2_pow), 5)}
        out["exponential"] = {"rate": round(float(b2[0]), 4), "r2": round(float(r2_exp), 5)}
        out["preferred"] = "power-law" if r2_pow > r2_exp else "exponential"
        out["theory_1D_exponent"] = -1.5
        # 与 1D 精确首返律逐点对照（只在 D=1 时有意义）
        if D == 1:
            even = (vals.astype(int) % 2 == 0)
            vv, cc = vals[even].astype(int), cnt[even]
            if len(vv) >= 5:
                pred = np.array([_first_return_1d(int(v) // 2) for v in vv])
                pred = pred / pred.sum() * cc.sum()
                # 卡方式偏差（相对）
                rel = np.abs(cc - pred) / np.maximum(pred, 1e-12)
                out["vs_exact_first_return_1D"] = {
                    "n_points": int(len(vv)), "median_rel_dev": round(float(np.median(rel)), 4),
                    "note": "与无限直线首返律 C(2m,m)/((2m-1)4^m) 的逐点相对偏差"}
    return out


def main(argv):
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "loop": "闭合→写记录→毁灭→重播种 r 条→继续（Z0① 字面：闭合不停）",
           "criticality": "每条闭合移除 1 条、重播 r 条 ⟹ r*=1 是精确的临界点（计数论证，非拟合）",
           "runs": []}

    # ---------- ① 相图：r 扫描 ----------
    print("=" * 92)
    print("① 相图：重播种强度 r 的活动量长期行为（Γ = ℤ² 环面 16×16，初始 200 条走者）")
    A = torus_adj(2, 16)
    steps = 400
    rows = []
    for r in (0.0, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0):
        res = run_loop(A, 200, steps, r=r, seed=1)
        act = np.array(res["activity"], float)
        tail = act[len(act) // 2:]
        growth = (np.log(max(tail[-1], 1e-9)) - np.log(max(tail[0], 1e-9))) / max(1, len(tail))
        rows.append({"r": r, "A_start": int(act[0]), "A_mid": int(act[len(act) // 2]),
                     "A_end": int(act[-1]), "log_growth_per_step": round(float(growth), 6),
                     "closed_total": int(res["closed_cum"][-1])})
        print(f"   r={r:<4}: A(0)={rows[-1]['A_start']:>7} A(mid)={rows[-1]['A_mid']:>9} "
              f"A(end)={rows[-1]['A_end']:>10}  对数增长率/步={rows[-1]['log_growth_per_step']:+.5f}  "
              f"闭合总数={rows[-1]['closed_total']}")
    out["phase"] = rows
    rec("LOOP", "重播种相图", "长期活动量 A(end) vs r",
        [x["r"] for x in rows], [x["A_end"] for x in rows], {"detail": rows},
        note="r<1 熄灭 / r=1 守恒 / r>1 爆炸；r*=1 精确")

    # ---------- ② 稳态结构：远足长度分布（D=1,2,3 与原生 𝒢_T）----------
    print("=" * 92)
    print("② 稳态远足长度分布 ρ(L)（r=1，临界稳态）")
    out["steady"] = []
    for D, L in ((1, 64), (2, 16), (3, 8)):
        A = torus_adj(D, L)
        res = run_loop(A, 400, 20000, r=1.0, seed=7)
        st = excursion_structure(res["excursion_lengths"], D=D, Lmax=max(8, (L ** D) // 2))
        act = np.array(res["activity"], float)
        entry = {"graph": f"torus_D{D}_L{L}", "D": D, "N": L ** D,
                 "activity_min": int(act.min()), "activity_max": int(act.max()),
                 "activity_const": bool(act.min() == act.max()),
                 "records": len(res["excursion_lengths"]),
                 "record_rate_per_step": round(len(res["excursion_lengths"]) / len(act), 4),
                 "structure": st}
        out["steady"].append(entry)
        print(f"   D={D} L={L:>3} N={L**D:>5}: 活动量 [{int(act.min())}, {int(act.max())}] "
              f"恒定={entry['activity_const']}  闭合 {entry['records']} 次 "
              f"({entry['record_rate_per_step']}/步)")
        if st:
            print(f"        ρ(L): 幂律指数={st['power_law']['exponent']} (R²={st['power_law']['r2']}) "
                  f"vs 指数率={st['exponential']['rate']} (R²={st['exponential']['r2']}) "
                  f"⟹ {st['preferred']}")
            print(f"        平均远足长度（未截断）={st['mean_len_all']} vs Kac 引理预言 N={L**D}  "
                  f"最长={st['max_len']}（拟合上限 {st['fit_Lmax']}）" + (f"   对 1D 精确首返律中位相对偏差="
                  f"{st.get('vs_exact_first_return_1D',{}).get('median_rel_dev')}" if D == 1 else ""))
        rec("LOOP", f"torus D={D} L={L}", "远足长度分布 ρ(L) vs L",
            st["lengths"] if st else [], st["counts"] if st else [],
            {"D": D, "N": L ** D, "records": entry["records"],
             "power_law": st["power_law"] if st else None,
             "exponential": st["exponential"] if st else None},
            note="极点（粒子）⟺ ρ 有 δ；支割 ⟺ 连续幂律/指数")

    # ---------- ③ 原生图 𝒢_T 上的环 + 稳态 L2 区域数 ----------
    print("=" * 92)
    print("③ 原生 Γ = 𝒢_T 上的自持环 + 稳态下自发形成几个 L2 区域")
    import l0_closure as L0
    import observable_sweep as OS
    import l2_excursions as LX
    from zcl import Engine
    eng = Engine()
    for T in (10, 12, 14, 16):
        reps, _ = L0.enumerate_necklaces(T, eng, chunk_bits=26, verbose=False)
        A = OS.build_sparse(T, reps, eng)
        res = run_loop(A, 300, 3000, r=1.0, seed=11, collect_support=True)
        act = np.array(res["activity"], float)
        # 用稳态前半段的远足碰撞图定区域
        sup = res["supports"]
        regions = None
        if sup and len(sup) > 20:
            half = sup[:len(sup) // 2]
            ex = [[(ln, s) for ln, s in half]]
            st = LX.region_stats(ex, A.shape[0], 1)
            regions = {"n_regions": len(st["regions"]),
                       "sizes_top": [r["n_excursions"] for r in st["regions"][:5]],
                       "claims_for_display": st["n_excursions"]}
        entry = {"T": T, "N": int(A.shape[0]), "activity_const": bool(act.min() == act.max()),
                 "activity_min": int(act.min()), "activity_max": int(act.max()),
                 "records": len(res["excursion_lengths"]),
                 "record_rate_per_step": round(len(res["excursion_lengths"]) / max(1, len(act)), 4),
                 "regions": regions,
                 "structure": excursion_structure(res["excursion_lengths"])}
        out["steady"].append({"graph": f"G_{T}", **entry})
        print(f"   T={T:>2} N={entry['N']:>6}: 活动量 [{entry['activity_min']},{entry['activity_max']}] "
              f"闭合 {entry['records']} 次 ({entry['record_rate_per_step']}/步)"
              + (f"  L2 区域数={regions['n_regions']} 最大区域={regions['sizes_top'][:3]}"
                 if regions else ""))
        del A

    with open(os.path.join(OUT, "z0_loop.json"), "w") as f:
        json.dump({**out, "records": RECORDS}, f, ensure_ascii=False, indent=1)
    print("=" * 92)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_loop.json")
    print("结论读法：r*=1 是**精确**临界点（每条闭合 −1、+r）。Z0① 的字面读法（闭合不停 = 继续走）")
    print("          恰好落在 r=1 ⟹ **不调参自持**；Z2 的全分支读法 r=deg v>1 ⟹ 爆炸相。")


if __name__ == "__main__":
    main(sys.argv[1:])
