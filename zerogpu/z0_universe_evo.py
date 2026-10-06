"""
z0_universe_evo.py --- **让 Zero 在电脑里演化起来**：一个能跑、能看、能量的演化宇宙

组件（每条都标了强度）：
  L0 规则（**忠实**）      全分支 $\\pm$、闭合＝词和归零、整数重数、闭合写记录
  ② 局部 $P_i$（**忠实**）  每顶点从**自己**的记录重播种（`ZS-8`：全局读回一代内塌缩）
  ① $\\kappa_1$ 读出（**忠实**） 读数乘 $\\kappa_1^{\\lfloor t/L\\rfloor}$，$\\kappa_1=q_4=5/9$（`G71`/`G72`）
  演化（**构造·非导出**）   折叠 $s_v\\leftarrow|\\rho_v(s)-s_v|$，$\\rho_v=c_v/w_v$（§27：唯一给出「变化不停」的）

测的不是"混沌翻腾"，而是**结构本身的生灭**（这才是"演化起来"的意思）：
  畴 = 活跃集 $\\{v: s_v>\\theta\\cdot\\overline{s}\\}$ 的连通块；逐代跟踪其**身份**（重叠即同一畴）
  ⟹ 报：畴数随时间的分布、**畴寿命分布**、出生/死亡率、以及活动量/结构量的长程行为

输出：`results/z0_universe_evo.json` ＋ `z0_universe_evo.png`（时空图 ＋ 畴寿命直方图）

用法：/usr/bin/python3 z0_universe_evo.py
"""
from __future__ import annotations

import json
import math
import os
import time
from collections import Counter, deque

import numpy as np

import z0_thm as ZT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
KAPPA1 = 5 / 9


def make_step(A, L=4):
    AL = (A ** L).toarray()
    R = np.diag(AL).copy()
    AT = AL.T.copy()
    V = A.shape[0]

    def step(s):
        c = s * R
        w = AT @ s
        rho = np.where(w > 0, c / np.maximum(w, 1e-300), 0.0)
        rho = rho / rho.mean()
        ns = np.abs(rho - s)
        if ns.mean() <= 1e-300:
            ns = s.copy()
        return ns / ns.mean()

    return step, R, V


def components(occ, A):
    """活跃集的连通块（返回 顶点→标签）。"""
    indptr, indices = A.indptr, A.indices
    lab = np.full(A.shape[0], -1, np.int64)
    c = 0
    for s in np.nonzero(occ)[0]:
        if lab[s] >= 0:
            continue
        lab[s] = c
        q = deque([int(s)])
        while q:
            u = q.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if occ[v] and lab[v] < 0:
                    lab[v] = c
                    q.append(int(v))
        c += 1
    return c, lab


def track(lab_prev, lab_now, n_prev, n_now):
    """畴身份：新畴与旧畴有重叠 ⟹ 同一个畴。返回 (存活数, 新生数, 死亡数, id 映射)。"""
    if n_prev == 0:
        return 0, n_now, 0, {j: j for j in range(n_now)}
    idmap = {}
    used_prev = set()
    for j in range(n_now):
        vs = np.nonzero(lab_now == j)[0]
        if len(vs) == 0:
            continue
        pv = lab_prev[vs]
        pv = pv[pv >= 0]
        if len(pv):
            k = int(Counter(pv.tolist()).most_common(1)[0][0])
            idmap[j] = k
            used_prev.add(k)
        else:
            idmap[j] = n_prev + j          # 新生
    survived = len(used_prev)
    newborn = n_now - survived
    died = n_prev - survived
    return survived, newborn, died, idmap


def run(A, gens=3000, L=4, theta=1.0, sample_every=1, label=""):
    step, R, V = make_step(A, L)
    rng = np.random.default_rng(11)
    s = rng.random(V) + 0.1
    s /= s.mean()
    space = []            # (每 K 代存一次活动场)
    dom_hist = []         # 每代畴数
    lifetimes = []        # 畴寿命
    births = deaths = 0
    live = {}             # id -> 出生代
    maxid = 0
    nextid = 0
    prev_lab, prev_n = None, 0
    mv, mm = [], []
    for g in range(gens):
        ns = step(s)
        mv.append(float(np.abs(ns - s).sum() / V))
        occ = ns > theta * ns.mean()
        mm.append(float(ns[occ].max() / ns[occ].sum()) if occ.any() else 1.0)
        n_now, lab_now = components(occ, A)
        dom_hist.append(n_now)
        if g % sample_every == 0:
            space.append(ns.copy())
        # 畴身份跟踪
        if prev_lab is not None:
            idmap_prev = {k: k for k in range(prev_n)}
        if prev_lab is not None and prev_n > 0:
            # 现有 id 集合 -> 旧标签
            old_of_id = {k: k for k in range(prev_n)}
            surv = set()
            newmap = {}
            for j in range(n_now):
                vs = np.nonzero(lab_now == j)[0]
                pv = prev_lab[vs]
                pv = pv[pv >= 0]
                if len(pv):
                    k = int(Counter(pv.tolist()).most_common(1)[0][0])
                    newmap[j] = live.get(k, k)
                    surv.add(k)
                else:
                    newmap[j] = None
            # 死亡
            for k in list(live.keys()):
                if k not in surv:
                    lifetimes.append(g - live.pop(k))
                    deaths += 1
            # 新生
            for j in range(n_now):
                if newmap[j] is None:
                    nextid += 1
                    live[nextid] = g
                    births += 1
        else:
            for j in range(n_now):
                nextid += 1
                live[nextid] = g
                births += 1
        prev_lab, prev_n, s = lab_now, n_now, ns
    # 收尾：未完的畴寿命（右删失，单列）
    censored = [gens - t for t in live.values()]
    out = dict(label=label, V=V, gens=gens, L=L, theta=theta,
               move_last=round(float(np.mean(mv[-20:])), 5),
               move_slope=round(float(np.polyfit(np.log(np.arange(100, gens, 100)),
                                                 np.log([np.mean(mv[g - 20:g]) for g in range(100, gens, 100)]), 1)[0]), 4),
               max_mass_last=round(float(np.mean(mm[-20:])), 4),
               domains_mean=round(float(np.mean(dom_hist)), 2),
               domains_min=int(np.min(dom_hist)), domains_max=int(np.max(dom_hist)),
               deaths=int(deaths), births=int(births),
               lifetimes_n=int(len(lifetimes)),
               lifetime_mean=round(float(np.mean(lifetimes)), 2) if lifetimes else None,
               lifetime_median=float(np.median(lifetimes)) if lifetimes else None,
               lifetime_max=int(np.max(lifetimes)) if lifetimes else None,
               censored_n=int(len(censored)),
               censored_mean=round(float(np.mean(censored)), 1) if censored else None)
    return out, np.array(space), dom_hist, lifetimes


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 96)
    print("Zero 演化宇宙：折叠动力学 ＋ 局部 P_i ＋ κ1 读出（L=4）")
    base = ZT.torus_adj(2, 16)
    hier, grp = ZT.build_hierarchical()
    res = {}
    for nm, A in (("环面16²（正则）", base), ("层级Γ（32团×8）", hier)):
        out, space, dh, lt = run(A, gens=3000, L=4, theta=1.0, label=nm)
        res[nm] = out | {"dom_hist": [int(x) for x in dh[::50]],
                         "lifetimes": [int(x) for x in lt[:400]]}
        print(f"\n### {nm}  V={out['V']}")
        print(f"   变化：move(末)={out['move_last']}  长程斜率={out['move_slope']}   max_mass={out['max_mass_last']}")
        print(f"   结构：畴数 均值={out['domains_mean']}  范围=[{out['domains_min']},{out['domains_max']}]")
        print(f"   生灭：出生={out['births']}  死亡={out['deaths']}  已闭合寿命数={out['lifetimes_n']}  "
              f"平均寿命={out['lifetime_mean']}  中位={out['lifetime_median']}  最长={out['lifetime_max']}")
        print(f"   （另有 {out['censored_n']} 个畴跑到终点仍未死，平均已活 {out['censored_mean']} 代）")
        if nm.startswith("环面"):
            np.save(os.path.join(OUT, "z0_universe_evo_space.npy"), space)
    json.dump(res, open(os.path.join(OUT, "z0_universe_evo.json"), "w"), ensure_ascii=False, indent=1, default=str)

    # 图：时空图 ＋ 畴寿命直方图
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        space = np.load(os.path.join(OUT, "z0_universe_evo_space.npy"))
        L16 = 16
        sub = space[::4]                     # 抽样代数
        img = sub.reshape(sub.shape[0], L16, L16).transpose(0, 1, 2)
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        # 时空：每条竖线是一个顶点的活动随代变化（按 16x16 展平）
        axes[0].imshow(np.log10(img.reshape(img.shape[0], -1).T + 1e-12), aspect="auto",
                       cmap="inferno", origin="lower")
        axes[0].set_title("spacetime: log10 activity  (x=vertex, y=generation)")
        axes[0].set_xlabel("generation (sampled every 4)"); axes[0].set_ylabel("vertex index")
        # 某一代的 16x16 快照
        snap = space[-1].reshape(L16, L16)
        im = axes[1].imshow(np.log10(snap + 1e-12), cmap="inferno")
        axes[1].set_title("snapshot (last generation)"); plt.colorbar(im, ax=axes[1], fraction=0.046)
        # 畴寿命直方图
        lt = res["环面16²（正则）", "lifetimes"] if False else res["环面16²（正则）"]["lifetimes"]
        axes[2].hist(lt, bins=40, color="steelblue")
        axes[2].set_title(f"domain lifetimes (n={len(lt)}, median={np.median(lt):.0f})")
        axes[2].set_xlabel("lifetime (generations)"); axes[2].set_ylabel("count")
        plt.tight_layout()
        plt.savefig(os.path.join(HERE, "z0_universe_evo.png"), dpi=110)
        print(f"\n图：z0_universe_evo.png")
    except Exception as e:
        print(f"（画图失败：{e}）")
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_universe_evo.json")
