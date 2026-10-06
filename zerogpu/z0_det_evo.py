"""
z0_det_evo.py --- **忠实版演化宇宙**：把 `z0_evo.py` 的"加权抽样"换成**确定性选择**

症结（本轮查明）：
  `z0_evo.py` 能长出**畴 ＋ 稳态生灭**（`EVO_LAYERS.md`），但它用 `rng.random` 做**加权抽样** ⟹ **违反 Z0③**；
  而它之所以自持，是因为**走者数守恒**（每个走者每步只走一条边）。
$$
\\textbf{Z0③}\\Rightarrow\\textbf{全分支}\\Rightarrow\\textbf{走者数不守恒};(\\text{无外加饱和})\\Rightarrow\\text{指数增长}.
\\qquad
\\textbf{守恒}\\Rightarrow\\textbf{每步只能选一条边}\\Rightarrow\\textbf{需要一条"选择规则"}.
$$
  概率是一条选择规则，但被 Z0③ 排除。**本文件试另一条：用模相位做确定性选择。**

规则（全部确定性、无概率、无自由参数）：
  · 状态：$n[v][b]$ ＝ 顶点 $v$ 上平衡为 $b$ 的走者数；**走者总数守恒**（$\\sum n$ 常数）
  · **分支选择（新）**：$d_v=\\operatorname{sign}\\cos(\\varphi\\,g_v)$，$\\varphi=\\log 7.2$（§19 的 $L{=}4$ 转动数），
    $g_v$ ＝ 该顶点自己的**闭合代数**（只在自己闭合时 $+1$）⟹ **逐顶点异质、确定性、无参数**
  · 闭合：$b+d_v=0$ ⟹ 写记录 $R_v{+}{=}1$，走者**重开**（Z0① 闭合不停）：$b\\leftarrow d_v$
  · **寿命回收（Z4 的守恒版）**：$|b+d_v|>L$（未闭合就走满寿命）⟹ 进终端账本并**原地重开** $b\\leftarrow d_v$
    —— 这是唯一能让"守恒"与"$|b|$ 有界"并存的办法（语料的 Z4 是带走不还，那样计数不守恒）
  · 传输：$v$ 的走者按邻居权重 $w_{v\\to u}=1+\\min(R_u,2)$（`D222` 最高两层）用**最大余数法**整数分配
    —— 这就是"加权抽样"的**确定性替身**（无随机、且精确守恒）

测（与 `EVO_LAYERS.md` 同一套，便于对照）：区域数、**生灭速率**、区域**寿命**、活动空间自相关、
记录速率（自发性老化）、以及长程 `move` 斜率。

用法：/usr/bin/python3 z0_det_evo.py      输出：results/z0_det_evo.json
"""
from __future__ import annotations

import json
import math
import os
import time
from collections import Counter, deque

import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
PHI = math.log(7.2)
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def largest_remainder(counts, weights):
    """整数分配 + 最大余数法：把 sum(counts) 个走者按 weights 分给各目标（确定性、精确守恒）。"""
    total = int(counts.sum())
    if total == 0:
        return np.zeros_like(weights, dtype=np.int64)
    w = np.maximum(weights, 0.0)
    if w.sum() <= 0:
        w = np.ones_like(w)
    quota = total * w / w.sum()
    base = np.floor(quota).astype(np.int64)
    rem = total - int(base.sum())
    if rem > 0:
        order = np.argsort(-(quota - base), kind="stable")   # 余数大的先拿（稳定排序 ⟹ 确定性）
        base[order[:rem]] += 1
    return base


def comps(occ, A):
    indptr, indices = A.indptr, A.indices
    V = A.shape[0]
    lab = np.full(V, -1, np.int64)
    c = 0
    for s0 in np.nonzero(occ)[0]:
        if lab[s0] >= 0:
            continue
        lab[s0] = c
        q = deque([int(s0)])
        while q:
            u = q.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if occ[v] and lab[v] < 0:
                    lab[v] = c
                    q.append(int(v))
        c += 1
    return c, lab


def run(A, mode="phase", steps=6000, L=4, n_walk=None, walk=1.0, seed=0, window=200):
    V = A.shape[0]
    rng = np.random.default_rng(seed)
    if n_walk is None:
        n_walk = V
    # 初始：每个顶点 n 个走者，平衡 ±1 各半
    n = np.zeros((V, 2 * L + 1), dtype=np.int64)          # 平衡 b ∈ [-L, L]，索引 b+L
    for v in range(V):
        n[v, 1 + L] = n_walk // V // 2 + 1                # b=+1
        n[v, -1 + L] = n_walk // V // 2 + 1               # b=-1
    term = 0                                              # 终端账本（寿命回收计数）
    R = np.zeros(V)
    g = np.zeros(V, np.int64)
    total0 = int(n.sum())
    # 邻接（列索引）
    indptr, indices = A.indptr, A.indices
    rows = []
    mv, doms, births, deaths, rec_rates = [], [], [], [], []
    live, nid, lifes, pl = {}, 0, [], None
    prev_R = R.copy()
    for t in range(steps):
        # ① 分支选择（确定性；同时依赖**余额符号**与**逐顶点相位**，否则同顶点走者会锁步）
        bsgn = np.sign(np.arange(-L, L + 1))          # 每一档余额的符号
        bsgn[bsgn == 0] = 1.0
        phase = np.cos(PHI * g)
        if mode == "phase":
            # 相长 ⟹ 朝 0 走（闭合、写记录）；相消 ⟹ 远离 0（继续探索）
            d = (-np.outer(np.sign(phase), bsgn)).astype(np.int64)          # [V, 2L+1]
            d[d == 0] = 1
        elif mode == "alt":
            d = np.where(t % 2 == 0, -np.sign(bsgn), np.sign(bsgn)).astype(np.int64)[None, :].repeat(V, 0)
        else:
            d = np.ones((V, 2 * L + 1), np.int64)
        newn = np.zeros_like(n)
        moved = 0
        # ② 逐顶点：传输（最大余数法）＋ 分支
        for v in range(V):
            nbrs = indices[indptr[v]:indptr[v + 1]]
            wts = 1.0 + np.minimum(R[nbrs], 2.0)                 # D222 最高两层
            per_b = n[v]
            if per_b.sum() == 0:
                continue
            alloc = largest_remainder(per_b, wts)                # 每个邻居分到几个走者（不分平衡）
            # 该顶点的走者统一按 d_v 选分支；闭合者重开
            for k, u in enumerate(nbrs):
                m = int(alloc[k])
                if m == 0:
                    continue
                # 把 m 个走者按原平衡分布摊开（保持平衡结构）
                src = per_b
                if src.sum() > 0:
                    take = largest_remainder(src, src.astype(float))
                    # 简化：按比例取（确定性）
                    take = np.floor(m * src / src.sum()).astype(np.int64)
                    short = m - int(take.sum())
                    if short > 0:
                        order = np.argsort(-src, kind="stable")[:short]
                        take[order] += 1
                for bi in np.nonzero(take)[0]:
                    cnt = int(take[bi])
                    b = bi - L
                    nb = b + int(d[v, bi])
                    if nb == 0:                                   # 闭合
                        R[v] += cnt
                        g[v] += 1
                        nb = int(d[v, bi])                        # 重开（Z0①）
                    elif abs(nb) > L:                             # 寿命用尽（Z4）
                        term += cnt
                        nb = int(d[v, bi])                        # 回收重开（守恒版 Z4）
                    newn[u, nb + L] += cnt
        n = newn
        tot = int(n.sum())
        if tot != total0:                                        # 守恒断言（见下）
            pass
        # ③ 观测
        act = n.sum(axis=1).astype(float)
        occ = act > act.mean()
        nd, lab = comps(occ, A)
        doms.append(nd)
        mv.append(float(np.abs(n - (n if t == 0 else n_prev)).sum()) if t > 0 else 0.0)
        n_prev = n.copy()
        if pl is not None:
            surv = set()
            for j in range(nd):
                pv = pl[lab == j]
                pv = pv[pv >= 0]
                if len(pv):
                    surv.add(int(np.bincount(pv).argmax()))
            for k in list(live):
                if k not in surv:
                    lifes.append(t - live.pop(k))
            newborn = 0
            for j in range(nd):
                pv = pl[lab == j]
                pv = pv[pv >= 0]
                if not len(pv):
                    nid += 1
                    live[nid] = t
                    newborn += 1
            births.append(newborn)
            deaths.append(len(surv) and 0 or 0)
        else:
            for j in range(nd):
                nid += 1
                live[nid] = t
        pl = lab
        rec_rates.append(float(R.sum() - prev_R.sum()))
        prev_R = R.copy()
    mv = np.array(mv)
    rrate = np.array(rec_rates)
    half = steps // 2
    out = dict(mode=mode, V=V, steps=steps, total_conserved=int(n.sum()) == total0, terminal=term,
               total_now=int(n.sum()),
               domains_mean=float(np.mean(doms)), domains_min=int(np.min(doms)), domains_max=int(np.max(doms)),
               lifetime_n=len(lifes),
               lifetime_mean=float(np.mean(lifes)) if lifes else 0.0,
               lifetime_med=float(np.median(lifes)) if lifes else 0.0,
               lifetime_p95=float(np.percentile(lifes, 95)) if lifes else 0.0,
               lifetime_max=int(max(lifes)) if lifes else 0,
               birth_per_win=float(np.mean(births)) if births else 0.0,
               move_last=float(np.mean(mv[-100:])),
               rec_rate_first=float(np.mean(rrate[:half])), rec_rate_last=float(np.mean(rrate[half:])),
               R_total=float(R.sum()))
    return out


if __name__ == "__main__":
    t0 = time.time()
    A = torus_adj(2, 16)
    print("=" * 100)
    print("忠实版演化宇宙：确定性相位选择 ＋ 最大余数法传输（无概率、无自由参数、走者守恒）")
    for mode in ("phase", "alt", "allplus"):
        r = run(A, mode=mode, steps=4000)
        print(f"\n### 模式 {mode}")
        print(f"   走者守恒={r['total_conserved']}（{r['total_now']}）  区域数 均={r['domains_mean']:.1f} "
              f"范围=[{r['domains_min']},{r['domains_max']}]")
        print(f"   区域寿命：n={r['lifetime_n']} 均={r['lifetime_mean']:.2f} 中位={r['lifetime_med']:.0f} "
              f"95%={r['lifetime_p95']:.0f} 最长={r['lifetime_max']}")
        print(f"   每窗新生={r['birth_per_win']:.2f}   move(末)={r['move_last']:.3f}   "
              f"记录速率 {r['rec_rate_first']:.3f}→{r['rec_rate_last']:.3f}（老化？）   ΣR={r['R_total']:.0f} 终端={r['terminal']}")
        RES[mode] = r
    ok_cons = all(RES[m]["total_conserved"] for m in RES)
    check("走者数精确守恒（三个模式）", ok_cons)
    ph = RES["phase"]
    check("相位选择下：区域生灭在稳态下持续（寿命中位 < 步数，且区域数 > 1）",
          ph["domains_mean"] > 1.5 and ph["lifetime_n"] > 0,
          f"区域均 {ph['domains_mean']:.1f}，寿命数 {ph['lifetime_n']}，中位 {ph['lifetime_med']:.0f}")
    check("相位选择下出现**长寿畴**（95% 寿命 ≥ 5 代）", ph["lifetime_p95"] >= 5,
          f"95%={ph['lifetime_p95']:.0f} 最长={ph['lifetime_max']}")
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_det_evo.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_det_evo.json")
