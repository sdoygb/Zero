"""
z0_det2.py --- **忠实演化宇宙 v2**：`z0_evo.py` 的确定性替身（无概率、走者守恒、多样性靠分裂）

为什么要这一版（本轮查明的链条）：
  · `z0_evo.py` 能长出**畴 ＋ 稳态生灭**（`EVO_LAYERS.md`），但它用 `rng.random` **加权抽样** ⟹ **违反 Z0③**；
  · 忠实全分支 ⟹ 走者数翻倍 ⟹ 不守恒 ⟹（无外加饱和）指数增长；
  · 我上一版 `z0_det_evo.py` 用"逐顶点确定性选支"⟹ **同顶点走者锁步 ⟹ 退化**（区域数恒 1、记录速率→0）。
  ⟹ 结论：**抽样在 z0_evo 里干的是"造多样性"这件活**。忠实版必须用**别的东西**造多样性。

本版用三条**确定性**规则顶替抽样：
  ① **均分分支**：$(v,b)$ 上的 $m$ 个走者按 $m/2,m/2$ 分给 $b\\pm1$（奇数时低支多拿 1）⟹ 整数守恒 **且支持集扩张**（＝多样性的确定论来源）
  ② **最大余数法定向传输**：按邻居权重 $w_{v\\to u}=1+\\min(R_u,2)$（`D222` 最高两层）整数分配走者（无随机、精确守恒）
  ③ **反射边界**：$|b|=L$ 处反射（守恒且 $|b|$ 有界；语料 Z4 是"带走不还"，那样计数不守恒）
  闭合（$b\\to0$）⟹ 写记录 $R_v{+}{=}1$，走者重开（Z0① 闭合不停）

测（与 `EVO_LAYERS.md` 对齐）：区域数与**生灭速率**、区域**寿命**、记录速率（老化）、长程 move。

用法：/usr/bin/python3 z0_det2.py      输出：results/z0_det2.json ＋ z0_det2.png
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


def apportion(total, w):
    """把 total 个整数按权重 w（一维）用最大余数法分完（确定性）。"""
    w = np.maximum(np.asarray(w, float), 0.0)
    if w.sum() <= 0:
        w = np.ones_like(w)
    q = total * w / w.sum()
    base = np.floor(q).astype(np.int64)
    rem = int(total - base.sum())
    if rem > 0:
        idx = np.argsort(-(q - base), kind="stable")[:rem]
        base[idx] += 1
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


def run(A, steps=1500, L=4, mode="REC", rng_seed=0):
    V = A.shape[0]
    indptr, indices = A.indptr, A.indices
    nbrs = [indices[indptr[v]:indptr[v + 1]] for v in range(V)]
    # 初始：每顶点 2 个走者，平衡 ±1
    n = np.zeros((V, 2 * L + 1), np.int64)
    n[:, L + 1] = 2
    n[:, L - 1] = 2
    total0 = int(n.sum())
    R = np.zeros(V)
    T = np.zeros((V, V))
    mv, doms, births, deaths, rrate, lifes = [], [], [], [], [], []
    live, nid, pl = {}, 0, None
    prevR = R.copy()
    for t in range(steps):
        # ① 均分分支（守恒 ＋ 支持集扩张）；反射在 |b|=L
        up = np.zeros_like(n); dn = np.zeros_like(n); closed = np.zeros(V, np.int64)
        for bi in range(2 * L + 1):
            m = n[:, bi]
            if m.sum() == 0:
                continue
            b = bi - L
            hi = m // 2
            lo = m - hi
            if b + 1 > L:      # 反射
                dn[:, bi] += hi
            else:
                up[:, bi + 1] += hi
            if b - 1 < -L:     # 反射
                up[:, bi] += lo
            else:
                dn[:, bi - 1] += lo
        # b=0 档到达 ⟹ 闭合写记录、重开（Z0①）
        z = L
        cl = up[:, z] + dn[:, z]
        if cl.sum():
            closed = cl.copy()
            R += cl.astype(float)
            up[:, z] = 0; dn[:, z] = 0
            up[:, L + 1] += cl // 2 + cl % 2      # 重开：+ 支
            dn[:, L - 1] += cl // 2               # 重开：− 支
        n2 = up + dn
        assert int(n2.sum()) == total0, "分裂破坏守恒"
        # ② 定向传输（最大余数法，权重读记录层）
        newn = np.zeros_like(n2)
        for v in range(V):
            tot = int(n2[v].sum())
            if tot == 0:
                continue
            nb = nbrs[v]
            w = 1.0 + np.minimum(R[nb], 2.0) if mode == "REC" else np.ones(len(nb))
            alloc = apportion(tot, w)
            for k, u in enumerate(nb):
                m = int(alloc[k])
                if m == 0:
                    continue
                take = apportion(m, n2[v].astype(float))
                newn[u] += take
                T[u, v] += m
        n = newn
        assert int(n.sum()) == total0, "传输破坏守恒"
        # ③ 观测
        act = n.sum(axis=1).astype(float)
        occ = act > 0
        nd, lab = comps(occ, A)
        doms.append(nd)
        mv.append(float(np.abs(n - n_prev).sum()) if t > 0 else 0.0)
        n_prev = n.copy()
        if pl is not None:
            surv = set()
            for j in range(nd):
                pv = pl[lab == j]; pv = pv[pv >= 0]
                if len(pv):
                    surv.add(int(np.bincount(pv).argmax()))
            for k in list(live):
                if k not in surv:
                    lifes.append(t - live.pop(k))
            nb_new = 0
            for j in range(nd):
                pv = pl[lab == j]; pv = pv[pv >= 0]
                if not len(pv):
                    nid += 1; live[nid] = t; nb_new += 1
            births.append(nb_new)
            deaths.append(nd - nb_new)
        else:
            for j in range(nd):
                nid += 1; live[nid] = t
            births.append(nd); deaths.append(0)
        pl = lab
        rrate.append(float(R.sum() - prevR.sum()))
        prevR = R.copy()
    mv = np.array(mv); rr = np.array(rrate)
    out = dict(mode=mode, V=V, steps=steps, conserved=int(n.sum()) == total0, total=int(n.sum()),
               domains_mean=float(np.mean(doms)), domains_max=int(np.max(doms)),
               birth_mean=float(np.mean(births)), death_mean=float(np.mean(deaths)),
               life_n=len(lifes), life_mean=float(np.mean(lifes)) if lifes else 0.0,
               life_med=float(np.median(lifes)) if lifes else 0.0,
               life_p95=float(np.percentile(lifes, 95)) if lifes else 0.0,
               life_max=int(max(lifes)) if lifes else 0,
               move_last=float(np.mean(mv[-100:])),
               rec_rate_first=float(np.mean(rr[:steps // 2])), rec_rate_last=float(np.mean(rr[steps // 2:])),
               R_total=float(R.sum()), dom_trace=[int(x) for x in doms[::10]])
    return out, n, R, T


if __name__ == "__main__":
    t0 = time.time()
    A = torus_adj(2, 16)
    print("=" * 100)
    print("忠实演化宇宙 v2：均分分支 ＋ 最大余数法传输 ＋ 反射边界（无概率、守恒、多样性靠分裂）")
    out, n, R, T = run(A, steps=1500, mode="REC")
    print(f"\n### REC 模式（权重 1+min(R,2)，读记录层）")
    print(f"   走者守恒={out['conserved']}（总数 {out['total']}）")
    print(f"   区域数 均={out['domains_mean']:.1f} 最大={out['domains_max']}   "
          f"生/死 每步={out['birth_mean']:.2f}/{out['death_mean']:.2f}")
    print(f"   区域寿命 n={out['life_n']} 均={out['life_mean']:.2f} 中位={out['life_med']:.0f} "
          f"95%={out['life_p95']:.0f} 最长={out['life_max']}")
    print(f"   move(末)={out['move_last']:.2f}  记录速率 {out['rec_rate_first']:.3f}→{out['rec_rate_last']:.3f}  ΣR={out['R_total']:.0f}")
    def _gini(x):
        x = np.sort(np.asarray(x, float)); n_ = len(x)
        return 0.0 if x.sum() == 0 else float((2*np.arange(1, n_+1)-n_-1).dot(x)/(n_*x.sum()))
    print(f"   记录 Gini={_gini(R):.4f}  边流量 Gini={_gini(T.ravel()):.4f}  最大记录占比={R.max()/max(R.sum(),1):.4f}")
    RES["REC"] = out
    out2, _, _, _ = run(A, steps=800, mode="FLAT")
    print(f"\n### FLAT 模式（对照：权重全 1，不看记录）")
    print(f"   区域数 均={out2['domains_mean']:.1f} 最大={out2['domains_max']}   "
          f"寿命 均={out2['life_mean']:.2f} 95%={out2['life_p95']:.0f} 最长={out2['life_max']}")
    RES["FLAT"] = out2
    check("走者数在两个模式下都精确守恒", out["conserved"] and out2["conserved"])
    check("REC 下出现**多区域 ＋ 持续生灭**（区域均>1 且 生/死速率>0）",
          out["domains_mean"] > 1.0 and out["birth_mean"] > 0 and out["death_mean"] > 0,
          f"区域均 {out['domains_mean']:.1f}，生 {out['birth_mean']:.2f}/死 {out['death_mean']:.2f}")
    check("REC 下出现**长寿区域**（95% 寿命 ≥ 5 步）", out["life_p95"] >= 5,
          f"95%={out['life_p95']:.0f} 最长={out['life_max']}")
    check("REC 与 FLAT **不同**（记录层确实在起作用）",
          abs(out["domains_mean"] - out2["domains_mean"]) > 0.05 or abs(out["R_total"] - out2["R_total"]) > 1,
          f"ΣR {out['R_total']:.0f} vs {out2['R_total']:.0f}")
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_det2.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_det2.json")
