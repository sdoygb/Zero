"""
z0_life.py --- **"Zero 演化起来"的相图**：变化、记忆、结构持久性，能不能同时拿到？

背景：§27 的折叠给了「变化不停」（混沌、非周期、正 Lyapunov），但 §35 之前查到的两个疑点：
  ① 细尺度畴**寿命恒为 1 代**（重排，不是成形/存活/解体）
  ② 宏观（块平均）场自相关 $\tau_{1/e}$ 也只有 **1–6 代**（无时间记忆）
本文件把"结构持久性"当第三个可观测量，扫**点火率** $f$（= 每代被更新的顶点比例）：

  同步（`sync`，$f{=}1$）／按 **Kac 钟**的局部异步（`async_kac`，$f$ 由 $\tau_v=2|E|/\deg v$ 的尺度定）／
  事件驱动（`async_event`，$f$ 由"该顶点自己超过自身滑动均值"定）
  —— 全部确定性、无概率（Z0③）。

判据（三个都要）：
  **变化** `move`(末) $>0.05$ ／ **记忆** 全局活动 ACF $\tau_{1/e}\ge3$ ／ **结构持久** 畴寿命中位 $>3$

用法：/usr/bin/python3 z0_life.py      输出：results/z0_life.json
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
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


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


def run(A, mode="sync", scale=0.25, gens=4000, L=4, theta=1.0, seed=11):
    AL = (A ** L).toarray()
    R = np.diag(AL).copy()
    AT = AL.T.copy()
    V = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel()
    E2 = float(deg.sum())
    tau = np.maximum(1, np.round(scale * E2 / np.maximum(deg, 1))).astype(np.int64)
    rng = np.random.default_rng(seed)
    s = rng.random(V) + 0.1
    s /= s.mean()
    clock = np.zeros(V, np.int64)
    runmean = s.copy()
    mv, lifes, live, nid, pl = [], [], {}, 0, None
    fired = []
    for g in range(gens):
        c = s * R
        w = AT @ s
        rho = np.where(w > 0, c / np.maximum(w, 1e-300), 0.0)
        rho /= rho.mean()
        d = np.abs(rho - s)
        prev = s.copy()
        if mode == "sync":
            s = d
            fired.append(V)
        elif mode == "async_kac":
            clock += 1
            fire = clock >= tau
            s = np.where(fire, d, s)
            clock[fire] = 0
            fired.append(int(fire.sum()))
        elif mode == "async_event":
            ev = c > runmean
            runmean = 0.99 * runmean + 0.01 * c
            s = np.where(ev, d, s)
            fired.append(int(ev.sum()))
        if s.mean() <= 1e-300:
            s = prev
        s = s / s.mean()
        mv.append(float(np.abs(s - prev).sum() / V))
        occ = s > theta * s.mean()
        n, lab = comps(occ, A)
        surv = set()
        if pl is not None:
            for j in range(n):
                pv = pl[lab == j]
                pv = pv[pv >= 0]
                if len(pv):
                    surv.add(int(np.bincount(pv).argmax()))
            for k in list(live):
                if k not in surv:
                    lifes.append(g - live.pop(k))
            for j in range(n):
                pv = pl[lab == j]
                pv = pv[pv >= 0]
                if not len(pv):
                    nid += 1
                    live[nid] = g
        else:
            for j in range(n):
                nid += 1
                live[nid] = g
        pl = lab
    mv = np.array(mv)
    acf = [float(((mv[:-k] - mv.mean()) * (mv[k:] - mv.mean())).mean() / max(mv.var(), 1e-30))
           for k in range(1, 200)]
    tau_acf = next((k for k, a in enumerate(acf, 1) if a < 1 / math.e), None)
    sl = float(np.polyfit(np.log(np.arange(200, gens, 200)),
                          np.log([max(mv[q - 20:q].mean(), 1e-12) for q in range(200, gens, 200)]), 1)[0])
    return dict(mode=mode, scale=scale, f=round(float(np.mean(fired)) / V, 4),
                move=round(float(mv[-20:].mean()), 4), slope=round(sl, 4), acf_tau=tau_acf,
                n_life=int(len(lifes)),
                life_mean=round(float(np.mean(lifes)), 2) if lifes else 0.0,
                life_med=float(np.median(lifes)) if lifes else 0.0,
                life_p95=float(np.percentile(lifes, 95)) if lifes else 0.0,
                life_max=int(max(lifes)) if lifes else 0)


if __name__ == "__main__":
    t0 = time.time()
    A, grp = ZT.build_hierarchical()
    print("=" * 100)
    print("「Zero 演化起来」的相图：变化／记忆／结构持久性（层级 Γ，L=4，4000 代）")
    print(f"{'模式':>12} {'scale':>6} {'f':>6} {'move':>8} {'斜率':>8} {'ACF τ':>6} "
          f"{'寿命均':>7} {'中位':>5} {'95%':>6} {'最长':>6}  判定")
    rows = []
    for mode, scale in (("sync", 1.0), ("async_kac", 0.25), ("async_kac", 0.10), ("async_kac", 0.05),
                        ("async_kac", 0.02), ("async_kac", 0.01), ("async_event", 1.0)):
        r = run(A, mode=mode, scale=scale)
        rows.append(r)
        ok = (r["move"] > 0.05) and (r["acf_tau"] or 0) >= 3 and r["life_med"] > 3
        v = "**三样齐**" if ok else ("长寿但冻" if r["life_mean"] > 3 else
                                     ("变化但短命" if r["move"] > 0.05 else "死"))
        print(f"{mode:>12} {scale:6.2f} {r['f']:6.3f} {r['move']:8.4f} {r['slope']:+8.4f} "
              f"{str(r['acf_tau']):>6} {r['life_mean']:7.2f} {r['life_med']:5.0f} "
              f"{r['life_p95']:6.0f} {r['life_max']:6d}  {v}")
    RES["phase"] = rows
    both = [r for r in rows if r["move"] > 0.05 and r["life_med"] > 3]
    check("相图里**没有**任何一档同时「变化(move>0.05)」与「结构持久(寿命中位>3)」",
          len(both) == 0, f"同时满足的档 = {[(r['mode'], r['scale']) for r in both]}" if both else "全部二选一")
    check("存在「长寿但冻」档（有结构、无变化）", any(r["life_mean"] > 3 and r["move"] < 0.05 for r in rows),
          f"{[(r['mode'], r['scale'], r['life_mean'], r['move']) for r in rows if r['life_mean'] > 3]}")
    check("存在「变化但短命」档（有变化、无结构持久）",
          any(r["move"] > 0.05 and r["life_med"] <= 3 for r in rows), "sync 与高点火率档")
    print(f"\n   ⟹ **结论**：{len(both)} 档三样齐；相图是二分的 —— 要变化就得放弃结构持久性，"
          f"要结构持久性就得冻住。（**已被 §36–§40 超越**：真解是零和层上的无偏好核。）")
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_life.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_life.json")
