"""
z0_evolve.py --- 攻「演化」这道墙：把 §14.1 剩下**未被否证**的两条候选 ＋ §19.8 的杠杆做成对照实验

背景（Z0_CORE §7/§8/§14/§16）：
  · 忠实框架里「稳定」与「有空间结构」互斥：要么爆炸（全局同步＋`all`/`top2`），要么稳定但种子均匀（`cap`）。
  · §8 用 **Kac 寿命** τ_v=2|E|/deg(v) 做局部异步毁灭 ⟹ 非正则 Γ 上活动 Gini 0→0.43、
    corr(活动,deg)→−0.44，**但活跃簇只有 1–3 个**（平滑调制，不是成畴），且 τ_v 只依赖 deg（**静态**）。
  · §14.1 穷尽否证"靠类竞争成畴"之后，只剩两条没试：
      ① Γ 的**层次结构**：让回返概率 (A^L)_vv 本身成块；
      ② Z4 的寿命 τ_v 依赖**局部历史**（而非仅 deg v）。
  · §19.8（上一轮）：模相位对动力学**不承重**（`amp` 是纯读出）；把读出接回种子是唯一
    在不动 L、不动 q_4 的前提下让状态扇区继承旋转的地方。

本机做三个对照（都不引入新自由参数）：
  E1 候选②：把 τ_v 从"Kac 静态"换成"**局部历史**"——① 自观测周期（上次两次点火的实际间隔）
           ② 记录增益匹配（本次循环攒到的记录 ≥ 上次循环攒到的）。看是否出现
           "有界 ＋ 持续变化（不冻结）＋ 多个持续簇"。
  E2 候选①：造**层级 Γ**（嵌套模块：小团 → 中模块 → 大模块，三层稀疏度），
           先量 (A^L)_vv 是否真的**按层分块**，再跑 §8 的局部异步引擎，看簇数是否 >
           平坦基线（环面＋枢纽）的 1–3 个。
  E3 §19.8：在 L=4 词级引擎上把**干涉**接进 Z4 的清空选择（相消时清掉弱模的载体），
           种子仍 = 闭合分支数（**保持 r*=1 临界**，既不爆炸也不死亡）⟹ 看状态周期是否
           从 2 变成"无周期"（被无理旋转驱动）。

用法：/usr/bin/python3 z0_evolve.py        （约 1–2 分钟）
输出：results/z0_evolve.json
"""
from __future__ import annotations

import json
import os
import time
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations

import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}


# ================================================================ 图
def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def add_hubs(A, frac=0.10, extra=4, seed=3):
    rng = np.random.default_rng(seed)
    V = A.shape[0]
    A = A.tolil()
    for u in rng.choice(V, size=max(1, int(frac * V)), replace=False):
        for _ in range(extra):
            v = int(rng.integers(0, V))
            if v != u:
                A[u, v] = 1.0
                A[v, u] = 1.0
    return csr_matrix(A)


def build_hierarchical(n_mod=32, m=8, group=4, seed=0, p_in=0.0, p_mid=1, p_top=1):
    """
    **层级 Γ**（Ravasz–Barabási 式嵌套模块）：V = n_mod*m = 256
      level 0（模块内）：团（每对相连）⟹ 度高、回返概率大
      level 1（组内模块之间）：每对模块加 p_mid 条边 ⟹ 稀疏
      level 2（组与组之间）：每对组加 p_top 条边 ⟹ 更稀疏
    返回 A、以及每个顶点所属的 (module, group)。
    """
    V = n_mod * m
    A = np.zeros((V, V))
    for c in range(n_mod):
        idx = np.arange(c * m, (c + 1) * m)
        for i in idx:
            for j in idx:
                if i < j:
                    A[i, j] = A[j, i] = 1.0
    rng = np.random.default_rng(seed)
    n_grp = n_mod // group
    for g in range(n_grp):
        mods = [g * group + k for k in range(group)]
        for a in range(len(mods)):
            for b in range(a + 1, len(mods)):
                for _ in range(p_mid):
                    u = int(rng.integers(mods[a] * m, (mods[a] + 1) * m))
                    v = int(rng.integers(mods[b] * m, (mods[b] + 1) * m))
                    A[u, v] = A[v, u] = 1.0
    for g in range(n_grp):
        h = (g + 1) % n_grp
        for _ in range(p_top):
            u = int(rng.integers(g * group * m, (g + 1) * group * m))
            v = int(rng.integers(h * group * m, (h + 1) * group * m))
            A[u, v] = A[v, u] = 1.0
    mod_of = np.repeat(np.arange(n_mod), m)
    grp_of = np.repeat(np.arange(n_grp), group * m)
    return csr_matrix(A), mod_of, grp_of


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def clusters(occ, A):
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
    return c


# ================================================================ E1/E2 引擎：局部异步毁灭 ＋ 可换时钟
def run_local(A, clock_rule="kac", K=2, steps=3000, sample_every=30, tau_scale=0.25,
              seed=0, trace=False):
    """
    与 z0_local.py 同一引擎（全分支 N←NA、归一化、闭合写记录、每顶点自己的时钟、
    局部清空+从该顶点记录层重播种 cap K），**只把时钟换成可插拔的**：
      kac      τ_v = round(tau_scale·2|E|/deg(v))                    （§8 基线，静态）
      history  τ_v ← 该顶点**上次两次点火的实测间隔**（自观测周期，局部历史）
      gain     点火条件换成"本循环攒到的记录 ≥ 上一循环攒到的记录"（记录增益匹配）
    """
    V = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel()
    E2 = float(deg.sum())
    tau = np.maximum(1, np.round(tau_scale * E2 / np.maximum(deg, 1))).astype(np.int64)
    N = np.eye(V) * (1.0 / V)
    clock = np.zeros(V, np.int64)
    record = np.zeros(V)
    term = np.zeros(V)
    last_fire = np.zeros(V, np.int64)
    gain_ref = np.zeros(V)          # 上一循环攒到的记录（gain 规则用）
    rec_at_fire = np.zeros(V)
    boot = np.zeros(V, bool)
    rows = []
    for t in range(steps):
        N = N @ A
        mx = N.max()
        if mx > 0:
            N = N / mx
        cl = np.diag(N).copy()
        if cl.sum() > 0:
            record += cl
            N[np.arange(V), np.arange(V)] = 0.0
        clock += 1
        if clock_rule == "gain":
            gain = record - rec_at_fire
            fire = np.nonzero((~boot) | (gain >= gain_ref))[0]
            boot[:] = True
        else:
            fire = np.nonzero(clock >= tau)[0]
        n_fire = len(fire)
        for v in fire:
            col = N[:, v].copy()
            if col.sum() > 0:
                term[v] += col.sum()
                term -= col
                N[:, v] = 0.0
            N[v, v] += min(record[v], K)
            if clock_rule == "history":
                if last_fire[v] > 0:
                    tau[v] = max(1, int(t + 1 - last_fire[v]))     # ★ 自观测周期
            elif clock_rule == "gain":
                gain_ref[v] = record[v] - rec_at_fire[v]
                rec_at_fire[v] = record[v]
            last_fire[v] = t + 1
            clock[v] = 0
        if (t + 1) % sample_every == 0:
            act = N.sum(axis=1)
            occ = act > act.mean() * 0.5
            rows.append({"t": t + 1,
                         "act_gini": round(gini(act), 5),
                         "act_maxfrac": round(float(act.max() / max(act.sum(), 1e-30)), 5),
                         "record_gini": round(gini(record), 5),
                         "clusters": int(clusters(occ, A)),
                         "n_fire": n_fire,
                         "tau_mean": float(tau.mean()), "tau_std": float(tau.std())})
    return rows, tau


def drift(rows, key):
    """后半段的漂移（末值 − 半程值）与后半段的相对起伏 ⟹ 判"是否到稳态/是否还在变"。"""
    v = np.array([r[key] for r in rows], float)
    h = len(v) // 2
    d = float(v[-1] - v[h])
    rel = float(v[h:].std() / max(abs(v[h:].mean()), 1e-30))
    return d, rel


# ================================================================ E3 词级引擎：干涉接进 Z4 的清空
def qL_exact(L):
    seen = {}
    for ones in combinations(range(L), L // 2):
        w = [1] * L
        for i in ones:
            w[i] = -1
        c = min(tuple(w[i:] + w[:i]) for i in range(L))
        seen[c] = seen.get(c, 0) + 1
    return sum(Fraction(o, math_comb(L, L // 2)) ** 2 * n for o, n in Counter(seen.values()).items())


def math_comb(n, k):
    from math import comb
    return comb(n, k)


def build_tower(L):
    from math import comb, log
    q = qL_exact(L)
    blk = Counter()
    for k in range(L + 1):
        for s in range(-k, k + 1, 2):
            rem = L - k
            if (rem - s) % 2:
                continue
            a = (rem - s) // 2
            if 0 <= a <= rem:
                blk[abs(s)] += comb(rem, a) * (float(q) ** abs(s))
    hs = sorted(blk)
    ws = np.array([blk[h] for h in hs], float)
    ws /= ws.sum()
    freq = np.array([-np.log(w / 2) for w in ws])
    return q, hs, ws, freq


def word_engine(L=4, T=20000, mode="plain"):
    """
    mode="plain"     ：§19 的现引擎（种子 = 闭合分支数）⟹ 已知：6 步后精确 2-周期
    mode="phase_clear"：**构造**（非导出）：闭合步若两模**相消**（Re[a1·conj(a2)]<0，规范不变），
                        就在该步清掉**弱模载体**（偶数拍的 (++),(--)）——Z4 的清空由干涉选择；
                        种子仍 = 闭合分支数 ⟹ r*=1 临界不变（既不爆炸也不死亡）。
    """
    q, hs, ws, freq = build_tower(L)
    L2 = Counter({(): 1})
    g = 0
    traj, comp, acts = [], [], []
    cleared = 0
    for _ in range(T):
        newL2 = Counter()
        nclo = 0
        amp = np.zeros(len(hs), dtype=complex)
        for w, mult in L2.items():
            if mult <= 0:
                continue
            for d in (+1, -1):
                nw = w + (d,)
                if sum(nw) == 0 and len(nw) >= 2:
                    nclo += 1
                    h = max(abs(sum(nw[: k + 1])) for k in range(len(nw)))
                    hh = hs.index(h) if h in hs else len(hs) - 1
                    amp[hh] += mult * ws[hh] * np.exp(-1j * freq[hh] * g)
                elif len(nw) < L:
                    newL2[nw] += mult
        if nclo > 0:
            newL2[()] += nclo
            if mode == "phase_clear":
                nz = [v for v in amp if abs(v) > 0]
                if len(nz) >= 2 and (nz[0] * np.conj(nz[1])).real < 0:
                    # 相消 ⟹ 清掉弱模（h=2）的载体 (++),(--)
                    for w in ((1, 1), (-1, -1)):
                        if newL2.get(w, 0) > 0:
                            del newL2[w]
                    cleared += 1
            g += 1
        L2 = Counter({w: m for w, m in newL2.items() if m > 0})
        traj.append(tuple(sorted(L2.items())))
        comp.append((sum(m for w, m in L2.items() if len(w) == 2),
                     sum(m for w, m in L2.items() if len(w) == 1)))
        acts.append(sum(L2.values()))
    return dict(traj=traj, comp=np.array(comp), acts=np.array(acts), cleared=cleared, g=g)


def state_period(traj, cap=400000):
    seen = {}
    for i, s in enumerate(traj):
        if s in seen:
            return seen[s] + 1, i - seen[s]
        seen[s] = i
        if i > cap:
            return None, None
    return None, None


def honest_period(traj, min_win=200):
    """诚实判周期：找到第一对重复状态 (i,j) 后**验证未来是否逐位一致**。
    只看'第一次重复'会把'状态投影重复但相位不同'误判成周期（§20 的 E3 就是这样）。"""
    seen = {}
    for k, s in enumerate(traj):
        if s in seen:
            i = seen[s]
            win = min(min_win, len(traj) - k)
            same = all(traj[i + d] == traj[k + d] for d in range(win))
            return i + 1, k - i, same, win
        seen[s] = k
    return None, None, None, 0


# ================================================================ 主程序
def part_E1():
    print("=" * 100)
    print("E1（候选②）把 τ_v 从 Kac 静态换成**局部历史** —— Γ = 环面16²＋10% 枢纽×4（与 §8 同基线）")
    base = torus_adj(2, 16)
    A = add_hubs(base, frac=0.10, extra=4, seed=3)
    rows_out = {}
    print(f"{'时钟':>10} {'末活动Gini':>10} {'后半漂移':>9} {'后半起伏':>9} {'末记录Gini':>10} "
          f"{'末簇数':>7} {'簇数范围':>9} {'τ末均值':>9} {'τ末std':>9}")
    for rule in ("kac", "history", "gain"):
        rows, tau = run_local(A, clock_rule=rule, K=2, steps=3000, sample_every=30)
        d, rel = drift(rows, "act_gini")
        cl = [r["clusters"] for r in rows]
        last = rows[-1]
        print(f"{rule:>10} {last['act_gini']:10.4f} {d:9.4f} {rel:9.4f} {last['record_gini']:10.4f} "
              f"{last['clusters']:7d} {f'{min(cl)}-{max(cl)}':>9} {last['tau_mean']:9.2f} {last['tau_std']:9.2f}")
        rows_out[rule] = dict(rows=rows, tau_last_mean=last["tau_mean"], tau_last_std=last["tau_std"],
                              gini_last=last["act_gini"], clusters_last=last["clusters"],
                              clusters_range=[int(min(cl)), int(max(cl))], drift=d, rel_std=rel)
    RES["E1_history_clocks"] = rows_out
    return rows_out


def part_E2():
    print("=" * 100)
    print("E2（候选①）层级 Γ：先量 (A^L)_vv 是否**按层分块**，再看局部异步引擎的簇数")
    A_h, mod_of, grp_of = build_hierarchical()
    base = torus_adj(2, 16)
    A_f = add_hubs(base, frac=0.10, extra=4, seed=3)
    L = 4
    tab = []
    for name, A in (("层级Γ(32团×8)", A_h), ("平坦基线(环面+枢纽)", A_f)):
        deg = np.asarray(A.sum(1)).ravel()
        AL = (A ** L).toarray() if A.shape[0] <= 300 else None
        dg = np.diag(AL)
        # 按层分组统计回返概率
        if name.startswith("层级"):
            grp_mean = [float(dg[grp_of == g].mean()) for g in range(grp_of.max() + 1)]
            grp_std = [float(dg[grp_of == g].std()) for g in range(grp_of.max() + 1)]
        else:
            grp_mean, grp_std = [float(dg.mean())], [float(dg.std())]
        between = float(np.std(grp_mean) / max(abs(np.mean(grp_mean)), 1e-30))
        print(f"  {name:20s} V={A.shape[0]} 度[{deg.min():.0f},{deg.max():.0f}] "
              f"(A^{L})_vv ∈ [{dg.min():.0f},{dg.max():.0f}] 均值={dg.mean():.1f} "
              f"层间/层内离散度比={between:.4f}")
        rows, tau = run_local(A, clock_rule="kac", K=2, steps=3000, sample_every=30, tau_scale=0.25)
        cl = [r["clusters"] for r in rows]
        d, rel = drift(rows, "act_gini")
        last = rows[-1]
        print(f"    → 末活动Gini={last['act_gini']:.4f} 后半漂移={d:+.4f} 后半起伏={rel:.4f} "
              f"末簇数={last['clusters']} 簇数范围={min(cl)}-{max(cl)} 记录Gini={last['record_gini']:.4f}")
        tab.append(dict(name=name, V=int(A.shape[0]), deg_min=int(deg.min()), deg_max=int(deg.max()),
                        diagL_mean=float(dg.mean()), diagL_min=float(dg.min()), diagL_max=float(dg.max()),
                        level_between_disp=between, gini_last=last["act_gini"], drift=d, rel_std=rel,
                        clusters_last=last["clusters"], clusters_range=[int(min(cl)), int(max(cl))],
                        record_gini_last=last["record_gini"], rows=rows))
    RES["E2_hierarchical"] = tab
    return tab


def part_E3():
    print("=" * 100)
    print("E3（§19.8 杠杆）把干涉接进 Z4 的清空选择（构造·非导出）；种子仍 = 闭合分支数（r*=1 不变）")
    for mode in ("plain", "phase_clear"):
        r = word_engine(L=4, T=20000, mode=mode)
        tr, per = state_period(r["traj"])
        htr, hper, hsame, hwin = honest_period(r["traj"])
        acts = r["acts"]
        comp = r["comp"]
        nz = comp[comp.sum(axis=1) > 0]
        # 组分（长度1载体 vs 长度2载体）的分布
        vals = Counter(map(tuple, nz.tolist()))
        print(f"  mode={mode:12s} 朴素(瞬态/周期)={str((tr, per)):>12} "
              f"诚实(首次重复/间隔/未来一致?)={(htr, hper, hsame)} "
              f"活动: 末={acts[-1]} 峰={acts.max()} 清空={r['cleared']:5d} 组分取值数={len(vals):3d}")
        RES.setdefault("E3_interference_clear", {})[mode] = dict(
            transient=tr, period=per,
            honest_first_repeat=htr, honest_gap=hper, honest_future_matches=hsame,
            act_last=int(acts[-1]), act_max=int(acts.max()),
            cleared=int(r["cleared"]), n_distinct_composition=len(vals),
            top_compositions=[[list(k), int(v)] for k, v in vals.most_common(5)])
    return RES["E3_interference_clear"]


if __name__ == "__main__":
    t0 = time.time()
    part_E1()
    part_E2()
    part_E3()
    json.dump(RES, open(os.path.join(OUT, "z0_evolve.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print("=" * 100)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_evolve.json")
