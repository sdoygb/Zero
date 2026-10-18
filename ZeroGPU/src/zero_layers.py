#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_layers.py —— L0/L1/L2 的准确表达与“量之间的相互影响”
=========================================================
上一轮 `zero_evo.py` 把多个扇区放进一个状态；这一轮专做两件事：

  A. **层账**：L0（舞台图，0 层）／L1（高度层 H=L/2，静态）／L2（3 频率槽）。
     其中 L2 的“3 层”是 **3 个独立相邻比**（4 顶点占据的归一化自由度），
     不是 4 个顶点——这是 Z129 §5.1 的口径必须说清的一处。
  B. **相互影响**：给每个量一个因果探针（脉冲响应矩阵），并测 L1 种子
     → L2 槽、宏观扇区 → L2 槽的响应。

核心发现（先说）：
  1. **闭合口径分叉**。`persistent`（Z020 的 (起点,当前位置) 记法）与
     `fresh`（词层口径：每代从当前位置起一条新的长度 L 词）给出的 L1/L2
     末态不同：
       · persistent K4：record=seed 比 [0.5517,0.2759,0.1034,0.069]，
         3 槽 [0.4055,0.6931,0.9808]，净环流 0.540（复现 Z129）；
       · fresh K4：record=[0.25,0.25,0.25,0.25]，3 槽全 0，净环流 0。
     词层口径的 L1 高度分布与静态 b_h **逐位相符**（cycle16：L=4 给
     [2/3,1/3]；L=8 给 [0.2286,0.5429,0.2000,0.0286]）；persistent 口径
     的高度分布被跨代长走者偏置（L=4 给 [0.374,0.626]）。
  2. **L1 高度层是静态读出**：在现有公理/引擎里它没有自身更新律；唯一
     的 L1 动力学是记录累积（边③）。把顶点 i 硬映射到高度层 i mod H 的
     做法（旧 `l1_height_layer_engine.py`）没有台账依据。
  3. **L2 槽的相互影响**：persistent 下种子→槽是锁死的（比例继承）；
     fresh 下种子被图混合，槽退化为 0；打开宏观回写（几何/相位/容量）后
     fresh 的槽可以被重新激发。脉冲响应矩阵给出各量的因果强度。

输出：results/zero_layers.json、zero_layers_l1.png、zero_layers_l2.png、
      zero_layers_influence.png

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_layers.py
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo, Params, lap_matrix  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
OUT: dict = {}
CK: list[dict] = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════
#  L0 / L1 / L2 的静态账
# ══════════════════════════════════════════════════════════════════
def l0_summary(adj) -> dict:
    V = len(adj)
    E = sum(len(a) for a in adj) // 2
    deg = [len(a) for a in adj]
    L = lap_matrix(adj)
    w = np.linalg.eigvalsh(L)
    nz = w[w > 1e-9]
    g = zs.girth(adj)
    return {"V": V, "E": E, "deg_min": min(deg), "deg_max": max(deg),
            "girth": (int(g) if math.isfinite(g) else None),
            "lam_min": float(nz.min()) if nz.size else 0.0,
            "lam_max": float(nz.max()) if nz.size else 0.0}


def l1_static(L: int) -> dict:
    H = L // 2
    dist = zs.layer_distribution(L)
    b = [int(dist.get(h, 0)) for h in range(1, H + 1)]
    return {"L": L, "H": H, "b_h": b, "total": int(sum(b)),
            "balanced": int(math.comb(L, L // 2)),
            "all_even": all(x % 2 == 0 for x in b),
            "top_is_2": (b[-1] == 2 if b else False)}


def height_counts(L: int) -> dict:
    b = Counter()
    for w in itertools.product([1, -1], repeat=L):
        if sum(w) != 0:
            continue
        s = mx = 0
        for x in w:
            s += x
            mx = max(mx, abs(s))
        b[mx] += 1
    return dict(sorted(b.items()))


# ══════════════════════════════════════════════════════════════════
#  L2：4 顶点占据 → 3 相邻比（3 频率槽）
# ══════════════════════════════════════════════════════════════════
def l2_metrics(record, W=None, eps=1e-300) -> dict:
    r = np.asarray(record, dtype=float)
    p = r / r.sum() if r.sum() > 0 else np.ones_like(r) / len(r)
    ratios = [float(p[i] / max(eps, p[i + 1])) for i in range(len(p) - 1)]
    freqs = [float(abs(math.log(x))) for x in ratios]
    out = {"p": [float(x) for x in p], "ratios": ratios, "freqs": freqs,
           "n_slots": len(freqs), "slot_spread": float(max(freqs) - min(freqs))}
    if W is not None:
        W = np.array(W, dtype=float)
        np.fill_diagonal(W, 0.0)
        S = float(np.abs(W).sum())
        A = (W - W.T) / 2.0
        netout = [float(A[i].sum()) for i in range(W.shape[0])]
        out.update({"W_norm": (W / S).tolist() if S > 0 else W.tolist(),
                    "circ": float(np.abs(A).sum() / S) if S > 0 else 0.0,
                    "netout": netout})
    return out


# ══════════════════════════════════════════════════════════════════
#  高度层：fresh（词层口径）与 persistent（Z020 口径）
# ══════════════════════════════════════════════════════════════════
def _step_height(F, V, L, kind, Hmax):
    G = np.zeros_like(F)
    for u in range(V):
        for si in range(2 * L + 1):
            s = si - L
            row = F[u, si]
            if not row.any():
                continue
            for h in range(Hmax + 1):
                c = row[h]
                if not c:
                    continue
                for ds in (1, -1):
                    s2 = s + ds
                    if abs(s2) > L:
                        continue
                    if kind == "path":
                        w = u + s2
                        if w < 0 or w >= V:
                            continue
                    h2 = max(h, abs(s2))
                    G[u, s2 + L, h2] += c
    return G


def height_fresh(V, L, gens, seed0, kind="cycle", T=4):
    """词层口径：每代重新起词；返回闭合高度直方图（h=0..H）。"""
    H = L // 2
    Hmax = L
    s0 = np.array(seed0, dtype=float)
    rec_h = np.zeros(Hmax + 1)
    rec_u = np.zeros(V)
    layers = []
    for g in range(gens):
        F = np.zeros((V, 2 * L + 1, Hmax + 1))
        F[:, L, 0] += s0
        for _ in range(L):
            F = _step_height(F, V, L, kind, Hmax)
        c = F[:, L, :].copy()
        end = np.zeros(V)
        for u in range(V):
            for si in range(2 * L + 1):
                s = si - L
                row = F[u, si]
                if not row.any():
                    continue
                w = (u + s) % V if kind == "cycle" else u + s
                if 0 <= w < V:
                    end[w] += row.sum()
        rec_h += c.sum(axis=0)
        rec_u += c.sum(axis=1)
        clo_u = c.sum(axis=1)
        surv = end - clo_u
        layers.append(clo_u.copy())
        if len(layers) > 2:
            layers.pop(0)
        seedn = (layers[-1] + layers[-2]) if len(layers) >= 2 else layers[-1]
        if seedn.sum() == 0:
            seedn = np.ones(V)
        if T and (g + 1) % T == 0:
            surv = np.zeros(V)
        s0 = surv + seedn
        m = max(s0.sum(), rec_h.sum(), rec_u.sum(), F.sum())
        if m > 1e200:
            s0 /= m
            rec_h /= m
            rec_u /= m
            layers = [lay / m for lay in layers]
    return rec_h[:H + 1], rec_u


def height_persistent(V, L, gens, seed0, kind="cycle", T=4):
    """Z020 口径：走者保留原起点；高度是整条长走者的高度（会偏置）。"""
    H = L // 2
    Hmax = L
    F = np.zeros((V, 2 * L + 1, Hmax + 1))
    seed = np.array(seed0, dtype=float)
    rec_h = np.zeros(Hmax + 1)
    layers = []
    for g in range(gens):
        F[:, L, 0] += seed
        for _ in range(L):
            F = _step_height(F, V, L, kind, Hmax)
        c = F[:, L, :].copy()
        rec_h += c.sum(axis=0)
        clo_u = c.sum(axis=1)
        F[:, L, :] = 0.0
        layers.append(clo_u.copy())
        if len(layers) > 2:
            layers.pop(0)
        seed = (layers[-1] + layers[-2]) if len(layers) >= 2 else layers[-1]
        if seed.sum() == 0:
            seed = np.ones(V)
        m = max(seed.sum(), F.sum(), rec_h.sum())
        if m > 1e200:
            F /= m
            seed /= m
            rec_h /= m
            layers = [lay / m for lay in layers]
        if T and (g + 1) % T == 0:
            F[:] = 0.0
    return rec_h[:H + 1]


def exp_l1_area():
    """L1 的第二个再分坐标：面积 A(w)=Σ_k 前缀和_k（Z048）。"""
    rows = []
    for L in (4, 8, 12, 16):
        words = zs.balanced_words(L)
        areas = [zs.word_area(w) for w in words]
        rows.append({"L": L, "n": len(words),
                     "zero": sum(1 for a in areas if a == 0),
                     "min": min(areas), "max": max(areas),
                     "hist": {str(k): int(v) for k, v in
                              sorted(Counter(areas).items())}})
    return rows


# ══════════════════════════════════════════════════════════════════
#  实验 1：闭合口径分叉（K4 / path16）
# ══════════════════════════════════════════════════════════════════
def exp_semantics(gens_k4=60000, gens_path=40000):
    out = {}
    K4 = zs.k4()
    for sem in ("persistent", "fresh"):
        p = Params(micro=True, micro_semantics=sem, L=4, T=4,
                   retention="top2", wave_on=False, poisson_on=False)
        e = ZeroEvo(K4, p, seed=np.array([8.0, 4.0, 1.5, 1.0]))
        e.run(gens_k4, every=1000)
        e.step()                      # 让 W 非零（非毁灭代）
        rec = e.record / e.record.sum()
        m = l2_metrics(e.record, e.W)
        out[f"K4_{sem}"] = {
            "record": m["p"], "freqs": m["freqs"], "circ": m["circ"],
            "netout": m["netout"], "lock_seed": float(np.linalg.norm(
                e.seed / e.seed.sum() - rec)),
            "lock_closure": float(np.linalg.norm(
                e.gen_closures / max(1e-30, e.gen_closures.sum()) - rec)),
            "closure": l2_metrics(e.gen_closures)["p"],
            "closure_freqs": l2_metrics(e.gen_closures)["freqs"],
            "W_norm": m["W_norm"],
        }
        print(f"  K4 {sem:10s} p={np.round(rec,4)} slots={np.round(m['freqs'],4)} "
              f"circ={m['circ']:.4f}")

    path = zs.path(16)
    for sem in ("persistent", "fresh"):
        p = Params(micro=True, micro_semantics=sem, L=4, T=4,
                   retention="top2", wave_on=False, poisson_on=False)
        e = ZeroEvo(path, p, seed=np.array([2.0 ** (15 - i) for i in range(16)]))
        e.run(gens_path, every=500)
        rec = e.record / e.record.sum()
        support = int((rec > 1e-3).sum())
        top = np.argsort(-rec)[:4]
        out[f"path16_{sem}"] = {
            "record": [float(x) for x in rec], "support": support,
            "top4": [int(i) for i in top],
            "top_ratio": float(rec[top[0]] / max(1e-30, rec[top[1]])),
        }
        print(f"  path16 {sem:10s} support={support} top4={top} "
              f"ratio={rec[top[0]]/max(1e-30,rec[top[1]]):.3f}")
    return out


# ══════════════════════════════════════════════════════════════════
#  实验 2：L1 高度层（fresh vs persistent vs 静态 b_h）
# ══════════════════════════════════════════════════════════════════
def exp_l1(gens=200):
    out = {"static": [], "height": []}
    for L in (4, 8, 16):
        out["static"].append(l1_static(L))
    for kind, V in (("cycle", 16), ("path", 16)):
        for L in (4, 8):
            b = np.array(l1_static(L)["b_h"], dtype=float)
            b = b / b.sum()
            rh_f, _ = height_fresh(V, L, gens, np.ones(V), kind)
            rh_p = height_persistent(V, L, gens, np.ones(V), kind)
            hf = rh_f[1:] / rh_f[1:].sum()
            hp = rh_p[1:] / rh_p[1:].sum()
            out["height"].append({
                "kind": kind, "L": L, "H": L // 2,
                "static": [float(x) for x in b],
                "fresh": [float(x) for x in hf],
                "persistent": [float(x) for x in hp],
                "fresh_max_dev": float(np.abs(hf - b).max()),
                "persistent_max_dev": float(np.abs(hp - b).max()),
            })
            print(f"  L1 {kind:5s} L={L}: static={np.round(b,4)} "
                  f"fresh={np.round(hf,4)} persistent={np.round(hp,4)}")
    return out


# ══════════════════════════════════════════════════════════════════
#  实验 3：L2 槽的相互影响（seed / 宏观扇区）
# ══════════════════════════════════════════════════════════════════
def _slot_seed_response(base_seed, sem="persistent", T=4, gens=5000):
    K4 = zs.k4()
    p = Params(micro=True, micro_semantics=sem, L=4, T=T,
               retention="top2", wave_on=False, poisson_on=False)
    e = ZeroEvo(K4, p, seed=np.array(base_seed, dtype=float))
    e.run(gens, every=500)
    m = l2_metrics(e.record, e.W)
    return m


def _slot_trace(sem, params, gens, seed):
    K4 = zs.k4()
    e = ZeroEvo(K4, params, seed=np.array(seed, dtype=float))
    tr = []
    for _ in range(gens):
        e.step()
        tr.append(l2_metrics(e.gen_closures)["freqs"])
    return tr


def exp_l2(seed=(8.0, 4.0, 1.5, 1.0)):
    out = {"seed_response": [], "macro": []}
    base_p = _slot_seed_response(seed, "persistent")
    base_f = _slot_seed_response(seed, "fresh")
    out["baseline"] = {"persistent": base_p, "fresh": base_f}
    # persistent：逐顶点 ±10% 种子扰动 → 3 槽变化
    for i in range(4):
        for sgn in (+1, -1):
            s = list(seed)
            s[i] *= (1.0 + sgn * 0.10)
            m = _slot_seed_response(s, "persistent")
            out["seed_response"].append({
                "vertex": i, "sign": sgn,
                "d_freqs": [m["freqs"][j] - base_p["freqs"][j] for j in range(3)],
                "d_p": [m["p"][j] - base_p["p"][j] for j in range(4)],
            })
    # fresh：打开宏观回写，看槽能否被重新激发
    K4 = zs.k4()
    base = dict(micro=True, micro_semantics="fresh", L=4, T=4,
                retention="top2", eta_n=0.2, eta_h=0.2, eta_psi=0.2,
                mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
                capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
                Gamma_h=0.2, kappa_h=0.5, alpha=0.0)
    cases = {
        "fresh_no_macro": dict(eta_n=0.0, eta_h=0.0, eta_psi=0.0,
                               g_src=0.0, kappa_h=0.0),
        "fresh_geom": dict(),
        "fresh_phase": dict(beta=0.5, g_wave=0.5, kappa_psi=0.5,
                            Gamma_w=0.5, wave_on=True),
        "fresh_capacity": dict(capacity="hard", K_hard=1.0, reset_amp=0.05),
    }
    rng = np.random.default_rng(7)
    traces = {}
    for name, kw in cases.items():
        p = Params(**{**base, **kw})
        x0 = rng.normal(size=4)
        x0 -= x0.mean()
        x0 *= 0.05 / max(1e-12, np.abs(x0).max())
        e = ZeroEvo(K4, p, seed=np.array(seed), n0=1.0 + x0)
        e.run(10000, every=500)
        m = l2_metrics(e.record, e.W)
        out["macro"].append({"case": name, "record": m["p"],
                             "freqs": m["freqs"], "circ": m["circ"],
                             "lock": float(np.linalg.norm(
                                 e.gen_closures / max(1e-30, e.gen_closures.sum())
                                 - e.record / e.record.sum())),
                             "n_std": float(e.n.std()), "events": len(e.events)})
        traces[name] = _slot_trace("fresh", p, 400, seed)
        print(f"  L2 {name:18s} p={np.round(m['p'],4)} "
              f"slots={np.round(m['freqs'],4)} circ={m['circ']:.4f}")
    # persistent baseline slot trace
    pp = Params(micro=True, micro_semantics="persistent", L=4, T=4,
                retention="top2", wave_on=False, poisson_on=False)
    traces["persistent"] = _slot_trace("persistent", pp, 400, seed)
    out["slot_trace"] = traces
    return out


# ══════════════════════════════════════════════════════════════════
#  实验 4：宏观量的脉冲响应矩阵（cycle16）
# ══════════════════════════════════════════════════════════════════
OBS = ["n_std", "q", "h_energy", "psi_norm", "psi_coh",
       "pr", "n_mean", "fiedler", "phi_energy"]


def _observables(e):
    s = e.snapshot()
    return np.array([s[k] for k in OBS], dtype=float)


def _one_phase(e_phase, steps, eps, probes):
    base = e_phase.clone()
    base_obs = []
    for _ in range(steps):
        base.step()
        base_obs.append(_observables(base))
    base_obs = np.asarray(base_obs)
    bstd = base_obs.std(axis=0)
    bstd[bstd == 0] = 1.0
    mats, raws = [], []
    for name, fn in probes.items():
        q = e_phase.clone()
        fn(q)
        obs = []
        for _ in range(steps):
            q.step()
            obs.append(_observables(q))
        obs = np.asarray(obs)
        d = np.abs(obs - base_obs).max(axis=0) / eps
        raws.append(d)
        mats.append(d / bstd)
    return np.asarray(mats), np.asarray(raws), bstd


def exp_influence(steps=5000, eps=1e-4):
    adj = zs.cycle(16)
    V = len(adj)
    w, Q = np.linalg.eigh(lap_matrix(adj))
    f = Q[:, 1] / np.linalg.norm(Q[:, 1])
    p = Params(micro=False, mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
               capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.7,
               Gamma_h=0.2, kappa_h=0.5, beta=0.5, g_wave=0.5,
               kappa_psi=0.5, Gamma_w=0.5, c_w=0.5, omega0=0.2,
               wave_on=True, poisson_on=False)
    e = ZeroEvo(adj, p, n0=1.0 + 0.3 * f,
                psi0=np.exp(2j * np.pi * np.arange(V) / V))
    e.run(100000, every=1000)
    from zero_evo_long import lyap_local
    lam = lyap_local(e, nsteps=3000, eps=1e-8)
    root = e.clone()

    probes = {
        "n_mode": lambda x: setattr(x, "n", x.n + eps * f),
        "n_point": lambda x: x.n.__setitem__(0, x.n[0] + eps),
        "h_point": lambda x: x.h.__setitem__(0, x.h[0] + eps),
        "hv_point": lambda x: x.hv.__setitem__(0, x.hv[0] + eps),
        "psi_phase": lambda x: setattr(x, "psi", x.psi * np.exp(1j * eps)),
        "psi_amp": lambda x: setattr(x, "psi", x.psi * (1.0 + eps)),
    }
    offsets = (0, 800, 1600, 2400)
    mats, raws, stds = [], [], []
    for off in offsets:
        ep = root.clone()
        for _ in range(off):
            ep.step()
        m, raw, bstd = _one_phase(ep, steps, eps, probes)
        mats.append(m)
        raws.append(raw)
        stds.append(bstd)
    M = np.mean(mats, axis=0)
    S = np.std(mats, axis=0)
    R = np.mean(raws, axis=0)
    bstd = np.mean(stds, axis=0)
    mat = []
    for i, name in enumerate(probes):
        mat.append({"probe": name,
                    "response": [float(x) for x in R[i]],
                    "normalized": [float(x) for x in M[i]],
                    "normalized_std": [float(x) for x in S[i]]})
        print(f"  influence {name:10s} " +
              " ".join(f"{OBS[j]}:{M[i,j]:.2f}" for j in range(len(OBS))))
    return {"steps": steps, "eps": eps, "observables": OBS,
            "matrix": mat, "lyapunov": float(lam),
            "phase_offsets": list(offsets),
            "baseline_std": [float(x) for x in bstd],
            "baseline_last": root.snapshot()}


# ══════════════════════════════════════════════════════════════════
#  画图
# ══════════════════════════════════════════════════════════════════
def make_plots(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # L1：每个 L 的高度分布，三口径并排
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    for row, kind in enumerate(("cycle", "path")):
        for col, L in enumerate((4, 8)):
            ax = axes[row, col]
            r = [x for x in res["l1"]["height"]
                 if x["kind"] == kind and x["L"] == L][0]
            H = r["H"]
            x = np.arange(H)
            for i, key in enumerate(("static", "fresh", "persistent")):
                ax.bar(x + (i - 1) * 0.26, r[key], width=0.24, label=key)
            ax.set_xticks(x)
            ax.set_xticklabels([f"h={h}" for h in range(1, H + 1)])
            ax.set_ylabel("weight")
            ax.set_title(f"L1 {kind}, L={L}")
            if row == 0 and col == 0:
                ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "zero_layers_l1.png"), dpi=150)
    plt.close(fig)

    # L2
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    ax = axes[0]
    labels = ["persistent", "fresh"]
    for i, k in enumerate(labels):
        ax.bar(np.arange(3) + i * 0.35, res["l2"]["baseline"][k]["freqs"],
               width=0.32, label=k)
    ax.set_xticks(np.arange(3) + 0.17)
    ax.set_xticklabels(["f(0/1)", "f(1/2)", "f(2/3)"])
    ax.set_ylabel("$|\\log(p_i/p_{i+1})|$")
    ax.set_title("(a) L2 slots: persistent vs fresh")
    ax.legend(fontsize=7)

    ax = axes[1]
    for key in ("persistent", "fresh_no_macro", "fresh_geom", "fresh_phase"):
        tr = np.asarray(res["l2"]["slot_trace"][key], dtype=float)
        ax.plot(np.mean(tr, axis=1), lw=1.0, label=key)
    ax.set_xlabel("generation")
    ax.set_ylabel("mean of 3 slot frequencies")
    ax.set_title("(b) slot dynamics (instantaneous closures)")
    ax.legend(fontsize=7)

    ax = axes[2]
    names = [r["case"] for r in res["l2"]["macro"]]
    for j in range(3):
        ax.plot(range(len(names)), [r["freqs"][j] for r in res["l2"]["macro"]],
                "o-", label=f"slot {j}")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("slot frequency")
    ax.set_title("(c) fresh slots under macro couplings")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "zero_layers_l2.png"), dpi=150)
    plt.close(fig)

    # influence
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    M = np.array([r["normalized"] for r in res["influence"]["matrix"]])
    im = ax.imshow(np.log10(1.0 + M), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(OBS)))
    ax.set_xticklabels(OBS, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels([r["probe"] for r in res["influence"]["matrix"]],
                       fontsize=8)
    ax.set_title("(d) impulse response: $\\log_{10}(1+\\Delta/\\sigma)$")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "zero_layers_influence.png"), dpi=150)
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════
def main():
    os.makedirs(RES, exist_ok=True)
    print("=" * 78)
    print("L0：舞台图（0 层）")
    print("=" * 78)
    OUT["l0"] = {}
    for name, adj in (("K4", zs.k4()), ("path16", zs.path(16)),
                      ("cycle16", zs.cycle(16)), ("Heawood", zs.heawood())):
        OUT["l0"][name] = l0_summary(adj)
        print(f"  {name:8s} {OUT['l0'][name]}")

    print("=" * 78)
    print("L1：高度层（H=L/2，静态）")
    print("=" * 78)
    OUT["l1"] = exp_l1(200)
    OUT["l1_area"] = exp_l1_area()
    print("  面积零计数：", [(r["L"], r["zero"]) for r in OUT["l1_area"]])
    for r in OUT["l1"]["static"]:
        print(f"  L={r['L']:2d} H={r['H']} b_h={r['b_h']} "
              f"总={r['total']}=C(L,L/2)={r['balanced']}")

    print("=" * 78)
    print("口径分叉：persistent vs fresh（K4 / path16）")
    print("=" * 78)
    OUT["semantics"] = exp_semantics()

    print("=" * 78)
    print("L2 槽的相互影响")
    print("=" * 78)
    OUT["l2"] = exp_l2()

    print("=" * 78)
    print("宏观脉冲响应（cycle16，c_h=1.7，波扇区开）")
    print("=" * 78)
    OUT["influence"] = exp_influence()

    print("=" * 78)
    print("层账判据")
    print("=" * 78)
    l1h = OUT["l1"]["height"]
    ck("L1 fresh 在 cycle 上逐位等于 b_h",
       all(r["fresh_max_dev"] < 1e-12 for r in l1h if r["kind"] == "cycle"),
       "max dev = 0")
    ck("L1 fresh 在 path 上偏差 < 0.01",
       all(r["fresh_max_dev"] < 0.01 for r in l1h if r["kind"] == "path"),
       "max = " + f"{max(r['fresh_max_dev'] for r in l1h if r['kind']=='path'):.4f}")
    ck("L1 persistent 偏置 > 0.1",
       all(r["persistent_max_dev"] > 0.1 for r in l1h),
       "min = " + f"{min(r['persistent_max_dev'] for r in l1h):.4f}")
    k4p = OUT["semantics"]["K4_persistent"]
    k4f = OUT["semantics"]["K4_fresh"]
    ck("K4 persistent record=seed（锁死）", k4p["lock_seed"] < 1e-12,
       f"lock = {k4p['lock_seed']:.2e}")
    ck("K4 fresh record 均匀", max(abs(x - 0.25) for x in k4f["record"]) < 1e-12,
       f"max dev = {max(abs(x-0.25) for x in k4f['record']):.2e}")
    ck("path16 persistent 支撑 2 / 比 2.0",
       OUT["semantics"]["path16_persistent"]["support"] == 2
       and abs(OUT["semantics"]["path16_persistent"]["top_ratio"] - 2.0) < 1e-6,
       f"support={OUT['semantics']['path16_persistent']['support']}")
    ck("path16 fresh 支撑 16",
       OUT["semantics"]["path16_fresh"]["support"] == 16,
       f"support={OUT['semantics']['path16_fresh']['support']}")
    ds = [abs(x) for r in OUT["l2"]["seed_response"] for x in r["d_freqs"] if abs(x) > 1e-12]
    exp_vals = (0.095310, 0.105361)   # log(1.1) 与 |log(0.9)|
    dev = [min(abs(x - e) for e in exp_vals) for x in ds]
    ck("种子→槽响应 = ±log(1.1) / |log(0.9)|", bool(ds) and max(dev) < 1e-4,
       f"max dev = {max(dev):.2e}" if ds else "none")
    ck("L1 面积零计数 = 2/8/58/526（Z048）",
       [r["zero"] for r in OUT["l1_area"]] == [2, 8, 58, 526],
       str([r["zero"] for r in OUT["l1_area"]]))
    ck("fresh 基线 3 槽全 0",
       max(OUT["l2"]["baseline"]["fresh"]["freqs"]) < 1e-12,
       f"max = {max(OUT['l2']['baseline']['fresh']['freqs']):.2e}")
    OUT["checks"] = CK
    OUT["n_fail"] = sum(1 for c in CK if not c["pass"])
    print(f"合计 {len(CK)} 条：PASS {len(CK)-OUT['n_fail']} / FAIL {OUT['n_fail']}")

    make_plots(OUT)
    OUT["meta"] = {"engine": "zero_evo.py", "date": "2026-10-18"}
    fp = os.path.join(RES, "zero_layers.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1)
    print(f"\n写入 {fp}")
    print("图：zero_layers_l1.png / zero_layers_l2.png / zero_layers_influence.png")


if __name__ == "__main__":
    main()
