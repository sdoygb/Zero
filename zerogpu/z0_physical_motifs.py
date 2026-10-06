"""
z0_physical_motifs.py --- **物理 Γ 的结构基元**：哪一块给出"亚稳畴"？（任务 C 的可量化判据）

背景：§37 证了亚稳畴＝零和主方程 $\\mathcal L_\\kappa$ 的**局域慢模**，且 $\\Gamma$ 的层级决定寿命。
但物理 $\\Gamma$（$\\mathcal G_T$ $V{=}246$／层级图 $V{=}256$）的零和盒无法枚举。
本文件改测**它的结构基元**：固定 $V{=}8,\\ L{=}2$（态空间恒定 $38165$），只换 $\\Gamma$ —— 受控比较。

层级 Γ 的构造是「模块（团）＋模块间弱连」，所以基元正是：
  · $K_8$ ＝**一个模块**（团）        · 两个 $K_4$ 弱连 ＝**一对模块**
  · 两个 $K_3$＋悬点 ＝更松的层级      · 环 8 ＝正则对照
  · 星 $K_{1,7}$ ＝**只非均匀、不分块**的对照   · $Q_3$ ＝中间情形

判据：$\\tau_2=1/\\lvert\\lambda_2\\rvert$（最慢非平凡模的寿命）与最局域模的参与率 $/N$。

用法：/usr/bin/python3 z0_physical_motifs.py      输出：results/z0_physical_motifs.json
"""
from __future__ import annotations

import json
import os
import time
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
V, LB = 8, 2
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return 0.0 if x.sum() == 0 else float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def fams():
    K = lambda ids: [(a, b) for i, a in enumerate(ids) for b in ids[i + 1:]]
    return {
        "环8（正则对照）": [(i, (i + 1) % 8) for i in range(8)],
        "K8（＝层级的**模块**）": K(range(8)),
        "两个K4弱连（＝**模块对**）": K([0, 1, 2, 3]) + K([4, 5, 6, 7]) + [(3, 4)],
        "两个K3＋悬点（更松层级）": K([0, 1, 2]) + K([3, 4, 5]) + [(2, 3), (5, 6), (6, 7)],
        "星K1,7（只非均匀、不分块）": [(0, i) for i in range(1, 8)],
        "立方体Q3（中间）": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
                            (0, 4), (1, 5), (2, 6), (3, 7)],
    }


def build_L(eg):
    st = [n for n in product(range(-LB, LB + 1), repeat=V) if sum(n) == 0]
    idx = {n: k for k, n in enumerate(st)}
    rows, cols, vals = [], [], []
    for n in st:
        k = idx[n]; d = 0
        for (i, j) in eg:
            for (a, b) in ((i, j), (j, i)):
                if n[a] - 1 < -LB or n[b] + 1 > LB:
                    continue
                nn = list(n); nn[a] -= 1; nn[b] += 1
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0); d += 1
        rows.append(k); cols.append(k); vals.append(-d)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), st


def analyse(eg, k=8):
    Lm, st = build_L(eg)
    N = Lm.shape[0]
    ev, evec = eigsh(Lm.tocsc(), k=k, which="SM")
    o = np.argsort(ev.real)
    vals, vecs = ev.real[o], evec[:, o].real
    deg = np.asarray(-Lm.diagonal()).ravel()
    nz = [v for v in vals if abs(v) > 1e-9]
    tau2 = max(1 / abs(v) for v in nz)
    prs, taus = [], []
    for j in range(len(vals)):
        if abs(vals[j]) < 1e-9:
            continue
        vv = vecs[:, j]
        prs.append(float((vv ** 2).sum() ** 2 / max((vv ** 4).sum(), 1e-30)) / N)
        taus.append(float(1 / abs(vals[j])))
    jm = int(np.argmin(prs))
    return dict(N=N, edges=len(eg), pi_gini=round(gini(deg), 4), tau2=round(tau2, 2),
                taus=sorted([round(x, 2) for x in taus], reverse=True)[:5],
                pr_min=round(min(prs), 6), pr_min_tau=round(taus[jm], 2),
                pr_mean=round(float(np.mean(prs)), 4))


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 100)
    print(f"物理 Γ 的结构基元（V={V}, L={LB}，态空间恒定 N；只换 Γ）")
    print(f"{'Γ':>26} {'边':>4} {'N':>7} {'π Gini':>8} {'λ2':>9} {'τ2':>7} {'最慢5τ':>26} {'最局域PR/N':>11}")
    out = {}
    for nm, eg in fams().items():
        r = analyse(eg)
        out[nm] = r
        print(f"{nm:>26} {r['edges']:4d} {r['N']:7d} {r['pi_gini']:8.4f} {-1/r['tau2']:9.4f} "
              f"{r['tau2']:7.2f} {str(r['taus']):>26} {r['pr_min']:11.6f} (τ={r['pr_min_tau']})")
    RES["motifs"] = out
    k8 = out["K8（＝层级的**模块**）"]
    pair = out["两个K4弱连（＝**模块对**）"]
    hier2 = out["两个K3＋悬点（更松层级）"]
    star = out["星K1,7（只非均匀、不分块）"]
    ring = out["环8（正则对照）"]
    check("**模块本身（K8）最快**：τ2 < 1（密团内部没有亚稳结构）", k8["tau2"] < 1.0, f"τ2={k8['tau2']}")
    check("**模块对最慢**：τ2 > 5（层级结构的来源在**模块之间**，不在模块内部）",
          pair["tau2"] > 5 and hier2["tau2"] > 5, f"模块对 {pair['tau2']}／更松层级 {hier2['tau2']}")
    check("**只非均匀、不分块（星图）不慢**：τ2 < 5 ⟹ 判据是「分块」不是「非均匀」",
          star["tau2"] < 5.0, f"星 τ2={star['tau2']}（π Gini 反而是最大 {star['pi_gini']}）")
    check("**模块对含极强局域模**（PR/N < 0.01）—— 这就是最小尺度的亚稳畴",
          pair["pr_min"] < 0.01, f"PR/N={pair['pr_min']}（≈{pair['pr_min']*pair['N']:.0f} 个态）")
    print(f"\n   ⟹ **任务 C（`L2-LATTICE-ORIGIN`）的可量化判据**：原生 Γ 要有**两层分块**（密块＋弱连），"
          f"使得 $\\mathcal L_\\kappa$ 出现 $\\tau_2\\gg1$ 的**局域**慢模；'非均匀'不够（星图不慢）。")
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_physical_motifs.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_physical_motifs.json")
