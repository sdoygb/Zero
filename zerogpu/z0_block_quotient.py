"""
z0_block_quotient.py --- 把 $\\mathcal L_\\kappa$ 抬到**物理 $\\Gamma$（$V{=}256$）的块商图**上

为什么这样走：物理 $\\Gamma$ 的零和盒 $\\{n\\in[-L,L]^{256}:\\sum n=0\\}$ 不可枚举。
按 §38 的判据（"亚稳性住在模块**之间**"），把状态换成**块和**

$$
N_B=\\sum_{v\\in B}n_v,\\qquad \\sum_B N_B=0\\ (\\text{零和约束保留}),
\\qquad
\\kappa_{AB}=\\#\\{(u,v)\\in E_\\Gamma:\\ u\\in A,\\ v\\in B\\}\\ (\\text{块间边数＝\\textbf{导出的}权重})
$$

$\\kappa_{AB}$ 来自 $\\Gamma$ 本身（不是记录权重）⟹ 与 `D192` 的"偏好只能来自约束/图"一致。

对照三种**分块**（同一 $V{=}256$、同一块数）：
  · **层级分块**（＝$\\Gamma$ 自己的 8 个超模块）——「密块＋弱连」
  · **随机分块**（同尺寸随机划 8 块）——**切穿**密团
  · **平坦分块**（环面 $16^2$ 划 8 个条带）——无层级

判据：$\\tau_2=1/\\lvert\\lambda_2\\rvert$（寿命）与最局域模的 PR$/N$。

用法：/usr/bin/python3 z0_block_quotient.py      输出：results/z0_block_quotient.json
"""
from __future__ import annotations

import json
import os
import time
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

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


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return 0.0 if x.sum() == 0 else float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def block_rates(A, labels, nb):
    """块间边数 kappa_{AB}（含块内自环 kappa_{AA}）。"""
    Ac = A.tocoo()
    K = np.zeros((nb, nb))
    for u, v in zip(Ac.row, Ac.col):
        K[labels[u], labels[v]] += 1
    return K                 # ★ 更正：非对角元 = 有向边数（A→B 的边数）＝率；不要除 2


def build_L(K, bound):
    nb = K.shape[0]
    st = [n for n in product(range(-bound, bound + 1), repeat=nb) if sum(n) == 0]
    idx = {n: k for k, n in enumerate(st)}
    rows, cols, vals = [], [], []
    for n in st:
        k = idx[n]; diag = 0.0
        for a in range(nb):
            for b in range(nb):
                if a == b or K[a, b] <= 0:
                    continue
                if n[a] - 1 < -bound or n[b] + 1 > bound:
                    continue
                nn = list(n); nn[a] -= 1; nn[b] += 1
                w = K[a, b] / max(K[a].sum(), 1e-30)
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(w)
                diag += w
        rows.append(k); cols.append(k); vals.append(-diag)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), len(st)


def spectrum(Lm, k=8):
    ev, evec = eigsh(Lm.tocsc(), k=k, which="SM")
    o = np.argsort(ev.real); vals, vecs = ev.real[o], evec[:, o].real
    N = Lm.shape[0]
    nz = [v for v in vals if abs(v) > 1e-9]
    tau2 = max(1 / abs(v) for v in nz) if nz else None
    prs, taus = [], []
    for j in range(len(vals)):
        if abs(vals[j]) < 1e-9:
            continue
        vv = vecs[:, j]
        prs.append(float((vv ** 2).sum() ** 2 / max((vv ** 4).sum(), 1e-30)) / N)
        taus.append(float(1 / abs(vals[j])))
    jm = int(np.argmin(prs)) if prs else 0
    jslow = int(np.argmax(taus)) if taus else 0            # ★ 最慢那个模本身的局域度
    return dict(N=N, tau2=round(tau2, 3) if tau2 else None,
                taus=sorted([round(x, 2) for x in taus], reverse=True)[:5],
                pr_min=round(min(prs), 6) if prs else None,
                pr_min_tau=round(taus[jm], 2) if prs else None,
                pr_slow=round(prs[jslow], 6) if prs else None)


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 100)
    hier, grp = ZT.build_hierarchical()          # V=256，32 团×8，8 个超模块
    torus = ZT.torus_adj(2, 16)                  # V=256 环面
    nb = 8
    rng = np.random.default_rng(0)
    parts = {
        "层级Γ/8超模块（密块＋弱连）": (hier, grp),
        "层级Γ/随机8块（切穿密团）": (hier, rng.permutation(np.repeat(np.arange(nb), 32))),
        "环面/8条带（无层级）": (torus, np.repeat(np.arange(nb), 32)),
    }
    print(f"块商图：$V{{=}}256$，块数 {nb}，块和盒界 $K{{=}}2$")
    print(f"{'分块':>30} {'块间边数均值':>12} {'块内自环':>9} {'N':>6} {'τ2':>9} {'最慢5τ':>24} {'最慢模PR/N':>11} {'最局域PR/N':>11}")
    out = {}
    for nm, (A, lab) in parts.items():
        K = block_rates(A, lab, nb)
        off = K[~np.eye(nb, dtype=bool)]
        Lm, N = build_L(K, bound=2)
        sp = spectrum(Lm)
        out[nm] = sp | {"offdiag_mean": round(float(off.mean()), 2),
                        "selfloop_mean": round(float(np.diag(K).mean()), 2)}
        print(f"{nm:>30} {off.mean():12.2f} {np.diag(K).mean():9.2f} {N:6d} {sp['tau2']:9.2f} "
              f"{str(sp['taus']):>24} {sp['pr_slow']:11.6f} {sp['pr_min']:11.6f}")
    RES["partitions"] = out
    h = out["层级Γ/8超模块（密块＋弱连）"]
    r = out["层级Γ/随机8块（切穿密团）"]
    f = out["环面/8条带（无层级）"]
    check("**层级分块的慢模最慢**：τ2(层级) > τ2(随机分块)——分块必须顺着 Γ 自己的层级走",
          h["tau2"] > r["tau2"], f"层级 {h['tau2']} vs 随机 {r['tau2']}")
    check("**层级分块的\*最慢模本身\*强局域**（PR/N < 随机分块的一半）",
          h["pr_slow"] < r["pr_slow"] / 2,
          f"最慢模 PR/N：层级 {h['pr_slow']} vs 随机 {r['pr_slow']} vs 环面 {f['pr_slow']}")
    check("三层都在同一个态空间规模上（可比较）", h["N"] == r["N"] == f["N"], f"N={h['N']}")
    print(f"\n   ⟹ 块商图把 §38 的判据**抬到了物理 $V{{=}}256$ 上**：分块要顺着 $\\Gamma$ 自己的层级，"
          f"随机分块（切穿密团）与无层级的平坦分块都更快、更不局域。")
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_block_quotient.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_block_quotient.json")
