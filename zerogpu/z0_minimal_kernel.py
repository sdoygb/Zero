"""
z0_minimal_kernel.py --- 旧理论库的机制：**D192 最小核（零和图上的无偏好转移）＋ D191 转移模式**

来源（都在 `../modular-equilibrium/derivations/`，`lh` 语料**没有继承**）：
  · [`D192_zero_sum_conversion_dynamics`]：零和约束 ⟹ 允许的补偿转化 $n\\to n+e_j-e_i$；
    "分化＋零和＋**无偏好乱动**" ⟹ **最小转移核** $\\kappa_{ij}(n)=1$（**不加任何额外权重**）；
    零和图 $G_0$ 上的对称随机游走 ⟹ **度加权平稳分布** $\\pi\\propto\\deg$（非均匀、无自由参数）。
  · [`D191_dynamic_zero_modes_without_storage`]：零层 $\\mathcal Z=\\{Q=0\\}$ 上移动；
    每步按 $\\Delta B=B(y)-B(x)$（$B=b_1$ 圈数）分成 **建造／保持／消亡** 三集合
    ——「分类属于**转移**，不属于储存」。

**方法学更正（这正是我前几轮缺的）**：我之前把"畴"定义成**活动连通块**并追踪其**身份** ⟹ 寿命恒 1 代
（对象在重排）。D191 说结构应读在**转移**上：$\\Delta B$ 的统计——**建造/消亡是否自动平衡**。

本文件（全确定性、无概率、无自由参数）：
  ① 造零和图 $G_0$（状态 $=\\{n\\in[-L,L]^m:\\sum n=0\\}$，边＝沿 $\\Gamma$ 的补偿转化）
  ② 算 $\\deg$ 分布与平稳测度 $\\pi\\propto\\deg$（`G62` 要的非均匀 $\\pi$：**不需要额外输入**）
  ③ 算**全部边**的 $\\Delta B$ 统计：$P(\\text{建造})/P(\\text{保持})/P(\\text{消亡})$，
     并看是否自动平衡（＝`H10` 说的 SOC 式临界，但**无需调参**）

用法：/usr/bin/python3 z0_minimal_kernel.py      输出：results/z0_minimal_kernel.json
"""
from __future__ import annotations

import json
import math
import os
import time
from collections import Counter, defaultdict
from itertools import product

import numpy as np

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
    x = np.sort(np.asarray(x, float))
    n = len(x)
    return 0.0 if n == 0 or x.sum() == 0 else float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def b1(n, edges):
    """配置 n 的支撑子图的独立圈数 b1 = |E|-|V|+C。"""
    sup = {i for i, v in enumerate(n) if v != 0}
    if not sup:
        return 0
    E = [(u, v) for (u, v) in edges if u in sup and v in sup]
    # 连通块
    adj = defaultdict(list)
    for u, v in E:
        adj[u].append(v)
        adj[v].append(u)
    seen, C = set(), 0
    for s in sup:
        if s in seen:
            continue
        C += 1
        st = [s]
        seen.add(s)
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
    return len(E) - len(sup) + C


def build(m=6, L=3, ring=True):
    """零和图 G0：状态 = 零和盒 B_L ⊂ Z^m；边 = 沿 Γ 的补偿转化 n → n+e_j-e_i。"""
    edges_g = [(i, (i + 1) % m) for i in range(m)] if ring else \
              [(i, j) for i in range(m) for j in range(i + 1, m)]
    states = [n for n in product(range(-L, L + 1), repeat=m) if sum(n) == 0]
    idx = {n: k for k, n in enumerate(states)}
    deg = np.zeros(len(states), np.int64)
    E0 = []
    for n in states:
        k = idx[n]
        for (i, j) in edges_g:
            for (a, b) in ((i, j), (j, i)):
                if n[a] - 1 < -L or n[b] + 1 > L:
                    continue
                nn = list(n)
                nn[a] -= 1
                nn[b] += 1
                E0.append((k, idx[tuple(nn)]))
                deg[k] += 1
    return states, idx, deg, E0, edges_g


def gamma_families(m=6):
    """几个小 Γ（都小到可枚举零和盒）：正则环／带弦（非均匀）／两块弱连（层级）／完全图。"""
    ring = [(i, (i + 1) % m) for i in range(m)]
    fam = {"ring6(正则)": ring}
    fam["ring6+弦(非均匀)"] = ring + [(0, 3)]
    fam["两块弱连(层级)"] = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 5)]
    fam["K6(完全)"] = [(i, j) for i in range(m) for j in range(i + 1, m)]
    return fam


def analyse_graph(edges_g, m, L):
    states = [n for n in product(range(-L, L + 1), repeat=m) if sum(n) == 0]
    idx = {n: k for k, n in enumerate(states)}
    deg = np.zeros(len(states), np.int64)
    E0 = []
    for n in states:
        k = idx[n]
        for (i, j) in edges_g:
            for (a, b) in ((i, j), (j, i)):
                if n[a] - 1 < -L or n[b] + 1 > L:
                    continue
                nn = list(n); nn[a] -= 1; nn[b] += 1
                E0.append((k, idx[tuple(nn)])); deg[k] += 1
    Bs = np.array([b1(n, edges_g) for n in states], np.int64)
    act = np.array([sum(abs(x) for x in n) for n in states], np.int64)      # 活动量 Σ|n_v|
    dBs = np.array([Bs[b] - Bs[a] for (a, b) in E0], np.int64)
    dAs = np.array([act[b] - act[a] for (a, b) in E0], np.int64)
    tot = len(dBs)
    bp = float((dBs > 0).sum()) / tot; kp = float((dBs == 0).sum()) / tot; dp = float((dBs < 0).sum()) / tot
    pi = deg / deg.sum()
    src = np.array([a for a, b in E0])
    drift_B = float((pi[src] / deg[src] * dBs).sum())      # = (1/|E|)Σ dBs（因 π∝deg）
    drift_A = float((pi[src] / deg[src] * dAs).sum())
    return dict(nodes=len(states), edges=tot, deg_mean=float(deg.mean()), deg_max=int(deg.max()),
                pi_gini=round(gini(deg), 5), B_mean=float(Bs.mean()), B_max=int(Bs.max()),
                p_build=round(bp, 5), p_keep=round(kp, 5), p_die=round(dp, 5),
                ratio=round(bp / max(dp, 1e-12), 6), drift_B=round(drift_B, 12),
                drift_A=round(drift_A, 8), act_mean=float(act.mean()))


def build_L(edges_g, m, L):
    """无偏好核的**生成元** L（L_{b,a}=kappa=1，对角=-度）：主方程 dP/dt = L P 的确定性对象。"""
    from scipy.sparse import csr_matrix
    states = [n for n in product(range(-L, L + 1), repeat=m) if sum(n) == 0]
    idx = {n: k for k, n in enumerate(states)}
    rows, cols, vals = [], [], []
    for n in states:
        k = idx[n]; d = 0
        for (i, j) in edges_g:
            for (a, b) in ((i, j), (j, i)):
                if n[a] - 1 < -L or n[b] + 1 > L:
                    continue
                nn = list(n); nn[a] -= 1; nn[b] += 1
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0)
                d += 1
        rows.append(k); cols.append(k); vals.append(-d)
    return csr_matrix((vals, (rows, cols)), shape=(len(states), len(states)))


def spectrum(edges_g, m, L, k=8):
    from scipy.sparse.linalg import eigsh
    Lm = build_L(edges_g, m, L); N = Lm.shape[0]
    ev, evec = eigsh(Lm.tocsc(), k=k, which="SM")
    o = np.argsort(ev.real); vals = ev.real[o]; vecs = evec[:, o].real
    nz = [v for v in vals if abs(v) > 1e-9]
    tau = [float(1 / abs(v)) for v in nz]
    # 逐个模算参与率，找**最局域**的那个（亚稳结构 = 局域的慢模）
    prs, taus_all = [], []
    for j in range(len(vals)):
        if abs(vals[j]) < 1e-9:
            continue
        vv = vecs[:, j]
        prs.append(float((vv ** 2).sum() ** 2 / max((vv ** 4).sum(), 1e-30)) / N)
        taus_all.append(float(1 / abs(vals[j])))
    jmin = int(np.argmin(prs)) if prs else 0
    return dict(N=N, slowest_tau=max(tau) if tau else None, taus=sorted(tau, reverse=True)[:5],
                pr=round(min(prs), 5) if prs else None,
                pr_tau=round(taus_all[jmin], 3) if prs else None,
                pr_mean=round(float(np.mean(prs)), 5) if prs else None,
                lam=[round(float(x), 5) for x in vals])


def main():
    t0 = time.time()
    print("=" * 104)
    print("§36 处方的 Γ 族扫描：无偏好核（Z0③ 全分支）在零和层上 —— 平衡是否普遍？活动量有无漂移？")
    print(f"{'Γ':>18} {'节点':>7} {'边':>8} {'π Gini':>8} {'建造':>7} {'保持':>7} {'消亡':>7} "
          f"{'建:消':>8} {'⟨ΔB⟩':>10} {'⟨Δ(活动)⟩':>11}")
    out = {}
    for nm, eg in gamma_families(6).items():
        r = analyse_graph(eg, 6, 3)
        out[nm] = r
        print(f"{nm:>18} {r['nodes']:7d} {r['edges']:8d} {r['pi_gini']:8.4f} {r['p_build']:7.4f} "
              f"{r['p_keep']:7.4f} {r['p_die']:7.4f} {r['ratio']:8.4f} {r['drift_B']:+10.2e} {r['drift_A']:+11.2e}")
    RES["gamma_scan"] = out
    bk = list(out.values())
    check("**建造/消亡精确平衡对四个 Γ 族都成立**（|⟨ΔB⟩|<1e-9）",
          all(abs(r["drift_B"]) < 1e-9 for r in bk), "普遍性 = 反对称性的推论，与 Γ 无关")
    check("π∝deg 在四个 Γ 族上都非均匀（Gini>0.02）",
          all(r["pi_gini"] > 0.02 for r in bk), f"Gini 范围 [{min(r['pi_gini'] for r in bk)}, {max(r['pi_gini'] for r in bk)}]")
    check("**活动量 Σ|n_v| 也零漂移**（说明稳态不只圈数临界，活动量也临界）",
          all(abs(r["drift_A"]) < 1e-6 for r in bk),
          f"⟨ΔA⟩ 范围 [{min(r['drift_A'] for r in bk):+.2e}, {max(r['drift_A'] for r in bk):+.2e}]")
    # ---- 谱：亚稳结构住在生成元的慢模里（这才是"畴的生灭/寿命"）
    print(f"\n{'Γ':>18} {'N':>6} {'λ2(最慢)':>10} {'τ2(步)':>8} {'最慢5个 τ':>28} {'慢模参与率/N':>12}")
    spec = {}
    for nm, eg in gamma_families(6).items():
        sp = spectrum(eg, 6, 3)
        spec[nm] = sp
        print(f"{nm:>18} {sp['N']:6d} {sp['lam'][-2]:10.4f} {sp['slowest_tau']:8.2f} "
              f"{str([round(t,2) for t in sp['taus']]):>28} {sp['pr']:12.4f} (τ={sp['pr_tau']})")
    RES["spectrum"] = spec
    blk = spec["两块弱连(层级)"]
    ring = spec["ring6(正则)"]
    k6 = spec["K6(完全)"]
    check("**分层 Γ 的慢模最慢**（τ2 最大）—— Γ 的层级造出长寿命结构",
          blk["slowest_tau"] > ring["slowest_tau"] > k6["slowest_tau"],
          f"层级 {blk['slowest_tau']:.2f} > 环 {ring['slowest_tau']:.2f} > K6 {k6['slowest_tau']:.2f}")
    check("**分层 Γ 存在强局域慢模**（最小参与率 ≪ 正则环）—— 这就是「亚稳畴」",
          blk["pr"] < ring["pr"] / 2,
          f"层级 {blk['pr']:.4f}(τ={blk['pr_tau']}) vs 环 {ring['pr']:.4f}(τ={ring['pr_tau']}) "
          f"vs K6 {k6['pr']:.4f}(τ={k6['pr_tau']})")
    print(f"\n   建造/消亡速率随 Γ 变（{min(r['p_build'] for r in bk):.4f}–{max(r['p_build'] for r in bk):.4f}），"
          f"但**平衡与漂移是 Γ 无关的** ⟹ 结构差异住在**速率**里，不住在平衡里。")
    print("=" * 104)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_minimal_kernel.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_minimal_kernel.json")


if __name__ == "__main__":
    main()
