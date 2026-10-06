"""
z0_chaos.py --- **第一个"每次毁灭重长都不一样"的构造**（Z0_CORE §15 的那个目标）

链条：
  ① §26：把重播种从**闭合计数**换成**局部闭合概率** $\\rho_v=c_v/w_v$（$w_v$=到达走量）⟹ **治好凝聚**（32 畴）
  ② 但 $\\rho$ 的线性化谱半径 $=1$（临界）⟹ 变化按 $g^{-1/2}$ 衰减（§26.2）
  ③ 本轮：在 $\\rho$ 上做 §15 点名的**非单调折叠**（stretch-and-fold）
$$
\\boxed{\\ s_v\\ \\longleftarrow\\ \\bigl|\\,\\rho_v(s)-s_v\\,\\bigr|\\ }
\\qquad
\\rho_v=\\frac{s_vR_v}{\\sum_u s_u(A^L)_{uv}},\\quad R_v=(A^L)_{vv}
$$
  两项都是导出的：$\\rho$ 是闭合概率；$|\\rho-s|$ 是"新闭合概率与现种子的**差**"。
  **无参数、确定性、无概率（Z0③）、无吸收态**（差值的绝对值只在 $\\rho=s$ 处为 0，而那是零测集）。

测三件事（§15 的目标就是这个）：
  A **变化不停**：$\\texttt{move}(g)$ 的长程斜率 ≈ 0（对照 §26 的 $-0.5$、§25 的 $-36.7$）
  B **非周期**：4000 代内无重复状态
  C **初值敏感（混沌）**：两个相差 $10^{-8}$ 的初值，距离指数放大 ⟹ Lyapunov $>0$

用法：/usr/bin/python3 z0_chaos.py      输出：results/z0_chaos.json
"""
from __future__ import annotations

import json
import math
import os
import time

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


def main():
    t0 = time.time()
    A, grp = ZT.build_hierarchical()
    step, R, V = make_step(A, 4)
    print("=" * 100)
    print(f"Γ：层级图 V={V}，L=4；映射 s ← |ρ(s) − s|（无参数、确定性、无概率）")

    # A 长程：变化是否停
    rng = np.random.default_rng(11)
    s = rng.random(V) + 0.1
    s /= s.mean()
    mv, mm = [], []
    GENS = 20000
    for _ in range(GENS):
        ns = step(s)
        mv.append(float(np.abs(ns - s).sum() / V))
        occ = ns > 0.5 * ns.mean()
        mm.append(float(ns[occ].max() / ns[occ].sum()))
        s = ns
    pts = [(g, float(np.mean(mv[g - 20:g]))) for g in (100, 400, 1000, 2000, 5000, 10000, 20000)]
    slope = float(np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0])
    print(f"\nA 长程 {GENS} 代：")
    for g, v in pts:
        print(f"    g={g:>6}: move={v:.4f}   max_mass={np.mean(mm[g-20:g]):.4f}")
    print(f"    双对数斜率 = {slope:+.4f}（对照 §26 的 −0.501、§25 mult 的 −36.7）")
    check("A 变化不停（|斜率| < 0.1）", abs(slope) < 0.1, f"斜率={slope:+.4f}")
    RES["A_sustained"] = dict(points=[[int(g), round(v, 6)] for g, v in pts],
                              slope=round(slope, 4),
                              max_mass_track=[[int(g), round(float(np.mean(mm[g-20:g])), 4)]
                                              for g in (100, 1000, 10000, 20000)])

    # B 非周期
    s = rng.random(V) + 0.1
    s /= s.mean()
    T = 4000
    seen = {}
    first_rep = None
    for k in range(T):
        key = tuple(np.round(s, 7))
        if key in seen:
            first_rep = (seen[key] + 1, k - seen[key])
            break
        seen[key] = k
        s = step(s)
    print(f"\nB 非周期：4000 代内{'**无**重复状态' if first_rep is None else f'首次重复={first_rep}'}")
    check("B 4000 代内无重复状态（非周期）", first_rep is None, f"{first_rep}")
    RES["B_aperiodic"] = dict(first_repeat=first_rep, window=T)

    # C 初值敏感
    def traj(gens, eps=0.0, seed=11):
        r = np.random.default_rng(seed)
        x = r.random(V) + 0.1
        x /= x.mean()
        if eps:
            x = np.abs(x + eps * r.standard_normal(V))
            x /= x.mean()
        out = []
        for _ in range(gens):
            x = step(x)
            out.append(x)
        return out
    a = traj(140)
    b = traj(140, eps=1e-8)
    d = [float(np.abs(x - y).sum() / V) for x, y in zip(a, b)]
    sat = next((g for g in range(1, len(d)) if d[g] > 1e-2), None)
    lam = math.log(d[sat - 1] / d[0]) / (sat - 1) if sat and d[0] > 0 else float("nan")
    print(f"\nC 初值敏感：d(1)={d[0]:.2e} → d(10)={d[9]:.2e} → d(40)={d[39]:.2e} → d(80)={d[79]:.2e}"
          f"  饱和于 g={sat}")
    print(f"    估计 Lyapunov 指数 ≈ {lam:+.4f} /代")
    check("C 初值差异被放大（Lyapunov > 0）", lam > 0, f"λ≈{lam:+.4f}")
    RES["C_lyapunov"] = dict(d=[float(x) for x in d[:20]], sat_gen=sat, lyapunov=round(lam, 4))

    # D Γ 判别量：R_v=(A^L)_vv 的展布 vs 混沌
    print("\nD Γ 判别量（§27.2）：R 展布小 ⟹ 混沌；展布大 ⟹ 精确冻结")
    base = ZT.torus_adj(2, 16)
    fam = [("环面16²(正则)", base),
           ("环面+枢纽10%", ZT.add_hubs(base, frac=0.10, extra=4, seed=3)),
           ("环面+枢纽5%x12", ZT.add_hubs(base, frac=0.05, extra=12, seed=3))]
    for sd in (0, 1, 2):
        h, _ = ZT.build_hierarchical(seed=sd)
        fam.append((f"层级Γ(seed={sd})", h))
    print(f"   {'Γ':>16} {'R展布':>8} {'move(末)':>9} {'max_mass':>9} {'Lyapunov':>10}  判定")
    famrows = []
    for nm, G in fam:
        st, RG, VG = make_step(G, 4)
        r = np.random.default_rng(11)
        x = r.random(VG) + 0.1
        x /= x.mean()
        mvs, mms = [], []
        for _ in range(4000):
            nx = st(x)
            mvs.append(float(np.abs(nx - x).sum() / VG))
            occ = nx > 0.5 * nx.mean()
            mms.append(float(nx[occ].max() / nx[occ].sum()))
            x = nx

        def tr(e):
            rr = np.random.default_rng(11)
            y = np.abs(rr.random(VG) + 0.1 + e * rr.standard_normal(VG))
            y /= y.mean()
            o = []
            for _ in range(140):
                y = st(y)
                o.append(y)
            return o
        a2, b2 = tr(0.0), tr(1e-8)
        dd = [float(np.abs(u - v).sum() / VG) for u, v in zip(a2, b2)]
        sat2 = next((g for g in range(1, len(dd)) if dd[g] > 1e-2), None)
        lam2 = math.log(dd[sat2 - 1] / dd[0]) / (sat2 - 1) if sat2 and dd[0] > 0 else float("nan")
        spread = float(RG.max() / RG.min())
        v = "**混沌**" if (np.mean(mvs[-20:]) > 0.05 and lam2 == lam2 and lam2 > 0) else ("**精确冻结**" if np.mean(mvs[-20:]) < 1e-6 else "其它")
        famrows.append(dict(name=nm, spread=round(spread, 3), move=round(float(np.mean(mvs[-20:])), 4),
                            max_mass=round(float(np.mean(mms[-20:])), 4),
                            lyapunov=None if lam2 != lam2 else round(lam2, 4), verdict=v))
        print(f"   {nm:>16} {spread:8.3f} {np.mean(mvs[-20:]):9.4f} {np.mean(mms[-20:]):9.4f} "
              f"{lam2:+10.4f}  {v}")
    RES["D_gamma_family"] = famrows
    ch = [r_ for r_ in famrows if "混沌" in r_["verdict"]]
    fz = [r_ for r_ in famrows if "冻结" in r_["verdict"]]
    check("D 判别量：混沌组的 R 展布 < 冻结组", bool(ch) and bool(fz)
          and max(r_["spread"] for r_ in ch) < min(r_["spread"] for r_ in fz),
          f"混沌展布≤{max((r_['spread'] for r_ in ch), default=None)} < 冻结展布≥{min((r_['spread'] for r_ in fz), default=None)}")

    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_chaos.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_chaos.json")


if __name__ == "__main__":
    main()
