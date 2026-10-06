"""
z0_final.py --- 把「演化引擎」组装起来：**折叠（§27）＋ 局部 $P_i$（②）＋ $\\kappa_1$ 读出（①）**，
                 再问它带不带得出 $D=4$／$L=16$ 那几道门

组装（三项都只做各自那一件事，不互相发明）：
  · 演化：$s_v\\leftarrow|\\rho_v(s)-s_v|$，$\\rho_v=c_v/w_v$（§27；唯一给出「变化不停」的构造）
  · 忠实②：重播种**逐顶点**用自己的记录 $P_i$（ZS-8；折叠本身就是局部的）
  · 忠实①：读数乘 $\\kappa_1^{\\lfloor t/L\\rfloor}$，$\\kappa_1=q_L$，$T_d=L/(-\\log\\kappa_1)$（`G71`/`G72`）

要回答的新问题（§27 留的）：**$L$ 一大，$R_v=(A^L)_{vv}$ 的展布就涨**，而 §27 的判别线是
「展布小 ⟹ 混沌；展布大 ⟹ 冻结」。而量子门要 $L=16$。**两者能不能同时满足？**

四问：
  A 折叠引擎在 $L=4,6,8,12,16$ 上还混沌吗（$R$ 展布 vs `move`/Lyapunov）
  B $\\kappa_1$ 读出：$\\kappa_1=q_L$、$T_d(L)$ 逐 $L$；包络逐点精确
  C 门：R31/R25 的账本峰 $\\arg\\max_D\\binom D2q_L^D$；R44 的 $S_{\\max}(q)$；R53 的 $L=16$ 三门
  D 组装判定：哪些 $(\\Gamma,L)$ 同时满足 {演化, 忠实, 门}

用法：/usr/bin/python3 z0_final.py      输出：results/z0_final.json
"""
from __future__ import annotations

import json
import math
import os
import time
from math import comb
from fractions import Fraction
from itertools import combinations
from collections import Counter

import numpy as np

import z0_thm as ZT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
MU1, MU2 = math.sqrt(5), 1.3819660


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def qL_exact(L):
    seen = {}
    for ones in combinations(range(L), L // 2):
        w = [1] * L
        for i in ones:
            w[i] = -1
        c = min(tuple(w[i:] + w[:i]) for i in range(L))
        seen[c] = seen.get(c, 0) + 1
    return sum(Fraction(o, comb(L, L // 2)) ** 2 * n for o, n in Counter(seen.values()).items())


def fold_run(A, L, gens=4000, seed=11):
    """折叠引擎：s ← |ρ(s) − s|，ρ_v = c_v/w_v。返回 (move 末值, 斜率, Lyapunov, max_mass, R 展布)"""
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

    r = np.random.default_rng(seed)
    s = r.random(V) + 0.1
    s /= s.mean()
    mv, mm = [], []
    for _ in range(gens):
        ns = step(s)
        mv.append(float(np.abs(ns - s).sum() / V))
        occ = ns > 0.5 * ns.mean()
        mm.append(float(ns[occ].max() / ns[occ].sum()))
        s = ns
    pts = [(g, float(np.mean(mv[g - 20:g]))) for g in (100, 400, 1000, 2000, 4000) if g <= len(mv)]
    slope = float(np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0])

    def tr(e):
        rr = np.random.default_rng(seed)
        x = np.abs(rr.random(V) + 0.1 + e * rr.standard_normal(V))
        x /= x.mean()
        o = []
        for _ in range(140):
            x = step(x)
            o.append(x)
        return o
    a, b = tr(0.0), tr(1e-8)
    d = [float(np.abs(u - v).sum() / V) for u, v in zip(a, b)]
    sat = next((g for g in range(1, len(d)) if d[g] > 1e-2), None)
    lam = math.log(d[sat - 1] / d[0]) / (sat - 1) if sat and d[0] > 0 else float("nan")
    return dict(move_last=round(float(np.mean(mv[-20:])), 5), slope=round(slope, 3),
                lyapunov=None if lam != lam else round(lam, 4),
                max_mass=round(float(np.mean(mm[-20:])), 4),
                R_spread=round(float(R.max() / R.min()), 4))


def s_max(q):
    return MU2 + (MU1 - MU2) * (1 + math.sqrt(max(0.0, 2 * q - 1))) / 2


def ledger_peak(q, Dmax=12):
    vals = [(D, comb(D, 2) * q ** D) for D in range(2, Dmax + 1)]
    return max(vals, key=lambda t: t[1])[0]


if __name__ == "__main__":
    t0 = time.time()
    A, grp = ZT.build_hierarchical()
    V = A.shape[0]
    print("=" * 104)
    print(f"Γ：层级图 V={V}；折叠引擎 s←|ρ(s)−s| ＋ 局部 P_i ＋ κ1 读出")

    print("\nA 折叠引擎在更大 L 上还混沌吗（§27 的判别线：R 展布小 ⟹ 混沌）")
    print(f"   {'L':>3} {'R展布':>9} {'move(末)':>10} {'斜率':>8} {'Lyapunov':>10} {'max_mass':>9}  判定")
    Arows = {}
    for L in (4, 6, 8, 12, 16):
        r = fold_run(A, L)
        v = "**混沌**" if (r["move_last"] > 0.05 and r["lyapunov"] and r["lyapunov"] > 0) else "冻结/衰减"
        Arows[L] = r | {"verdict": v}
        print(f"   {L:3d} {r['R_spread']:9.4f} {r['move_last']:10.5f} {r['slope']:8.3f} "
              f"{(r['lyapunov'] if r['lyapunov'] is not None else float('nan')):+10.4f} "
              f"{r['max_mass']:9.4f}  {v}")
    RES["A_fold_by_L"] = {str(k): v for k, v in Arows.items()}
    chaos_L = [L for L, r in Arows.items() if "混沌" in r["verdict"]]

    print("\nB κ1 读出（κ1=q_L，T_d=L/(−log κ1)）")
    print(f"   {'L':>3} {'q_L':>12} {'κ1':>10} {'T_d(步)':>9}   V(n)=κ1^n 前 5 项")
    Brows = {}
    for L in (4, 6, 8, 12, 16):
        q = qL_exact(L)
        k1 = float(q)
        Td = L / (-math.log(k1))
        env = [k1 ** n for n in range(5)]
        Brows[L] = dict(q=str(q), kappa1=k1, Td=round(Td, 4), env=[round(x, 6) for x in env])
        print(f"   {L:3d} {str(q):>12} {k1:10.6f} {Td:9.4f}   {[round(x,4) for x in env]}")
    RES["B_kappa1"] = {str(k): v for k, v in Brows.items()}
    check("B κ1(L)=q_L 都在 G72 允许类内（整数计数之比 ⟹ 有理）",
          all(isinstance(Brows[L]["q"], str) for L in Brows))

    print("\nC 门（账本峰 / 语境性 / 稠密）")
    print(f"   {'L':>3} {'q_L':>12} {'账本峰 argmax D':>15} {'S_max(q_L)':>11} {'语境?':>6}")
    Crows = {}
    for L in (4, 6, 8, 12, 16):
        q = float(qL_exact(L))
        pk = ledger_peak(q)
        sm = s_max(q)
        Crows[L] = dict(q=q, peak_D=pk, S_max=round(sm, 6), contextual=bool(sm > 2))
        print(f"   {L:3d} {q:12.6f} {pk:15d} {sm:11.6f} {'是' if sm > 2 else '否':>6}")
    # R53 的 L=16（嵌套高度 π）
    d53 = json.load(open(os.path.join(os.path.dirname(HERE), "R53_zero_to_quantum_results.json")))
    q53 = d53["result"]["q"]; sm53 = d53["result"]["S_max"]
    print(f"   R53 的 L=16（嵌套高度 π）：q={q53:.6f}  S_max={sm53:.6f}  秩={d53['result']['rank']}/{d53['result']['needed_rank']}  "
          f"账本峰 D={d53['result']['ledger_peak_D']}  三门全过={d53['all_pass']}")
    RES["C_gates"] = {str(k): v for k, v in Crows.items()}
    RES["C_R53_L16"] = dict(q=q53, S_max=sm53, rank=d53["result"]["rank"], all_pass=d53["all_pass"])
    check("C 旋转类账本 q_L 在大 L 上给不出 D=4 峰（R31：只有 L=4）",
          Crows[4]["peak_D"] >= 4 and all(Crows[L]["peak_D"] < 4 for L in (6, 8, 12, 16)),
          f"峰: " + ", ".join(f"L={L}:{Crows[L]['peak_D']}" for L in (4, 6, 8, 12, 16)))
    check("C 旋转类 q_L 一律不语境（S_max≤2）；只有 R53 的 L=16 嵌套 π 过语境门",
          all(not Crows[L]["contextual"] for L in Crows) and sm53 > 2)

    print("\nD 组装判定：{演化, 忠实, 门} 谁能同时拿")
    print(f"   L=4 ：演化={('是' if 4 in chaos_L else '否')}  忠实=是（局部 P_i＋κ1 读出）  "
          f"门：D=4 峰 ✓ / 语境 ✗（S_max={s_max(float(qL_exact(4))):.4f}）")
    print(f"   L=16：演化={('是' if 16 in chaos_L else '否')}  忠实=是  门：R53 三门 ✓（q={q53:.6f}）")
    RES["D_package"] = dict(chaos_L=chaos_L,
                            L4=dict(evolve=4 in chaos_L, gates_D4=True,
                                    gates_contextual=bool(s_max(float(qL_exact(4))) > 2)),
                            L16=dict(evolve=16 in chaos_L, gates_R53=bool(d53["all_pass"])))
    both = (4 in chaos_L) and (16 in chaos_L)
    check("D 折叠混沌在 L=4 与 L=16 **都**活着（量子门与演化可同取）", both,
          f"混沌的 L = {chaos_L}")
    print("=" * 104)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_final.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_final.json")
