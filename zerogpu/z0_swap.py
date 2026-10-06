"""
z0_swap.py --- 治「衰减」：让相位在**两种归一化之间切换**，并让局部时钟 $g_v$ 真正异质

诊断（§26）：`prob`（$\\rho_v=c_v/w_v$）已经治好**凝聚**（32 畴），但变化按 $g^{-1/2}$ 衰减
—— 那是**线性化谱半径 $=1$（临界／边缘）**的特征，不是普通收敛。

本文件试两条无参数设计：

**A `swap`：相位在两种归一化之间切换**（都不产生吸收态、都 $\\in[0,1]$）
$$
\\text{相长}:\\ s_v\\leftarrow\\rho_v=\\frac{c_v}{\\sum_u s_u(A^L)_{uv}}\\ (\\text{按\\textbf{到达}量归一的闭合概率});
\\qquad
\\text{相消}:\\ s_v\\leftarrow\\sigma_v=\\frac{c_v}{\\sum_u s_v(A^L)_{vu}}\\ (\\text{按\\textbf{出发}量归一})
$$
两个映射的不动点不同 ⟹ 被无理旋转切换时，轨道"追两个不动点"而**不落定**（若成立，衰减就停）。

**B `swap_local`：让 $g_v$ 异质**。§26.3 查明 $g_v\\equiv n$（每顶点每代都闭合）⟹ 门退化成全局。
这里改成**只在"该顶点自己的事件"上 $+1$**：$c_v$ 超过它自己的滑动均值时才 $+1$ ⟹ 安静顶点的钟慢。

对照：`prob`（无门）、`prob_phase`（保留式门）。

判据（§26 的口径）：`max_mass`、`Gini`、畴数、以及**长程双对数斜率**（斜率近 0 = 衰减停住）。

用法：/usr/bin/python3 z0_swap.py      输出：results/z0_swap.json
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
PHI = math.log(7.2)


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def run(A, arm, s0, gens=2000, L=4, phi=PHI, sample=(100, 400, 1000, 2000)):
    V = A.shape[0]
    AL = (A ** L).toarray()
    R = np.diag(AL).copy()
    AT = AL.T.copy()
    s = np.asarray(s0, float).copy(); s /= s.mean()
    g = np.zeros(V, np.int64)
    gG = 0
    run_mean = s.copy()
    mv = []
    prev = None
    for n in range(gens):
        c = s * R
        w_in = AT @ s                      # 到达量
        w_out = AL @ s                     # 出发量
        rho = np.where(w_in > 0, c / np.maximum(w_in, 1e-300), 0.0)
        sig = np.where(w_out > 0, c / np.maximum(w_out, 1e-300), 0.0)
        # 事件：c_v 超过自己的滑动均值 ⟹ 才是"该顶点的事件"
        ev = c > run_mean
        run_mean = 0.99 * run_mean + 0.01 * c
        if arm == "prob":
            ns = rho
        elif arm == "prob_phase":
            g += (c > 0).astype(np.int64)
            ns = np.where(np.cos(phi * g) > 0, rho, s)
        elif arm == "swap":
            g += (c > 0).astype(np.int64)
            ns = np.where(np.cos(phi * g) > 0, rho, sig)
        elif arm == "swap_local":
            g += ev.astype(np.int64)       # ★ 只在"自己的事件"上 +1 ⟹ 异质
            ns = np.where(np.cos(phi * g) > 0, rho, sig)
        elif arm == "swap_global":
            gG += 1
            ns = rho if math.cos(phi * gG) > 0 else sig
        if ns.mean() <= 0:
            ns = s.copy()
        ns = ns / ns.mean()
        if n > 0:
            mv.append(float(np.abs(ns - prev).sum() / V))
        prev = ns
        s = ns
    pts = []
    for G in sample:
        if G <= len(mv):
            pts.append((G, float(np.mean(mv[max(0, G - 20):G]))))
    slope = float(np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0]) if len(pts) > 1 else float("nan")
    occ = s > 0.5 * s.mean()
    return dict(pts=[(g_, round(v, 9)) for g_, v in pts], slope=round(slope, 3),
                max_mass=round(float(s[occ].max() / s[occ].sum()), 5) if occ.any() else 1.0,
                gini=round(ZT.gini(s), 4), domains=int(ZT.domains_lab(occ, A)[0]),
                g_distinct=int(len(set(g.tolist()))))


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 100)
    A, grp = ZT.build_hierarchical()
    V = A.shape[0]
    rng = np.random.default_rng(11)
    s0 = rng.random(V) + 0.1
    print(f"Γ：V={V}，L=4，φ=log 7.2；判据=长程双对数斜率（越接近 0 = 衰减越停）")
    print(f"\n{'臂':>13} {'斜率':>8} {'max_mass':>9} {'Gini':>7} {'畴数':>5} {'g_v取值数':>9}   分段 move")
    out = {}
    for arm in ("prob", "prob_phase", "swap", "swap_local", "swap_global"):
        r = run(A, arm, s0)
        out[arm] = r
        seg = "  ".join(f"{g_}:{v:.1e}" for g_, v in r["pts"])
        print(f"{arm:>13} {r['slope']:8.3f} {r['max_mass']:9.4f} {r['gini']:7.4f} "
              f"{r['domains']:5d} {r['g_distinct']:9d}   {seg}")
    RES["arms"] = out
    best = min(out.items(), key=lambda kv: abs(kv[1]["slope"]))
    print(f"\n  衰减最慢的臂：**{best[0]}**（斜率 {best[1]['slope']}）")
    check("`swap`（相位在两种归一化间切换）比 `prob` 衰减更慢",
          abs(out["swap"]["slope"]) < abs(out["prob"]["slope"]),
          f"prob={out['prob']['slope']} swap={out['swap']['slope']}")
    check("`swap_local` 的 g_v 真正异质（取值数 > 1）",
          out["swap_local"]["g_distinct"] > 1,
          f"g_v 取值数={out['swap_local']['g_distinct']}（对照 prob_phase={out['prob_phase']['g_distinct']}）")
    check("没有任何臂做到「斜率≈0」（衰减停住）", all(abs(v["slope"]) > 0.05 for v in out.values()),
          "全部仍在衰减" if all(abs(v["slope"]) > 0.05 for v in out.values()) else "有臂停住了")
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL, slowest=best[0])
    json.dump(RES, open(os.path.join(OUT, "z0_swap.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_swap.json")
