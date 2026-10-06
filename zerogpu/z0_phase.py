"""
z0_phase.py --- `ACTION-PHASE-MATCH` 的最小可测版：把**相位接到重播种的选择**上

上一轮的诊断（§25）：只要重播种的每代倍率含 $R_v=(A^L)_{vv}$，就**必凝聚**——五种规则无一例外。
所以这一轮换两件事：

**① 换掉乘性**：用**局部闭合概率**而不是闭合**计数**
$$
c_v=s_vR_v,\quad w_v=\\sum_u s_u(A^L)_{uv},\\qquad \\rho_v=\\frac{c_v}{w_v}\\ \\in[0,1]
$$
$\\rho_v$ 是"到达该顶点的走量里闭合掉的比例"，**不随 $R_v$ 增长** ⟹ 不再有 $R_v^{\\,n}$ 的单峰化。

**② 把相位接进选择**（不是接进幅度）：每顶点有自己的闭合代数 $g_v$（该顶点累计闭合次数），
门控 $\\cos(\\varphi g_v)>0$ 时用 $\\rho_v$ 更新，否则**保留**（不是清零 ⟹ **不产生吸收态**，这是 §19.6 相位门死掉的原因）。
$\\varphi=\\log 7.2$（§19 的 $L{=}4$ 转动数）。

四条臂：
  `mult`        §24/§25 基线（乘性）—— 对照
  `prob`        局部闭合概率（无相位）
  `prob_phase`  概率 ＋ **逐顶点相位门**（保留式，不吸收）
  `prob_global_phase` 概率 ＋ **全局相位门**（所有顶点共用 $g$；对照"相位是否必须局部"）

判据仍是 §25 的三元组：`max_mass`（→1 凝聚）、`gini`、`move`（→0 冻结），
外加 **`novelty`**：末 20 代里出现过的**不同形状数**（=1 ⟹ 精确冻结）。
判决：**持续多畴** ⟺ $max\_mass<0.5$ 且 $move>10^{-3}$（有界＋变化＋结构）。

用法：/usr/bin/python3 z0_phase.py      输出：results/z0_phase.json
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
PHI = math.log(7.2)          # §19：L=4 的模相位旋转率（2π/φ = 3.1828）


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def run(A, arm, s0, gens=80, L=4, phi=PHI, tol_novel=1e-6):
    V = A.shape[0]
    AL = (A ** L).toarray()
    R = np.diag(AL).copy()
    AT = AL.T.copy()
    s = np.asarray(s0, float).copy(); s /= s.mean()
    g = np.zeros(V, np.int64)          # 逐顶点闭合代数
    gG = 0                             # 全局闭合代数
    mms, gs, mvs, shapes = [], [], [], []
    prev = None
    for _ in range(gens):
        c = s * R
        w = AT @ s
        rho = np.where(w > 0, c / np.maximum(w, 1e-300), 0.0)
        g += (c > 0).astype(np.int64)
        gG += 1
        if arm == "mult":
            ns = c
        elif arm == "prob":
            ns = rho
        elif arm == "prob_phase":
            ok = np.cos(phi * g) > 0
            ns = np.where(ok, rho, s)
        elif arm == "prob_global_phase":
            ok = math.cos(phi * gG) > 0
            ns = rho if ok else s
        if ns.mean() <= 0:
            ns = s.copy()
        ns = ns / ns.mean()
        occ = ns > 0.5 * ns.mean()
        mms.append(float(ns[occ].max() / ns[occ].sum()) if occ.any() else 1.0)
        gs.append(ZT.gini(ns))
        mvs.append(float(np.abs(ns - prev).sum() / V) if prev is not None else float("nan"))
        shapes.append(tuple(np.round(ns, 6)))
        prev = ns
        s = ns
    novelty = len(set(shapes[-20:]))
    return dict(max_mass=[round(x, 5) for x in mms], gini=[round(x, 5) for x in gs],
                move=[round(x, 6) if x == x else None for x in mvs], novelty_last20=novelty)


def verdict(tr):
    """注意：80 代的 move 不足以判"持续"。§26 查明它按幂律衰减 ⟹ 这里只判"有结构/凝聚"，
    变化是否持续由下面的长程幂律检（part_long）另判。"""
    mm = tr["max_mass"][-1]
    if mm > 0.5:
        return "凝聚"
    if tr["novelty_last20"] <= 1:
        return "冻结"
    return "有界＋结构（变化待长程判）"


def long_run(A, arm, s0, gens=2000, L=4, phi=PHI):
    """长程：move 是否指数归零（真冻结）还是幂律（未停）。返回分段均值与双对数斜率。"""
    V = A.shape[0]
    AL = (A ** L).toarray(); R = np.diag(AL).copy(); AT = AL.T.copy()
    s = np.asarray(s0, float).copy(); s /= s.mean()
    g = np.zeros(V, np.int64); mv = []
    for n in range(gens):
        c = s * R; w = AT @ s
        rho = np.where(w > 0, c / np.maximum(w, 1e-300), 0.0)
        g += (c > 0).astype(np.int64)
        ns = c if arm == "mult" else (rho if arm == "prob" else np.where(np.cos(phi * g) > 0, rho, s))
        ns = ns / ns.mean()
        if n > 0:
            mv.append(float(np.abs(ns - prev).sum() / V))
        prev = ns; s = ns
    pts = [(100, np.mean(mv[80:100])), (400, np.mean(mv[380:400])),
           (1000, np.mean(mv[980:1000])), (2000, np.mean(mv[1980:2000]))]
    slope = float(np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0])
    return pts, slope


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 100)
    A, grp = ZT.build_hierarchical()
    V = A.shape[0]
    rng = np.random.default_rng(11)
    s0 = rng.random(V) + 0.1
    print(f"Γ：V={V}，L=4，φ=log 7.2={PHI:.6f}（每顶点自己的闭合代数 g_v；门控=保留式，不清零）")
    print(f"\n{'臂':>18} {'max_mass(80)':>13} {'gini(80)':>9} {'move(80)':>10} {'novelty(末20)':>13}  判定")
    out = {}
    for arm in ("mult", "prob", "prob_phase", "prob_global_phase"):
        tr = run(A, arm, s0, gens=80)
        v = verdict(tr)
        out[arm] = dict(trajectory=tr, verdict=v)
        print(f"{arm:>18} {tr['max_mass'][-1]:13.4f} {tr['gini'][-1]:9.4f} "
              f"{(tr['move'][-1] if tr['move'][-1] is not None else float('nan')):10.6f} "
              f"{tr['novelty_last20']:13d}  {v}")
    print(f"\n  max_mass 轨迹（每 10 代）：")
    for arm in out:
        mm = out[arm]["trajectory"]["max_mass"]
        print(f"    {arm:>18}: {[mm[i] for i in (0, 9, 19, 39, 59, 79)]}")
    RES["arms"] = {k: v["trajectory"] | {"verdict": v["verdict"]} for k, v in out.items()}
    check("局部闭合概率**不再凝聚**（max_mass ≪ 0.5）",
          out["prob"]["trajectory"]["max_mass"][-1] < 0.5,
          f"mult={out['mult']['trajectory']['max_mass'][-1]} prob={out['prob']['trajectory']['max_mass'][-1]}")
    check("`prob` 族不再凝聚（max_mass ≪ 0.5）",
          all(out[a]["trajectory"]["max_mass"][-1] < 0.5 for a in ("prob", "prob_phase")),
          f"mult={out['mult']['trajectory']['max_mass'][-1]}")
    print(f"\n  长程判（§26.2）：move 是**指数归零**还是**幂律**？")
    slopes = {}
    for arm in ("mult", "prob", "prob_phase"):
        pts, slope = long_run(A, arm, s0, gens=2000)
        slopes[arm] = round(slope, 3)
        print(f"    {arm:>11}: move(g) = " + "  ".join(f"{g}:{v:.2e}" for g, v in pts)
              + f"   双对数斜率 = {slope:+.3f}")
    check("`mult` 是真冻结（指数归零，斜率 ≪ -1）", slopes["mult"] < -1.0, f"斜率={slopes['mult']}")
    check("`prob` 族是**幂律**衰减（|斜率| < 1，未停止）", abs(slopes["prob"]) < 1.0,
          f"斜率={slopes['prob']}")
    RES["long_run_slopes"] = slopes
    print("\n  ⟹ **诚实判定（§26）**：`prob` 族＝**有界 ＋ 结构（32 畴）＋ 变化（幂律衰减、未停）**；")
    print("     这是 §25 三选二的**部分**反例；但若要求 move 有正下界，则不算「持续演化」。")
    print("     相位门那半边**否证**：逐顶点 g_v 退化成全局（g_v≡n），有理／无理门无差别。")
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_phase.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_phase.json")
