"""
z0_four.py --- 按①→④把四处**不忠实**逐个接上，并做**增量归因**：哪一步真的改变了结局？

四步（顺序即用户给的顺序）：
  ① 终端压制 $\\kappa_1^{\\lfloor t/L\\rfloor}$（`G71`；$\\kappa_1=q_4=5/9$，$T_d=6.805$ 步）
     —— 两种读法都测：(a) **只作用于读出**（`G71` 原文：$V(t)$ 是干涉可见度）；
        (b) 作用到重播种上（把相干衰减喂进动力学）⟹ 预期**衰亡**。
  ② 重播种改**局部 $P_i$**（每顶点自己的记录）而不是全局读回（`ZS-8`：全局读回把 6 种结果塌成 1 种）
  ③ 让 **L1 承重**：重播种读**保留的记录层**（`D222` 最近两层），而不是瞬时闭合量
  ④ 饱和换**导出口径**：`G32` 的正交极小投影 $2(T+1)$，而不是塞进去的 `cap`

配置（累加）：S0 基线 → S1 ＋①(b) → S2 ＋② → S3 ＋③ → S4 ＋④
判据（§25/§26/§27 同一套）：`max_mass`（→1 凝聚）、`Gini`、畴数、`move` 斜率（≈0 = 变化不停，负 = 衰减/冻结）

用法：/usr/bin/python3 z0_four.py      输出：results/z0_four.json
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
KAPPA1 = 5 / 9          # G71/G72：kappa1 = q_4


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def run(A, sw, gens=1200, L=4, cap_invented=2, T=4):
    V = A.shape[0]
    AL = (A ** L).toarray()
    R = np.diag(AL).copy()
    rng = np.random.default_rng(11)
    s = rng.random(V) + 0.1
    s /= s.mean()
    layers = []
    mv = []
    prev = None
    for n in range(gens):
        c = s * R                                   # 该代闭合量（局部）
        layers.append(c.copy())
        if len(layers) > 2:
            layers.pop(0)
        # ③ L1 承重：从保留层重播种（D222 最近两层），否则用瞬时闭合量
        src = sum(layers) if sw.get("L1") else c
        # ④ 饱和：导出口径 2(T+1)（G32） vs 塞进去的 cap
        if sw.get("cap_invented"):
            src = np.minimum(src, cap_invented)
        if sw.get("cap_derived"):
            src = np.minimum(src, 2 * (T + 1))
        # ② 局部 vs 全局
        ns = np.full(V, float(src.mean())) if sw.get("global_read") else src.copy()
        # ①(b) 相干衰减喂进动力学
        if sw.get("kappa_dyn"):
            ns = ns * (KAPPA1 ** n)
        if ns.mean() <= 1e-300:
            ns = s.copy()
        ns = ns / ns.mean()
        if n > 0:
            mv.append(float(np.abs(ns - prev).sum() / V))
        prev = ns
        s = ns
    occ = s > 0.5 * s.mean()
    pts = [(g, float(np.mean(mv[g - 20:g]))) for g in (100, 400, 800, 1200) if g <= len(mv)]
    slope = float(np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0]) if len(pts) > 1 else float("nan")
    return dict(max_mass=round(float(s[occ].max() / s[occ].sum()), 5) if occ.any() else 1.0,
                gini=round(ZT.gini(s), 4), domains=int(ZT.domains_lab(occ, A)[0]),
                move_last=round(float(np.mean(mv[-20:])), 6), slope=round(slope, 3),
                pts=[[int(g), round(v, 7)] for g, v in pts])


def verdict(r):
    if r["max_mass"] > 0.5:
        return "凝聚"
    if r["move_last"] < 1e-6:
        return "冻结（精确）"
    if r["slope"] < -0.1:
        return "衰减（未停但趋零）"
    return "**变化不停**"


if __name__ == "__main__":
    t0 = time.time()
    A, grp = ZT.build_hierarchical()
    V = A.shape[0]
    print("=" * 104)
    print(f"Γ：层级图 V={V}，L=4；κ1=q_4={KAPPA1:.4f}，T_d=L/(-log κ1)={4/(-math.log(KAPPA1)):.3f} 步")
    print("\n① 的两种读法（先说清）：")
    print("   (a) 只作用于**读出**（G71 原文：V(t) 是干涉可见度）⟹ 不改动力学，只让读数按 κ1^{⌊t/L⌋} 衰减")
    print("   (b) 把相干衰减**喂进重播种** ⟹ 那就是给动力学加了一个几何阻尼（下表的 S1）")
    configs = [
        # 轨道 A：现引擎的**全局读回**语义（ZS-8 的 global_P）
        ("A0 全局读回（现引擎那一类）",      dict(global_read=True)),
        ("A1 ＋①(b) κ1 阻尼",              dict(global_read=True, kappa_dyn=True)),
        ("A2 ＋③ L1 承重",                 dict(global_read=True, kappa_dyn=True, L1=True)),
        ("A3 ＋④ 导出饱和",                dict(global_read=True, kappa_dyn=True, L1=True, cap_derived=True)),
        # 轨道 B：②局部 P_i 之后，再加 ①③④
        ("B0 ②局部 P_i（无幅）",            dict()),
        ("B1 ＋①(b) κ1 阻尼",              dict(kappa_dyn=True)),
        ("B2 ＋③ L1 承重",                 dict(kappa_dyn=True, L1=True)),
        ("B3 ＋④ 导出饱和 2(T+1)",          dict(kappa_dyn=True, L1=True, cap_derived=True)),
        ("B3′ ④ 换成塞进去的 cap2（对照）",  dict(kappa_dyn=True, L1=True, cap_invented=2)),
    ]
    print(f"\n{'配置':>30} {'max_mass':>9} {'Gini':>7} {'畴数':>5} {'move(末)':>10} {'斜率':>8}  判定")
    out = {}
    for nm, sw in configs:
        r = run(A, sw)
        out[nm] = r | {"verdict": verdict(r)}
        print(f"{nm:>30} {r['max_mass']:9.4f} {r['gini']:7.4f} {r['domains']:5d} "
              f"{r['move_last']:10.6f} {r['slope']:8.3f}  {verdict(r)}")
    RES["configs"] = out
    # 增量归因
    print("\n增量归因（每一步相对上一步改变了什么）：")
    keys = [c[0] for c in configs]
    for i in range(1, len(keys)):
        if keys[i][0] != keys[i - 1][0]:
            print(f"   —— 换轨道：{keys[i-1][:14]} → {keys[i][:14]} ——")
        a, b = out[keys[i - 1]], out[keys[i]]
        d = {k: round(b[k] - a[k], 5) for k in ("max_mass", "gini", "move_last", "slope") if isinstance(a[k], float)}
        print(f"   {keys[i-1][:14]} → {keys[i][:14]}: Δ={d}  判定 {a['verdict']} → {b['verdict']}")
    print("\n①(a) 纯读出读法的核验：")
    L = 4
    A1c, A2c = 2.52853, 0.35119
    phi = math.log(7.2)
    steps = np.arange(0, 24)
    g = steps / 2.0
    env = KAPPA1 ** np.floor(steps / L)
    P = (A1c ** 2 + A2c ** 2 + 2 * A1c * A2c * np.cos(phi * g)) * env
    ef = int(np.argmax(P < P[0] / math.e))
    print(f"   P(g) 前 10 步 = {np.round(P[:10], 3).tolist()}   e 折 ≈ {ef} 步（理论 T_d={4/(-math.log(KAPPA1)):.2f}）")
    Vn = KAPPA1 ** np.floor(np.arange(0, 24) / L)
    check("①(a) 包络 V(n)=κ1^{⌊t/L⌋} 逐点精确、T_d=L/(-log κ1)=6.805（公式级断言）",
          abs(float(Vn[4]) - KAPPA1) < 1e-15 and abs(4 / (-math.log(KAPPA1)) - 6.8052) < 1e-3,
          f"V(4)={Vn[4]:.6f}=κ1；T_d={4/(-math.log(KAPPA1)):.4f}；P 自身表观 e 折={ef} 步（含 cos 调制，另记）")
    best = [k for k, v in out.items() if v["verdict"] == "**变化不停**"]
    check("四步里**没有任何一步**做到「变化不停」", len(best) == 0,
          f"做到的是：{best}" if best else "全部落到 凝聚/冻结/衰减")
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL, evolving=best)
    print("=" * 104)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    json.dump(RES, open(os.path.join(OUT, "z0_four.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_four.json")
