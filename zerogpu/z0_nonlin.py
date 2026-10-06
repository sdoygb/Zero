"""
z0_nonlin.py --- 重播种非线性的判决性对照：**有没有任何无参数规则给出"持续多畴"？**

动机（§16 与 §24.2）：
  引擎里所有现行重播种都对 $R_v=(A^L)_{vv}$ **乘性** ⟹ $s_v\\propto R_v^{\\,n}$ ⟹ 必然凝聚（§24.2 已量到速率）。
  §16 说：唯一天然的非线性是**类计数**（合并更多世界线不产生更多不同的类，故对重数是次线性/饱和的）；
  §15 另列了两种折叠候选：**秩截断**（不连续）与**单调饱和**。三者都没在 Γ 引擎上测过。

本文件四条臂（全部确定性、无参数、无概率）：
  `mult`       基线：$s\\leftarrow s\\cdot R$（§24.2 的 `diag`）
  `classcount` 类计数：$s_v\\leftarrow$ 该顶点这一代**实现出的不同闭合类个数**（§16 的唯一天然非线性）
  `rank`       秩保留：$s_v\\leftarrow \\mathrm{rank}(c_v)$（§15 候选 2，不连续）
  `saturate`   单调饱和：$s_v\\leftarrow c_v/(1+c_v)$（§15 候选 1）

判据（三个数就够）：
  `max_mass(g)` 最大畴质量（→1 = 凝聚）
  `gini(g)`     分布不均度
  `mv(g)`       相邻代形状的 L1 变化率（→0 = **冻结**）
  判定：`凝聚`（max_mass>0.5）/ `冻结`（mv<1e-3 且 max_mass<0.5）/ **`持续多畴`**（两者都不满足且 max_mass<0.5）

用法：/usr/bin/python3 z0_nonlin.py      输出：results/z0_nonlin.json
"""
from __future__ import annotations

import json
import math
import os
import time
from collections import defaultdict

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


def class_count_by_vertex(A, L=4):
    """每个顶点能实现出的**不同闭合类**个数（旋转类）。
    类 = 本原 Dyck 词的旋转类；顶点可实现的类 = 存在该长度的闭合走。
    对 min degree>=2 的图，短类人人可及 ⟹ 该数**与顶点无关**（这正是 §16 的"饱和"）。"""
    def prims(L):
        out = []
        def rec(w):
            if w and sum(w) == 0:
                out.append(tuple(w)); return
            if len(w) == L:
                return
            for d in (+1, -1):
                rec(w + [d])
        for d in (+1, -1):
            rec([d])
        return out
    words = prims(L)
    classes = {min(tuple(w[i:] + w[:i]) for i in range(len(w))) for w in words}
    V = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel()
    AL = (A ** L).toarray()
    cnt = np.zeros(V)
    for v in range(V):
        # 能走回自身的长度 ℓ 上，类的可实现性：长度 ℓ 的闭合走存在 ⟺ (A^ℓ)_vv>0
        realizable = 0
        for cw in classes:
            ell = len(cw)
            if ell <= L and (A ** ell).toarray()[v, v] > 0:
                realizable += 1
        cnt[v] = realizable
    return cnt, len(classes), words


def run_arm(A, R, arm, s0, gens=60, L=4, cls=None):
    V = A.shape[0]
    s = np.asarray(s0, float).copy()
    s = s / s.mean()
    mms, gs, mvs = [], [], []
    prev = None
    for g in range(gens):
        c = s * R                                   # 该代闭合量
        if arm == "mult":
            ns = c
        elif arm == "classcount":
            ns = cls.copy()                         # 类个数（饱和量）
        elif arm == "rank":
            order = np.argsort(np.argsort(c)).astype(float) + 1.0
            ns = order
        elif arm == "saturate":
            ns = c / (1.0 + c)
        elif arm == "dev":                          # §15 说的"非单调折叠"（V 形，无参数）
            ns = np.abs(c - np.median(c))
        ns = ns / ns.mean()
        occ = ns > 0.5 * ns.mean()
        mm = float(ns[occ].max() / ns[occ].sum()) if occ.any() else 1.0
        mms.append(mm); gs.append(ZT.gini(ns))
        mvs.append(float(np.abs(ns - prev).sum() / V) if prev is not None else float("nan"))
        prev = ns
        s = ns
    return dict(max_mass=[round(x, 5) for x in mms], gini=[round(x, 5) for x in gs],
                move=[round(x, 5) if x == x else None for x in mvs])


def verdict(tr):
    mm60, mv60 = tr["max_mass"][-1], tr["move"][-1]
    if mm60 > 0.5:
        return "凝聚"
    if mv60 is not None and mv60 < 1e-3:
        return "冻结"
    return "持续多畴"


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 100)
    A, grp = ZT.build_hierarchical()
    V = A.shape[0]
    L = 4
    R = np.diag((A ** L).toarray()).copy()
    cls, nclass, words = class_count_by_vertex(A, L)
    print(f"Γ：V={V}  闭合词 {len(words)} 个 → **旋转类 {nclass} 个**（L={L}）")
    print(f"类计数逐顶点：min={cls.min():.0f} max={cls.max():.0f} 不同取值={len(set(cls.tolist()))}"
          f"  ⟹ {'**与顶点无关（§16 的饱和）**' if len(set(cls.tolist()))==1 else '有顶点依赖'}")
    rng = np.random.default_rng(7)
    s0 = rng.random(V) + 0.1
    print(f"\n{'臂':>11} {'max_mass(g=60)':>15} {'gini(60)':>9} {'move(60)':>9}   判定")
    out = {}
    for arm in ("mult", "classcount", "rank", "saturate", "dev"):
        tr = run_arm(A, R, arm, s0, gens=60, L=L, cls=cls)
        v = verdict(tr)
        out[arm] = dict(trajectory=tr, verdict=v)
        print(f"{arm:>11} {tr['max_mass'][-1]:15.4f} {tr['gini'][-1]:9.4f} "
              f"{(tr['move'][-1] if tr['move'][-1] is not None else float('nan')):9.5f}   {v}")
    print(f"\n  轨迹抽样（max_mass, 每 10 代）：")
    for arm in out:
        mm = out[arm]["trajectory"]["max_mass"]
        print(f"    {arm:>11}: {[mm[i] for i in (0, 9, 19, 29, 39, 49, 59)]}")
    RES["arms"] = {k: {"verdict": v["verdict"], "max_mass": v["trajectory"]["max_mass"],
                       "gini": v["trajectory"]["gini"], "move": v["trajectory"]["move"]}
                   for k, v in out.items()}
    RES["class_count"] = dict(n_classes=nclass, per_vertex_min=float(cls.min()),
                              per_vertex_max=float(cls.max()),
                              vertex_independent=len(set(cls.tolist())) == 1)
    check("类计数确实与顶点无关（§16 的饱和机制在 Γ 上成立）",
          len(set(cls.tolist())) == 1, f"类数={nclass}")
    pers = [k for k, v in out.items() if v["verdict"] == "持续多畴"]
    check("五条无参数规则里没有任何一条给出「持续多畴」", len(pers) == 0,
          f"持续多畴的臂={pers}" if pers else "全部落到 凝聚/冻结")
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_nonlin.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_nonlin.json")
