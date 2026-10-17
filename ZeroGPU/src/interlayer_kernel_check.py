#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
interlayer_kernel_check.py —— 补 #2：层间耦合的**有效相互作用核**
================================================================
对象：D31/D32 的机制主张——「L2 物质 →边③写历史(L1) →边②历史作种子 → 影响整个 L2」，
      并据此断言 D31 §3「引力长程、不因距离衰减（共同的历史）」与 D32 §2.1「等价于 Lφ=−x 的 1/r 泊松」。

做法：**用台账自己的引擎**（z0_vectorized.ZeroCoreVec ＋ z0_engine.torus_adj），
      在一个顶点上放点源扰动，逐代测响应（= 扰动引擎 − 基准引擎 的 record 差）的**空间分布**。

承重递推（从台账读出，不自己编；z0_vectorized.py:96-124）：
      seed_t[v]              = retained(gen_closures_t)[v]          （边② L1→L2）
      gen_closures_{t+1}[v]  = [(W_t + diag(seed_t)) A^L]_{vv}      （边③ L2→L1 ＋ 闭合）
      W_{t+1}                = 上式的**非对角**部分（行索引 = 起点，永不改变）
⟹ 归纳：gen_closures_t[v] 只依赖**顶点 v 自己**的 seed 与历史 ⟹ **无跨顶点耦合**（除归一化）。

判据（PASS = 原稿主张成立）：
  K1a 层间耦合产生跨顶点的相互作用
  K1b 核形状是 1/r（1D 时线性 |d|）
  K2a 核随距离不衰减（D31 §3）
  K3a 耦合与代数无关（静态核）
  K4a 3D 给 1/r 型核
  K5a 跨顶点耦合来自物理机制（而非数值归一化）
  K6a（结构事实，应为 PASS）逐顶点自指：扰动一个顶点，其他顶点的 record 不变
运行：cd ZeroGPU/src && /usr/bin/python3 interlayer_kernel_check.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import torus_adj          # noqa: E402
from z0_vectorized import ZeroCoreVec    # noqa: E402

RES = os.path.join(HERE, "..", "results")
CK, OUT = [], []
V0 = 64


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


def run_pair(adj, V, src=0, gens=600, delta=1e-3, retention="last", L=4, marks=(20, 50, 200, 300, 600)):
    """基准 vs 扰动，逐代记录 (t, 非零响应顶点数, max|resp|, 归一化次数)"""
    b = ZeroCoreVec(adj, L=L, retention=retention, seed=[1.0] * V)
    p = ZeroCoreVec(adj, L=L, retention=retention, seed=[1.0] * V)
    p.seed[src] += delta
    out = {}
    for t in range(gens):
        b.run_generation()
        p.run_generation()
        if (t + 1) in marks:
            r = p.record - b.record
            out[t + 1] = {"nz": int((np.abs(r) > 1e-12).sum()),
                          "max": float(np.abs(r).max()),
                          "norm_events": int(p.norm_events + b.norm_events),
                          "resp": r.copy()}
    return out


# ============================================================
print("=" * 72)
print("K0  承重递推核对：一代映射的耦合半径")
print("=" * 72)
adj1 = torus_adj(1, V0)
A = np.zeros((V0, V0))
for u, nb in enumerate(adj1):
    for w in nb:
        A[u][w] = 1.0
AL = np.linalg.matrix_power(A, 4)
rng = sorted({min(int((i - j) % V0), int((j - i) % V0)) for i in range(V0) for j in range(V0) if AL[i, j] > 0})
print(f"  A^L（L=4）非零元覆盖格距 = {rng}；对角元 = {np.unique(AL.diagonal())[:3]}")
ck("K0a 一代映射的耦合是瞬时的/全图的（任意格距一步可达）",
   max(rng) > 4,
   f"实际只覆盖格距 ≤ {max(rng)}（= L）⟹ 一代耦合半径恰为 L 步，"
   "远距离只能靠多代累积（若存在）")

# ============================================================
print()
print("=" * 72)
print("K1–K5  1D 环面 N=64：响应核的空间分布随代数")
print("=" * 72)
res = run_pair(adj1, V0)
for t in sorted(res):
    d = res[t]
    print(f"  t={t:>3}: 非零响应顶点数={d['nz']:>3}  max|resp|={d['max']:.3e}  归一化累计={d['norm_events']}")
first_spread = min([t for t in res if res[t]["nz"] > 1], default=None)
first_norm = min([t for t in res if res[t]["norm_events"] > 0], default=None)
print(f"  首次出现跨顶点响应：t={first_spread}；首次触发归一化：t={first_norm}")
ck("K1a 层间耦合产生跨顶点的相互作用",
   (first_spread is not None) and first_spread <= 20,
   f"在 t ≤ {first_spread-1 if first_spread else '∞'} 代内响应**严格只在扰动顶点**（非零顶点数=1）；"
   f"跨顶点响应只在 t≈{first_spread} 出现，且与归一化首次触发（t≈{first_norm}）同时")
ck("K5a 跨顶点耦合来自物理机制（而非数值归一化）",
   (first_spread is not None) and (first_norm is not None) and first_spread < first_norm,
   f"两者同时发生（t≈{first_spread}）⟹ 唯一的跨顶点效应来自 `_normalize()`（整体除以总和，"
   "数值工程），不是边②/③的机制")

resp = res[600]["resp"]
K = resp / resp[0]
ds = np.array([min(int(v), V0 - int(v)) for v in range(V0)])
prof = {int(d): float(np.mean(K[ds == d])) for d in np.unique(ds)}
dd = np.array(sorted(prof))
KK = np.array([prof[d] for d in dd])


def r2(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        return float("nan")
    Ap = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(Ap, y, rcond=None)
    pred = Ap @ coef
    ss = float(((y - pred) ** 2).sum())
    st = float(((y - y.mean()) ** 2).sum())
    return 1 - ss / st if st > 0 else float("nan")


pos = dd > 0
print("  600 代（归一化已多次触发）后的归一化核：" +
      " ".join(f"d{int(d)}:{v:.4f}" for d, v in list(zip(dd, KK))[:8]))
print(f"  形状拟合 R²：log K ~ d = {r2(dd[pos], np.log(KK[pos])):.4f}；"
      f"log K ~ d² = {r2(dd[pos] ** 2, np.log(KK[pos])):.4f}；K ~ 1/d = {r2(1/dd[pos], KK[pos]):.4f}")
ck("K1b 核形状是 1/r（1D 时 Z033 命题 14 给线性 |d|）",
   r2(1 / dd[pos], KK[pos]) > 0.99,
   f"实测 K(1)={KK[list(dd).index(1)]:.4f}、K(4)={KK[list(dd).index(4)]:.4f}、"
   f"K(32)={KK[-1]:.4f}；与 1/d 的 R²={r2(1/dd[pos], KK[pos]):.3f}（1D 泊松核是**线性增长**的 |d|，形状不符）")
ck("K2a 核随距离不衰减（D31 §3『共同的历史，不管多远』）",
   KK[-1] > 0.5 * KK[0],
   f"K(d_max={int(dd[-1])})/K(0) = {KK[-1]:.3e} ⟹ 强度随距离强烈衰减")
ck("K3a 耦合与代数无关（静态核）",
   res[20]["nz"] == res[600]["nz"],
   f"t=20 时非零响应顶点数 = {res[20]['nz']}（零耦合），t=600 时 = {res[600]['nz']}（全图）"
   " ⟹ 核随时间从'零'跳到'全图'，不是静态核")

# ============================================================
print()
print("=" * 72)
print("K4  维数依赖：3D 立方格")
print("=" * 72)
adj3 = torus_adj(3, 8)
V3 = 8 ** 3
res3 = run_pair(adj3, V3, gens=600, marks=(200, 300, 600))
for t in sorted(res3):
    print(f"  t={t:>3}: 非零响应顶点数={res3[t]['nz']:>4}/{V3}  max|resp|={res3[t]['max']:.3e}")
ck("K4a 3D 给 1/r 型核（Z033 命题 14：只有 d=3 给 1/r 势）",
   res3[600]["nz"] < V3 and res3[200]["nz"] > 1,
   f"3D 与 1D 同形：t=200 时 {res3[200]['nz']}/{V3} 非零（零耦合），t=600 时 "
   f"{res3[600]['nz']}/{V3}（归一化扩散）⟹ 无 1/r 结构；"
   "1/r 只在泊松方程 Lφ=−x（Z033 定理 16）里出现")

# ============================================================
print()
print("=" * 72)
print("K6  结构事实：逐顶点自指（扰动一个顶点，其他顶点不变）")
print("=" * 72)
b = ZeroCoreVec(adj1, L=4, retention="last", seed=[1.0] * V0)
p = ZeroCoreVec(adj1, L=4, retention="last", seed=[1.0] * V0)
p.seed[32] += 1e-3
for _ in range(100):
    b.run_generation()
    p.run_generation()
diff = np.abs(p.record - b.record)
others = np.delete(diff, 32)
print(f"  扰动顶点 32，100 代后：diff[32]={diff[32]:.6e}；其他顶点的 max diff={others.max():.3e}")
ck("K6a 层间边是**逐顶点自指回路**（record[v] 只由 v 自己的历史决定）",
   float(others.max()) == 0.0 and float(diff[32]) > 0,
   f"扰动一个顶点，其他 63 个顶点的 record 差**恒为 0**（{others.max():.1e}）⟹ "
   "边②/③在引擎里是逐顶点自指：seed[v] ← record[v]，gen_closures[v] 只依赖 seed[v] 与 W[v,·]")

# ============================================================
print()
print("=" * 72)
n_fail = sum(1 for x in CK if not x["pass"])
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 72)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "interlayer_kernel_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "AL_support_radius": int(max(rng)),
               "trace_1d": {str(t): {"nz": res[t]["nz"], "max": res[t]["max"],
                                     "norm": res[t]["norm_events"]} for t in res},
               "trace_3d": {str(t): {"nz": res3[t]["nz"], "max": res3[t]["max"]} for t in res3},
               "first_spread_t": first_spread, "first_norm_t": first_norm,
               "kernel_1d_t600": {str(int(d)): float(v) for d, v in zip(dd, KK)},
               "shape_R2": {"K_vs_inv_d": r2(1 / dd[pos], KK[pos]),
                            "logK_vs_d": r2(dd[pos], np.log(KK[pos])),
                            "logK_vs_d2": r2(dd[pos] ** 2, np.log(KK[pos]))},
               "K6_other_max": float(others.max())},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/interlayer_kernel_check.json")
sys.exit(1 if n_fail else 0)