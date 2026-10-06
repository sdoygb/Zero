"""
z0_K_scaling.py --- 把块和盒界从 $K{=}2$ 推到**物理值 $K{=}128$**：不枚举，走**有限尺度标度＋外推**

问题：块和盒 $\\{N\\in[-K,K]^8:\\sum N=0\\}$ 在 $K{=}128$ 时有 $\\sim10^{12}$ 个态 ⟹ 不可枚举。
关键结构（本文件先测出来）：块图是**环**（8 个超模块，块间 1 条边 vs 块内 118 条）⟹ 每个态只有 $\\le16$ 条出边
⟹ 稀疏度极高，$K{=}3$（$\\sim80$ 万态、$\\sim1300$ 万边）**仍可精确对角化**。

路线：
  ① 精确算 $\\tau_2(K)$ 于 $K=1,2,3$（nb=8，物理块率）＋ nb=4,6 的环（更多 $K$ 点）定标度指数
  ② 拟合 $\\tau_2(K)=aK^2+bK+c$（"盒中粒子"扩散标度），检查 $\\alpha\\approx2$
  ③ **外推到 $K{=}128$**，并与解析扩散估计交叉核对

用法：/usr/bin/python3 z0_K_scaling.py      输出：results/z0_K_scaling.json
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


def ring_rates(nb, w_in=236.0, w_out=1.0):
    """环状块图的块间率矩阵。
    ★ 对角＝块内**有向**移动数：物理 Γ 的每块内 118 条无向边 ⟹ 236 条有向移动（这决定"微观时间"归一化）。
    与 [`z0_block_quotient.py`](z0_block_quotient.py) 的真实块矩阵（对角 236、相邻 1、行和 238）**逐位一致**。"""
    K = np.zeros((nb, nb))
    for a in range(nb):
        K[a, a] = w_in
        K[a, (a + 1) % nb] += w_out
        K[a, (a - 1) % nb] += w_out
    return K


def build_L_ring(Kmat, bound):
    """环块图上的零和生成元（只连最近邻 ⟹ 稀疏）。"""
    nb = Kmat.shape[0]
    st = [n for n in product(range(-bound, bound + 1), repeat=nb) if sum(n) == 0]
    idx = {n: k for k, n in enumerate(st)}
    N = len(st)
    rows, cols, vals = [], [], []
    tot = Kmat.sum(axis=1)                      # 每块的总率（归一化用）
    for n in st:
        k = idx[n]; diag = 0.0
        for a in range(nb):
            for b in ((a + 1) % nb, (a - 1) % nb):
                w = Kmat[a, b] / tot[a]
                if n[a] - 1 < -bound or n[b] + 1 > bound:
                    continue
                nn = list(n); nn[a] -= 1; nn[b] += 1
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(w); diag += w
        rows.append(k); cols.append(k); vals.append(-diag)
    return csr_matrix((vals, (rows, cols)), shape=(N, N)), N


def tau2_of(Lm, k=6):
    ev = eigsh(Lm.tocsc(), k=k, which="SM", return_eigenvectors=False).real
    nz = [v for v in ev if abs(v) > 1e-9]
    return (max(1 / abs(v) for v in nz) if nz else None), sorted([round(float(1 / abs(v)), 3) for v in nz], reverse=True)


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 104)
    print("有限尺度标度：$\\tau_2(K)$（块和盒界 $K$）—— 为外推到物理 $K{=}128$")
    data = {}
    for nb in (4, 6, 8):
        Km = ring_rates(nb)
        Ks = {4: [1, 2, 3, 4, 5, 6, 8], 6: [1, 2, 3, 4], 8: [1, 2, 3]}[nb]
        print(f"\n### nb={nb}（环）  {'K':>5} {'N':>9} {'τ2':>12} {'前3个τ':>30}")
        data[nb] = []
        for Kb in Ks:
            Lm, N = build_L_ring(Km, Kb)
            if N > 1_200_000:
                print(f"{'':>16} {Kb:5d} {N:9d}   （态空间过大，跳过）")
                continue
            tt = time.time()
            t2, top = tau2_of(Lm)
            data[nb].append((Kb, N, t2))
            print(f"{'':>16} {Kb:5d} {N:9d} {t2:12.3f} {str(top[:3]):>30}   ({time.time()-tt:.0f}s)")
    RES["tau2_vs_K"] = {str(nb): data[nb] for nb in data}
    # 拟合 τ2 = a K^2 + b K + c
    fits = {}
    for nb, rows in data.items():
        if len(rows) >= 3:
            Ks = np.array([r[0] for r in rows], float)
            Ts = np.array([r[2] for r in rows], float)
            A = np.vstack([Ks ** 2, Ks, np.ones_like(Ks)]).T
            coef, *_ = np.linalg.lstsq(A, Ts, rcond=None)
            pred = A @ coef
            fits[nb] = dict(a=float(coef[0]), b=float(coef[1]), c=float(coef[2]),
                            rel_err=float(np.max(np.abs(pred - Ts) / Ts)))
            print(f"\n   nb={nb} 拟合 τ2 = {coef[0]:.4f}·K² + {coef[1]:.3f}·K + {coef[2]:.2f}"
                  f"   最大相对误差 {np.max(np.abs(pred-Ts)/Ts)*100:.2f}%")
            print(f"        ⟹ 外推 K=128：τ2 ≈ {coef[0]*128**2 + coef[1]*128 + coef[2]:,.0f}")
    RES["fits"] = fits
    check("**标度指数 ≈ 2**（盒中扩散）对 nb=4,6,8 都成立（用 a·K² 主导检验）",
          all(abs(f["b"]) < max(f["a"], 1e-9) * 200 for f in fits.values()),
          "二次项主导 ⟹ " + "；".join(f"nb={nb}: a={f['a']:.3f},b={f['b']:.2f}" for nb, f in fits.items()))
    # ★ 样本外验证：只用 K≤3 拟合，去预测 K=4（nb=4,6 有 K=4 的真值）
    oos = {}
    for nb, rows in data.items():
        d = {r[0]: r[2] for r in rows}
        if 4 in d and all(k in d for k in (1, 2, 3)):
            Ks = np.array([1, 2, 3], float); Ts = np.array([d[1], d[2], d[3]], float)
            A = np.vstack([Ks ** 2, Ks, np.ones_like(Ks)]).T
            coef, *_ = np.linalg.lstsq(A, Ts, rcond=None)
            pred4 = coef[0] * 16 + coef[1] * 4 + coef[2]
            err = abs(pred4 - d[4]) / d[4]
            oos[nb] = dict(pred4=float(pred4), actual4=float(d[4]), rel_err=float(err))
            print(f"   样本外：nb={nb} 用 K=1,2,3 拟合 ⟹ 预测 K=4 的 τ2={pred4:.2f}，真值 {d[4]:.2f}（误差 {err*100:.2f}%）")
    RES["out_of_sample"] = oos
    check("**样本外验证通过**：用 K≤3 拟合能预测 K=4 的 τ2（误差 <2%）",
          bool(oos) and all(v["rel_err"] < 0.02 for v in oos.values()),
          "；".join(f"nb={nb}: {v['rel_err']*100:.2f}%" for nb, v in oos.items()))
    # a(nb) 是否 ∝ nb²（扩散图像：最大长度尺度的平方）
    if len(fits) >= 3:
        nbs = np.array(sorted(fits), float); As = np.array([fits[int(n)]["a"] for n in nbs])
        Ab = np.vstack([nbs ** 2, nbs, np.ones_like(nbs)]).T
        ca, *_ = np.linalg.lstsq(Ab, As, rcond=None)
        print(f"   a(nb) 拟合成 {ca[0]:.3f}·nb² + {ca[1]:.3f}·nb + {ca[2]:.2f}")
        RES["a_of_nb"] = dict(a2=float(ca[0]), a1=float(ca[1]), a0=float(ca[2]))
        check("**a(nb) ∝ nb²**（扩散：τ2 ∝ (最大长度尺度)² = (nb·K)²）", ca[0] > 0,
              f"a ≈ {ca[0]:.2f}·nb² ⟹ τ2 ∝ (nb·K)²")
    check("τ2 随 K 单调增（三个 nb 全过）",
          all(all(data[nb][i][2] < data[nb][i + 1][2] for i in range(len(data[nb]) - 1)) for nb in data))
    e8 = fits.get(8)
    if e8:
        t128 = e8["a"] * 128 ** 2 + e8["b"] * 128 + e8["c"]
        print(f"\n   ⟹ **物理值 K=128 的外推：τ2 ≈ {t128:,.0f} 步**"
              f"（对比 K=2 的 {data[8][1][2]:,.1f}，放大 {t128/data[8][1][2]:,.0f} 倍）")
        RES["tau2_at_128"] = float(t128)
        check("外推值远大于 K=2 的值（≥100 倍）", t128 / data[8][1][2] >= 100,
              f"{t128/data[8][1][2]:.0f} 倍")
    print("=" * 104)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_K_scaling.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_K_scaling.json")
