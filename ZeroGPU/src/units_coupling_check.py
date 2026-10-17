#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
units_coupling_check.py —— 单位/耦合账：κ 是约定、G 的量纲缺口
==============================================================
处理 D18–D24 审计报告 §7 的缺口 2（K 与 s₀）与缺口 5（质量/能量标度）。

C1 κ（边导纳整体值）是能量单位的约定：L→κL ⟹ φ→φ/κ，且**比值**可观测
C2 谱只被 κ 整体缩放 ⟹ κ 不进入任何无量纲结论
C3 D23 §6.3 的「K=1/(4π)」与 §6.2 的 4πG=1/(Ks₀) 合起来只是 G=1/s₀ 的恒等式
   （K 的辨识：K = 1/(4πG_disc)，即 G 在键长单位下 = A/(4π)，A 由舞台的 σ 给出）
C4 量纲秩：G 需 (L,M,T) 三个单位；离散模型只提供 ℓ_b 与 τ 两个 ⟹ 缺秩 1
C5 台账量清单（全是无量纲比值）⟹ 形状可算、绝对标度不可算
C6 Z121 的质量只给比值（√m 比），整体 A 是锚 ⟹ 与 G 同病
运行：cd ZeroGPU/src && /usr/bin/python3 units_coupling_check.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zero_sim as zs  # noqa: E402

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


# ============================================================
print("=" * 72)
print("C1/C2  κ（边导纳的整体值）是单位约定")
print("=" * 72)
N = 6
idx = {(i, j, k): (i * N + j) * N + k for i in range(N) for j in range(N) for k in range(N)}
n = len(idx)
L0 = np.zeros((n, n))
for (i, j, k), u in idx.items():
    for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        q = ((i + d[0]) % N, (j + d[1]) % N, (k + d[2]) % N)
        v = idx[q]
        L0[u, u] += 1
        L0[v, v] += 1
        L0[u, v] -= 1
        L0[v, u] -= 1
rho = np.zeros(n)
rho[0] = 1.0
rho -= rho.mean()
phis, ratios = {}, []
src2 = rho.copy()
src2[idx[(2, 2, 2)]] += 1.0                      # 第二个源（用来看比值）
for kappa in (0.25, 1.0, 4.0):
    Lk = kappa * L0
    phi1 = -np.linalg.lstsq(Lk, rho, rcond=None)[0]
    phi2 = -np.linalg.lstsq(Lk, src2, rcond=None)[0]
    phis[kappa] = float(phi1[idx[(1, 0, 0)]])
    ratios.append(float(phi1[idx[(1, 0, 0)]] / phi2[idx[(1, 0, 0)]]))
print(f"  κ=0.25/1/4 时 φ(1,0,0) = {[round(phis[k], 4) for k in (0.25, 1.0, 4.0)]}")
print(f"  两个源在同一测点的响应**比值** = {[round(r, 6) for r in ratios]}")
ck("C1a φ ∝ 1/κ（整体值只缩放）",
   abs(phis[0.25] / phis[1.0] - 4) < 0.02 and abs(phis[4.0] / phis[1.0] - 0.25) < 0.02,
   f"κ 缩小 4 倍 ⟹ φ 放大 {phis[0.25]/phis[1.0]:.3f} 倍")
ck("C1b 可观测量（同舞台两源的响应比）与 κ 无关",
   max(ratios) - min(ratios) < 1e-9, f"比值 {ratios} 在 κ 变化 16 倍下不变 ⟹ κ 不可测（Z002 引理 4 只给均匀性）")
lam0 = np.sort(np.linalg.eigvalsh(L0))
lamk = np.sort(np.linalg.eigvalsh(4.0 * L0))
ck("C2a 谱只被 κ 整体缩放（不改变任何无量纲谱量）",
   np.allclose(lamk, 4.0 * lam0, atol=1e-9),
   f"λ(K₄ 6³ 环面) 前 3 个/后 1 个：{np.round(lam0[:3],4).tolist()} / {lam0[-1]:.3f}；"
   "κ 只改整体 ⟹ 谱隙比、维数、围长等无量纲量都不含 κ")

# ============================================================
print()
print("=" * 72)
print("C3  D23 的『K=1/(4π)』是什么")
print("=" * 72)
A_cubic = 1.0                                     # stage_ledger_check.py B1：单位立方舞台 A=1
G_disc = A_cubic / (4 * np.pi)
print(f"  单位立方舞台（键长 1）：势幅度 A=1 ⟹ G_disc = A/(4π) = {G_disc:.5f}")
print(f"  D23 §6.3 报『测量的 K = 0.08014』，理论 1/(4π) = {1/(4*np.pi):.5f}，偏差 0.70%")
print(f"  D23 §6.2 的关系 4πG = 1/(K s₀)（s₀=1）⟹ K ≡ 1/(4πG)：与上式**同一个数**")
ck("C3a D23 的『K=1/(4π)』是独立于 G 的纤维弹性系数",
   False,
   f"1/(4π)={1/(4*np.pi):.5f} 与 G_disc={G_disc:.5f} 是彼此倒数（K≡1/(4πG)）；"
   "所谓『数值验证 K』= 复读『G_disc=A/4π』，A=1 又是单位立方舞台的约定 ⟹ 循环")
note("K 的正确辨识", "K = 1/(4πG_disc) = 1/A", "A 由舞台键几何（σ、ℓ_b）决定，见 stage_ledger_check.py B1")

# ============================================================
print()
print("=" * 72)
print("C4  量纲秩：G 需要三个单位，离散模型只给两个")
print("=" * 72)
# 单位向量（L, M, T）
U_len = np.array([1.0, 0, 0])          # ℓ_b
U_time = np.array([0, 0, 1.0])         # τ
U_c = U_len - U_time                   # c = ℓ/τ
G_dim = np.array([3.0, -1.0, -2.0])    # G = [L³ M⁻¹ T⁻²]
M_un = np.array([U_len, U_time, U_c])
rank = int(np.linalg.matrix_rank(M_un))
print(f"  可用单位：ℓ_b={U_len}、τ={U_time}、c={U_c}；矩阵秩 = {rank}")
print(f"  G 的量纲向量 = {G_dim}（L³M⁻¹T⁻²）")
ck("C4a 离散模型（ℓ_b、τ、c）能定出 G 的量纲组合",
   rank == 3,
   f"秩 {rank} < 3 ⟹ 缺 1 个独立单位（质量）；G 的表达式里必有一个自由指数，"
   "这与 D21 §5.2 的 G_macro = G_disc·a³/τ² 缺 kg⁻¹ 是同一件事")
note("缺失的那一个自由度", "质量单位（Z025 定理 3：Zero 无物质扇区）", "长度与时间单位可由 ℓ_b、τ 给，质量单位台账里没有")

# ============================================================
print()
print("=" * 72)
print("C5  台账量清单：全部是无量纲比值")
print("=" * 72)
LEDGER = {"log2（Z087 稳态 2:1）": np.log(2), "log(3/2)、log(8/3)（Z131 三频率槽）": np.log(1.5),
          "层权 2/3 : 1/3（Z079）": 2 / 3, "κ₁ = 5/9（Z079）": 5 / 9,
          "2cos(2πk/8)=±√2、0（Z121）": np.sqrt(2), "Koide 2/3（Z121 质量比）": 2 / 3,
          "δ₀ = π/12 − 2/9（Z123）": np.pi / 12 - 2 / 9, "Heawood |V|=14、6 边、度 3（Z021）": 14,
          "seed 比 [8,4,1.5,1]（Z131）": 8 / 1.0}
for k, v in LEDGER.items():
    print(f"    {k:<34} = {v:.6f}（纯数）")
ck("C5a 台账量里有能定绝对标度的有量纲量",
   False,
   "清单里全部是无量纲比值（频率比、层权、计数、角）；Z135 已把『能量/波长标度』登记为具名输入② ⟹ "
   "任何有量纲常数（G、c、ħ、Λ、m）的**数值**都需要至少一个锚")
note("Z121 的对照", "质量只给比值（√m 比 [1,14.38,58.97]）", "整体因子 A 未定，δ₀ 判为 per-界原初数据 ⟹ 与 G 同病")

# ============================================================
print()
print("=" * 72)
n_fail = sum(1 for c in CK if not c["pass"])
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 72)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "units_coupling_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT, "phi_by_kappa": {str(k): v for k, v in phis.items()},
               "ratio_by_kappa": ratios, "G_disc_unit_cubic": G_disc,
               "dim_rank": rank, "ledger": {k: float(v) for k, v in LEDGER.items()}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/units_coupling_check.json")
sys.exit(1 if n_fail else 0)
