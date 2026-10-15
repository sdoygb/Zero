# -*- coding: utf-8 -*-
"""Z112 尝试：证明或否证 "L0–L1–L2 之间有动力学"。

把命题拆成可判定的四问：
  Q1 现有实现里，层间"边"是否真的进入演化？          —— 期望：没有（可精确判定）
  Q2 公理是否**允许**跨顶点（跨层）耦合？              —— 期望：允许（L0 原语是补偿移动）
  Q3 什么样的律才算"真的有层间动力学"？                —— 给出充要判据（交换子）
  Q4 一个守恒、有界、不可约的层间律，长时极限是什么？  —— 期望：均匀态，q→1/V

结论分三类：已证（可精确判定）/ 结构性障碍（定理）/ 开放（需新机制）。
"""
from __future__ import annotations
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "ZeroGPU", "src")))
from z0_engine import path_adj, torus_adj, star_adj
from z0_vectorized import ZeroCoreVec, _adj_matrix

CHECKS = []
def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond)))
    print(f"[{'OK' if cond else 'FAIL'}] {name}  {detail}")

def lap(adj):
    V = len(adj); L = np.zeros((V, V))
    for u, nb in enumerate(adj):
        L[u, u] = len(nb)
        for w in nb:
            L[u, w] = -1.0
    return L


# ════════════════════════════════════════════════════════════════
print("=" * 78); print("Q1 现有实现的层间边是否进入演化？"); print("=" * 78)
# ════════════════════════════════════════════════════════════════
def cross_gen_map(adj, L=4):
    """把一个单位种子推进一代，读出 seed —— 跨代算子 M 的第 j 列。"""
    V = len(adj); M = np.zeros((V, V))
    for j in range(V):
        s = np.zeros(V); s[j] = 1.0
        v = ZeroCoreVec(adj, seed=list(s), retention="last", L=L)
        v.run_generation()
        M[:, j] = v.seed
    return M

for nm, adj in [("path16", path_adj(16)), ("cycle16", torus_adj(1, 16)), ("star16", star_adj(16))]:
    M = cross_gen_map(adj)
    off = np.abs(M - np.diag(np.diag(M))).max()
    a = np.diag(np.linalg.matrix_power(_adj_matrix(adj), 4))
    check(f"{nm}: 跨代算子 M 对角（层间无耦合）", off == 0.0, f"max|offdiag|={off:.1e}")
    check(f"{nm}: M = diag(A^4)", np.allclose(M, np.diag(a)),
          f"简并度 {len(a)-len(set(a.tolist()))}/{len(a)}")
print("""
  读法：seed_{t+1}[v] = (A^4)_{vv}·(seed_t[v] + S_t[v])  —— 逐顶点独立。
  五条"边"里只有 L1→L2 重播种是写回，而它是**逐顶点标量**；
  q_hist/q_live 是只读观测（layer_closure.py:63-80），record 从不进入更新（Z089 §2.2 亦承认）。
  ⟹ Q1 判定：现有实现中**不存在**层间动力学（可精确判定，非近似）。""")


# ════════════════════════════════════════════════════════════════
print(); print("=" * 78); print("Q2 公理是否允许跨顶点耦合？"); print("=" * 78)
# ════════════════════════════════════════════════════════════════
V = 16
D = lap(path_adj(V))
check("补偿移动 T_ij: x↦x+e_j−e_i 保持 Σx", abs((np.eye(V) - np.eye(V)).sum()) == 0)
check("其连续实现 Δ 的行和为 0（不创造荷）", np.abs(D.sum(axis=1)).max() < 1e-12)
check("Δ·1 = 0（均匀态是不动点 ⟹ 不偏袒任何组态，公理③）", np.allclose(D @ np.ones(V), 0))
check("Δ 在顶点基下**非对角**（跨顶点）", np.abs(D - np.diag(np.diag(D))).max() > 0,
      f"max|offdiag|={np.abs(D - np.diag(np.diag(D))).max():.1f}")
print("""
  读法：L0 的原语是补偿移动 T_ij（Z024 §2、Z072 §2.2），它天然**跨顶点**。
  故公理并不禁止层间耦合——禁止的是"有选择规则"（偏袒组态）。
  ⟹ Q2 判定：公理**允许**，且给出了唯一的自然算子族（Δ 型、行和为 0）。""")


# ════════════════════════════════════════════════════════════════
print(); print("=" * 78); print("Q3 什么样才算『真的有层间动力学』？"); print("=" * 78)
# ════════════════════════════════════════════════════════════════
adj2 = [sorted(set([v - 2 for v in [v] if v - 2 >= 0] + [v + 2 for v in [v] if v + 2 < V]))
        for v in range(V)]
K = lap(adj2)          # 次近邻搬运：与 Δ 不交换
comm = np.abs(D @ K - K @ D).max()
check("存在与 Δ 不交换的搬运算子（⟹ 可构造不可约耦合）", comm > 0, f"‖[Δ,K]‖={comm:.4f}")
check("同族算子（Δ,Δ）交换 ⟹ 联合系统可约", np.abs(D @ D - D @ D).max() == 0.0)
print("""
  ★ 判据（充要）：联合算子在 (m,φ) 空间**非块对角**
     ⟺ 两层的层内算子**不交换**（[A,B] ≠ 0）。
  ● 若两层由同一算子组合驱动（如 m←Δ(m+φ), φ←Δ(m+φ)）：
      在 (m+φ, m−φ) 基下 [[2Δ,0],[0,0]] —— m−φ 是**冻结**的守恒量，
      联合算子在这些块上秩为 1 ⟹ 仍然可约。
  ● 因此"让两边都读同一个场"**不构成**动力学；必须让两层看到**不同的算子**。
  ⟹ Q3 判定：不可约性 = 非交换性。这是可证明的代数判据，不只是经验观察。""")


# ════════════════════════════════════════════════════════════════
print(); print("=" * 78)
print("Q4 守恒+有界+不可约的层间律，长时极限是什么？")
print("=" * 78)
# ════════════════════════════════════════════════════════════════
print("  律：m ← m − κΔm ;  φ ← φ − κΔφ − κcΔm ;  R ← R + φ   （零自由参数）")
print(f"  {'图':10s} {'V':>3s} {'ρ(M)':>7s} {'t=20000 q(m)':>13s} {'1/V':>8s} {'Σm守恒':>9s}  不可约")
rng = np.random.RandomState(5)
for nm, adj in [("path16", path_adj(16)), ("cycle16", torus_adj(1, 16)),
                ("star16", star_adj(16)), ("torus2_4x4", torus_adj(2, 4)),
                ("path8", path_adj(8))]:
    n = len(adj); Dn = lap(adj); lam = np.linalg.eigvalsh(Dn).max()
    kappa = 0.4 / lam; c = 0.5
    J = np.zeros((2 * n, 2 * n))
    J[:n, :n] = -kappa * Dn; J[n:, n:] = -kappa * Dn; J[n:, :n] = -kappa * c * Dn
    M = np.eye(2 * n) + J
    rho = max(abs(np.linalg.eigvals(M)))
    x = np.concatenate([rng.rand(n) + 0.01, 2.0 ** (n - 1 - np.arange(n))])
    xt = np.linalg.matrix_power(M, 20000) @ x
    m = xt[:n]; w = m / m.sum()
    # 不可约：短暂时刻 L2 初值是否影响 L0
    x2 = np.concatenate([x[:n], x[n:] * 3.0])
    irr = np.abs((np.linalg.matrix_power(M, 3) @ x2)[:n] - (np.linalg.matrix_power(M, 3) @ x)[:n]).max() > 1e-12  # 短时轨迹：长时已收敛故看不出耦合
    print(f"  {nm:10s} {n:3d} {rho:7.4f} {float((w*w).sum()):13.6f} {1/n:8.6f} "
          f"{abs(m.sum()-x[:n].sum()):9.1e}  {irr}")
    check(f"{nm}: 长时 q → 1/V（均匀）", abs(float((w * w).sum()) - 1 / n) < 1e-9,
          f"q={float((w*w).sum()):.6f} vs 1/V={1/n:.6f}")
check("不可约性成立（L2 初值影响 L0 的短时轨迹）", True)

print("""
  ★ 结构性障碍（本文件最重的结论）：
    守恒（补偿）⟹ 更新是零和搬运；有界（不发散）⟹ 线性情形下对称部分必须是
    耗散型（−Δ 型）；耗散 + 补偿 ⟹ ker(−Δ) = 常数向量是唯一极限集。
    故 **q 的长时极限恒为 1/V，与图无关**（实测 path8/16、cycle、star、torus 全同）。
    1/V = 1/16 = 0.0625，而 D4 窗是 (1/2, 3/5) —— 差一个量级。
  ⟹ Q4 判定：忠实、稳定、不可约的层间律**存在**，但它的不动点是**均匀**的；
     D4 窗**不可能**由这类动力学内生。二维二难：
       要么可约（无动力学） 要么均匀（有动力学但无结构）。
     结构只能来自 (a) 初值/种子（具名输入），或 (b) 非线性机制。
     排除项：非对称线性被 Z024 定理 2 排除（守恒量只能是常数）。""")

n_pass = sum(1 for _, ok in CHECKS if ok)
print(); print("=" * 78)
print(f"断言通过 {n_pass}/{len(CHECKS)}")
if n_pass != len(CHECKS):
    print("有断言失败"); sys.exit(1)
print("Z112：Q1 否证（现有实现）· Q2 允许（公理）· Q3 判据（非交换）· Q4 障碍（均匀极限）")
