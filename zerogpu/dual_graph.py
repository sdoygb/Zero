"""dual_graph.py --- L0 的**双图对照**：语料里 "L0 内容图" 的两个不同实现。

  (i)  E1 词图  ：节点 = 长度 L 的全部平衡 ±1 **词**（C(L,L/2) 个）
                  边   = **线性**相邻异号对换 i=0..L-2（不含绕回对 (L-1,0)）
                  出处 ：E1_NN_check.py:59-78（`moves()` 用 range(len(w)-1)）
  (ii) 闭链图   ：节点 = 长度 T 的平衡 **项链**（旋转类）
                  边   = **循环**相邻异号对换 (i, (i+1) mod T)，**简单图**去重
                  出处 ：simulations/zero_sum_geometry_probe.py:102-114

LAYER_LEDGER §3.1 把「Γ_L 恒有 1 个零模；λ1≈π²/L²；λmax→4L/3」列为 **L0** 事实，
但没说清是哪张图。本节把两者放在同一批可观测量下实测。
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, json, os
from zcl import Engine
import l0_closure as L0
from analyze_l0 import e1_word_graph, cyc_necklace_graph


def spec(n, E):
    A = sp.coo_matrix((np.ones(len(E)), (E[:, 0], E[:, 1])), shape=(n, n)).tocsr()
    A = A.maximum(A.T)
    deg = np.asarray(A.sum(1)).ravel()
    Lap = (sp.diags(deg) - A).tocsr()
    if n < 8:
        ev = np.sort(np.linalg.eigvalsh(Lap.toarray()))
        nz = int((np.abs(ev) < 1e-5).sum())
        return nz, float(ev[nz]), float(ev[-1]), float(deg.mean())
    lmax = float(spla.eigsh(Lap, k=1, which="LA", return_eigenvectors=False)[0])
    lo = np.sort(spla.eigsh(Lap, k=3, which="SA", return_eigenvectors=False))
    nz = int((np.abs(lo) < 1e-5).sum())
    return nz, float(lo[nz]), lmax, float(deg.mean())


eng = Engine()
res = {}
print("=" * 78)
print("§A1  E1 词图（平衡词 + **线性**对换，无绕回）  —— C(L,L/2) 个节点")
print(f"{'L':>3} {'N':>8} {'边':>8} {'均度':>8} {'零模':>5} {'λ1':>10} {'λ1·L²':>9} {'λmax':>9} {'λmax/L':>8}")
e1 = []
for L in range(4, 17, 2):
    n, r, c = e1_word_graph(L)
    E = np.unique(np.sort(np.stack([r, c], 1), 1), axis=0)
    nz, l1, lmax, dm = spec(n, E)
    e1.append(dict(L=L, N=n, edges=len(E), mean_deg=round(dm, 4), zero=nz, l1=l1,
                   l1_L2=round(l1 * L * L, 4), lmax=round(lmax, 4), lmax_L=round(lmax / L, 4)))
    print(f"{L:>3} {n:>8,} {len(E):>8,} {dm:>8.4f} {nz:>5} {l1:>10.6f} "
          f"{l1*L*L:>9.4f} {lmax:>9.4f} {lmax/L:>8.4f}")
res["e1_word_graph"] = e1

print("\n" + "=" * 78)
print("§A2  闭链图（项链 + **循环**对换，简单图）  —— K(T) 个节点")
print(f"{'T':>3} {'N':>8} {'边':>9} {'均度':>8} {'零模':>5} {'λ1':>10} {'λ1·T²':>9} {'λmax':>9} {'λmax/T':>8}")
cy = []
for T in range(4, 19, 2):
    reps, _ = L0.enumerate_necklaces(T, eng, verbose=False)
    n, r, c = cyc_necklace_graph(T, reps, eng)
    E = np.unique(np.sort(np.stack([r, c], 1), 1), axis=0)
    nz, l1, lmax, dm = spec(n, E)
    cy.append(dict(T=T, N=n, edges=len(E), mean_deg=round(dm, 4), zero=nz, l1=l1,
                   l1_T2=round(l1 * T * T, 4), lmax=round(lmax, 4), lmax_T=round(lmax / T, 4)))
    print(f"{T:>3} {n:>8,} {len(E):>9,} {dm:>8.4f} {nz:>5} {l1:>10.6f} "
          f"{l1*T*T:>9.4f} {lmax:>9.4f} {lmax/T:>8.4f}")
res["cyc_necklace_graph"] = cy
os.makedirs("results", exist_ok=True)
json.dump(res, open("results/dual_graph.json", "w"), indent=2)
print("\n判读：均度 L/2、λ1·L²→9.838、λmax/L→4/3 三条**只属于 §A1 的词图**；")
print("      §A2 的闭链图均度 ≈ (T+1)/2、λ1·T² ≈ 87-91（无 π² 律）、λmax/T → ~1.95。")
print("saved results/dual_graph.json")
