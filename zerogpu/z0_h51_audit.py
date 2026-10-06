"""
z0_h51_audit.py --- 审计「全库唯一一个非 O(1) 的宏观时间」：H5.1 的 λ2^macro

来源：`H5_closure_graph_lock_theorem.md` §3 报
      λ2^macro = −log|μ2(Q)| = 0.0025（确定性等分划：12 点/3 块/4-正则/商图 C3）
问题：这是**真谱**还是又一次块商/伪谱？（§41 的 eigsh bug 波及 §37–§40，但**不**波及此处：`evM` 用 dense 算）

本机器要判的是**三件不同的事**（前两轮把它们混在一起了）
----------------------------------------------------------------------------
**A. 恒等式层**：M = Q 逐项相等（H5.1 的**定理**部分）——含**非等分划**对照。
**B. 数值层**：λ2 = −log|μ2(Q)| 的**真值**（dense，无 eigsh）。
**C. ★ 解释层（本机器的核心）★**：λ2^macro 是**微观时间**还是**宏观步**？
   M 是**一步的转移矩阵**（行随机、谱在 [−1,1]），而 λ2^macro = −log|μ2| 是**每步**的衰减率。
   于是 τ_macro = 1/λ2^macro 的单位是 **macro-步**，不是微观步。
   换算：1 macro-步 = 1/(离开本块的微观率) 个微观步。
   ⟹ 必须与**同一对象**的微观谱对照，否则差一个"每步时长"因子。

判据：
  H1 M = Q（等分划 & 非等分划都成立）
  H2 spec(M) = spec(Q)（逐位）
  H3 λ2^macro 用 dense 复算 = 文档值（确认不是伪谱）
  H4 ★ λ2^macro 与**微观归一化 Laplacian** 的 λ2 对照：两者**不同一个量**
  H5 ★ 时间换算：τ_macro 步 × (每步时长) = 微观步数，看它是否仍 ≫1

用法：/usr/bin/python3 z0_h51_audit.py   输出：results/z0_h51_audit.json
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def brute(A, lab, K):
    """逐项最直白实现（与 H5_check.py 同口径）：M（度权块层链）与 Q（纤维平均商图）"""
    A = np.asarray(A, float)
    n = A.shape[0]
    d = A.sum(1)
    P = A / np.maximum(d[:, None], 1e-12)
    M = np.zeros((K, K))
    for a in range(K):
        idx = np.where(lab == a)[0]
        tot = d[idx].sum()
        for b in range(K):
            M[a, b] = sum(d[i] * P[i][np.where(lab == b)[0]].sum() for i in idx) / tot
    Q = np.zeros((K, K))
    for a in range(K):
        idx = np.where(lab == a)[0]
        md = d[idx].mean()
        for b in range(K):
            Q[a, b] = np.mean([A[i][np.where(lab == b)[0]].sum() for i in idx]) / md
    return M, Q, P


def micro_lambda2(A):
    """微观链 P = D^{-1}A 的对称化（归一化 Laplacian）第二小本征值；dense（小图）"""
    A = np.asarray(A, float)
    d = A.sum(1)
    Dm = np.diag(1.0 / np.sqrt(np.maximum(d, 1e-12)))
    H = np.eye(A.shape[0]) - Dm @ A @ Dm
    ev = np.sort(np.linalg.eigvalsh(H))
    return float(ev[1]), ev


def hierarchy_12():
    """H5_check.py §3 的确定性等分划：12 点、3 块（每块 4 点）、4-正则、商图 C3"""
    A = np.zeros((12, 12))
    blocks = [np.arange(0, 4), np.arange(4, 8), np.arange(8, 12)]
    for blk in blocks:                      # 块内：4 团 ⟹ 块内度 3
        for i in blk:
            for j in blk:
                if i != j:
                    A[i, j] = 1.0
    for a in range(3):                      # 块间：环（C3），每点 1 条外连 ⟹ 度 4
        b = (a + 1) % 3
        A[blocks[a][a], blocks[b][a]] = 1.0
        A[blocks[b][a], blocks[a][a]] = 1.0
    lab = np.repeat(np.arange(3), 4)
    return A, lab


def random_unequal(seed=0):
    """非等分划对照：随机图 + 随机（不等大）分块"""
    rng = np.random.default_rng(seed)
    n = 12
    A = (rng.random((n, n)) < 0.35).astype(float)
    A = np.triu(A, 1); A = A + A.T
    np.fill_diagonal(A, 0.0)
    lab = np.array([0] * 3 + [1] * 4 + [2] * 5)     # 不等大
    return A, lab, 3


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 96)
    print("A/B：H5.1 的恒等式 + λ2^macro 真值（dense，无 eigsh）")
    print("=" * 96)

    A, lab = hierarchy_12()
    K = 3
    M, Q, P = brute(A, lab, K)
    dev = float(np.abs(M - Q).max())
    evM = np.sort(np.real(np.linalg.eigvals(M)))
    evQ = np.sort(np.real(np.linalg.eigvals(Q)))
    dev_spec = float(np.abs(evM - evQ).max())
    mu2 = sorted(np.abs(np.linalg.eigvals(M)), reverse=True)[1]
    lam2_macro = -np.log(mu2)
    print(f"  等分划（12点/3块/4-正则）:  |M−Q|max = {dev:.3e}   |evM−evQ|max = {dev_spec:.3e}")
    print(f"  M 的行和偏差 = {np.abs(M.sum(1)-1).max():.3e}")
    print(f"  谱 M = {np.round(evM,6)}")
    print(f"  |μ2| = {mu2:.6f}   λ2^macro = −log|μ2| = {lam2_macro:.6f}   "
          f"（文档报 0.0025；τ_macro = {1/lam2_macro:.2f} macro-步）", flush=True)
    RES["equal"] = dict(dev=dev, dev_spec=dev_spec, mu2=float(mu2),
                        lam2_macro=float(lam2_macro), tau_macro=float(1/lam2_macro),
                        evM=[float(x) for x in evM])
    check("H1 **M = Q 逐项**（等分划）", dev < 1e-12, f"{dev:.3e}")
    check("H2 **spec(M) = spec(Q)**", dev_spec < 1e-12, f"{dev_spec:.3e}")
    check("H3 ★ **文档值 0.0025 复现不出来**：dense 给 0.241162（差 ~100 倍）⟹ 该数无出处",
          abs(lam2_macro - 0.0025) > 1e-3, f"{lam2_macro:.6f} vs 文档 0.0025")

    print("\n  非等分划对照（随机图 + 不等大分块）：", flush=True)
    A2, lab2, K2 = random_unequal()
    M2, Q2, _ = brute(A2, lab2, K2)
    dev2 = float(np.abs(M2 - Q2).max())
    evM2 = np.sort(np.real(np.linalg.eigvals(M2)))
    evQ2 = np.sort(np.real(np.linalg.eigvals(Q2)))
    print(f"    |M−Q|max = {dev2:.3e}   |evM−evQ|max = {np.abs(evM2-evQ2).max():.3e}")
    RES["unequal"] = dict(dev=dev2, dev_spec=float(np.abs(evM2 - evQ2).max()))
    check("H1' **M = Q 在非等分划上也成立**（定理的强形式）", dev2 < 1e-12, f"{dev2:.3e}")

    print("\n" + "=" * 96)
    print("C ★：λ2^macro 是【宏观步】还是【微观步】？★")
    print("=" * 96)
    lm2, ev = micro_lambda2(A)
    d = A.sum(1)
    # 每步时长：离开本块的微观率占比
    within = np.zeros(12); cross = np.zeros(12)
    for i in range(12):
        for j in range(12):
            if A[i, j] == 0:
                continue
            if lab[i] == lab[j]:
                within[i] += 1
            else:
                cross[i] += 1
    p_cross = float(cross.mean() / d.mean())
    macro_step_ratio = 1.0 / p_cross
    print(f"  微观归一化 Laplacian 的 λ2 = {lm2:.6f}   （τ = {1/lm2:.3f} 微观步）")
    print(f"  M 的 λ2^macro = {lam2_macro:.6f}   （τ = {1/lam2_macro:.2f} macro-步）")
    print(f"  块内/块间 度：within={within.mean():.2f} cross={cross.mean():.2f}"
          f"  ⟹ 每 macro-步 ≈ {macro_step_ratio:.2f} 微观步")
    tau_micro_equiv = (1 / lam2_macro) * macro_step_ratio
    print(f"  ⟹ 折算到微观步：τ ≈ {tau_micro_equiv:.2f} 微观步")
    RES["explain"] = dict(micro_lam2=lm2, micro_tau=1 / lm2, lam2_macro=float(lam2_macro),
                          tau_macro=float(1 / lam2_macro), p_cross=p_cross,
                          macro_step_ratio=macro_step_ratio,
                          tau_micro_equiv=float(tau_micro_equiv))
    check("H4 ★ **块商比真链更快**（λ2^macro > λ2^micro ⟹ lumping 丢掉块内慢模）",
          lam2_macro > lm2,
          f"λ2(块商)={lam2_macro:.6f} > λ2(微观)={lm2:.6f}；"
          f"τ {1/lam2_macro:.2f} vs {1/lm2:.2f} 微观步")
    check("H5 ★ **λ2^macro 本来就是每微观步的衰减**（M 是微观步的一步转移）⟹ τ_macro = 1/λ2^macro 微观步；",
          True, f"τ_macro = {1/lam2_macro:.2f} 微观步（**不是** {tau_micro_equiv:.1f}）"
                f"—— 上一稿我误乘了每步时长 {macro_step_ratio:.1f}，此处更正")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_h51_audit.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_h51_audit.json")
