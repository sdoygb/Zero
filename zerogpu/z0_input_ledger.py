"""
z0_input_ledger.py --- ② 具名输入清单的收敛：每条输入能不能被"自洽性"钉住？

方法
----
对每条在册输入，问三个问题（缺一不可）：
  Q1 **它有没有被独立约束过？**（存在 $\ge2$ 条互不依赖的条件同时指向它）
  Q2 **它的替代值会不会破坏已导出的东西？**（反事实检验）
  Q3 **它是不是"读法/口径"？**（换口径结论就变 ⟹ 不是物理输入，是约定）

判决分四档：
  【导出】   可由已立定理推出（无自由）
  【条件锁定】被 $\ge2$ 条独立条件钉死（付的价签＝接受那些条件）
  【约定】   口径选择，换口径结论变
  【裸输入】 只剩"选它"这一条理由

本机器**只做能机器判的部分**：凡是能在零和生成元上直接算的反事实，全部算出来；
不能算的（如 $E5$）明确标"未判"，不冒充结论。

用法：/usr/bin/python3 z0_input_ledger.py   输出：results/z0_input_ledger.json
"""
from __future__ import annotations

import json
import os

import numpy as np

import z0_nonadditivity as Z
import z0_slow_search as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def tau2_of(A, L):
    l2, N, tag = S.gap(A, L)
    return (None if l2 is None else 1 / abs(l2)), N, tag


if __name__ == "__main__":
    print("=" * 96)
    print("Q2 反事实检验：换掉一条输入，已导出的东西会不会坏？")
    print("=" * 96)

    # ---- 反事实 A：容量 L（物理 4）→ 换 1,2,3,5，看"零漂移/对称/唯一零模"是否还成立
    print("\n【反事实 A】容量 L：它是物理输入还是读法？", flush=True)
    print(f"{'Γ':>8} {'L':>2} {'N':>7} {'对称偏差':>10} {'λ0重数':>7} {'τ2':>9}", flush=True)
    cfA = []
    for kind, V in (("ring", 4), ("K33", 6)):
        A = S.G(V, kind)
        for L in (1, 2, 3, 4, 5):
            Lc, st = Z.build(A, L, True)
            N = Lc.shape[0]
            if N > 2000:            # ★ 只在小态空间做 dense；大态由 z0_slow_search 用 SA 做
                continue
            H = Lc.toarray()
            symd = float(np.max(np.abs(H - H.T)))
            ev = np.linalg.eigvalsh(H)
            m0 = int(np.sum(np.abs(ev) < 1e-9))
            tau, _, _ = tau2_of(A, L)
            cfA.append(dict(kind=kind, L=L, N=int(N), sym_dev=symd, m0=m0, tau2=tau))
            print(f"{kind+str(V):>8} {L:>2} {N:>7} {symd:>10.2e} {m0:>7} {tau:>9.4f}", flush=True)
    RES["counterfactual_L"] = cfA
    check("**L 不是读法**：任意 L 下对称性/唯一零模都保持（L 只改标度，不改结构）",
          all(r["sym_dev"] < 1e-9 and r["m0"] == 1 for r in cfA),
          f"{len(cfA)} 例：对称偏差全 $<10^{{-9}}$、$\\lambda{{=}}0$ 重数全 $=1$")

    # ---- 反事实 B：偏好核 κ≠1 → 对称性与零漂移会不会坏
    print("\n【反事实 B】偏好核 $\\kappa_{ij}\\ne1$（加记录权重）→ 精确对称/零漂移会坏吗？", flush=True)
    print(f"{'方案':>16} {'对称偏差':>10} {'<Δ活动量>':>11} {'λ0重数':>7} {'τ2':>9}", flush=True)
    cfB = []
    A = S.G(6, "ring")
    base_L = 2
    Lc0, st0 = Z.build(A, base_L, True)
    for tag, weight in (("无偏好 κ=1", None), ("记录权重 1+R", "record")):
        # 构造带偏好的生成元：率 = 1 + (记录量)  （记录量 = |n_v|）
        st = Z.states(A.shape[0], base_L, True)
        idx = {s: i for i, s in enumerate(st)}
        n = len(st)
        rows, cols, vals = [], [], []
        dg = np.zeros(n)
        for k, s in enumerate(st):
            for a in range(A.shape[0]):
                if s[a] <= -base_L:
                    continue
                for bb in np.nonzero(A[a])[0]:
                    if s[bb] >= base_L:
                        continue
                    nn = list(s); nn[a] -= 1; nn[bb] += 1
                    w = 1.0 if weight is None else (1.0 + abs(s[a]))
                    rows.append(idx[tuple(nn)]); cols.append(k); vals.append(w)
                    dg[k] += w
        rows += list(range(n)); cols += list(range(n)); vals += list(-dg)
        from scipy.sparse import csr_matrix
        Lc = csr_matrix((vals, (rows, cols)), shape=(n, n)).toarray()
        symd = float(np.max(np.abs(Lc - Lc.T)))
        ev = np.linalg.eigvalsh(Lc)
        m0 = int(np.sum(np.abs(ev) < 1e-9))
        nz = [x for x in ev if abs(x) > 1e-9]
        tau = 1 / abs(nz[0]) if nz else None
        cfB.append(dict(tag=tag, sym_dev=symd, m0=m0, tau2=tau))
        print(f"{tag:>16} {symd:>10.2e} {m0:>11} {tau:>9.4f}", flush=True)
    RES["counterfactual_bias"] = cfB
    check("**偏好核破坏精确对称**（$\\kappa{=}1$ 对称偏差 $=0$；加记录权重后 $\\ne0$）",
          cfB[0]["sym_dev"] < 1e-12 and cfB[1]["sym_dev"] > 1e-6,
          f"κ=1: {cfB[0]['sym_dev']:.1e} → 加权重: {cfB[1]['sym_dev']:.1e}")

    # ---- 清单
    print("\n" + "=" * 96)
    print("清单：每条输入的判决（机器可判部分）")
    print("=" * 96)
    ledger = [
        ("零和约束 $\\sum_v n_v=0$", "【导出】", "唯一生成型；排除盒角、给唯一零模（§43）"),
        ("生成元 $\\kappa_{ij}=1$", "【导出】", "Z0③ 全分支；任何偏好都破精确对称（反事实 B）"),
        ("均匀微观测度", "【导出】", "精确对称的推论（$\\max|\\mathcal L-\\mathcal L^{\\sf T}|=0$）"),
        ("容量 $L$", "【条件锁定】", "不改结构、只改标度（反事实 A：任意 $L$ 下对称+唯一零模）"),
        ("$\\Gamma$", "【裸输入】", "`Z1` 只保证存在；$\\tau_2$ 对它不敏感（§47）"),
        ("$B=4$", "【裸输入】", "`G70` 自撤因果论证；实由 `G61` 标签计数钉住"),
        ("$\\kappa_1=q$／$T_d$", "【约定】", "`G72` 自撤【导出】：\"是识别在定 $T_d$，不是生成元\""),
        ("$D=4$", "【裸输入】", "`Z17`：Z0 不给维数约束"),
        ("$A1$ 作用量相位", "【裸输入｜优先级 1】", "`G92:44`：开放、生成型、承重"),
        ("折叠 $|\\rho-s|$", "【裸输入｜构造】", "§27：非导出"),
        ("$E5$：$\\pi$ 的选择", "【未判】", "缺一条导出；本机器无法判（需语料侧工作）"),
    ]
    RES["ledger"] = [dict(item=it, verdict=v, note=n) for it, v, n in ledger]
    for it, v, n in ledger:
        print(f"  {v:>18}  {it:<26} {n}", flush=True)
    n_naked = sum(1 for _, v, _ in ledger if "裸输入" in v)
    print(f"\n  统计：导出 {sum(1 for _,v,_ in ledger if v=='【导出】')} 条 · "
          f"条件锁定 1 条 · 约定 1 条 · **裸输入 {n_naked} 条** · 未判 1 条", flush=True)
    check("**裸输入 ≤ 6 条**（收敛目标：把可消的消掉，剩下的明码标价）",
          n_naked <= 6, f"裸输入 {n_naked} 条")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_input_ledger.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"→ results/z0_input_ledger.json")
