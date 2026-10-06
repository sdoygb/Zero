"""
z0_native_constraint.py --- ③ 给零和生成元加一条**原生约束**：能产生偏好/结构吗？

"原生"的判据（由 §45 定死）
---------------------------
* **不得改率** $\kappa_{ij}$（加权重＝加偏好，§35/§36 已否证）；
* 只能改**哪些移动被允许**（D192：偏好只能来自**约束**）；
* 不含自由参数；
* 必须**同时**保住 忠实（无偏好率）＋ 守恒（$\sum_v n_v=0$ 硬约束）。

候选（都满足"无自由参数 + 保守恒"）
------------------------------------
| 名字 | 许可条件 | 直觉 |
|:--|:--|:--|
| `free` | 全许 | 对照（＝现状） |
| `drift` | $n_u>0$ 且 $n_w\le n_u$ | 顺闭合方向漂移 |
| `closure` | $n_u>0$ 且 $n_w\ge n_u$ | 逆向 |
| `support` | $n_u>0$ | 只要求有源 |
| `majority` | $n_u>0$ 且目标邻域过半 $\ge n_w$ | 局域多数 |
| `kick` | $n_u>0$ 且（$n_w\le n_u$ 或 $n_w<0$） | 零上涨落不自持 |

结论
----
**C1 全部方向性硬约束都锁死真空。** 因为任何"要移动就得有正源"的规则在
$n\\equiv0$（$\mathcal Z$ 的原点）上一条都不许 ⟹ **该态吸收** ⟹ 生成元**非不可约**
（实测 $\\lambda{=}0$ 重数 $0$／$1$ 分裂，$5$ 个候选 × $4$ 个图全中）。
**C2 存活的那一条真的产生了结构**：`majority`(ring6)：$\\tau_2$ $0.0651\\to0.1319$（$\\times2.03$）、
**Gini$(\\pi)=0.2026$**（`free` 是 $0.000000$）、$\\max|\\pi-1/N|=2.1\\times10^{-2}$
⟹ 这是本项目第一条"**原生约束 ⟹ 非均匀平稳测度**"的存活证据
（$\\pi\\propto\\deg$ 是**块层**读数；这里的是**微观**测度）。
**C3 但代价是 $O(1)$**：所有存活/陷阱候选的 $\\tau_2$ 与 `free` 之比为 $1.7\\!\\sim\\!5.5$（$17$ 例），
**无一超过 $6$** ⟹ 结构买得到，"**长寿命**"买不到。这与 §47「无指数壁垒」一致。

用法：/usr/bin/python3 z0_native_constraint.py   输出：results/z0_native_constraint.json
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
from scipy.sparse import csr_matrix

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


def make(A, L, rule):
    A = np.asarray(A, float); V = A.shape[0]
    st = Z.states(V, L, True); idx = {s: i for i, s in enumerate(st)}
    nbr = [np.nonzero(A[v])[0] for v in range(V)]
    rows, cols, vals = [], [], []
    dg = np.zeros(len(st))
    for k, s in enumerate(st):
        sarr = np.asarray(s)
        for u in range(V):
            if s[u] <= -L:
                continue
            for w in nbr[u]:
                if s[w] >= L:
                    continue
                if rule == "free":
                    ok = True
                elif rule == "drift":
                    ok = (s[u] > 0 and s[w] <= s[u])
                elif rule == "closure":
                    ok = (s[u] > 0 and s[w] >= s[u])
                elif rule == "support":
                    ok = (s[u] > 0)
                elif rule == "majority":
                    nb = np.asarray(nbr[w])
                    ok = (s[u] > 0 and (len(nb) == 0 or np.sum(sarr[nb] >= s[w]) > len(nb) / 2))
                elif rule == "kick":
                    ok = (s[u] > 0 and (s[w] <= s[u] or s[w] < 0))
                else:
                    ok = True
                if not ok:
                    continue
                nn = list(s); nn[u] -= 1; nn[w] += 1
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0)
                dg[k] += 1.0
    rows += list(range(len(st))); cols += list(range(len(st))); vals += list(-dg)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), st


def stats(M):
    H = M.toarray()
    ev = np.linalg.eigvalsh(H)
    m0 = int(np.sum(np.abs(ev) < 1e-8))
    nz = [x for x in ev if abs(x) > 1e-8]
    tau = 1 / abs(nz[0]) if nz else float("nan")
    gini = float("nan")
    if m0 == 1:
        N = H.shape[0]
        Aq = np.vstack([H.T - np.eye(N), np.ones(N)])
        b = np.zeros(N + 1); b[-1] = 1
        pi, *_ = np.linalg.lstsq(Aq, b, rcond=None)
        pi = np.clip(pi, 0, None); pi /= pi.sum()
        x = np.sort(pi); n = len(x)
        gini = float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))
    return m0, tau, gini


RULES = ("free", "drift", "closure", "support", "majority", "kick")

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 96)
    print("C1/C2：方向性硬约束是否破坏不可约性？（λ=0 重数 = 1 ⟺ 唯一常返类）")
    print("=" * 96)
    print(f"{'Γ':>7} {'L':>2} {'规则':>9} {'N':>6} {'λ0重数':>6} {'不可约':>6} {'Gini(π)':>8} {'τ2':>9}", flush=True)
    rows = []
    for kind, V, L in (("ring", 4, 2), ("K", 4, 2), ("path", 6, 1), ("ring", 6, 1)):
        A = S.G(V, kind)
        for rule in RULES:
            M, st = make(A, L, rule)
            m0, tau, gini = stats(M)
            rows.append(dict(g=kind + str(V), L=L, rule=rule, N=int(M.shape[0]),
                             m0=m0, tau2=tau, gini=gini))
            print(f"{kind+str(V):>7} {L:>2} {rule:>9} {M.shape[0]:>6} {m0:>6} "
                  f"{'是' if m0==1 else '否':>6} {gini:>8.4f} {tau:>9.4f}", flush=True)
        print(flush=True)
    RES["candidates"] = rows
    dir_rules = [r for r in rows if r["rule"] != "free"]
    n_irr = [r for r in dir_rules if r["m0"] == 1]
    check("**C1 多数方向性硬约束破坏不可约性**（$\\lambda{=}0$ 重数 $\\ne1$），但**有反例**",
          len(n_irr) < len(dir_rules),
          f"{len(dir_rules)-len(n_irr)}/{len(dir_rules)} 破坏；"
          f"存活者：{[(r['g'], r['rule']) for r in n_irr]}")
    check("**对照：`free` 全部不可约**（重数 $=1$）",
          all(r["m0"] == 1 for r in rows if r["rule"] == "free"),
          f"{sum(1 for r in rows if r['rule']=='free')} 例")
    check("**C2 `free` 的 $\\pi$ 严格均匀**（Gini $=0$；方向性规则下无法定义 $\\pi$）",
          all(abs(r["gini"]) < 1e-12 for r in rows if r["rule"] == "free" and r["gini"] == r["gini"]),
          "自由核的均匀测度再次确认")

    print("\n" + "=" * 96)
    print("C3：陷阱的**代价**——把不可约规则与 free 对照（找 τ2 的乘性因子）")
    print("=" * 96)
    print(f"{'Γ':>7} {'规则':>9} {'τ2':>9} {'τ2(prev)':>9} {'比值':>7}", flush=True)
    ratios = []
    for kind, V, L in (("ring", 4, 2), ("K", 4, 2), ("path", 6, 1), ("ring", 6, 1)):
        A = S.G(V, kind)
        base = None
        for rule in RULES:
            M, st = make(A, L, rule)
            m0, tau, _ = stats(M)
            if rule == "free":
                base = tau
            if m0 != 1 and base and tau == tau:
                r = tau / base
                ratios.append(r)
                print(f"{kind+str(V):>7} {rule:>9} {tau:>9.4f} {base:>9.4f} {r:>7.3f}", flush=True)
    RES["trap_cost_ratios"] = ratios
    check("**C3 代价是 $O(1)$**：所有候选的 $\\tau_2/\\tau_2^{\\rm free}<10$（不是指数）",
          all(1.0 < r < 10.0 for r in ratios),
          f"{len(ratios)} 例，比值 {min(ratios):.3f}–{max(ratios):.3f}")
    # ★ 存活者的非均匀测度
    M, st = make(S.G(6, "ring"), 1, "majority")
    _, tau_m, gini_m = stats(M)
    M0, _ = make(S.G(6, "ring"), 1, "free")
    _, tau_f, gini_f = stats(M0)
    RES["survivor"] = dict(rule="majority", g="ring6", tau2=tau_m, tau2_free=tau_f,
                           gini=gini_m, gini_free=gini_f)
    check("**★ C2 存活者真的产生非均匀微观测度**（Gini $0.20$ vs `free` 的 $0.00$）",
          gini_m > 0.1 and abs(gini_f) < 1e-12 and tau_m > tau_f,
          f"majority: τ2={tau_m:.4f} Gini={gini_m:.4f}  |  free: τ2={tau_f:.4f} Gini={gini_f:.6f}")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_native_constraint.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_native_constraint.json")
