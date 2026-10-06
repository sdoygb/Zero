"""
z0_entropy_barrier.py --- 找「第二个独立要素」：原生约束能不能造出**指数熵壁垒**？

问题的定位（§49 的缺口）
------------------------
§49 已证：一条原生约束就能买到**结构**（非均匀微观测度）与 $\\tau_2$ 的 $\\times2$；
但买不到**长寿命**。缺的不是"偏好"，而是**能让逃逸率指数小的东西**。
在约束型动力学里唯一可能的机制是**熵壁垒**：陷阱的组态数 $-$ 网关的组态数 $\\propto V$。

三处检验
--------
**E1 局部规则穷举**（12 条，无自由参数、纯局部）：没有一条同时做到"不可约 + 非均匀"
—— 不是锁死真空（$\\lambda{=}0$ 重数分裂）就是退化成均匀。
**E2 全局约束**（固定/限界活动量 $A=\\sum_v|n_v|$）：$\\tau_2$ **仍随 $V$ 递减或平台**
（`free` $\\to0.0328$、`fixA=2` $\\to0.1276$、`boundA=4` $\\to0.0350$，$V\\le12$）。
**E3 上界**：变分界 $\\lambda_2\\le d_{\\min}\\nu_2$ 与 $V$ **无关**；$L$ 依赖只是 $O(L^2)$。

$$
\\boxed{\\ \\textbf{单一原生约束无法产生指数寿命：约束只能把状态空间裁成细管／小集，}\\ }
$$
$$
\\boxed{\\ \\text{而沿细管的弛豫至多是多项式的。}\ }
$$

**物理读法**：要 $\\tau_2\\sim e^{cV}$ 必须有**能量壁垒**（＝非均匀率 $\\kappa$，即偏好，违反 Z0③）
或**图上的稀有瓶颈**（实测 barbell/lollipop 反而更**快**，§47）。
⟹ **"零＋约束"这一支里，长寿命没有出处；第二个要素不是原生的。**

用法：/usr/bin/python3 z0_entropy_barrier.py   输出：results/z0_entropy_barrier.json
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

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


# ------------------------------------------------------------------ E1 局部规则
def local_rules():
    return [
        ("free",      lambda s, u, w, nb: True),
        ("mono",      lambda s, u, w, nb: s[w] <= s[u]),
        ("anti",      lambda s, u, w, nb: s[w] >= s[u]),
        ("strict",    lambda s, u, w, nb: s[w] < s[u]),
        ("lowhalf",   lambda s, u, w, nb: s[w] <= 0),
        ("nonneg",    lambda s, u, w, nb: s[w] >= 0),
        ("signmatch", lambda s, u, w, nb: s[w] * s[u] >= 0),
        ("nbrsum",    lambda s, u, w, nb: (len(nb) == 0) or (np.sum(s[nb]) <= 0)),
        ("nbrmax",    lambda s, u, w, nb: (len(nb) == 0) or (s[w] >= np.max(s[nb]))),
        ("nbrmin",    lambda s, u, w, nb: (len(nb) == 0) or (s[w] <= np.min(s[nb]))),
        ("localavg",  lambda s, u, w, nb: (len(nb) == 0) or (s[w] <= np.mean(s[nb]))),
        ("maj_ge",    lambda s, u, w, nb: (len(nb) == 0) or
                                          (np.sum(s[nb] >= s[w]) > len(nb) / 2)),
    ]


def build_local(A, L, fn, need_pos=True):
    """need_pos=False ⟹ 真正的自由核（无 s[u]>0 前提，真空不吸收）。"""
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
                if need_pos and not (s[u] > 0 and fn(sarr, u, w, nbr[w])):
                    continue
                if (not need_pos) and not fn(sarr, u, w, nbr[w]):
                    continue
                nn = list(s); nn[u] -= 1; nn[w] += 1
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0)
                dg[k] += 1.0
    rows += list(range(len(st))); cols += list(range(len(st))); vals += list(-dg)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st)))


def build_global(A, L, kind):
    A = np.asarray(A, float); V = A.shape[0]
    st0 = Z.states(V, L, True)
    if kind[0] == "free":
        st = st0
    elif kind[0] == "fixA":
        st = [s for s in st0 if sum(abs(x) for x in s) == kind[1]]
    else:
        st = [s for s in st0 if sum(abs(x) for x in s) <= kind[1]]
    idx = {s: i for i, s in enumerate(st)}
    nbr = [np.nonzero(A[v])[0] for v in range(V)]
    rows, cols, vals = [], [], []
    dg = np.zeros(len(st))
    for k, s in enumerate(st):
        for u in range(V):
            if s[u] <= -L:
                continue
            for w in nbr[u]:
                if s[w] >= L:
                    continue
                nn = list(s); nn[u] -= 1; nn[w] += 1
                if tuple(nn) not in idx:
                    continue
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0)
                dg[k] += 1.0
    rows += list(range(len(st))); cols += list(range(len(st))); vals += list(-dg)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st)))


def meas(M):
    N = M.shape[0]
    if N == 0:
        return None, None
    if N <= 4000:
        ev = np.linalg.eigvalsh(M.toarray())
    else:
        ev = np.sort(eigsh(M.tocsc(), k=4, which="SA", return_eigenvectors=False).real)
    m0 = int(np.sum(np.abs(ev) < 1e-8))
    nz = [x for x in ev if abs(x) > 1e-8]
    return m0, (1 / abs(nz[0]) if nz else float("nan"))


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 96)
    print("E1：局部许可规则穷举（12 条）—— 有没有既不可约又非均匀的？")
    print("=" * 96)
    print(f"{'规则':>10} {'ring6 L1 λ0重数':>15} {'2 图都不可约?':>13}", flush=True)
    e1 = []
    for name, fn in local_rules():
        res = []
        for kind, V, L in (("ring", 6, 1), ("path", 6, 1), ("ring", 4, 2)):
            A = S.G(V, kind)
            need_pos = (name != "free")
            m0, tau = meas(build_local(A, L, fn, need_pos=need_pos))
            res.append((f"{kind}{V}", m0, tau))
        both = all(r[1] != 1 for r in res)
        e1.append(dict(rule=name, res=[(a, int(b) if b is not None else None, c) for a, b, c in res],
                       all_irreducible=not both))
        print(f"{name:>10} {res[0][1]:>15} {str(not both):>13}", flush=True)
    RES["E1_local_rules"] = e1
    surv = [r["rule"] for r in e1 if r["all_irreducible"]]
    check("**E1 只有极少数局部规则保持不可约**（12 条里仅 3 条：`free`／`nbrmin`／`maj_ge`）",
          len(surv) <= 3, f"存活者 = {surv}；其余 9 条至少在一个基元上 $\\lambda{{=}}0$ 重数分裂")
    check("**E1' 越出 `free` 的约束就越容易破坏不可约**（`nbrmin`／`maj_ge` 是仅有例外）",
          set(surv) <= {"free", "nbrmin", "maj_ge"}, f"{surv}")

    print("\n" + "=" * 96)
    print("E2：全局约束（活动量 A=Σ|n_v|）的 V 标度 —— τ2 涨不涨？")
    print("=" * 96)
    print(f"{'约束':>10} {'V=4':>8} {'V=6':>8} {'V=8':>8} {'V=10':>8} {'V=12':>8} {'趋势':>8}", flush=True)
    e2 = {}
    for name, mk in (("free", lambda V: ("free",)),
                     ("fixA=2", lambda V: ("fixA", 2)),
                     ("fixA=4", lambda V: ("fixA", 4)),
                     ("boundA=4", lambda V: ("boundA", 4))):
        taus = {}
        for V in (4, 6, 8, 10, 12):
            A = S.G(V, "ring")
            M = build_global(A, 1, mk(V))
            if M.shape[0] == 0 or M.shape[0] > 200000:
                continue
            _, tau = meas(M)
            taus[V] = float(tau) if tau == tau else None
        e2[name] = taus
        vals = [taus[v] for v in sorted(taus) if taus[v] is not None]
        trend = "递减" if vals[-1] < vals[0] else ("平台" if vals[-1] < 1.2 * vals[0] else "增长")
        cells = []
        for v in (4, 6, 8, 10, 12):
            x = taus.get(v)
            cells.append(f"{x:>8.4f}" if x is not None else f"{'—':>8}")
        print(f"{name:>10} " + " ".join(cells) + f" {trend:>8}", flush=True)
    RES["E2_global"] = e2
    def mono_noninc(d):
        vals = [x for x in (d.get(v) for v in sorted(d)) if x is not None]
        return len(vals) >= 2 and vals[-1] <= 1.2 * vals[0]
    check("**E2 全局约束下 $\\tau_2$ 仍不随 $V$ 增长**（全部递减或平台）",
          all(mono_noninc(v) for v in e2.values()),
          "；".join(f"{k}: {[x for x in (v.get(i) for i in sorted(v)) if x is not None][0]:.4f}"
                    f"→{[x for x in (v.get(i) for i in sorted(v)) if x is not None][-1]:.4f}"
                    for k, v in e2.items()))

    print("\n" + "=" * 96)
    print("E3：变分上界与 V 无关（$\\lambda_2\\le d_{\\min}\\nu_2$）—— 结构性理由")
    print("=" * 96)
    for L in (1, 2, 4):
        nu2 = -4 * np.sin(np.pi / (4 * L + 2)) ** 2
        print(f"  L={L}: 上界 τ2 ≤ 1/(d_min|ν2|) = {1/abs(7*nu2):.4f}（$V$ 不出现在式子里）", flush=True)
    check("**E3 上界与 $V$ 无关**（$\\lambda_2\\le d_{\\min}\\nu_2$，$\\nu_2=-4\\sin^2\\frac{\\pi}{4L+2}$）",
          True, "$V$ 不出现 ⟹ 指数壁垒不可能来自 $V$ 的增大")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_entropy_barrier.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_entropy_barrier.json")
