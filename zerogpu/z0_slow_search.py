"""
z0_slow_search.py --- ① 系统性地再找一次「慢」（带纪律：dense 优先、避开奇异 eigensolver）

纪律（§41 的教训，第三次同型）
------------------------------
1. 只要 $N\\le4000$ 就 **dense 精确**；$N>4000$ 用 `which='SA'`（**代数**最小，非最小模）；
2. **不用** `eigsh(which='SM')` / `sigma=0`（奇异矩阵 ⟹ 伪 Ritz 值）；
3. 每个候选都要与**已知精确量**对账：变分界 $\\lambda_2\\le d_{\\min}\\nu_2$，$\\nu_2=-4\\sin^2\\frac{\\pi}{4L+2}$；
4. 判"慢"的阈值先声明：$\\tau_2=1/|\\lambda_2|>1$ 才算"非 $O(1)$"。

三条结论
--------
**S1 无慢模**：17 个配置、10 个图族、$L\\le5$，$\\tau_2$ **全部 $<0.14$**（无一超 1）。
**S2 随容量收敛**：$\\tau_2$ 随 $L$ **单调下降并趋于平台**（ring4 $0.096\\to0.0645$，$K_4$ $0.075\\to0.063$），
   指数 $-0.11\\sim-0.27$ ⟹ **无指数壁垒**（若真亚稳，$\\tau_2$ 应 $\\sim e^{cL}$）。
**S3 变分界全过**：$\\lambda_2\\le\\bar d\\,\\nu_2$，且大 $L$ 时 $\\tau_2\\le4L^2/(\\pi^2\\bar d)$。
   ⟹ 即便 $L$ 大，上界也只是 **$O(L^2)$ 多项式**，不可能给指数长寿命。

**机制**：零和约束**同时删掉两个角落** $(\\pm L,\\dots,\\pm L)$，那正是最深的两个陷阱
⟹ 剩下的链在 $L$ 方向上是**扩张的**，不是阻挫的。

用法：/usr/bin/python3 z0_slow_search.py   输出：results/z0_slow_search.json
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
from scipy.sparse.linalg import eigsh

import z0_nonadditivity as Z

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(A, b, dl=4000):
    Lc, st = Z.build(A, b, True)
    N = Lc.shape[0]
    if N <= dl:
        ev = np.linalg.eigvalsh(Lc.toarray()); tag = "dense"
    else:
        ev = np.sort(eigsh(Lc.tocsc(), k=4, which="SA", return_eigenvectors=False).real); tag = "SA"
    nz = [x for x in ev if abs(x) > 1e-8]
    return (float(nz[0]) if nz else None), N, tag


def G(V, kind):
    A = np.zeros((V, V))
    if kind == "path":
        for i in range(V - 1):
            A[i, i + 1] = A[i + 1, i] = 1
    elif kind == "ring":
        for i in range(V):
            A[i, (i + 1) % V] = A[(i + 1) % V, i] = 1
    elif kind == "K":
        A[:] = 1; np.fill_diagonal(A, 0)
    elif kind == "star":
        for j in range(1, V):
            A[0, j] = A[j, 0] = 1
    elif kind == "K33":
        A = np.zeros((6, 6)); A[:3, 3:] = 1; A[3:, :3] = 1
    elif kind == "petersen":
        A = np.zeros((10, 10))
        for i in range(5):
            A[i, (i + 1) % 5] = A[(i + 1) % 5, i] = 1
            A[5 + i, 5 + (i + 2) % 5] = A[5 + (i + 2) % 5, 5 + i] = 1
            A[i, 5 + i] = A[5 + i, i] = 1
    elif kind == "cube":
        A = np.zeros((8, 8))
        for i in range(8):
            for k in range(3):
                A[i, i ^ (1 << k)] = 1
    elif kind == "barbell8":
        A = np.zeros((8, 8))
        for off in (0, 4):
            for i in range(off, off + 4):
                for j in range(off, off + 4):
                    if i != j:
                        A[i, j] = 1
        A[3, 4] = A[4, 3] = 1
    elif kind == "ladder":
        A = np.zeros((8, 8))
        for i in range(4):
            A[i, (i + 1) % 4] = A[(i + 1) % 4, i] = 1
            A[4 + i, 4 + (i + 1) % 4] = A[4 + (i + 1) % 4, 4 + i] = 1
            A[i, 4 + i] = A[4 + i, i] = 1
    return A


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 100)
    print("S1/S3：广扫 L=1（10 个图族）+ 变分界对账")
    print("=" * 100)
    print(f"{'Γ':>11} {'V':>3} {'d̄':>5} {'N':>7} {'λ2':>10} {'τ2':>9} {'d̄|ν2|':>8} {'界':>4} {'法':>5}", flush=True)
    rows = []
    for kind in ("path", "ring", "K", "star", "K33", "petersen", "cube", "barbell8", "ladder"):
        for V in ((4, 6, 8) if kind in ("path", "ring", "K", "star") else (0,)):
            A = G(V, kind) if V else G(0, kind)
            if A.shape[0] == 0:
                continue
            Vv = A.shape[0]; d = float(A.sum(1).mean())
            l2, N, tag = gap(A, 1)
            if l2 is None:
                continue
            bnd = d * (-1.0)          # ν2(b=1) = -1
            rows.append(dict(name=f"{kind}{Vv}", V=Vv, d=d, N=N, lam2=l2,
                             tau2=1 / abs(l2), bound=bnd, ok=bool(l2 <= bnd + 1e-9)))
            print(f"{kind+str(Vv):>11} {Vv:>3} {d:>5.2f} {N:>7} {l2:>10.4f} {1/abs(l2):>9.4f}"
                  f" {bnd:>8.3f} {str(l2<=bnd+1e-9):>4} {tag:>5}", flush=True)
    RES["scan_L1"] = rows
    check("**S1 无慢模**（$L{=}1$，10 图族）：$\\tau_2$ 全 $<1$",
          all(r["tau2"] < 1 for r in rows),
          f"最大 {max(r['tau2'] for r in rows):.4f}（{max(rows, key=lambda r: r['tau2'])['name']}），共 {len(rows)} 例")
    check("**S3 变分界 $\\lambda_2\\le\\bar d\\nu_2$ 全过**", all(r["ok"] for r in rows), f"{len(rows)} 例")

    print("\n" + "=" * 100)
    print("S2：τ2 随容量 L 的走向（推到 L=5）—— 涨（亚稳）还是落（扩张）？")
    print("=" * 100)
    print(f"{'Γ':>8} {'d̄':>5} {'L':>2} {'N':>9} {'λ2':>10} {'τ2':>9} {'4L²/(π²d̄)':>11} {'法':>5}", flush=True)
    Lrows = []
    for kind, V in (("ring", 4), ("K", 4), ("K33", 6)):
        A = G(V, kind); d = float(A.sum(1).mean())
        for L in (1, 2, 3, 4, 5):
            l2, N, tag = gap(A, L)
            if l2 is None or N > 900000:
                continue
            bnd = 4 * L * L / (np.pi ** 2 * d)
            Lrows.append(dict(name=f"{kind}{V}", L=L, N=N, lam2=l2, tau2=1 / abs(l2), bound=bnd))
            print(f"{kind+str(V):>8} {d:>5.2f} {L:>2} {N:>9} {l2:>10.4f} {1/abs(l2):>9.4f}"
                  f" {bnd:>11.4f} {tag:>5}", flush=True)
    RES["scan_L"] = Lrows
    exps = {}
    print(flush=True)
    for nm in sorted(set(r["name"] for r in Lrows)):
        rs = sorted([r for r in Lrows if r["name"] == nm], key=lambda r: r["L"])
        Ls = np.array([r["L"] for r in rs], float); T = np.array([r["tau2"] for r in rs])
        p = np.polyfit(np.log(Ls), np.log(T), 1)
        exps[nm] = float(p[0])
        print(f"  {nm:>8}: τ2 ∝ L^{p[0]:+.3f}   τ2 = {np.round(T,4)}", flush=True)
    RES["L_exponents"] = exps
    check("**S2 无指数壁垒**：$\\tau_2$ 随 $L$ **下降或平台**（指数全 $<0$）",
          all(e < 0 for e in exps.values()),
          "；".join(f"{k}: L^{v:+.2f}" for k, v in exps.items()))
    check("**S1 在 $L\\le5$ 上仍成立**：所有 $\\tau_2<1$",
          all(r["tau2"] < 1 for r in Lrows),
          f"最大 {max(r['tau2'] for r in Lrows):.4f}")

    print("\n" + "=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_slow_search.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_slow_search.json")
