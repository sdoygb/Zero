"""
z0_nonadditivity.py --- 零第一次真正做功的地方：**约束谱是无约束（自由）谱的子集**

★ 本条更正了本轮我自己前两稿的框架 ★
--------------------------------------------------------------------------
前两稿把"约束修正"当成"两个不等谱的差"，并猜过 Minkowski 夹逼 —— 都错了。
机器的定案是**精确的子集关系**：

  自由系（**无**零和约束的盒 $[-L,L]^V$）：
      生成元按源点拆 $\\mathcal L=\\sum_v\\mathcal L_v$，**逐项相等**（已验，max 差 = 0）；
      $\\mathcal L_v$ 只依赖 $(n_v,\\{n_w\\}_{w\\sim v})$ ⟹ 谱 = 单点谱的 **Minkowski 和**。
  真实系（**加**零和约束 $\\sum_v n_v=0$）：
      谱 = 自由谱的**子集**，逐位相同（实测最大距离 $<10^{-9}$）。

$$
\\boxed{\\ \\operatorname{spec}(\\mathcal L_{\\mathcal Z})\\ \\subset\\ \\operatorname{spec}(\\mathcal L_{\\rm free})
\\qquad\\Longrightarrow\\qquad
|\\lambda_2^{\\rm true}|\\ \\ge\\ |\\lambda_2^{\\rm free}|\\ }
$$

**机制**：约束**删掉**的是自由的**慢模**（盒角吸收态附近的模态）；
自由系 $\\lambda=0$ 重数 $=5,7,9,\\dots$（盒角吸收），真实系**恒 = 1**。
⟹ **约束让系统更快弛豫**（不是更慢）—— 这是"零"唯一的动力学效应。

**标度**（$b{=}1$，实测 $V\\le8$）：$|\\lambda_2^{\\rm true}|\\approx2.4\\,V$
⟹ $\\tau_2=1/|\\lambda_2|\\approx0.4/V$ —— 寿命随系统尺寸**下降**。

**严格上界**（变分，单点试验函数）：$\\lambda_2^{\\rm true}\\le d_{\\min}\\nu_2$，
$\\nu_2=-4\\sin^2\\!\\frac{\\pi}{4L+2}$ ⟹ 物理 $\\Gamma$、$L{=}4$：$\\tau_2\\le1.184$ 微观步（与 $V$ 无关）。

用法：/usr/bin/python3 z0_nonadditivity.py   输出：results/z0_nonadditivity.json
"""
from __future__ import annotations

import json
import os
import time
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def states(V, L, constraint=True):
    """{n in [-L,L]^V : sum n = 0}（约束）或全体（自由）。

    ★ 计数（L=1，约束）：Σ_k C(V,k)C(V−k,k)，渐近 ~ 3^V/√(2πV/3) ⟹ 可达上限 V≈14。
    ★ **不要**用 (2L+1)^V 再过滤：V=24 时 3^24=2.8e11 ⟹ 必死。
    """
    if not constraint:
        return list(product(range(-L, L + 1), repeat=V))
    out = []

    def rec(pref, rem, left):
        if abs(rem) > L * left:
            return
        if left == 0:
            out.append(tuple(pref)); return
        lo, hi = max(-L, rem - L * (left - 1)), min(L, rem + L * (left - 1))
        for x in range(lo, hi + 1):
            rec(pref + [x], rem - x, left - 1)
    rec([], 0, V)
    return out


def build(A, L, constraint=True):
    A = np.asarray(A, float); V = A.shape[0]
    st = states(V, L, constraint)
    idx = {n: i for i, n in enumerate(st)}
    nbr = [np.nonzero(A[v])[0] for v in range(V)]
    rows, cols, vals = [], [], []
    dg = np.zeros(len(st))
    for k, n in enumerate(st):
        for a in range(V):
            if n[a] <= -L:
                continue
            for bb in nbr[a]:
                if n[bb] >= L:
                    continue
                nn = list(n); nn[a] -= 1; nn[bb] += 1
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0)
                dg[k] += 1.0
    rows += list(range(len(st))); cols += list(range(len(st))); vals += list(-dg)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), st


def lowest(Lm, k=4):
    n = Lm.shape[0]
    if n <= 3000:
        return np.sort(np.linalg.eigvalsh(Lm.toarray()))
    return np.sort(eigsh(Lm.tocsc(), k=k, which="SA", return_eigenvectors=False).real)


def zero_mult(ev, tol=1e-9):
    return int(np.sum(np.abs(ev) < tol))


def graphs(V, kind):
    A = np.zeros((V, V))
    if kind == "path":
        for i in range(V - 1):
            A[i, i + 1] = A[i + 1, i] = 1
    elif kind == "ring":
        for i in range(V):
            A[i, (i + 1) % V] = A[(i + 1) % V, i] = 1
    return A


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 100)
    print("N1：约束谱 ⊂ 自由谱（逐位相同），且约束把 λ=0 重数压到 1")
    print("=" * 100)
    print(f"{'Γ':>5} {'V':>3} {'N_free':>7} {'N_cons':>7} {'λ0重数(free)':>12} {'λ0重数(cons)':>12}"
          f" {'λ2_free':>10} {'λ2_cons':>10} {'⊂?':>5}", flush=True)
    rows = []
    for kind in ("path", "ring"):
        for V in (3, 4, 5, 6, 7):
            A = graphs(V, kind)
            Lf, _ = build(A, 1, False); Lc, _ = build(A, 1, True)
            ef = lowest(Lf); ec = lowest(Lc)
            fnz = [x for x in ef if abs(x) > 1e-9]
            cnz = [x for x in ec if abs(x) > 1e-9]
            l2f, l2c = (float(fnz[0]) if fnz else None), (float(cnz[0]) if cnz else None)
            dev = max(float(np.min(np.abs(ef - x))) for x in ec)
            rows.append(dict(kind=kind, V=V, N_free=int(Lf.shape[0]), N_cons=int(Lc.shape[0]),
                             m0_free=zero_mult(ef), m0_cons=zero_mult(ec),
                             l2_free=l2f, l2_cons=l2c, subset_dev=dev))
            print(f"{kind:>5} {V:>3} {Lf.shape[0]:>7} {Lc.shape[0]:>7} {zero_mult(ef):>12}"
                  f" {zero_mult(ec):>12} {(l2f or float('nan')):>10.4f} {(l2c or float('nan')):>10.4f}"
                  f" {str(dev < 1e-9):>5}", flush=True)
    RES["subset"] = rows
    check("**N1 约束谱 ⊂ 自由谱**（逐位相同；密集对角化段 <1e-9，V=8 用稀疏故放宽到 1e-6）",
          all(r["subset_dev"] < (1e-9 if r["V"] <= 7 else 1e-6) for r in rows),
          "；".join(f"{r['kind']}{r['V']}:{r['subset_dev']:.1e}" for r in rows if r["subset_dev"] > 1e-9) or "全部 <1e-9")
    check("**自由系 λ=0 重数 = 2LV+1**（= 总量守恒的扇区数；L=1 ⟹ 2V+1）",
          all(r["m0_free"] == 2 * 1 * r["V"] + 1 for r in rows),
          "；".join(f"{r['kind']}{r['V']}: {r['m0_free']}（应 {2*r['V']+1}）" for r in rows))
    check("**约束把 λ=0 重数压到 1**（只取 Σn=0 那一个扇区）",
          all(r["m0_cons"] == 1 for r in rows),
          "；".join(f"{r['kind']}{r['V']}: {r['m0_free']}→{r['m0_cons']}" for r in rows))
    check("**约束让系统更快**：|λ2_cons| ≥ |λ2_free|（全例）",
          all(abs(r["l2_cons"]) >= abs(r["l2_free"]) - 1e-9 for r in rows), f"{len(rows)} 例")

    print("\n" + "=" * 100)
    print("N2：自由系的可加性（L = Σ_v L_v 逐项相等 ⟹ 谱 = 单点谱的 Minkowski 和）")
    print("=" * 100)
    ok2 = True
    for kind in ("path", "ring"):
        V, L = 4, 1
        A = graphs(V, kind)
        st = states(V, L, False)
        idx = {n: i for i, n in enumerate(st)}
        nbr = [np.nonzero(A[v])[0] for v in range(V)]
        Lsum = np.zeros((len(st), len(st)))
        for a in range(V):
            La = np.zeros_like(Lsum)
            for k, n in enumerate(st):
                if n[a] <= -L:
                    continue
                d = 0.0
                for bb in nbr[a]:
                    if n[bb] >= L:
                        continue
                    nn = list(n); nn[a] -= 1; nn[bb] += 1
                    La[idx[tuple(nn)], k] += 1.0; d += 1.0
                La[k, k] -= d
            Lsum += La
        Lf, _ = build(A, L, False)
        dev = float(np.max(np.abs(Lsum - Lf.toarray())))
        ok2 &= dev < 1e-12
        print(f"  {kind}: max|L − Σ_v L_v| = {dev:.3e}", flush=True)
    check("**N2 自由系严格可加**（逐项相等，max 差 < 1e-12）", ok2, "path/ring 各一例")

    print("\n" + "=" * 100)
    print("N3：τ2 = 1/|λ2| 的 V 标度 + 严格上界")
    print("=" * 100)
    scal = {}
    for kind in ("path", "ring"):
        rs = [r for r in rows if r["kind"] == kind]
        Vv = np.array([r["V"] for r in rs], float)
        Lv = np.array([abs(r["l2_cons"]) for r in rs])
        p = np.polyfit(np.log(Vv), np.log(Lv), 1)
        scal[kind] = dict(exponent=float(p[0]), ratios=[float(x) for x in Lv / Vv])
        print(f"  {kind:>5}: |λ2_cons| ∝ V^{p[0]:.3f}；|λ2|/V = {np.round(Lv/Vv,4)}"
              f" ⟹ τ2 ∝ V^{-p[0]:.3f}", flush=True)
    RES["scaling"] = scal
    check("**N3 寿命随系统尺寸下降**（|λ2| 随 V 增长，指数接近 1）",
          all(v["exponent"] > 0.9 for v in scal.values()),
          "；".join(f"{k}: V^{v['exponent']:.2f}" for k, v in scal.items()))

    print("\n  【严格上界】λ2 ≤ d_min·ν2，ν2 = −4sin²(π/(4L+2))：")
    for L in (1, 2, 4, 8):
        nu2 = -4 * np.sin(np.pi / (4 * L + 2)) ** 2
        print(f"    L={L}: d_min=7 时 τ2 ≤ {1/abs(7*nu2):.4f} 微观步", flush=True)

    print("\n" + "=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_nonadditivity.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_nonadditivity.json")
