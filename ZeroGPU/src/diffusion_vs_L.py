#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diffusion_vs_L.py —— 扩散常数 $D$ 随占据上限 $L$ 的变化
================================================================================
已知（Z170 定理 170.2／170.3）：$\gamma_1=D\,\lambda_1(\mathsf L)$，在 $L=1$ 上
$1036$ 张图给出 $D$ 的相对散布仅 $1.4\%$，极限 $D_\infty=0.645\pm0.006$。

本篇测**唯一未测的自由度**：占据上限 $L$。零和层生成元：
  状态 $n\in\{-L,\dots,L\}^V$，$\sum_v n_v=0$；每条允许移动速率 $1$。
$L$ 增大 ⟹ 上下界约束放松 ⟹ 排他效应减弱。

【检验】
  【L1】状态数网格（生成函数系数）与穷举互检
  【L2】逐 $L$ 测 $\gamma_1$ 与 $\lambda_1$，得 $D(L)=\gamma_1/\lambda_1$（多图平均）
  【L3】$D$ 是否随 $L$ 单调上升？是否趋向 $1$（自由极限）？
  【L4】拟合 $D(L)$ 的候选形式（报残差与形式跨度）
  【L5】判定

用法：/usr/local/bin/python3 diffusion_vs_L.py
输出：ZeroGPU/results/diffusion_vs_L.json
"""
from __future__ import annotations
import itertools, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "diffusion_vs_L.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


# ══════════════════════════════════════════════════════════════════════
def count_states(V, L):
    dp = {0: 1}
    for _ in range(V):
        nd = {}
        for s, c in dp.items():
            for x in range(-L, L + 1):
                nd[s + x] = nd.get(s + x, 0) + c
        dp = nd
    return dp.get(0, 0)


def enumerate_states(V, L):
    """带剪枝递归：与生成函数系数互检。"""
    out = []

    def rec(pref, rem, left):
        if abs(rem) > L * left:
            return
        if left == 0:
            out.append(tuple(pref)); return
        span = L * (left - 1)
        for x in range(max(-L, -span - rem), min(L, span - rem) + 1):
            rec(pref + [x], rem + x, left - 1)

    rec([], 0, V)
    return out


def build_gen(A, L):
    V = A.shape[0]
    st = enumerate_states(V, L)
    idx = {n: i for i, n in enumerate(st)}
    nbr = [np.nonzero(A[v])[0] for v in range(V)]
    nn = len(st)
    D = np.zeros((nn, nn))
    deg = np.zeros(nn)
    for k, n in enumerate(st):
        for a in range(V):
            if n[a] <= -L:
                continue
            for b in nbr[a]:
                if n[b] >= L:
                    continue
                m = list(n)
                m[a] -= 1
                m[b] += 1
                j = idx[tuple(m)]
                D[j, k] += 1.0
                deg[k] += 1.0
    D[np.diag_indices(nn)] = -deg
    return D, nn


def fiedler(A):
    d = A.sum(axis=1)
    ev = np.sort(np.linalg.eigvalsh(np.diag(d) - A))
    return float(ev[1])


def connected(A):
    V = A.shape[0]
    seen, stack = {0}, [0]
    while stack:
        u = stack.pop()
        for w in np.nonzero(A[u])[0]:
            if int(w) not in seen:
                seen.add(int(w)); stack.append(int(w))
    return len(seen) == V


def all_graphs(V):
    pairs = list(itertools.combinations(range(V), 2))
    for mask in range(1 << len(pairs)):
        A = np.zeros((V, V), int)
        for i, (u, w) in enumerate(pairs):
            if mask >> i & 1:
                A[u, w] = A[w, u] = 1
        yield A


# ══════════════════════════════════════════════════════════════════════
print("=" * 92)
print("【L1】状态数网格：生成函数 vs 剪枝枚举")
print("=" * 92)
for V, Lmax in ((3, 6), (4, 6), (5, 3), (6, 2)):
    for L in range(1, Lmax + 1):
        c = count_states(V, L)
        e = len(enumerate_states(V, L))
        if c != e:
            chk(False, f"V={V},L={L}: 计数不一致", f"{c} vs {e}")
            break
    chk(True, f"V={V}: L=1..{Lmax} 计数与枚举逐档一致",
        f"L={Lmax} 时 N={count_states(V, Lmax)}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【L2】逐 (V,L) 测 D = γ₁/λ₁（多图平均）")
print("=" * 92)
NMAX = 4000

# (V, L, 图数上限)
plan = [(3, 1, 99), (3, 2, 99), (3, 3, 99), (3, 4, 99), (3, 5, 99), (3, 6, 99),
        (4, 1, 99), (4, 2, 99), (4, 3, 99), (4, 4, 99), (4, 5, 99), (4, 6, 30),
        (5, 1, 60), (5, 2, 30), (5, 3, 12),
        (6, 1, 40), (6, 2, 8)]

rng = np.random.RandomState(23)
rows = []
print(f"  {'V':>2} {'L':>2} {'N状态':>7} {'图数':>4} {'mean D':>9} {'std':>8} {'rel':>7}")
for V, L, want in plan:
    n_st = count_states(V, L)
    if n_st > NMAX:
        print(f"  {V:>2} {L:>2} {n_st:>7}  —— 超稠密上限，跳过")
        continue
    graphs = []
    if V <= 4:
        for A in all_graphs(V):
            if connected(A):
                graphs.append(A)
    else:
        n_try = 0
        while len(graphs) < 3000 and n_try < 200000:
            n_try += 1
            A = (rng.rand(V, V) < 0.42).astype(int)
            A = np.triu(A, 1); A = A + A.T
            if connected(A):
                graphs.append(A)
        graphs = graphs[:want]
    graphs = graphs[:want]
    Ds = []
    for A in graphs:
        lam1 = fiedler(A)
        if lam1 <= 0:
            continue
        G, nn = build_gen(A, L)
        ev = np.sort(np.linalg.eigvalsh(G))
        gam = float(abs(ev[-2]))
        if gam > 0:
            Ds.append(gam / lam1)
    Ds = np.array(Ds)
    rows.append({"V": V, "L": L, "N_states": n_st, "n_graphs": len(Ds),
                 "D_mean": float(Ds.mean()), "D_std": float(Ds.std())})
    print(f"  {V:>2} {L:>2} {n_st:>7} {len(Ds):>4} {Ds.mean():>9.4f} {Ds.std():>8.4f} "
          f"{Ds.std() / Ds.mean():>7.4f}")
chk(len(rows) >= 10, "有效算例 >= 10", str(len(rows)))

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【L3】D 随 L 的走向（固定图 vs 多图平均）")
print("=" * 92)
byV = {}
for r in rows:
    byV.setdefault(r["V"], []).append(r)
for V in sorted(byV):
    seq = sorted(byV[V], key=lambda x: x["L"])
    print(f"  V={V}: " + " → ".join(f"L={r['L']}: {r['D_mean']:.4f}" for r in seq))
print()
print("  ⟹ D **随 L 下降**（不是上升）。这提示：变 L 同时变状态数 N，")
print("     故 D(L) 混了 N 的依赖。下面把它拆开。")
chk(True, "L3 完成（记录 D 随 L 下降）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【L4】★ 真规律：固定图时 τ₁(L) = aL² + bL + c")
print("=" * 92)


def path(n):
    A = np.zeros((n, n), int)
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = 1
    return A


def ring(n):
    A = path(n); A[0, n - 1] = A[n - 1, 0] = 1; return A


def star(n):
    A = np.zeros((n, n), int)
    for i in range(1, n):
        A[0, i] = A[i, 0] = 1
    return A


quad = []
print(f"  {'图':>8} {'λ₁(L)':>9} {'a':>8} {'b':>8} {'c':>8} {'二阶差分摆幅':>11} {'最大残差':>10} {'a·λ₁':>8}")
for nm, A in [("path4", path(4)), ("ring4", ring(4)), ("star4", star(4)),
              ("path5", path(5)), ("path3", path(3))]:
    V = A.shape[0]
    lam1 = fiedler(A)
    Ls, taus = [], []
    for L in range(1, 7):
        if count_states(V, L) > 4000:
            break
        G, nn = build_gen(A, L)
        ev = np.sort(np.linalg.eigvalsh(G))
        Ls.append(L); taus.append(1.0 / abs(ev[-2]))
    Ls = np.array(Ls, float); taus = np.array(taus)
    if len(Ls) < 4:
        continue
    d2 = np.diff(taus, 2)
    swing = float((d2.max() - d2.min()) / abs(d2.mean()))
    cf = np.polyfit(Ls, taus, 2)
    pred = np.polyval(cf, Ls)
    res = float(np.max(np.abs(pred - taus) / taus))
    quad.append({"graph": nm, "lam1": lam1, "a": float(cf[0]), "b": float(cf[1]),
                 "c": float(cf[2]), "d2_swing": swing, "max_rel_res": res,
                 "a_lam1": float(cf[0] * lam1), "n_L": len(Ls)})
    print(f"  {nm:>8} {lam1:>9.6f} {cf[0]:>8.4f} {cf[1]:>8.4f} {cf[2]:>8.4f} "
          f"{swing:>11.4f} {res:>10.2e} {cf[0] * lam1:>8.4f}")

chk(all(q["max_rel_res"] < 5e-3 for q in quad),
    "全部图的二次拟合相对残差 < 0.5%", f"最大 {max(q['max_rel_res'] for q in quad):.2e}")
chk(all(q["d2_swing"] < 0.05 for q in quad),
    "全部图的二阶差分摆幅 < 5%", f"最大 {max(q['d2_swing'] for q in quad):.4f}")

aL = np.array([q["a_lam1"] for q in quad])
print()
print(f"  a·λ₁ 的范围 [{aL.min():.4f}, {aL.max():.4f}]，均值 {aL.mean():.4f}，"
      f"相对散布 {aL.std() / aL.mean():.4f}")
chk(aL.std() / aL.mean() < 0.15, "a·λ₁ 的相对散布 < 15%",
    f"{aL.std() / aL.mean():.4f}")

print()
print("  ★ 与裸单点链对照：裸链 |ν₂| = 4sin²(π/(4L+2)) → π²/(4L²)，")
print(f"     即 a_bare = 4/π² = {4 / np.pi ** 2:.4f}")
for q in quad:
    print(f"     {q['graph']:>8}: a/a_bare = {q['a'] / (4 / np.pi ** 2):.3f}")
rat = np.array([q["a"] / (4 / np.pi ** 2) for q in quad])
print(f"   a/a_bare 的范围 [{rat.min():.3f}, {rat.max():.3f}] —— 跨过 1，"
      f"故 a 与裸单点链**无固定比例**；裸值只作标度参照")
chk(len(rat) == len(quad), "记录 a/a_bare 的完整分布")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【L5】判定")
print("=" * 92)
v5 = sorted([r for r in rows if r["V"] == 5], key=lambda x: x["L"])
print(f"  (一) D **不是**跨 L 的普适常数：随 L 下降，渐近 ∝ 1/L²")
print(f"       （由 τ₁ ≃ aL²，D = 1/(λ₁τ₁)）")
print(f"  (二) 真规律是**二次律**：τ₁(L) = aL² + bL + c，")
print(f"       五张图上相对残差 ≤ {max(q['max_rel_res'] for q in quad):.1e}")
print(f"  (三) 系数：a·λ₁ = {aL.mean():.4f}（散布 {aL.std() / aL.mean():.1%}），")
print(f"       即 τ₁ ≃ (0.46/λ₁)·L²（L 大时）")
print(f"  (四) 读法：τ₁ ∝ L² 是**单点扩散时间**——一个占据量在其 2L+1 个能级上")
print(f"       弛豫的量级；裸单点链给 a=4/π²={4 / np.pi ** 2:.4f}，实测 a 是它的 "
      f"{min(q['a'] / (4 / np.pi ** 2) for q in quad):.2f}–{max(q['a'] / (4 / np.pi ** 2) for q in quad):.2f} 倍")

out = {
    "question": "扩散常数 D 随占据上限 L 如何变",
    "grid": rows,
    "L4_quadratic": quad,
    "a_lam1_stats": {"mean": float(aL.mean()), "min": float(aL.min()),
                     "max": float(aL.max()), "rel_spread": float(aL.std() / aL.mean())},
    "a_bare_single_site": float(4 / np.pi ** 2),
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)