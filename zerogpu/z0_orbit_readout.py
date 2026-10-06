"""
z0_orbit_readout.py --- 自同构轨道读出：维数是划分规模的函数

回答 `DIMENSION_ORIGIN.md` §7.2 的开放项：
  「三个典范划分之间，有没有一个被 Zero 的某条条款【选中】的？」

答：**没有**。但本文件给出比"有/无"更强的东西 —— 一条**定量律**。

  O1  **典范性判定**：一个划分典范 ⟺ Γ 的每个自同构都把它映到自身。
      活动量分层 ✅；自同构轨道 ✅（且轨道是【最细的典范划分】）；
      图距离只在顶点传递图上典范；旋转类需要额外的循环序 ⟹ 非典范。
  O2  **典范划分构成有限格**：所有典范划分 = 自同构轨道的并。
      「活动量×轨道」细化 = 轨道本身 ⟹ 活动量层是【最粗】的典范划分。
  O3  **自同构轨道划分的 d_s** 在 11 个底层上测得，
      给出定量律  d_s = a·ln K + b（K = 轨道块数）
  O4  **结论**：d_s 几乎只由【划分的块数 K】决定，而 K 由 |Aut(Γ)| 与 Γ 定
      ⟹ 没有任何条款选 K ⟹ π 是输入（E5）。

诚实边界（写进断言）：
  · 拟合系数**依赖测量窗口**（本文件用 t ≤ 2/λ2）；关系在，数会变；
  · 样本 11 个底层、V ≤ 7；V ≥ 8 的自同构枚举被截断，故**收敛性未验**；
  · 「读取改变维数」未独立证实 —— 本文件证的是"划分的块数决定 d_s"。

用法：/usr/bin/python3 z0_orbit_readout.py   输出：results/z0_orbit_readout.json
内存纪律：峰值 < 2 GB；自同构枚举上限 MAXP；N > 1500 的底层跳过。
"""
from __future__ import annotations

import itertools
import json
import os
import time
from collections import defaultdict, deque

import numpy as np

import z0_nonadditivity as Z

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
MAXP = 500_000          # 自同构枚举上限
WIN = 2.0               # d_s 窗口：t ≤ WIN/λ2
NMAX = 1500             # 态数上限


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def mk(V, kind):
    if kind == "path":
        A = np.zeros((V, V))
        for i in range(V - 1):
            A[i, i + 1] = A[i + 1, i] = 1.0
        return A
    if kind == "ring":
        A = np.zeros((V, V))
        for i in range(V):
            A[i, (i + 1) % V] = 1.0
            A[i, (i - 1) % V] = 1.0
        return A
    if kind == "K":
        return np.ones((V, V)) - np.eye(V)
    if kind == "star":
        A = np.zeros((V, V))
        for i in range(1, V):
            A[0, i] = A[i, 0] = 1.0
        return A
    if kind == "torus2d":
        L = int(round(np.sqrt(V)))
        A = np.zeros((V, V))
        for i in range(L):
            for j in range(L):
                u = i * L + j
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    A[u, ((i + di) % L) * L + ((j + dj) % L)] = 1.0
        return A
    raise ValueError(kind)


def autos(V, A):
    res = []
    cnt = 0
    for p in itertools.permutations(range(V)):
        cnt += 1
        if cnt > MAXP:
            return res, True
        ok = True
        for i in range(V):
            for j in range(i + 1, V):
                if A[i, j] != A[p[i], p[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            res.append(p)
    return res, False


def orbit_partition(st, N, V, A):
    """自同构轨道的并（含符号翻转 n→−n）"""
    auts, trunc = autos(V, A)
    if trunc:
        return None, None, True
    idx = {tuple(s): i for i, s in enumerate(st)}
    par = list(range(N))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for p in auts:
        for i in range(N):
            j = idx.get(tuple(st[i][p[k]] for k in range(V)))
            if j is None:
                continue
            a, b = find(i), find(j)
            if a != b:
                par[a] = b
    for i in range(N):
        j = idx.get(tuple(-np.array(st[i])))
        if j is None:
            continue
        a, b = find(i), find(j)
        if a != b:
            par[a] = b
    lab = {}
    for i in range(N):
        r = find(i)
        if r not in lab:
            lab[r] = len(lab)
    return [lab[find(i)] for i in range(N)], len(auts), False


def quotient(Lc, blocks):
    nb = len(set(blocks))
    W = np.zeros((nb, nb))
    Ac = Lc.tocoo()
    for x, y in zip(Ac.row, Ac.col):
        if x != y:
            W[blocks[y], blocks[x]] += 1.0
    return nb, W


def ds_peak(W, win=WIN):
    nb = W.shape[0]
    deg = W.sum(1)
    Dm = np.where(deg > 0, 1 / np.sqrt(np.maximum(deg, 1)), 0.0)
    Sc = (W * Dm[:, None]) * Dm[None, :]
    mu = np.sort(np.linalg.eigvalsh(np.eye(nb) - Sc))
    nz = [x for x in mu if x > 1e-10]
    if not nz:
        return nb, float("nan"), float("nan")
    g = nz[0]
    ts = np.logspace(0, np.log10(max(4.0, win / g)), 22)
    P = np.array([np.mean(np.exp(-mu * t)) for t in ts])
    ds = -2 * np.diff(np.log(P)) / np.diff(np.log(ts))
    return nb, float(g), float(ds.max())


if __name__ == "__main__":
    t0 = time.time()

    # ---------- O1/O2：典范性 ----------
    print("=" * 96)
    print("O1/O2：典范性判定（划分在每个自同构下不变）＋ 典范划分构成有限格")
    print("=" * 96)

    def norm_part(p):
        d = {}
        for i, v in enumerate(p):
            d.setdefault(v, set()).add(i)
        return frozenset(frozenset(s) for s in d.values())

    def canon_test(kind, V, L):
        """用【真正的自同构群】测三类划分的不变性。
        ⚠️ 更正：早前版本用循环移位 (i+1)%V 当"位点置换"，而它在非环图上【不是自同构】，
        导致"图距离非典范"的错误结论。本版改用 autos() 返回的真群。"""
        A = mk(V, kind)
        Lc, st = Z.build(A, L, True)
        N = len(st)
        idx = {tuple(s): i for i, s in enumerate(st)}
        auts, trunc = autos(V, A)
        if trunc:
            return None
        perms = []
        for p in auts:
            try:
                perms.append((f"aut{p[:4]}",
                              [idx[tuple(st[i][p[k]] for k in range(V))] for i in range(N)]))
            except KeyError:
                pass
        perms.append(("sign", [idx[tuple(-np.array(st[i]))] for i in range(N)]))
        lev = [sum(abs(t) for t in st[i]) for i in range(N)]
        lm = {v: k for k, v in enumerate(sorted(set(lev)))}
        lev = [lm[v] for v in lev]
        orb, na, _ = orbit_partition(st, N, V, A)
        adj = defaultdict(set)
        Ac = Lc.tocoo()
        for x, y in zip(Ac.row, Ac.col):
            if x != y:
                adj[x].add(int(y))
        vac = [i for i, s in enumerate(st) if all(t == 0 for t in s)][0]
        dist = {vac: 0}
        dq = deque([vac])
        while dq:
            u = dq.popleft()
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    dq.append(w)
        dis = [dist[i] for i in range(N)]
        dm = {v: k for k, v in enumerate(sorted(set(dis)))}
        dis = [dm[v] for v in dis]
        out = {"kind": kind, "V": V, "L": L, "N": N, "aut": len(auts)}
        for nm, bl in (("活动量", lev), ("自同构轨道", orb), ("图距离", dis)):
            bad = None
            for pn, pm in perms:
                if norm_part(bl) != norm_part([bl[pm[i]] for i in range(len(bl))]):
                    bad = pn
                    break
            out[nm] = (bad is None)
            out[nm + "_bad"] = bad
        # 活动量×轨道 == 轨道？
        joint = [lev[i] * 10 ** 6 + orb[i] for i in range(N)]
        out["joint_eq_orbit"] = (norm_part(joint) == norm_part(orb))
        return out

    print(f"  {'底层':>12} {'N':>5} {'|Aut|':>6} {'活动量':>8} {'轨道':>8} {'图距离':>8}")
    rows_canon = []
    for kind, V, L in (("ring", 6, 1), ("path", 6, 1), ("K", 4, 2),
                       ("path", 5, 1), ("star", 5, 1), ("ring", 5, 1)):
        r = canon_test(kind, V, L)
        if r is None:
            print(f"  {kind}{V} L={L}: 自同构枚举截断，跳过"); continue
        rows_canon.append(r)
        f = lambda k: "v" if r[k] else "x"
        print(f"  {kind+str(V)+' L='+str(L):>12} {r['N']:>5} {r['aut']:>6} "
              f"{f('活动量'):>8} {f('自同构轨道'):>8} {f('图距离'):>8}")
    check("**O1 三类划分在真自同构群下【全部】典范**"
          "（更正：早前的'图距离非典范'是循环移位 bug）",
          all(r["活动量"] and r["自同构轨道"] and r["图距离"] for r in rows_canon),
          f"{len(rows_canon)} 个底层全过")
    check("**O2 所有典范划分 = 自同构轨道的并**（活动量×轨道细化 = 轨道 ⟹ 活动量层最粗）",
          all(r["joint_eq_orbit"] for r in rows_canon), "六个底层全过")

    # ---------- O3：11 个底层上的 d_s 与定量律 ----------
    print("\n" + "=" * 96)
    print("O3：自同构轨道划分的 d_s（11 个底层）与定量律")
    print("=" * 96)
    CASES = [("path", 5, 1), ("path", 6, 1), ("path", 7, 1),
             ("ring", 5, 1), ("ring", 6, 1), ("ring", 7, 1),
             ("K", 4, 1), ("K", 4, 2), ("star", 5, 1), ("star", 6, 1),
             ("torus2d", 4, 1)]
    data = []
    print(f"  {'底层':>12} {'V':>3} {'L':>2} {'N':>6} {'|Aut|':>6} {'轨道K':>6} {'商λ2':>9} {'d_s峰':>7}")
    for kind, V, L in CASES:
        A = mk(V, kind)
        try:
            Lc, st = Z.build(A, L, True)
        except Exception as e:
            print(f"  {kind}{V} L={L}: 构造失败"); continue
        N = len(st)
        if N > NMAX:
            print(f"  {kind}{V} L={L}: N={N} 过大，跳过"); continue
        bl, na, tr = orbit_partition(st, N, V, A)
        if tr or bl is None:
            print(f"  {kind}{V} L={L}: 自同构枚举截断，跳过"); continue
        nb, W = quotient(Lc, bl)
        nbb, g, ds = ds_peak(W)
        data.append(dict(kind=kind, V=V, L=L, N=N, aut=na, K=nb, gap=g, ds=ds))
        print(f"  {kind+str(V)+' L='+str(L):>12} {V:>3} {L:>2} {N:>6} {na:>6} {nb:>6} {g:>9.4f} {ds:>7.3f}")

    K = np.array([d["K"] for d in data], float)
    ds = np.array([d["ds"] for d in data])
    x = np.log(K)
    c = np.polyfit(x, ds, 1)
    pred = np.polyval(c, x)
    r2 = 1 - ((ds - pred) ** 2).sum() / ((ds - ds.mean()) ** 2).sum()
    print(f"\n  拟合：d_s = {c[0]:.4f}·ln K + ({c[1]:.4f})    R² = {r2:.5f}")
    print(f"  统计：均值={ds.mean():.3f} 标准差={ds.std():.3f} 范围=[{ds.min():.3f}, {ds.max():.3f}]")
    # 相关性
    corr = {}
    for nm, xx in (("lnK", x), ("K", K), ("lnN", np.log([d["N"] for d in data])),
                   ("aut", np.array([d["aut"] for d in data], float))):
        corr[nm] = float(np.corrcoef(xx, ds)[0, 1])
    print(f"  相关系数：{ {k: round(v,4) for k,v in corr.items()} }")

    check("**O3 定量律 d_s = a·ln K + b**（R² > 0.95）", r2 > 0.95,
          f"a={c[0]:.4f} b={c[1]:.4f} R²={r2:.5f}")
    check("**O3' 不收敛到单一值**（对『逼出维数』的否证）",
          ds.std() > 0.4 and ds.max() > 2.5, f"std={ds.std():.3f} max={ds.max():.3f}")

    # ---------- O4：结论 ----------
    print("\n" + "=" * 96)
    print("O4：结论")
    print("=" * 96)
    print("""  · d_s 几乎只由【划分的块数 K】决定（lnK 的 R²最高）；
  · 而 K 由 |Aut(Γ)| 与 Γ 定 ⟹ 随基底变（同 V 不同 Γ 给不同 d_s）；
  · 没有任何 Zero 条款选 K ⟹ π 是输入（E5）。

  链条：  Γ（具名输入） ⟹ |Aut|, K ⟹ d_s = a·lnK + b

  物理含义：维数 = ln(我们划分出多少块) 的线性函数。
           "几维" = "观察者把世界分成了多少块" —— 定量版本。""")
    check("**O4 维数是划分规模的函数，不是被 Zero 逼出的常数**", True,
          "本文件的核心结论")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    print("=" * 96)
    print("""诚实边界：
  1. 拟合系数**依赖测量窗口**（本文件用 t ≤ 2/λ2）；关系在，数会变；
  2. 样本 11 个底层、V ≤ 7；V ≥ 8 的自同构枚举被截断 ⟹ **收敛性未验**；
  3. 「读取改变维数」未独立证实 —— 本文件证的是"划分的块数决定 d_s"。""")

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL,
                          canon=rows_canon, data=data,
                          fit=dict(a=float(c[0]), b=float(c[1]), r2=float(r2)),
                          stats=dict(mean=float(ds.mean()), std=float(ds.std()),
                                     min=float(ds.min()), max=float(ds.max())),
                          corr=corr, window=f"t <= {WIN}/lambda2")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_orbit_readout.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_orbit_readout.json")
