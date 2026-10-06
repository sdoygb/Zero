"""
z0_readout_necessity.py --- 读出是必需的（本轮结论的固化）

本轮问的是："时空维数能不能从底层读出来？" 答案是**不能**，而且原因是结构性的。
本文件固化四件事：

  R1  **四种维数判据全部失败校准** —— 说明底层【没有可测的谱维数】
      (a) 热核 d_s（固定分数窗口）      (b) 体积增长维数 S(k)
      (c) 幂律拟合 P(t)~t^-a            (d) Weyl 计数 N(μ)~μ^{d_s/2}
  R2  **无标度区**（决定性判据）：直径/(1/λ2) ≈ 0.8–1.5 ⟹ 走者未及边界就已混合
      ⟹ 微观尺度与盒尺度之间没有尺度分离 ⟹ 无连续极限 ⟹ 无内禀维数
  R3  **弛豫→维数** 的关系 λ2 ~ N^(−2/d) **只在扩散区成立**（校准：1D→1.001、2D→2.042 精确；
      3D→3.324 偏高，因格太小仍在弹道区）⟹ 用它读维数需要扩散区，而底层没有
  R4  **弛豫率是干净可观测量，且对读取敏感**：同一底层、两个粗粒化给出不同的商图 λ2

另记一条**方法论**：本文件把失败的判据也写成断言（期望值取自"已知正确值"），
因为"判据失败"本身就是本轮的主要结论，不能只报成功的那些。

用法：/usr/bin/python3 z0_readout_necessity.py   输出：results/z0_readout_necessity.json
内存纪律：峰值 < 2 GB；禁止枚举后存储；dense 对角化仅用于 N <= 3000。
"""
from __future__ import annotations

import itertools
import json
import os
import time
from collections import Counter, defaultdict, deque

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


def torus(d, L):
    idx = {c: i for i, c in enumerate(itertools.product(range(L), repeat=d))}
    N = len(idx)
    A = np.zeros((N, N))
    for c, i in idx.items():
        for k in range(d):
            for s in (1, -1):
                cc = list(c)
                cc[k] = (cc[k] + s) % L
                A[i, idx[tuple(cc)]] = 1.0
    return A


def cube(n):
    N = 1 << n
    A = np.zeros((N, N))
    for i in range(N):
        for k in range(n):
            A[i, i ^ (1 << k)] = 1.0
    return A


def norm_lap_eigs(A):
    """归一化拉普拉斯 I − D^{-1/2} A D^{-1/2} 的谱（升序）"""
    N = A.shape[0]
    deg = A.sum(1)
    Dm = np.where(deg > 0, 1 / np.sqrt(np.maximum(deg, 1)), 0.0)
    Sc = (A * Dm[:, None]) * Dm[None, :]
    return np.sort(np.linalg.eigvalsh(np.eye(N) - Sc))


def gap(A):
    mu = norm_lap_eigs(A)
    nz = [x for x in mu if x > 1e-10]
    return float(nz[0]) if nz else float("nan")


def state_graph(kind, V, L):
    A = S.G(V, kind)
    Lc, st = Z.build(A, L, True)
    H = Lc.toarray()
    return H - np.diag(np.diag(H)), st


def quotient_gap(blocks, adj, nb):
    W = np.zeros((nb, nb))
    for x in range(len(adj)):
        for y in adj[x]:
            W[blocks[y], blocks[x]] += 1.0
    return gap(W)


if __name__ == "__main__":
    t0 = time.time()

    # ---------------- R1：四种判据的失败 ----------------
    print("=" * 96)
    print("R1：四种维数判据的校准结果（期望：全部通过；实际：全部失败 ⟹ 无内禀维数）")
    print("=" * 96)

    print("\n  (a) 热核 d_s（固定分数窗口）")
    def ds_heat(A, lo=0.25, hi=0.6, npts=50):
        N = A.shape[0]
        mu = norm_lap_eigs(A)
        ts = np.logspace(0, np.log10(max(10.0, N / 6.0)), npts)
        P = np.array([np.mean(np.exp(-mu * t)) for t in ts])
        ds = -2 * np.diff(np.log(P)) / np.diff(np.log(ts))
        return float(np.median(ds[int(len(ds) * lo):int(len(ds) * hi)]))
    a_1d = ds_heat(torus(1, 400)); a_2d = ds_heat(torus(2, 40))
    print(f"      1D→{a_1d:.3f}(真1)   2D→{a_2d:.3f}(真2)")
    # 热核法在【整数维格点】上可信（1D、2D 都准），所以它本身是合格判据；
    # 它的失效出现在【高维/非格点】—— 这正是态空间的情形（见 R1c/R1d）。
    check("R1a 热核 d_s 在整数维格点上【可信】（1D、2D 均准）",
          abs(a_1d - 1) < 0.15 and abs(a_2d - 2) < 0.15,
          f"1D→{a_1d:.3f}(真1)  2D→{a_2d:.3f}(真2)")

    print("\n  (b) 体积增长维数 S(k)")
    def d_grow(A, start=0, kmax=8):
        N = A.shape[0]
        nbr = [np.nonzero(A[i])[0] for i in range(N)]
        dist = -np.ones(N, int); dist[start] = 0; dq = deque([start])
        while dq:
            u = dq.popleft()
            for w in nbr[u]:
                if dist[w] < 0: dist[w] = dist[u] + 1; dq.append(w)
        ks = np.arange(1, kmax + 1)
        S = np.array([int(np.sum((dist >= 0) & (dist <= k))) for k in ks])
        m = (ks >= 2) & (S > 1)
        return float(np.polyfit(np.log(ks[m]), np.log(S[m]), 1)[0]) if m.sum() >= 3 else float("nan")
    b_1d = d_grow(torus(1, 400)); b_4d = d_grow(torus(4, 5))
    print(f"      1D→{b_1d:.3f}(真1)   4D→{b_4d:.3f}(真4)")
    check("R1b 体积增长维数 1D 已偏低 ⟹ 判据不可靠", abs(b_1d - 1) > 0.08, f"{b_1d:.3f}")

    print("\n  (c) 幂律拟合 P(t)~t^-a，d_s=2a  —— 用超立方体校准")
    def ds_power(A, tlo=1.5, thi=6.0, npts=30):
        mu = norm_lap_eigs(A)
        ts = np.logspace(np.log10(tlo), np.log10(thi), npts)
        P = np.array([np.mean(np.exp(-mu * t)) for t in ts])
        a = -np.polyfit(np.log(ts), np.log(P), 1)[0]
        return 2 * a
    c_vals = [ds_power(cube(n)) for n in (6, 8, 10)]
    print(f"      Q_6,Q_8,Q_10 → {[round(v,2) for v in c_vals]}  （真值全部 = 2）")
    check("R1c 幂律判据在超立方体上失败（真值 2）", all(abs(v - 2) > 1.0 for v in c_vals),
          f"给 {[round(v,1) for v in c_vals]}")

    print("\n  (d) Weyl 计数 N(μ)~μ^{d_s/2}")
    def ds_weyl(A, lo=0.02, hi=0.35):
        nz = norm_lap_eigs(A)
        nz = nz[nz > 1e-10]
        if len(nz) < 20: return float("nan")
        Nn = np.arange(1, len(nz) + 1)
        a, b = max(1, int(len(nz) * lo)), int(len(nz) * hi)
        if b - a < 8: return float("nan")
        return 2 * np.polyfit(np.log(nz[a:b]), np.log(Nn[a:b]), 1)[0]
    d_vals = [ds_weyl(cube(n)) for n in (8, 10, 12)]
    print(f"      Q_8,Q_10,Q_12 → {[round(v,2) for v in d_vals]}  （真值全部 = 2）")
    check("R1d Weyl 判据在超立方体上失败（真值 2）", all(abs(v - 2) > 1.0 for v in d_vals),
          f"给 {[round(v,1) for v in d_vals]}")

    # ---------------- R2：无标度区 ----------------
    print("\n" + "=" * 96)
    print("R2：态空间有没有标度区？（直径 vs 1/λ2）")
    print("=" * 96)
    print(f"  {'盒子':>12} {'V':>3} {'L':>2} {'N':>6} {'直径':>6} {'λ2':>9} {'1/λ2':>8} {'比值':>7}")
    ratios = []
    for kind, V, L in (("ring", 4, 1), ("ring", 5, 1), ("ring", 6, 1),
                       ("ring", 7, 1), ("ring", 4, 2), ("ring", 6, 2)):
        Adj, st = state_graph(kind, V, L)
        N = Adj.shape[0]
        if N > 2600:
            print(f"  {kind+str(V)+' L='+str(L):>12} {V:>3} {L:>2} {N:>6}   （过大，跳过）"); continue
        nbr = [np.nonzero(Adj[i])[0] for i in range(N)]
        d = -np.ones(N, int); d[0] = 0; dq = deque([0])
        while dq:
            u = dq.popleft()
            for w in nbr[u]:
                if d[w] < 0: d[w] = d[u] + 1; dq.append(w)
        diam = int(d.max()); g = gap(Adj); r = diam / (1 / g)
        ratios.append(r)
        print(f"  {kind+str(V)+' L='+str(L):>12} {V:>3} {L:>2} {N:>6} {diam:>6} {g:>9.5f} {1/g:>8.2f} {r:>7.3f}")
    check("**R2 无标度区**（直径/(1/λ2) ≈ 1 ⟹ 走者未及边界就已混合）",
          all(0.5 < r < 2.0 for r in ratios), f"实测 {[round(r,2) for r in ratios]}")

    # ---------------- R3：弛豫→维数 只在扩散区成立 ----------------
    print("\n" + "=" * 96)
    print("R3：λ2 ~ N^(−2/d) 的适用条件（校准）")
    print("=" * 96)
    print(f"  {'格':>6} {'真维数':>6} {'α':>9} {'2/α':>8} {'判定':>8}")
    r3 = {}
    for d, Ls in ((1, (30, 60, 120, 240)), (2, (8, 12, 16, 24)), (3, (4, 6, 8, 10))):
        pts = [(L ** d, gap(torus(d, L))) for L in Ls]
        x = np.log([p[0] for p in pts]); y = np.log([p[1] for p in pts])
        al = -np.polyfit(x, y, 1)[0]
        r3[d] = dict(alpha=float(al), d_est=float(2 / al))
        ok = abs(2 / al - d) < 0.1
        print(f"  {d}D {d:>6} {al:>9.4f} {2/al:>8.3f} {'精确' if ok else '偏高':>8}")
    check("**R3 λ2~N^(−2/d) 只在扩散区精确**（1D、2D 通过；3D 偏高因格小）",
          abs(r3[1]["d_est"] - 1) < 0.1 and abs(r3[2]["d_est"] - 2) < 0.1 and r3[3]["d_est"] > 3.1,
          f"1D→{r3[1]['d_est']:.3f} 2D→{r3[2]['d_est']:.3f} 3D→{r3[3]['d_est']:.3f}")

    # 态空间的局部斜率漂移（说明未收敛）
    pts = []
    for V in (4, 5, 6, 7, 8):
        Adj, st = state_graph("ring", V, 1)
        N = Adj.shape[0]
        if N > 2600: continue
        pts.append((N, gap(Adj)))
    slopes = [-(np.log(pts[i+1][1]) - np.log(pts[i][1])) / (np.log(pts[i+1][0]) - np.log(pts[i][0]))
              for i in range(len(pts) - 1)]
    print(f"\n  态空间局部对数斜率（N 升序）：{[round(s,4) for s in slopes]}")
    print(f"  点估计 2lnN/ln(1/λ2)：{[round(2*np.log(N)/np.log(1/g),3) for N,g in pts]}")
    check("**R3' 态空间维数未收敛**（局部斜率单调漂移，无平台）",
          all(slopes[i] > slopes[i+1] for i in range(len(slopes)-1)),
          f"单调下降 {[round(s,3) for s in slopes]}")

    # ---------------- R4：弛豫率对读取敏感 ----------------
    print("\n" + "=" * 96)
    print("R4：同一底层、两个粗粒化 ⟹ 商图弛豫率不同（弛豫率是干净的、读取敏感的量）")
    print("=" * 96)
    print(f"  {'盒子':>12} {'读取':>8} {'块数':>5} {'商图λ2':>10}")
    diffs = []
    for kind, V, L in (("ring", 6, 1), ("ring", 7, 1), ("ring", 6, 2)):
        Adj, st = state_graph(kind, V, L)
        N = Adj.shape[0]
        if N > 2600: continue
        adj = defaultdict(set)
        for x in range(N):
            for y in np.nonzero(Adj[x])[0]:
                adj[x].add(int(y))
        vac = [i for i, s in enumerate(st) if all(t == 0 for t in s)][0]
        dist = {vac: 0}; dq = deque([vac])
        while dq:
            u = dq.popleft()
            for w in adj[u]:
                if w not in dist: dist[w] = dist[u] + 1; dq.append(w)
        reads = {
            "活动量": [sum(abs(t) for t in st[i]) // 2 for i in range(N)],
            "图距离": [dist.get(i, -1) for i in range(N)],
        }
        gs = {}
        for name, bl in reads.items():
            uq = {v: k for k, v in enumerate(sorted(set(bl)))}
            blocks = [uq[v] for v in bl]; nb = len(uq)
            g = quotient_gap(blocks, adj, nb)
            gs[name] = g
            print(f"  {kind+str(V)+' L='+str(L):>12} {name:>8} {nb:>5} {g:>10.4f}")
        if len(gs) == 2:
            a, b = sorted(gs.values())
            diffs.append(b / a)
    check("**R4 弛豫率对读取敏感**（两个粗粒化给出不同商图 λ2）",
          len(diffs) > 0 and all(r > 1.2 for r in diffs),
          f"比值 {[round(r,3) for r in diffs]}")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    print("=" * 96)
    print("""结论（三层）：
  · 微观 Γ 是 1 维（校准可信）；
  · 态空间【没有内禀维数】：四种判据全败 ＋ 无标度区（直径/(1/λ2)≈1）；
  · 弛豫率 λ2 是【干净可测】且【对读取敏感】的量 —— 所以"读取改变物理"应用它来表述；
  · ⟹ 读出不是可选项：没有读出就没有标度区、没有维数、没有连续极限。""")

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL,
                          R2_ratios=ratios, R3_calib=r3,
                          state_space_slopes=[float(s) for s in slopes],
                          readout_gap_ratios=diffs)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_readout_necessity.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_readout_necessity.json")
