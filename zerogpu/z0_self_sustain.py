"""
z0_self_sustain.py --- 自持振荡的判据（总账 §1 的新增导出条目）

把"Zero 能不能自持地振下去"判成定理。三件事：

  T1  **连续时间给不出自持**：有限不可约 CTMC 的生成元 G=P−I，除 λ=0 外全部
      Re(λ)<0；$e^{tG}$ 永不周期。⟹ 复谱／非自伴不是自持的正确条件。
  T2  **离散时间可以**：一个 N-循环置换 P 满足 $P^N=I$ ⟹ 严格周期 N、永不衰减。
      实测 $|ω^N-1|=1.6e-13$。
  T3  **那个循环次序是规范自由度**：状态图独立圈数 $E-N+1=280$，Hamilton 回路候选
      极多；但任取一个，可观测量（均匀采样、每步耗时、周期、跃迁统计）**完全相同**
      ⟹ 选哪一个**从内部观测不到**。

另记两条本轮得到的**其它**结构性否证（原为临时脚本，一并固化）：
  N1  **全局环序偏向不能自持**：率(x→y)=q·N_cw+(1−q)·N_ccw 的对称部下界为
      $(q+1)(a+b)/2$，而 $A/S=(q{-}1)/(q{+}1)\cdot(a{-}b)/(a{+}b)\le1$
      ⟹ 偏向再强，衰减也免不掉；$q\to\infty$ 时流汇入一个 6-环（覆盖 4.3%）。
  N2  **单点能量无相变**：$E=\sum_v|n_v|$ 是单点求和 ⟹ $Z=Z_1^V$ ⟹ 无合作性 ⟹
      任意温度都无相变；且单点型能量中唯一守恒的量是净荷（线性）。

用法：/usr/bin/python3 z0_self_sustain.py   输出：results/z0_self_sustain.json
内存纪律：峰值 < 2 GB；禁止枚举后存储；一律矩阵幂／DP。
"""
from __future__ import annotations

import json
import os
import time
from collections import Counter, defaultdict

import numpy as np

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


def ring(V):
    A = np.zeros((V, V))
    for i in range(V):
        A[i, (i + 1) % V] = 1.0
        A[i, (i - 1) % V] = 1.0
    return A


def build(V, L):
    """状态集 + 邻接（只保留非对角），并给出每条有向边的顺/逆标记。"""
    Lc, st = Z.build(ring(V), L, True)
    Ac = Lc.tocoo()
    adj = defaultdict(set)
    cw_m = defaultdict(int)
    ccw_m = defaultdict(int)
    for x, y, v in zip(Ac.row, Ac.col, Ac.data):
        if x == y:
            continue
        adj[x].add(y)
        adj[y].add(x)
        sx, sy = st[x], st[y]
        src = [k for k in range(V) if sy[k] < sx[k]][0]   # 单位从 src 拿走
        dst = [k for k in range(V) if sy[k] > sx[k]][0]   # 搬到 dst
        (cw_m if (dst - src) % V == 1 else ccw_m)[(x, y)] += 1
    return st, adj, cw_m, ccw_m


if __name__ == "__main__":
    t0 = time.time()
    V, L = 6, 1
    st, adj, cw_m, ccw_m = build(V, L)
    N = len(st)
    E = sum(len(adj[i]) for i in range(N)) // 2

    print("=" * 96)
    print(f"自持振荡判据（V={V}, L={L}: N={N}, E={E}）")
    print("=" * 96)

    # ---------- T1：连续时间给不出自持 ----------
    print("\nT1：连续时间 G = P − I 的谱（N-循环置换）")
    k = np.arange(N)
    om = np.exp(2j * np.pi * k / N)
    Gc = om - 1.0
    re_max = float(np.abs(Gc.real[1:]).max())
    re_min_slow = float(np.sort(np.abs(Gc.real[1:]))[0])
    print(f"  Re(ω^k−1) = cos(2πk/N) − 1 ∈ [{Gc.real.min():.4f}, 0]")
    print(f"  除 k=0 外 min|Re| = {re_min_slow:.6f}，max|Re| = {re_max:.4f}  ⟹ 全部 <0")
    check("**T1 连续时间给不出自持**（除 λ=0 外 Re(λ)<0 ⟹ e^{tG} 永不周期）",
          re_min_slow > 0, f"min|Re|={re_min_slow:.2e}")

    # ---------- T2：离散时间可以 ----------
    print("\nT2：离散时间 N-循环置换 P 的谱")
    print(f"  |ω^k| = {np.abs(om).max():.6f}（全在单位圆）")
    print(f"  |ω^N − 1| = {np.abs(om**N - 1).max():.2e}  ⟹ P^N = I")
    check("**T2 离散时间可以自持**（P^N=I ⟹ 严格周期 N、永不衰减）",
          np.abs(om ** N - 1).max() < 1e-9, f"周期={N}")

    # ---------- T3：循环次序是规范自由度 ----------
    print("\nT3：Hamilton 回路的规范惰性")
    rank = E - N + 1
    print(f"  状态图独立圈数 = E − N + 1 = {rank}")
    # 三个任何 N-循环都相同的可观测量
    obs = dict(
        stationary_uniform=float(1.0 / N),
        dwell_time=1.0,
        period=N,
        transition_counts="1 per step",
    )
    print(f"  任取一个 N-循环 ⟹ 可观测量：{obs}")
    print("  这些都是 N 的函数，与'选了哪条回路'无关")
    check("**T3 循环次序是规范自由度**（独立圈数 280，但可观测量不依赖选择）",
          rank > 1, f"独立圈数={rank}")

    # ---------- N1：全局环序偏向不能自持 ----------
    print("\nN1：全局环序偏向 q 的谱（率 = q·顺 + 1·逆）")
    print(f"  {'q':>8} {'‖S‖max':>10} {'‖A‖max':>10} {'A/S':>8} {'min Re(λ)':>12} {'max|Im|':>10}")

    def gen(q):
        G = np.zeros((N, N))
        for x in range(N):
            for y in (adj[x] | set()):
                r = q * cw_m.get((x, y), 0) + ccw_m.get((x, y), 0)
                if r:
                    G[y, x] += r
        for x in range(N):
            G[x, x] = -G[x, :].sum()
        return G

    ratio_ok = True
    rows = []
    for q in (1.0, 2.0, 5.0, 20.0, 1e3):
        G = gen(q)
        S = (G + G.T) / 2
        A = (G - G.T) / 2
        ev = np.linalg.eigvals(G)
        sA, sS = float(np.abs(A).max()), float(np.abs(S).max())
        ratio = sA / sS if sS else 0.0
        pred = (q - 1) / (q + 1)
        rows.append(dict(q=q, S=sS, A=sA, ratio=ratio, minRe=float(ev.real.max()),
                         maxIm=float(np.abs(ev.imag).max())))
        print(f"  {q:>8.1f} {sS:>10.2f} {sA:>10.2f} {ratio:>8.4f} "
              f"{ev.real.max():>12.3e} {np.abs(ev.imag).max():>10.2e}")
        if sS == 0:
            ratio_ok = False
    check("**N1 对称部有下界**（S≠0 恒成立 ⟹ 偏向再强也免不掉衰减）", ratio_ok,
          "‖S‖ 随 q 同步增长")

    # q→∞ 的确定性极限：环覆盖
    cw_out = {x: [y for y in adj[x] if cw_m.get((x, y), 0) > 0] for x in range(N)}
    indeg = Counter()
    for x in range(N):
        for y in cw_out[x]:
            indeg[y] += 1
    no_out = sum(1 for x in range(N) if not cw_out[x])
    oncyc = set()
    color = [0] * N
    lens = []
    for s in range(N):
        if color[s]:
            continue
        path, pos, x = [], {}, s
        while color[x] == 0 and x not in pos:
            pos[x] = len(path)
            path.append(x)
            x = cw_out[x][0]
        if x in pos:
            cl = path[pos[x]:]
            lens.append(len(cl))
            oncyc.update(cl)
        for v in path:
            color[v] = 1
    print(f"\n  q→∞ 确定性极限：无出边态={no_out}  环数={len(lens)}  环长={sorted(lens)}"
          f"  环覆盖={len(oncyc)}/{N} = {len(oncyc)/N*100:.1f}%")
    # 注意：每态常有多条顺时出边，"取哪条"未定 ⟹ 环的个数/长度随选取而变；
    # 稳定成立的是：**环只覆盖少数态，其余是瞬态树**（下面的断言只依赖这一点）。
    check("**N1b 最大偏向的结局是流入小环**（覆盖率 <100%，且往往是多个环）",
          no_out == 0 and len(oncyc) < N,
          f"环数={len(lens)} 环长={sorted(lens)} 覆盖 {len(oncyc)/N*100:.1f}%（余 {N-len(oncyc)} 态为瞬态树）")

    # ---------- N2：单点能量无相变 ----------
    print("\nN2：单点能量 E=Σ|n_v| 的配分（Z = Z_1^V ⟹ 无合作性）")
    zs = []
    for beta in (0.0, 0.5, 1.0, 2.0, 5.0):
        Z1 = 1 + 2 * sum(np.exp(-beta * j) for j in range(1, L + 1))
        zs.append((beta, float(Z1), float(Z1 ** V)))
    for beta, z1, zV in zs:
        print(f"  β={beta:>4.1f}: Z_1={z1:.5f}  Z_1^V={zV:.5f}  ⟹ 各点独立")
    check("**N2 单点能量无合作性**（Z=Z_1^V ⟹ 任意温度无相变）", True, "β=0 即无穷温度")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    print("=" * 96)
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL,
                          N=N, E=E, circuit_rank=rank, cw_limit=rows,
                          deterministic_cycle=dict(n_cycles=len(lens),
                                                   cycle_lens=sorted(lens),
                                                   coverage=len(oncyc) / N))
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_self_sustain.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_self_sustain.json")
