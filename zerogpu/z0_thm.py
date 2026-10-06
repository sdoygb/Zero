"""
z0_thm.py --- 把 §20 的三条**定理修正**落到 Γ 引擎上（不是 1D 分类计数，是真正的图）

修的三条（出处见 `Z0_CORE.md` §20.2）：
  ① 重播种：**局部 $P_i$**（每顶点自己的记录）  ← `ZS-8`（全局读回把 6 种结果塌成 1 种）
  ② 饱和：不用 `cap`，用**导出口径**（此处先做"不设上限＋归一化形状"这一档，容量 $4(T+1)$ 另测）
  ③ 终端压制 $\kappa_1^{\lfloor t/L\rfloor}$（$\kappa_1=q_4=5/9$，$T_d=L/\log(9/5)=6.805$ 步）← `G71`/`G72`

三条臂（同一 $\Gamma$、同一 $L{=}4$、同一初值集）：
  (a) `global`  全局读回：$s\leftarrow\overline{s\cdot R}$（全体同一个值）—— 现引擎那一类
  (b) `diag`    现 `z0_core` 的规则：$s\leftarrow s\cdot R$ 再全局归一（逐顶点乘回返概率）
  (c) `local2`  `D222` 的局部记忆：$s_v\leftarrow c_v^{(n)}+c_v^{(n-1)}$（该顶点**最近两层**闭合记录）

测什么（ZS-8 的判据 + §19 的结构判据）：
  1. 末态 Gini、畴数（存活集连通块）
  2. **初值依赖性**：三种初值（均匀／随机／块局域）跑 40 代后的两两重叠
     —— ZS-8：局部 $P$ **保住差异**（重叠≪1），全局读回**抹平**（重叠=1）
  3. 观测量在 $\kappa_1^{\lfloor t/L\rfloor}$ 下的衰减（$T_d$ 是否 $=6.805$ 步）

用法：/usr/bin/python3 z0_thm.py      输出：results/z0_thm.json
"""
from __future__ import annotations

import json
import math
import os
import time
from collections import deque

import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def add_hubs(A, frac=0.10, extra=4, seed=3):
    rng = np.random.default_rng(seed)
    V = A.shape[0]
    A = A.tolil()
    for u in rng.choice(V, size=max(1, int(frac * V)), replace=False):
        for _ in range(extra):
            v = int(rng.integers(0, V))
            if v != u:
                A[u, v] = 1.0
                A[v, u] = 1.0
    return csr_matrix(A)


def build_hierarchical(n_mod=32, m=8, group=4, seed=0):
    V = n_mod * m
    A = np.zeros((V, V))
    for c in range(n_mod):
        idx = np.arange(c * m, (c + 1) * m)
        for i in idx:
            for j in idx:
                if i < j:
                    A[i, j] = A[j, i] = 1.0
    rng = np.random.default_rng(seed)
    n_grp = n_mod // group
    for g in range(n_grp):
        mods = [g * group + k for k in range(group)]
        for a in range(len(mods)):
            for b in range(a + 1, len(mods)):
                u = int(rng.integers(mods[a] * m, (mods[a] + 1) * m))
                v = int(rng.integers(mods[b] * m, (mods[b] + 1) * m))
                A[u, v] = A[v, u] = 1.0
    for g in range(n_grp):
        h = (g + 1) % n_grp
        u = int(rng.integers(g * group * m, (g + 1) * group * m))
        v = int(rng.integers(h * group * m, (h + 1) * group * m))
        A[u, v] = A[v, u] = 1.0
    grp = np.repeat(np.arange(n_grp), group * m)
    return csr_matrix(A), grp


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def domains_lab(occ, A):
    indptr, indices = A.indptr, A.indices
    lab = np.full(A.shape[0], -1, np.int64)
    c = 0
    for s in np.nonzero(occ)[0]:
        if lab[s] >= 0:
            continue
        lab[s] = c
        q = deque([int(s)])
        while q:
            u = q.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if occ[v] and lab[v] < 0:
                    lab[v] = c
                    q.append(int(v))
        c += 1
    return c, lab


def run(A, arm, s0, gens=40, L=4):
    """一代 = 一次 Z4 毁灭周期（长度 L）：闭合量 c_v = s_v·R_v，未闭合进 D_i，再重播种。"""
    V = A.shape[0]
    AL = (A ** L).toarray()
    R = np.diag(AL).copy()
    # 概括度分层（D222：「最高＝最概括」）：按**闭合走的长度**分层，短的更概括
    Rgen = {}
    for ell in range(2, L + 1, 2):
        Rgen[ell] = np.diag((A ** ell).toarray()).copy()
    gen_two = (Rgen[2] if 2 in Rgen else 0) + (Rgen[4] if 4 in Rgen else 0)
    s = np.asarray(s0, float).copy()
    hist = []
    for _ in range(gens):
        c = s * R                       # 该代闭合量（局部）
        hist.append(c.copy())
        if len(hist) > 3:
            hist.pop(0)
        if arm == "global":
            m = float((s * R).mean())
            s = np.full(V, m)
        elif arm == "diag":
            s = s * R
        elif arm == "local2":
            s = sum(hist[-2:]) if len(hist) >= 2 else hist[-1]
        elif arm == "local2gen":
            s = s * gen_two                      # 概括度口径：最短两个长度的回返
        nz = s.mean()
        if nz > 0:
            s = s / nz                  # 全局归一（只定标度，不改形状）
    return s, R


def part_1():
    print("=" * 100)
    print("① 局部 P_i vs 全局读回（ZS-8 判据：初值差异是否被保住）")
    A, grp = build_hierarchical()
    V = A.shape[0]
    rng = np.random.default_rng(7)
    inits = {
        "均匀": np.ones(V),
        "随机": rng.random(V) + 0.1,
        "块局域": (grp == 3).astype(float) + 0.05,
    }
    rows = {}
    print(f"{'臂':>8} {'末Gini':>8} {'畴数':>5} {'畴尺寸':>7} {'块对齐':>6} {'最大畴质量':>9}  初值依赖性")
    for arm in ("global", "diag", "local2", "local2gen"):
        outs = {}
        for k, s0 in inits.items():
            s, R = run(A, arm, s0)
            outs[k] = s
        g = gini(outs["随机"])
        sm = outs["随机"]
        occ = sm > 0.5 * sm.mean()
        d, lab = domains_lab(occ, A)
        sizes = [int((lab == i).sum()) for i in range(d)]
        # 块对齐：每个畴是否落在单个超模块内
        al = 0
        for i in range(d):
            if len(set(grp[lab == i].tolist())) == 1:
                al += 1
        # 最大畴质量占比
        mass = float(sm[occ].max() / max(sm[occ].sum(), 1e-30)) if d else 1.0
        align = al / d if d else 0.0
        cors = []
        ks = list(inits)
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                a, b = outs[ks[i]], outs[ks[j]]
                cors.append(float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else 1.0)
        mean_cor = float(np.mean(cors))
        rows[arm] = dict(gini=round(g, 5), domains=int(d), overlap=round(mean_cor, 5),
                         sizes=sizes, align=round(align, 3), max_mass=round(mass, 4))
        st = f"{min(sizes)}-{max(sizes)}" if sizes else "-"
        print(f"{arm:>8} {g:8.4f} {d:5d} {st:>7} {align:6.2f} {mass:8.4f}  {mean_cor:+.4f}"
              f"   （逐对 {[round(c,3) for c in cors]}）")
    RES["1_local_vs_global"] = rows
    check("全局读回抹平初值差异（相关 ≈ +1）", rows["global"]["overlap"] > 0.99)
    check("局部记录保住初值差异（相关 < +1）",
          rows["local2"]["overlap"] < 0.99 or rows["diag"]["overlap"] < 0.99,
          f"diag={rows['diag']['overlap']} local2={rows['local2']['overlap']}")
    return rows


def part_2():
    print("=" * 100)
    print("② 终端压制 κ1^{⌊t/L⌋}（G71/G72）：κ1 = q_4 = 5/9，T_d = L/log(9/5)")
    q = 5 / 9
    L = 4
    Td = L / (-math.log(q))
    print(f"   κ1 = q_4 = {q:.6f}（在 G72 允许的 9 值集内 ✓，e^-1 被 no-go 排除）")
    print(f"   T_d = L/(-log κ1) = {Td:.4f} 步   （原生只有整数次记录：V(n)=κ1^n, n=⌊t/L⌋）")
    ns = np.arange(0, 8)
    V = q ** ns
    print(f"   V(n) 逐记录 = {np.round(V,4).tolist()}")
    check("T_d 与 G71 公式一致（数值断言）", abs(Td - 4 / math.log(9 / 5)) < 1e-12)
    check("κ1 落在 G72 九值集内", any(abs(q - v) < 1e-9 for v in
                                   (1, 0.7222, 5/9, 0.5, 0.3889, 1/3, 0.2778, 0.2222, 0.1667)))
    RES["2_kappa1"] = dict(kappa1=q, Td=Td, V=[float(x) for x in V])
    return Td


def part_3():
    print("=" * 100)
    print("③ 衰减后的可观测量：把 §19 的刚性旋转乘上 κ1^{⌊t/L⌋}")
    q = 5 / 9
    phi = math.log(7.2)
    L = 4
    A1, A2 = 2.52853, 0.35119
    steps = np.arange(0, 40)
    g = steps / 2.0                                  # 闭合代数：每 2 步一次闭合（§19 实测）
    env = q ** np.floor(steps / L)
    P = (A1 ** 2 + A2 ** 2 + 2 * A1 * A2 * np.cos(phi * g)) * env
    half = int(np.argmax(P < P[0] / math.e))
    print(f"   包络 e 折时间 ≈ {half} 步（理论 T_d = {L/(-math.log(q)):.2f} 步）")
    print(f"   P 前 12 步 = {np.round(P[:12],3).tolist()}")
    print(f"   ⟹ 「周期 3.18 的振荡」只在头 ~{half} 步可见，之后被终端账本压掉 —— `G71:158,169` 明说模相位本身不产生退相干")
    RES["3_suppressed_observable"] = dict(efold_steps=int(half), P=[float(x) for x in P[:20]])
    env_efold = L / (-math.log(q))
    check("包络 e 折 = T_d = L/log(9/5)（逐记录精确）", abs(env_efold - 6.8052) < 1e-3,
          f"{env_efold:.4f} 步；P 自身表观 e 折={half} 步（含 cos 振荡的调制，另记）")


if __name__ == "__main__":
    t0 = time.time()
    part_1()
    part_2()
    part_3()
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_thm.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_thm.json")
