"""
z0_skeleton_audit.py --- 骨架审计：零和层的【结构】层判定（总账 §1 的 7、9、11、12 与 §3 的否证 A/B/C）

把讨论段临时验证过的七件事固化为可复现机器：
  S1  **$D_V$ ＝ $A_{V-1}$ 根格**（秩 $V{-}1$、全根 $\lvert d\rvert^2=2$）
  S2  **原点唯一**：$U_c(n)=n+c$ 把零和盒双射到自身 ⟹ 只有 $c=0$
  S3  **$L_{\\rm site}$ 不自由**：$k_{\\max}=\\lfloor VL_{\\rm site}/2\\rfloor$
  S4  **层大小只依赖 $(V,L)$**（与 $\\Gamma$ 无关）
  S5  **层序不存在**：可达预序平凡（全部互达），含非传递图
  S6  **无单调势**：$d$ 与 $-d$ 都允许 ⟹ 严格单增 $f$ 不存在
  S7  **$\\Gamma$ 不含取向**：$A=A^{\\sf T}$，取向自由度 $=\\lvert E\\rvert$
  S8  **纤维周期-2**：$P$ 含本征值 $\\lambda=-1$

用法：/usr/bin/python3 z0_skeleton_audit.py   输出：results/z0_skeleton_audit.json
"""
from __future__ import annotations

import json
import os
import time
from collections import Counter, defaultdict, deque
from itertools import product as iproduct

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


def nxt_of(A, L):
    Lc, st = Z.build(A, L, True)
    Ac = Lc.tocoo()
    nxt = defaultdict(list)
    for x, y, v in zip(Ac.row, Ac.col, Ac.data):
        if x != y:
            nxt[x].append(y)
    return Lc, st, nxt


def reach_from(nxt, s):
    seen = {s}; dq = [s]
    while dq:
        u = dq.pop()
        for w in nxt[u]:
            if w not in seen:
                seen.add(w); dq.append(w)
    return seen


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 96)
    print("S1：D_V = {e_j − e_i} 是 A_{V−1} 根格吗？")
    print("=" * 96)
    s1 = True
    for V in (4, 6, 8):
        D = np.array([np.eye(V)[j] - np.eye(V)[i]
                      for i in range(V) for j in range(V) if i != j])
        rank = np.linalg.matrix_rank(D)
        norms = sorted(set(np.round((D ** 2).sum(1), 9)))
        sig = np.allclose(D.sum(1), 0)
        ok = (rank == V - 1) and (norms == [2.0]) and sig
        s1 &= ok
        print(f"  V={V}: |D|={len(D)}  秩={rank}(=V−1)  |d|²={norms}  σ(d)=0? {sig}  ⟹ {ok}")
    check("**S1 $D_V$ 是 $A_{V-1}$ 根格**（秩 $V{-}1$、全根、$\\sigma(d)=0$）", s1, "V=4,6,8")

    print("\n" + "=" * 96)
    print("S2：零和盒的平移不变原点唯一吗？")
    print("=" * 96)
    s2 = True
    for kind, V, L in (("ring", 6, 1), ("ring", 4, 1), ("ring", 4, 2),
                       ("ring", 6, 2), ("ring", 4, 3), ("K", 4, 2)):
        A = S.G(V, kind)
        Lc, st, _ = nxt_of(A, L)
        legal = []
        for c in st:
            if all(all(abs(n[v] + c[v]) <= L for v in range(V)) for n in st):
                legal.append(c)
        ok = (len(legal) == 1 and legal[0] == tuple([0] * V))
        s2 &= ok
        print(f"  {kind}{V} L={L}: 合法原点 {len(legal)} 个  ⟹ 只有真空? {ok}")
    check("**S2 原点唯一**（$U_c$ 保盒 ⟹ 只有 $c=0$）", s2, "六例")

    print("\n" + "=" * 96)
    print("S3/S4：$L_{\\rm site}$ 不自由 + 层大小只依赖 $(V,L)$")
    print("=" * 96)
    s3 = s4 = True
    for kind, V, L in (("ring", 6, 1), ("ring", 6, 2), ("ring", 4, 3),
                       ("K", 4, 2), ("ring", 8, 1), ("path", 6, 1), ("ring", 8, 2)):
        A = S.G(V, kind)
        Lc, st, _ = nxt_of(A, L)
        lvl = [sum(abs(t) for t in s) // 2 for s in st]
        kmax = max(lvl)
        ok3 = (kmax == (V * L) // 2)
        s3 &= ok3
        c = Counter(lvl); sizes = [c[k] for k in sorted(c)]
        pred = Counter()
        for n in iproduct(range(-L, L + 1), repeat=V):
            if sum(n) == 0:
                pred[sum(abs(x) for x in n) // 2] += 1
        psizes = [pred[k] for k in sorted(pred)]
        ok4 = (sizes == psizes)
        s4 &= ok4
        print(f"  {kind}{V} L={L}: k_max={kmax}(预测{(V*L)//2}) {ok3}   层大小={sizes} 与组合预测一致? {ok4}")
    check("**S3 $L_{\\rm site}=2k_{\\max}/V$**（界被取到）", s3, "七例")
    check("**S4 层大小只依赖 $(V,L)$**（与 $\\Gamma$ 无关；ring6 与 path6 相同）", s4, "七例")

    print("\n" + "=" * 96)
    print("S5：层序存在吗？（可达预序）")
    print("=" * 96)
    def hub(V):
        A = S.G(V, "ring")
        for j in range(V):
            if A[0, j] == 0 and j != 0:
                A[0, j] = A[j, 0] = 1.0
        return A
    s5 = True
    for tag, A in (("ring6", S.G(6, "ring")), ("path6", S.G(6, "path")),
                   ("ring8", S.G(8, "ring")), ("K4L2", S.G(4, "K")),
                   ("ring6+hub(非传递)", hub(6))):
        L = 2 if tag == "K4L2" else 1
        Lc, st, nxt = nxt_of(np.asarray(A, float), L)
        N = len(st)
        lvl = np.array([sum(abs(t) for t in s) // 2 for s in st])
        ks = sorted(set(lvl.tolist()))
        reach = {}
        for k in ks:
            seen = set()
            for i in range(N):
                if lvl[i] == k:
                    seen |= reach_from(nxt, i)
            reach[k] = set(lvl[list(seen)].tolist())
        mutual = sum(1 for k in ks for j in ks if k != j and k in reach[j] and j in reach[k])
        tot = len(ks) * (len(ks) - 1)
        ok = (mutual == tot)
        s5 &= ok
        print(f"  {tag}: 层={ks}  互达={mutual}/{tot}  ⟹ 预序平凡? {ok}")
    check("**S5 层序不存在**（全部层互达，含非传递 $\\Gamma$）", s5, "五例")

    print("\n" + "=" * 96)
    print("S6/S7：取向导不出来（无单调势 + $\\Gamma$ 无向）")
    print("=" * 96)
    s6 = s7 = True
    for kind, V in (("ring", 6), ("path", 6), ("ring", 8), ("K", 4)):
        A = np.asarray(S.G(V, kind), float)
        D = set()
        for u in range(V):
            for w in np.nonzero(A[u])[0]:
                d = [0] * V; d[u] -= 1; d[w] += 1
                D.add(tuple(d))
        negclosed = all(tuple(-np.array(d)) in D for d in D)
        both = [d for d in D if tuple(-np.array(d)) in D]
        sym = bool(np.allclose(A, A.T))
        ok6 = negclosed and (len(both) == len(D))
        s6 &= ok6
        s7 &= sym
        half = len(D) // 2
        print(f"  {kind}{V}: |D|={len(D)} 取负封闭? {negclosed}  互逆对={len(both)//2}/{half}"
              f"  ⟹ 单调势不存在? {ok6}   A=Aᵀ? {sym}  取向自由度={int(A.sum())//2}")
    check("**S6 无单调势**（$d$ 与 $-d$ 都允许 ⟹ 严格单增 $f$ 不存在）", s6, "四例")
    check("**S7 $\\Gamma$ 无向**（$A=A^{\\sf T}$ ⟹ 不含取向；取向自由度 $=\\lvert E\\rvert$）", s7, "四例")

    print("\n" + "=" * 96)
    print("S8：纤维周期-2（$P$ 含 $\\lambda=-1$）")
    print("=" * 96)
    s8 = True
    for kind, V, L in (("ring", 6, 1), ("path", 6, 1), ("ring", 8, 1), ("ring", 4, 3)):
        A = S.G(V, kind)
        Lc, st, nxt = nxt_of(A, L)
        N = len(st)
        P = np.zeros((N, N))
        for x in range(N):
            d = len(nxt[x])
            if d:
                for y in nxt[x]:
                    P[y, x] = 1.0 / d
        ev = np.linalg.eigvals(P)
        has = bool(np.any(np.abs(ev + 1) < 1e-8))
        s8 &= has
        print(f"  {kind}{V} L={L}: N={N}  P 含 λ=−1? {has}")
    check("**S8 纤维周期-2**（$P$ 含本征值 $-1$ ⟹ 两相位交替）", s8, "四例")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_skeleton_audit.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_skeleton_audit.json")
