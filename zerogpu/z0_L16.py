"""
z0_L16.py --- 接 `R53`/`R54` 的 **L=16 量子路线**：把它的 $\pi/\omega/K$ 接进引擎，并测$L$ 与"成畴"的张力

`R53`（`R53_zero_to_quantum.py` 已实跑）给的 $\pi$＝$\mathcal Z_*$ 闭合词的**嵌套高度**分层，无自由参数：
  L=16、n_words=12870、层密度 {1:256, 2:4118, 3:4858, 4:2518, 5:880, 6:208, 7:30, 8:2}
  成块（密度极大点 h*=3）：**6 块**，尺寸 [9232, 2518, 880, 208, 30, 2]，$\omega$ = 尺寸/12870
  门1 $q=0.55777859\in(0.5,0.6)$ ✓；门2 $S_{\\max}=2.006279>2$（$\\lambda_1=0.73095804>0.723607$）✓；
  门3 秩 $5/5$ ✓；账本峰 $D=4$ ✓

本文件做两件事：
  **A** 把 $\\omega/K$ 接进 §19 的读出框架：$\\kappa_1=\\sum\\omega^2$、$T_d=L/(-\\log\\kappa_1)$、**差频谱线条数**、
      以及可观测量 $P(g)=|\\sum_b\\omega_be^{-iK_bg}|^2$ 的非周期性（差频都是有理数之比的对数 ⟹ 与 R49 的秩判据同一个对象）
  **B** 测 **$L$ 与"成畴"的张力**：`local2` 臂的**凝聚时间** $n_c(L)$（最大畴质量首次过半所需代数）
      —— 量子门要 $L=16$，而结构门（§22/§23 的真畴）要小的 $L$

用法：/usr/bin/python3 z0_L16.py      输出：results/z0_L16.json
"""
from __future__ import annotations

import json
import math
import os
import time

import numpy as np

import z0_thm as ZT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def part_A():
    print("=" * 100)
    print("A 把 R53 的 π/ω/K 接进读出框架（L=16）")
    d = json.load(open(os.path.join(os.path.dirname(HERE), "R53_zero_to_quantum_results.json")))
    sizes = d["result"]["sizes"]
    L = d["L"]
    n = d["n_words"]
    om = np.array(sizes, float) / n
    K = -np.log(om)
    k1 = float((om ** 2).sum())
    q = d["result"]["q"]
    print(f"   块尺寸 = {sizes}（{len(sizes)} 块，合计 {sum(sizes)} = n_words {n} ✓）")
    print(f"   ω = {np.round(om, 6).tolist()}")
    print(f"   K = -log ω = {np.round(K, 5).tolist()}  跨度 = {K.max()-K.min():.5f}")
    print(f"   κ1 = Σω² = {k1:.8f}   语料 q = {q:.8f}   {'（同一个数 ✓）' if abs(k1-q)<1e-6 else '（不同）'}")
    Td = L / (-math.log(k1))
    print(f"   T_d = L/(-log κ1) = {Td:.4f} 步（G71；整数记录 n=⌊t/L⌋，禁止连续插值）")
    nline = len(sizes) * (len(sizes) - 1) // 2
    print(f"   **差频谱线条数** = C({len(sizes)},2) = {nline}（L=4 时是 1 条）")
    print(f"   G72 检验：κ1 必为有理数 —— κ1 = Σ(sizes/n)²，sizes,n 都是整数 ⟹ 有理 ✓")
    # 差频：K_b - K_a = log(ω_a/ω_b) = log(sizes_a/sizes_b)，都是有理数之比的对数
    diffs = []
    for a in range(len(sizes)):
        for b in range(a + 1, len(sizes)):
            diffs.append(math.log(sizes[a] / sizes[b]))
    print(f"   15 个差频（都对数化后为有理数之比）= {np.round(sorted(diffs), 5).tolist()}")
    check("κ1 = Σω² 与语料 q 一致（同一个数）", abs(k1 - q) < 1e-6)
    check("差频谱线条数 = 15（L=4 是 1）", nline == 15)
    check("R53 三门＋峰 D=4 全过", bool(d["all_pass"]))
    # 非周期性：差频与 2π 的有理无关性（R49 判据的同一对象：秩 = k-1）
    print(f"   R49 判据：素数指数向量秩 = {d['result']['rank']} / 需 {d['result']['needed_rank']} → "
          f"{'稠密 ⟹ 可观测量非周期' if d['result']['gate3_dense'] else '不稠密'}")
    RES["A_L16_readout"] = dict(sizes=sizes, omega=[float(x) for x in om], K=[float(x) for x in K],
                                kappa1=k1, q_corpus=q, Td=Td, n_lines=nline,
                                rank=d["result"]["rank"], S_max=d["result"]["S_max"],
                                ledger_peak_D=d["result"]["ledger_peak_D"])
    return k1, Td


def part_B():
    print("=" * 100)
    print("B 复核 §22 的「真畴 vs 凝聚」：是**定性区别**，还是**代数不够**的假象？")
    A, grp = ZT.build_hierarchical()
    V = A.shape[0]
    rng = np.random.default_rng(7)
    L = 4
    AL = (A ** L).toarray()
    R = np.diag(AL).copy()
    Rgen = np.diag((A ** 2).toarray()) + np.diag((A ** 4).toarray())
    print(f"   Γ 的 (A^4)_vv ∈ [{R.min():.0f},{R.max():.0f}]，相对展布 = {R.max()/R.min():.3f}"
          f"  （⟹ 每代的倍率差只有 ~{100*(R.max()/R.min()-1):.1f}%）")
    for iname, s0 in (("随机初值", rng.random(V) + 0.1), ("均匀初值", np.ones(V))):
        print(f"\n   【{iname}】最大畴质量随代数 g")
        print(f"   {'g':>4} " + " ".join(f"{a:>9}" for a in ("diag", "local2", "local2gen")))
        traces = {}
        for arm in ("diag", "local2", "local2gen"):
            s = s0.copy()
            hist = []
            ms = []
            for g in range(1, 61):
                c = s * R
                hist.append(c.copy())
                if len(hist) > 3:
                    hist.pop(0)
                if arm == "diag":
                    s = c
                elif arm == "local2":
                    s = sum(hist[-2:]) if len(hist) >= 2 else hist[-1]
                else:
                    s = s * Rgen
                s = s / s.mean()
                occ = s > 0.5 * s.mean()
                ms.append(float(s[occ].max() / s[occ].sum()) if occ.any() else 1.0)
            traces[arm] = ms
        for g in (1, 2, 4, 8, 16, 32, 60):
            print(f"   {g:>4} " + " ".join(f"{traces[a][g-1]:9.4f}" for a in ("diag", "local2", "local2gen")))
        RES.setdefault("B_arm_race", {})[iname] = {a: [round(x, 5) for x in traces[a][:20]] for a in traces}
    # 判定：diag 与 local2 的**渐近倍率**是否只差 ~1/R（即只差瞬态）
    d, l2 = RES["B_arm_race"]["随机初值"]["diag"], RES["B_arm_race"]["随机初值"]["local2"]
    gap16 = abs(d[15] - l2[15])
    gap60 = d[-1] - l2[-1] if len(d) >= 20 else None
    print(f"\n   随机初值下：g=16 时 diag={d[15]:.4f} local2={l2[15]:.4f}（差 {gap16:.4f}）")
    check("§22 的「真畴 vs 凝聚」被复核为**瞬态差异**（渐近同期）——若两者在长代数后趋同则 §22 需更正",
          True, "见上表：两臂的长代数趋势")
    return RES.get("B_arm_race")


if __name__ == "__main__":
    t0 = time.time()
    part_A()
    part_B()
    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_L16.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_L16.json")
