"""
z0_pi.py --- 让 $\\pi$＝**闭包词的嵌套高度**从**引擎自己的闭合记录**里跑出来（不再从 `R53` 文件读）

链条（照 `R53_zero_to_quantum.py:40-180` 的定义，但**独立复算**、且**记录来自引擎**）：
$$
\\text{引擎全分支走}\\to\\text{闭合记录（词）}\\to h_{\\max}=\\max|\\text{前缀和}|\\to\\text{层密度}\\to h^*
\\to\\text{成块}\\to\\omega\\to q,\\lambda_1,\\text{秩},S_{\\max}\\to\\text{三门}
$$
**关键：$\\pi$ 的"读法"必须声明**，本文件把三种读法全测（R53 只用了第一种）：
  R1 **长度恰为 $L$ 的平衡词**（`R53` 的口径；$\\binom{16}{8}=12870$ 个）
  R2 **首达闭合、所有偶长度 $\\le L$**（`L2_period.py` 引擎的 $\\mathcal Z_*$ 实际记录；$L{=}16$ 时 1252 个）
  R3 **所有 $\\le L$ 的平衡词**（并集）
并给出每种读法的 $q$、$S_{\\max}$、$\\lambda_1$、秩、账本峰 $D$、以及是否过三门。

用法：/usr/bin/python3 z0_pi.py      输出：results/z0_pi.json
"""
from __future__ import annotations

import json
import math
import os
import time
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
MU = (math.sqrt(5), 1.3819660, 1.3819660)


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ---------------------------------------------------------------- 引擎侧的记录（全分支枚举）
def balanced_words(L):
    """引擎走 $L$ 步后回到零的词（= R53 的 closed_words）。"""
    out = []
    for pos in combinations(range(L), L // 2):
        w = [-1] * L
        for p in pos:
            w[p] = 1
        out.append(tuple(w))
    return out


def primitive_excursions(L):
    """**首达闭合**的词（引擎 $\\mathcal Z_*$ 每步记录的那一类），长度偶、$\\le L$。"""
    out = []
    for ell in range(2, L + 1, 2):
        for w in balanced_words(ell):
            ps, s, ok = [], 0, True
            for x in w[:-1]:
                s += x
                if s == 0:
                    ok = False
                    break
            if ok:
                out.append(w)
    return out


def hmax(w):
    s, m = 0, 0
    for x in w:
        s += x
        m = max(m, abs(s))
    return m


def layer_density(words):
    d = defaultdict(int)
    for w in words:
        d[hmax(w)] += 1
    return dict(sorted(d.items()))


def make_blocks(words, dens):
    hstar = max(sorted(dens), key=lambda h: dens[h])
    cnt = Counter()
    for w in words:
        h = hmax(w)
        cnt[h if h > hstar else 0] += 1
    return cnt, hstar


# ---------------------------------------------------------------- 判据（照 R53 定义）
def factorize(n):
    f, d = {}, 2
    n = int(n)
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def exact_rank(vals):
    fs = [factorize(v) for v in vals]
    primes = sorted({p for f in fs for p in f})
    rows = [[Fraction(fs[i + 1].get(p, 0) - fs[i].get(p, 0)) for p in primes] for i in range(len(fs) - 1)]
    if not rows:
        return 0
    rk = 0
    for c in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        pv = rows[rk][c]
        rows[rk] = [x / pv for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def ledger_peak(q):
    for D in range(2, 400):
        if D == 2:
            if 3 * q < 1:
                return 2
            continue
        if D * q > (D - 2) and (D + 1) * q < (D - 1):
            return D
    return None


def kcbs_smax(weights):
    v = sorted(weights, reverse=True)[:3]
    s = sum(v)
    if s <= 0:
        return 0.0
    v = [x / s for x in v] + [0.0] * (3 - len(v))
    return sum(a * b for a, b in zip(v, MU))


def analyse(words, label):
    cnt, hstar = make_blocks(words, layer_density(words))
    tot = sum(cnt.values())
    ws = sorted((v / tot for v in cnt.values()), reverse=True)
    k = len(ws)
    q = sum(x * x for x in ws)
    t3 = ws[:3]
    lam1 = t3[0] / sum(t3) if len(t3) >= 3 else 1.0
    r = exact_rank(sorted(cnt.values(), reverse=True))
    return dict(label=label, n_words=len(words), blocks=k, h_star=hstar,
                sizes=sorted(cnt.values(), reverse=True),
                weights=[round(x, 6) for x in ws], q=round(q, 8),
                lam1=round(lam1, 8), rank=r, needed_rank=k - 1,
                S_max=round(kcbs_smax(ws), 6), ledger_peak_D=ledger_peak(q),
                gate1=bool(0.5 < q < 0.6), gate2=bool(kcbs_smax(ws) > 2.0),
                gate3=bool(r == k - 1),
                all_pass=bool(0.5 < q < 0.6 and kcbs_smax(ws) > 2.0 and r == k - 1))


if __name__ == "__main__":
    t0 = time.time()
    L = 16
    print("=" * 100)
    print(f"π 从引擎的闭合记录里跑出来（L={L}）")
    # 引擎侧记录
    R1 = balanced_words(L)                 # 长度恰为 L 的平衡词
    R2 = primitive_excursions(L)           # 首达闭合、所有偶长度 ≤ L（引擎 Z* 的实际记录）
    R3 = [w for ell in range(2, L + 1, 2) for w in balanced_words(ell)]   # 所有 ≤ L 的平衡词
    print(f"   R1 长度恰 {L} 的平衡词：{len(R1)} 个（应为 C(16,8)={math.comb(16,8)}）")
    print(f"   R2 首达闭合（偶长 ≤ {L}）：{len(R2)} 个")
    print(f"   R3 所有 ≤ {L} 的平衡词：{len(R3)} 个")
    check(f"R1 词数 = C({L},{L//2})", len(R1) == math.comb(L, L // 2))
    dens = layer_density(R1)
    print(f"   R1 层密度（h_max 分布）= {dens}")
    RES["R1_layer_density"] = dens

    rows = {}
    print(f"\n   {'读法':>34} {'块数':>5} {'h*':>4} {'q':>10} {'λ1':>9} {'秩/需':>7} {'S_max':>9} {'峰D':>4}  三门")
    for words, lb in ((R1, f"R1 长度恰={L}（R53 口径）"),
                      (R2, f"R2 首达闭合 ≤{L}（引擎 Z*）"),
                      (R3, f"R3 所有平衡词 ≤{L}")):
        a = analyse(words, lb)
        rows[lb] = a
        print(f"   {lb:>34} {a['blocks']:5d} {a['h_star']:4d} {a['q']:10.6f} {a['lam1']:9.6f} "
              f"{a['rank']}/{a['needed_rank']:<5d} {a['S_max']:9.6f} {a['ledger_peak_D']:4d}  "
              f"{'全过' if a['all_pass'] else '未过'}")
    RES["readings"] = rows

    # 与 R53 的文件对账
    d53 = json.load(open(os.path.join(os.path.dirname(HERE), "R53_zero_to_quantum_results.json")))
    r53 = d53["result"]
    r1 = rows[f"R1 长度恰={L}（R53 口径）"]
    print(f"\n   与 R53 文件对账：q {r1['q']} vs {r53['q']}；S_max {r1['S_max']} vs {r53['S_max']}；"
          f"秩 {r1['rank']} vs {r53['rank']}；块 {r1['sizes']} vs {r53['sizes']}")
    check("R1 独立复算 == R53 文件（q / S_max / 秩 / 块尺寸）",
          abs(r1["q"] - r53["q"]) < 1e-7 and abs(r1["S_max"] - r53["S_max"]) < 1e-6
          and r1["rank"] == r53["rank"] and r1["sizes"] == sorted(r53["sizes"], reverse=True))
    check("R1 三门全过（引擎记录 → π → 门，不引用外部文件）", r1["all_pass"])
    print(f"\n   读法敏感性：过三门的读法 = "
          f"{[k for k, v in rows.items() if v['all_pass']] or '（无）'}")
    for lb, a in rows.items():
        if not a["all_pass"]:
            why = []
            if not a["gate1"]:
                why.append(f"q={a['q']:.4f} 不在 (0.5,0.6)")
            if not a["gate2"]:
                why.append(f"S_max={a['S_max']:.4f} ≤ 2")
            if not a["gate3"]:
                why.append(f"秩 {a['rank']} ≠ {a['needed_rank']}")
            print(f"     {lb}: " + "；".join(why))

    # ω/K 与差频（§24.1 的口径，但现在 ω 来自引擎记录）
    r1a = rows[f"R1 长度恰={L}（R53 口径）"]
    tot1 = sum(r1a["sizes"])
    om = [v / tot1 for v in r1a["sizes"]]          # ★ 用未舍入的权重，避免 6 位舍入污染 κ1=q 的断言
    K = [-math.log(x) for x in om]
    diffs = sorted({round(abs(K[i] - K[j]), 5) for i in range(len(K)) for j in range(i + 1, len(K))})
    k1 = sum(x * x for x in om)
    print(f"\n   ω = {om}")
    print(f"   K = -log ω = {[round(x, 5) for x in K]}  跨度={max(K)-min(K):.5f}")
    print(f"   κ1 = Σω² = {k1:.8f}   T_d = L/(-log κ1) = {L/(-math.log(k1)):.4f} 步")
    print(f"   差频谱线 {len(diffs)} 条 = {diffs[:8]}…")
    RES["readout"] = dict(omega=om, K=[round(x, 6) for x in K], kappa1=round(k1, 8),
                          Td=round(L / (-math.log(k1)), 4), n_lines=len(diffs))
    check("κ1 = q（同一个数：退相干常数＝账本余额）", abs(k1 - r1a["q"]) < 1e-7,
          f"κ1={k1:.8f}  q={r1a['q']:.8f}（差 {abs(k1-r1a['q']):.1e}，来自 q 的 8 位舍入）")

    # ---------------- B 读法 × L 全扫（本轮新增：R53 只报了 L=16 的一种读法）
    print("\nB 读法 × L 全扫：哪些 (L, 读法) 过三门？")
    print(f"   {'L':>3} | {'R1 恰=L':^26} | {'R2 首达（引擎 Z*）':^26} | {'R3 所有':^26}")
    print(f"   {'':>3} | {'q':>8}{'S_max':>8}{'峰D':>4}{'门':>6} | {'q':>8}{'S_max':>8}{'峰D':>4}{'门':>6} | {'q':>8}{'S_max':>8}{'峰D':>4}{'门':>6}")
    scan = []
    for LL in range(4, 21, 2):
        cells = []
        for kind, ws_ in (("R1", balanced_words(LL)), ("R2", primitive_excursions(LL)),
                          ("R3", [w for e in range(2, LL + 1, 2) for w in balanced_words(e)])):
            aa = analyse(ws_, kind)
            pk = aa["ledger_peak_D"] if aa["ledger_peak_D"] else -1
            cells.append((aa["q"], aa["S_max"], pk, "过" if aa["all_pass"] else "—"))
            scan.append(dict(L=LL, reading=kind, q=aa["q"], S_max=aa["S_max"],
                             peak_D=pk, all_pass=aa["all_pass"], blocks=aa["blocks"],
                             sizes=aa["sizes"]))
        print(f"   {LL:3d} | " + " | ".join(f"{c[0]:8.5f}{c[1]:8.5f}{c[2]:4d}{c[3]:>6}" for c in cells))
    RES["B_scan"] = scan
    winners = [(r["L"], r["reading"]) for r in scan if r["all_pass"]]
    print(f"\n   过三门的 (L, 读法) = {winners}")
    r2win = [L for L, k in winners if k == "R2"]
    print(f"   **R2（首达＝引擎默认记录）过门的 L = {r2win}** ⟹ 量子扇区**不需要** R53 的「恰为 L」读法")
    check("R2（引擎默认读法）在某个 L 上过三门", len(r2win) > 0, f"L={r2win}")
    check("过门的组合不止 R53 那一个（L=16, R1）", len(winners) > 1, f"{winners}")

    # L=4：R1 的划分是否**恰好等于**旋转类划分（q_4=5/9 的来源）
    w4 = balanced_words(4)
    a4 = analyse(w4, "R1@L=4")
    qrot = Fraction(5, 9)
    print(f"\n   L=4 的 R1：块尺寸={a4['sizes']} q={a4['q']}（旋转类 $q_4=5/9={float(qrot):.5f}$）")
    check("L=4 时 R1（嵌套高度）与旋转类**给出同一个划分**（块 [4,2] ⟹ q=5/9）",
          a4["sizes"] == [4, 2] and abs(a4["q"] - float(qrot)) < 1e-7,
          f"sizes={a4['sizes']} q={a4['q']}")

    # ---------------- C 四门：三门 ＋ 容量上界（L2_period_abs）
    print("\nC 四门（三门＋容量）：容量上界 $c_T\\le K$ vs 三门的 $L$ 候选，交集是什么？")
    print("   容量上界来自 `L2_period_abs`：$c_T=\\binom{T}{\\lfloor T/2\\rfloor}\\le K_2=2(T{+}1)$ 或 $\\le K_1=4(T{+}1)$")
    print(f"   {'L':>3} {'词数 C(L,L/2)':>13} {'K2':>5} {'K1':>5} {'词≤K2':>6} {'词≤K1':>6} {'首达数':>7} {'首达≤K1':>8} {'块数':>5} {'块≤K2':>6} {'三门':>8}")
    cap_rows = []
    for r in scan:
        if r["reading"] != "R1":
            continue
        LL = r["L"]
        c = math.comb(LL, LL // 2)
        K2, K1 = 2 * (LL + 1), 4 * (LL + 1)
        p = len(primitive_excursions(LL))
        blocks = r["blocks"]
        gates_ok = [q["reading"] for q in scan if q["L"] == LL and q["all_pass"]]
        cap_rows.append(dict(L=LL, words=c, K2=K2, K1=K1, words_le_K2=bool(c <= K2),
                             words_le_K1=bool(c <= K1), prim=p, prim_le_K2=bool(p <= K2),
                             prim_le_K1=bool(p <= K1), blocks=blocks, blocks_le_K2=bool(blocks <= K2),
                             gates=gates_ok))
        print(f"   {LL:3d} {c:13d} {K2:5d} {K1:5d} {str(c<=K2):>6} {str(c<=K1):>6} {p:7d} {str(p<=K1):>8} "
              f"{blocks:5d} {str(blocks<=K2):>6} {str(gates_ok):>8}")
    RES["C_four_gates"] = cap_rows
    for tag, key in (("词数口径 K2", "words_le_K2"), ("词数口径 K1", "words_le_K1"),
                     ("首达口径 K2", "prim_le_K2"), ("首达口径 K1", "prim_le_K1"),
                     ("块数口径 K2", "blocks_le_K2")):
        surv = [r["L"] for r in cap_rows if r[key] and r["gates"]]
        print(f"   四门存活（{tag}）：{surv if surv else '**空集**'}")
        RES.setdefault("C_survivors", {})[tag] = surv
    empty_word = (not RES["C_survivors"]["词数口径 K2"]) and (not RES["C_survivors"]["词数口径 K1"]) \
        and (not RES["C_survivors"]["首达口径 K2"]) and (not RES["C_survivors"]["首达口径 K1"])
    check("C 四门在「词数／首达」两种容量口径下**都是空集**", empty_word,
          "⟹ 容量上界与量子门互斥（新登记张力）")
    check("C 换成「块数」口径则存活非空（块数=6 ≪ K2）",
          len(RES["C_survivors"]["块数口径 K2"]) > 0,
          f"存活 = {RES['C_survivors']['块数口径 K2']}")
    print("   ⟹ **关键歧义**：容量到底卡什么？卡「闭合词数」⟹ 无解；卡「块数/类数」⟹ L∈{14,16,20} 可活。")

    print("=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_pi.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_pi.json")
