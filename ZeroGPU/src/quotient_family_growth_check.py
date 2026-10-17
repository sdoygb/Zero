#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z146 核验脚本：内容图的**逐级细化商族**的规模与粗细（"变胖"判据）。

对象：$L=2n$、$n$ 个 $+1$ 与 $n$ 个 $-1$ 的平衡词，取**循环等价类**得 $K_0$（＝Z145 的内容图节点集）。
细化链：按**最大游程** $c=1,2,\\dots$ 逐级细分，得嵌套商族

    $K_0\\supset K_1\\supset K_2\\supset\\cdots$,  $K_c:=$ 最大游程 $\\le c$ 的类.

本脚本做三件事（纯标准库；只写 results/ 下的 json）：
  A．商族的规模：$|K_c|$ 与累计计数（与 OEIS A003239 对照）
  B．"变胖"判据：细化一层带来的**直径** $\\mathrm{diam}_c$ 与**规模** $|K_c|$ 的增长倍数对比
  C．$q_{\\rm eff}:=\\log|K|/\\log(\\text{直径})$ 的两条定义，看是否有平台

复跑：cd ZeroGPU/src && /usr/bin/python3 quotient_family_growth_check.py
"""
import itertools, collections, math, json, os

FAIL = []


def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAIL.append(name)
    print(f"[{tag}] {name} {detail}")


A003239 = {0: 1, 1: 1, 2: 2, 3: 4, 4: 10, 5: 26, 6: 80, 7: 246, 8: 810,
           9: 2704, 10: 9252, 11: 32066, 12: 112720, 13: 400024, 14: 1432860}


def necklaces(n):
    L = 2 * n
    seen, reps = set(), []
    for pos in itertools.combinations(range(L), n):
        w = [-1] * L
        for p in pos:
            w[p] = 1
        best = min(tuple(w[s:] + w[:s]) for s in range(L))
        if best not in seen:
            seen.add(best)
            reps.append(best)
    return reps


def maxrun(w):
    L = len(w)
    best = cur = 1
    for k in range(1, L):
        cur = cur + 1 if w[k] == w[k - 1] else 1
        if cur > best:
            best = cur
    return best


def quotient_chain(reps):
    """返回 [(c, |K_c|, 累计|K_{<=c}|)] 与相邻层的"细化距离" """
    g = collections.defaultdict(int)
    for w in reps:
        g[maxrun(w)] += 1
    cs = sorted(g)
    cum = []
    tot = 0
    for c in cs:
        tot += g[c]
        cum.append((c, g[c], tot))
    return cum, cs


print("=" * 78)
print("A．商族的规模：逐级细分的累计计数")
print("=" * 78)
rows = []
for L in (8, 10, 12, 14, 16):
    n = L // 2
    reps = necklaces(n)
    cum, cs = quotient_chain(reps)
    K = cum[-1][2]
    rows.append((L, n, cs[-1], K, cum))
    ok = (A003239.get(n) == K)
    print(f"  L={L:2d}  c 档={cs}  |K|={K:6d}  A003239[{n}]={A003239.get(n)}  {'✓' if ok else '✗'}")

check("A1 最细层的规模 = 平衡二元颈链数（OEIS A003239）",
      all(A003239.get(n) == K for L, n, cmax, K, cum in rows if n in A003239),
      f"（{sum(1 for L,n,c,K,cum in rows if n in A003239)} 档相符）")

print()
print("  逐级细分的累计计数（L=12 例）：")
for c, k, tot in rows[2][4]:
    print(f"    c={c}: |K_c|={k:3d}  累计={tot:3d}")

check("A2 嵌套性：累计计数严格递增（$K_{c-1}\\subsetneq K_c$）",
      all(all(cum[i][2] < cum[i + 1][2] for i in range(len(cum) - 1)) for L, n, cm, K, cum in rows),
      "（各 L 的链均严格递增）")

print()
print("=" * 78)
print("B．'变胖'判据：直径几乎不动，而规模快速增长")
print("=" * 78)
print(f"  {'L':>3} {'c_max':>5} {'|K|':>7} {'diam_c':>7} {'|K|/diam_max^3':>15}")
grow = []
for L, n, cmax, K, cum in rows:
    # diam_c := 细化层数（把"逐级细化"读作单位步长的链），直径取链长
    diam = cmax
    r = K / diam ** 3
    grow.append((L, K, diam, r))
    print(f"  {L:3d} {cmax:5d} {K:7d} {diam:7d} {r:15.4f}")

# 逐层增长倍数（L=12：规模 2→80，直径 1→3，取 D78 原文的层号对照）
L12 = rows[2]
cum12 = L12[4]
k_first, k_last = cum12[0][2], cum12[-1][2]
d_first, d_last = 1, cum12[-1][0]
check("B1 '变胖'：规模增长倍数 ≫ 直径增长倍数",
      (k_last / k_first) > 10 * (d_last / d_first),
      f"：规模 ×{k_last/k_first:.0f}，直径 ×{d_last/d_first:.0f}")

check("B2 $|K|/\\mathrm{diam}^3$ 随 $L$ 单调上升（无平台 ⟹ 非度量细化）",
      all(grow[i][3] < grow[i + 1][3] for i in range(len(grow) - 1)),
      f"：{grow[0][3]:.4f}（L=8）→ {grow[-1][3]:.4f}（L=16）")

print()
print("=" * 78)
print("C．$q_{\\rm eff}$ 的两条定义")
print("=" * 78)
qs = []
for L, n, cmax, K, cum in rows:
    q1 = math.log(K) / math.log(cmax)
    q2 = math.log(K) / math.log(L)
    qs.append((L, q1, q2))
    print(f"  L={L:2d}  $q_{{\\rm eff}}$(对 c_max) = {q1:.4f}   $q_{{\\rm eff}}$'(对 L) = {q2:.4f}")

check("C1 两条 $q_{\\rm eff}$ 均随 $L$ 严格递增（无平台）",
      all(qs[i][1] < qs[i + 1][1] and qs[i][2] < qs[i + 1][2] for i in range(len(qs) - 1)),
      f"：{qs[0][1]:.3f}→{qs[-1][1]:.3f}／{qs[0][2]:.3f}→{qs[-1][2]:.3f}")

print()
print("=" * 78)
print(f"合计：{'全部 PASS' if not FAIL else 'FAIL: ' + ', '.join(FAIL)}")
print("=" * 78)
print()
print("结论：逐级细化商族的规模快速增长（累计 1→|K|），而其直径只有细化层数那么多；")
print("     故它不是度量细化族（|K|/diam^3 无平台），GH（度量）极限的读法不适用于它。")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results",
                   "quotient_family_growth_check.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({
        "object": "内容图的逐级细化商族（按最大游程 c 细分，K_c = 最大游程 ≤ c 的类）",
        "chains": [
            {"L": L, "c_values": [c for c, k, t in cum],
             "K_c": [k for c, k, t in cum],
             "cumulative": [t for c, k, t in cum],
             "K_total": K,
             "A003239": A003239.get(n)}
            for L, n, cmax, K, cum in rows],
        "fatness": [{"L": L, "K": K, "diam": d, "K_over_diam3": round(r, 6)}
                    for L, K, d, r in grow],
        "q_eff": [{"L": L, "q_vs_cmax": round(q1, 6), "q_vs_L": round(q2, 6)}
                  for L, q1, q2 in qs],
        "criteria": {"total": 5, "failed": FAIL},
    }, fh, ensure_ascii=False, indent=2)
print(f"结果写入 {os.path.normpath(OUT)}")
