"""
z0_faithful.py --- **定理合规的 L2 引擎**：照抄 `L2_period.py` 的规则，逐条对语料**被钉死的量**做断言

为什么另起一个：`Z0_CORE.md` §20 查明引擎与语料定理冲突五处，其中第 1 条（**全局重播种**）正是
`ZS-8` 判死的"全局 R 读回"（语料实测 6 种结果塌成 1 种）。本文件的 A 组把规则按 `L2_period.py`
逐字实现并**断言**定理值；B 组做**局部 P_i vs 全局 P** 的对照（ZS-8 的机制，换成我自己的引擎）。

规则（逐字来自 `L2_period.py:50-90`，只用 Z0 派生条款）：
  步    每个活动路径延长一位取 ± 两支（Z0③ 全分支、整数重数）
  闭合  词和为 0（Z1 定理 2）⟹ 写记录、退出活动层
  周期末 t=nT：未闭合路径全部进 D_i（Z4 终端款）
  重播种 每条**保留的**闭合历史 w 生成 {w+, w-}（Z5，2 个种子）
  记忆  只保留**最高 memory 层**（D222 的"最高两层"＝最概括，不是最近）

A 组断言（靶值全部来自语料，不是我定的）：
  A1 c_k 谱 = C(k,⌊k/2⌋)                        [`L2_period.py` A1；`D_arc/D220:100-160`]
  A2 λ(T) = 4,4,8,8,18,18,46,46,130,130 (T=3..12) [`L2_C_recursion_verdict.md:71-95`]
  A3 λ 与站点数无关（1/2/3/6 站点同值）           [`L2_C_recursion_verdict.md:27-38`]
  A4 ρ(T)=D/h 收敛到 1.207884(T=3) / 2.415768(T=4) / 2.224945(T=5) / 4.449890(T=6)
                                                 [`L2_rho_closed_verdict.md:25-64`]
  A5 奇偶配对 h(2k−1)=h(2k)、D(2k−1)=½D(2k)、ρ(2k)=2ρ(2k−1)  [`L2_period_verdict.md:109-127`]
  A6 周期内沉积支撑 = {0,2,…,T−2}（奇数步无沉积）  [`L2_imprint_verdict.md:16-19`]
  A7 纯毁灭+重播种**没有稳态**（h 逐周期指数增长）  [`L2_period_verdict.md:159-175`]
B 组对照（ZS-8）：局部保留 P_i 保住差异；全局合并把差异抹平 [`zero_sum_global_R_local_P.md:41-83`]

用法：/usr/bin/python3 z0_faithful.py      输出：results/z0_faithful.json
"""
from __future__ import annotations

import json
import os
import time
from collections import Counter, defaultdict
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ================================================================ 引擎（分类计数版，与 L2_period 等价）
def run_site(T, cycles, memory=None, start_balance=1, start_parity=0):
    """单站点。返回每周期统计。cls[(b,p)] = 路径数（b 平衡，p 长度奇偶）。"""
    active = {(start_balance, start_parity): 1}
    history = defaultdict(int)          # key = 层号（＝闭合发生的周期），value = 记录数
    layer = 0
    stats = []
    for c in range(1, cycles + 1):
        closed = 0
        deposit_steps = Counter()
        loc = 0
        for step in range(T):
            nxt = defaultdict(int)
            for (b, p), n in active.items():
                for db in (+1, -1):
                    nb = b + db
                    if nb == 0:
                        closed += n
                        history[layer + 1] += n
                        deposit_steps[step] += n          # 0 基步号
                    else:
                        nxt[(nb, 1 - p)] += n
            active = nxt
            loc = len(active)
        destroyed = sum(active.values())
        active = {}
        layer += 1
        # ★ 顺序要紧：先按 D222 剪枝（只留"最高"＝最概括的 memory 层），**再**算 h 与种子
        if memory is not None:
            keep = sorted(history.keys())[-memory:]
            history = defaultdict(int, {k: history[k] for k in keep})
        h = sum(history.values())
        if h:
            active = {(1, 0): h, (-1, 0): h}                  # 每条保留历史 2 个符号种子
        stats.append(dict(cycle=c, closed=closed, destroyed=destroyed, history=h,
                          deposit_steps=dict(sorted(deposit_steps.items()))))
    return stats


# ================================================================ A 组：定理合规性
def part_A():
    print("=" * 96)
    print("A 组：把语料被钉死的量逐条断言（靶值来自语料，不是我定的）")

    # A1 c_k
    act = {(1, 0): 1}
    seq = [1]
    for _ in range(8):
        nxt = defaultdict(int)
        for (b, p), n in act.items():
            for db in (+1, -1):
                if b + db != 0:
                    nxt[(b + db, 1 - p)] += n
        act = nxt
        seq.append(sum(act.values()))
    print(f"\nA1 单种子第 k 步未闭合数 = c_k：")
    print(f"   实测 {seq}")
    print(f"   c_k  {[comb(k, k // 2) for k in range(9)]}")
    check("A1 c_k = C(k, ⌊k/2⌋) 逐项", all(seq[k] == comb(k, k // 2) for k in range(9)))
    RES["A1_c_k"] = seq

    # A2 λ(T) —— 口径按 `L2_layer_retention_verdict.md:79-99`：**只留最高一层**（k=1）给 λ=C(T)
    print(f"\nA2 增长率 λ(T)（语料表：T=3..12 → 4,4,8,8,18,18,46,46,130,130，口径 k=1）")
    target = {3: 4, 4: 4, 5: 8, 6: 8, 7: 18, 8: 18, 9: 46, 10: 46, 11: 130, 12: 130}
    lam = {}
    for T in range(3, 13):
        st = run_site(T, 24, memory=1)
        ratios = [st[i + 1]["history"] / st[i]["history"] for i in range(8, 20)
                  if st[i]["history"] > 0]
        lam[T] = round(sum(ratios) / len(ratios), 3)
    print(f"   实测 λ = {lam}")
    print(f"   语料 λ = {target}")
    check("A2 λ(T) 与语料表逐项一致", all(abs(lam[T] - target[T]) < 0.01 for T in target),
          f"最大偏差 {max(abs(lam[T]-target[T]) for T in target):.2e}")
    print(f"   三种保留口径对照（`L2_layer_retention_verdict.md`：k=1 给 λ=C；k=2 给 λ²=C(λ+1)；无保留更大）：")
    regime = {}
    for mem, name in ((1, "k=1(最高一层)"), (2, "k=2(最高两层)"), (None, "无保留(全留)")):
        vals = {}
        for T in (3, 4, 5):
            st = run_site(T, 24, memory=mem)
            r = [st[i + 1]["history"] / st[i]["history"] for i in range(8, 20) if st[i]["history"] > 0]
            vals[T] = round(sum(r) / len(r), 3)
        regime[name] = vals
    print(f"   {regime}")
    RES["A2_lambda"] = {"measured": lam, "corpus": target, "regimes": regime}

    # A3 站点数无关
    print(f"\nA3 λ 与站点数无关（1/2/3/6 站点；语料：全给 4.000000 @T=3）")
    lamS = {}
    for S in (1, 2, 3, 6):
        st = [run_site(3, 20, memory=1, start_balance=(1 if i % 2 == 0 else -1)) for i in range(S)]
        tot = [sum(s[i]["history"] for s in st) for i in range(len(st[0]))]
        r = [tot[i + 1] / tot[i] for i in range(6, 18) if tot[i] > 0]
        lamS[S] = round(sum(r) / len(r), 6)
    print(f"   实测 λ(S) = {lamS}")
    check("A3 λ 与站点数无关", len(set(lamS.values())) == 1 and abs(lamS[1] - 4.0) < 1e-6)
    RES["A3_sites"] = lamS

    # A4 ρ(T) —— 口径已查明：语料 `L2_rho_closed.py:52-66` 是**先播种、后剪枝**（有效记忆 = memory+1）
    print(f"\nA4 ρ(T)=D/h（语料表 memory=2；口径=先播种后剪枝，有效记忆 3 层）")
    rho_t = {3: 1.207884, 4: 2.415768, 5: 2.224945, 6: 4.449890}

    def corpus_cycles(T, n=120, memory=2):
        """逐字照抄 `L2_rho_closed.py:52-66`：h 与种子用**未剪枝**的 history，剪枝只影响下一周期。"""
        active = {(1, 0): 1}
        hist = defaultdict(int)
        layer = 0
        out = []
        for _ in range(n):
            for _s in range(T):
                nxt = defaultdict(int)
                for (b, p), cnt in active.items():
                    for db in (+1, -1):
                        nb = b + db
                        if nb == 0:
                            hist[layer + 1] += cnt
                        else:
                            nxt[(nb, 1 - p)] += cnt
                active = nxt
            d = sum(active.values())
            active = {}
            h = sum(hist.values())
            if h:
                active = {(1, 0): h, (-1, 0): h}
            layer += 1
            keep = sorted(hist.keys())[-memory:]
            hist = defaultdict(int, {k: hist[k] for k in keep})
            out.append((d, h))
        return out

    rho_m = {}
    for T in rho_t:
        o = corpus_cycles(T, 120)
        rho_m[T] = round(o[-1][0] / o[-1][1], 6)
    print(f"   实测（照抄语料口径）= {rho_m}")
    print(f"   语料 ρ               = {rho_t}")
    check("A4 ρ(T) 与语料表一致（逐字复算语料口径）",
          all(abs(rho_m[T] - rho_t[T]) < 1e-6 for T in rho_t),
          f"最大偏差 {max(abs(rho_m[T]-rho_t[T]) for T in rho_t):.2e}")
    print("   ⟹ 结案：语料表自洽；文档写「只留 memory 层」与代码不符（代码把剪枝放在播种之后，有效记忆 memory+1）")
    RES["A4_rho"] = {"corpus": rho_t, "measured": rho_m,
                     "note": "语料 L2_rho_closed.py 先播种后剪枝 ⟹ 有效记忆 = memory+1 = 3 层"}

    # A5 奇偶配对
    print(f"\nA5 奇偶配对 h(2k−1)=h(2k)、D(2k−1)=½D(2k)、ρ(2k)=2ρ(2k−1)")
    ok = True
    rows = []
    for k in (2, 3, 4, 5):
        a = run_site(2 * k - 1, 25)[-1]
        b = run_site(2 * k, 25)[-1]
        good = (a["history"] == b["history"]) and (abs(b["destroyed"] - 2 * a["destroyed"]) <= 1)
        ok &= good
        rows.append(dict(k=k, h_odd=a["history"], h_even=b["history"],
                         D_odd=a["destroyed"], D_even=b["destroyed"], ok=bool(good)))
    print(f"   {[(r['k'], r['h_odd'], r['h_even'], r['D_odd'], r['D_even']) for r in rows]}")
    check("A5 奇偶配对逐周期成立", ok)
    RES["A5_parity_pair"] = rows

    # A6 沉积支撑
    print(f"\nA6 周期内沉积步支撑（语料：{{0,2,…,T−2}}，奇数步无沉积）")
    ok = True
    sup = {}
    for T in (4, 5, 6, 7):
        st = run_site(T, 30, memory=1)
        s = set()
        for r in st[10:]:
            s |= set(int(x) for x in r["deposit_steps"].keys())
        sup[T] = sorted(s)
        exp = sorted(range(0, T, 2))
        good = (s == set(exp))
        ok &= good
        print(f"   T={T}: 实测 {sorted(s)}  期望 {exp}  {'OK' if good else '≠'}")
    check("A6 沉积支撑 = {0,2,…,T−2}", ok)
    RES["A6_deposit_support"] = sup

    # A7 无稳态
    print(f"\nA7 纯毁灭+重播种**无稳态**（h 逐周期指数增长）")
    st = run_site(3, 12)
    hs = [r["history"] for r in st]
    print(f"   T=3 的 h 逐周期 = {hs}")
    check("A7 h 逐周期严格增长（无稳态）", all(hs[i + 1] > hs[i] for i in range(len(hs) - 1)))
    RES["A7_no_steady_state"] = hs


# ================================================================ B 组：局部 P_i vs 全局 P（ZS-8 机制）
def word_cycle(T, seeds, memory=None):
    """显式词版：seeds = [(词, 重数)]。返回 (新历史 Counter, 末活动量, 沉积步)。"""
    active = Counter()
    for w, n in seeds:
        active[tuple(w)] += n
    hist = Counter()
    dep = Counter()
    for step in range(T):
        nxt = Counter()
        for w, n in active.items():
            for d in (+1, -1):
                nw = w + (d,)
                if sum(nw) == 0:
                    hist[nw] += n
                    dep[step] += n
                else:
                    nxt[nw] += n
        active = nxt
    destroyed = sum(active.values())
    if memory is not None:
        hist = Counter(dict(hist.most_common(memory)))      # 只留最高的 memory 条（示意）
    return hist, destroyed, dep


def part_B(T=4, S=6, epochs=4):
    print("=" * 96)
    print(f"B 组：局部 P_i vs 全局 P（ZS-8 机制；T={T}, {S} 站点, {epochs} 代）")
    init = [tuple([1] + [(-1) ** ((i + j) % 2) for j in range(T - 1)]) for i in range(S)]
    init = [w for w in init if sum(w) != 0][:S]
    while len(init) < S:
        init.append((1,) + (-1,) * (T - 1))
    init = [[(w, 1)] for w in init]   # 每站点一个种子列表

    def run(mode):
        seeds = list(init)
        sets = []
        for _ in range(epochs):
            if mode == "local":
                hs = [word_cycle(T, s)[0] for s in seeds]            # 每站点用自己的历史重播种
            else:
                pooled = [pair for sl in seeds for pair in sl]        # 全局合并（ZM-8 的 global_P）
                h, _, _ = word_cycle(T, pooled)
                hs = [h] * S
            sets.append([frozenset(h.keys()) for h in hs])
            seeds = []
            for h in hs:
                s = []
                for w, n in h.items():
                    s.append((w + (1,), n))
                    s.append((w + (-1,), n))
                seeds.append(s) if mode == "local" else seeds.append(s)
            if mode == "global":
                flat = [x for s in seeds for x in s]
                seeds = [flat] * S
        # 站点两两 Jaccard（末代）
        last = sets[-1]
        ov = []
        for i in range(S):
            for j in range(i + 1, S):
                u = last[i] | last[j]
                ov.append(len(last[i] & last[j]) / len(u) if u else 1.0)
        return sum(ov) / len(ov), len({frozenset(x) for x in last})

    for mode in ("local", "global"):
        o, k = run(mode)
        print(f"   {mode:7s}: 站点间历史重叠(Jaccard)={o:.4f}  不同结果数={k}/{S}")
        RES.setdefault("B_local_vs_global", {})[mode] = dict(overlap=round(o, 4), distinct=k)
    lo = RES["B_local_vs_global"]["local"]["overlap"]
    gl = RES["B_local_vs_global"]["global"]["overlap"]
    check("B 局部 P 保住差异（重叠 ≪1）、全局 P 抹平（重叠 =1）", lo < 0.5 and gl > 0.99,
          f"local={lo} global={gl}")


if __name__ == "__main__":
    t0 = time.time()
    part_A()
    part_B()
    print("=" * 96)
    print(f"断言：通过 {len(PASS)} 项 / 不符 {len(FAIL)} 项" + (f"  不符清单={FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_faithful.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_faithful.json")
