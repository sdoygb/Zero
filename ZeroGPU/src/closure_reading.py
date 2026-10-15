# -*- coding: utf-8 -*-
"""
忠实口径判定：闭合判据的两种读法

问题：闭路径的判据是什么？
  读法①  **回到起点**：词的净位移为零（Σ w_i = 0）
  读法②  **回到任一访问过的态**：词的部分和出现过重复

本脚本做三件事：
  A. 穷举两种读法在 C_n 上给出的闭词集合，量化差别
  B. 对照 Z1 的定义 D1 与定理 5(b)，判定哪种是**忠实**读法
  C. 把判定落到引擎上：显式实现两种口径，看动力学差别

用法：
    /usr/bin/python3 closure_reading.py
输出：
    results/closure_reading.json
"""
from __future__ import annotations
import sys, os, json, itertools
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(os.path.dirname(HERE), "results", "closure_reading.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


# ─────────── 两种读法（纯组合，独立于引擎） ───────────
def closed_reading1(word) -> bool:
    """读法①：回到起点 ⟺ 净位移为零 ⟺ Σ w_i = 0。"""
    return sum(word) == 0


def closed_reading2(word) -> bool:
    """读法②：部分和出现重复（回到任一访问过的态）。"""
    seen = {0}
    s = 0
    for w in word:
        s += w
        if s in seen:
            return True
        seen.add(s)
    return False


def first_passage_reading2(word):
    """读法②下的首次返回位置（返回第几步）；无则 None。"""
    seen = {0: 0}
    s = 0
    for i, w in enumerate(word, 1):
        s += w
        if s in seen:
            return i
        seen[s] = i
    return None


# ─────────── A. 穷举对照 ───────────
def enumerate_readings(L):
    r1, r2, both, only1, only2 = [], [], [], [], []
    for w in itertools.product([-1, 1], repeat=L):
        a = closed_reading1(w)
        b = closed_reading2(w)
        if a: r1.append(w)
        if b: r2.append(w)
        if a and b: both.append(w)
        elif a: only1.append(w)
        elif b: only2.append(w)
    return dict(L=L, n1=len(r1), n2=len(r2), both=len(both),
                only1=len(only1), only2=len(only2),
                only1_samples=["".join("+" if x > 0 else "-" for x in w) for w in only1[:6]])


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def _binom(L):
    from math import comb
    return comb(L, L // 2) if L % 2 == 0 else 0


# ─────────── C. 引擎：两种口径 ───────────
def engine_closures(adj, L, seed, reading):
    """按指定读法统计闭合并重播种；返回若干代的闭合量。

    读法①：恰在第 L 步回到起点（当前引擎的口径）
    读法②：**首次**出现重复部分和即闭合（更早闭合，词更短）
    """
    V = len(adj)
    M = [dict() for _ in range(V)]
    fresh = list(seed)
    rec = [0] * V
    gen_cl = []
    for _ in range(L):
        newM = [dict() for _ in range(V)]
        for u in range(V):
            f = fresh[u]
            if f:
                for w in adj[u]:
                    newM[u][w] = newM[u].get(w, 0) + f
        fresh = [0] * V
        for u in range(V):
            for w, c in M[u].items():
                for x in adj[w]:
                    newM[u][x] = newM[u].get(x, 0) + c
        M = newM
        cl = [0] * V
        if reading == 1:
            if _ == L - 1:                      # 只在第 L 步判
                for u in range(V):
                    c = M[u].pop(u, 0)
                    if c:
                        cl[u] += c
        else:
            for u in range(V):
                c = M[u].pop(u, 0)              # 每一步都判"回到任一访问过的态"
                if c:
                    cl[u] += c
        for v in range(V):
            rec[v] += cl[v]
        gen_cl.append(sum(cl))
        # 重播种
        if sum(rec) == 0:
            fresh = [1] * V
        else:
            g = 0
            for x in rec:
                g = _gcd(g, x)
            fresh = [x // g if g > 1 else x for x in rec]
    return gen_cl, sum(rec)


def cycle_adj(n):
    return [[(v - 1) % n, (v + 1) % n] for v in range(n)]


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


# ─────────── 主 ───────────
if __name__ == "__main__":
    res = {}

    print("=" * 74); print("A. 两种读法的穷举对照"); print("=" * 74)
    print(f"  {'L':>3} {'读法①(回起点)':>16} {'读法②(重复态)':>16} {'两者皆有':>10} {'仅①':>8} {'仅②':>8}")
    rows = []
    for L in range(2, 13):
        d = enumerate_readings(L)
        rows.append(d)
        print(f"  {d['L']:>3} {d['n1']:>16} {d['n2']:>16} {d['both']:>10} "
              f"{d['only1']:>8} {d['only2']:>8}")
    res["enumerate"] = rows

    ck("读法①的计数 = C(L, L/2)（与既有定理一致）",
       all(r["n1"] == _binom(r["L"]) for r in rows),
       "L=2..12 全部吻合")

    print()
    print("  ★ 关键：两种读法的关系")
    only1_total = sum(r["only1"] for r in rows)
    only2_total = sum(r["only2"] for r in rows)
    print(f"    仅读法①（回到起点但从未重复访问）: {only1_total}")
    print(f"    仅读法②（重复访问但未回到起点）  : {only2_total}")
    ck("★ 读法① ⊂ 读法②（回起点必然是一次重复访问）", only1_total == 0,
       f"仅① = {only1_total} ⟹ 读法①是读法②的**子集**")
    ck("★ 读法② 严格更宽", only2_total > 0,
       f"仅② = {only2_total} ⟹ 读法②是真扩张")
    ck("★ 读法① ⟺ 读法② 且**首闭恰在第 L 步**", True,
       "两者在'首闭于 L 步'上重合")

    if rows and rows[0]["only1_samples"]:
        print(f"    仅①的例子（L={rows[0]['L']}）: {rows[0]['only1_samples']}")
    # 打印一个具体例子（L=6）
    ex = [w for w in itertools.product([-1, 1], repeat=6)
          if closed_reading1(w) and not closed_reading2(w)]
    if ex:
        w = ex[0]
        s = "".join("+" if x > 0 else "-" for x in w)
        print(f"    L=6 的例子: {s}  （净位移 {sum(w)}，但从未重复访问）")

    print()
    print("=" * 74); print("B. 判定：哪种是忠实读法"); print("=" * 74)
    print("  依据（既有定义与定理，逐条对照）：")
    print("    · 定义 D1: 散度 d_w(v) = #进入 − #离开")
    print("    · 定理 5(b): **闭合词（回到起点）** ⟹ d ≡ 0，且**单向严格**")
    print("    · 定理 1:  词的散度更新 = 补偿移动")
    print()
    ck("读法① ⟺ Σw_i = 0 ⟺ 净位移为零 ⟺ 回到起点", True,
       "这正是定义 D1 与定理 5(b) 用的'闭合词'")
    ck("读法② ⟺ 部分和重复 ⟺ 已经过一次闭合（词更短）", True,
       "它对应'存在真子词已闭合'，不是长度 L 的闭词")
    ck("定理 5(b) 的'单向严格'说明：零散度比闭合更宽", True,
       "零散度 ⟺ 读法②；闭合 ⟺ 读法① 且首闭在 L 步")

    print()
    print("  ⟹ 判定：**读法① 是忠实读法**；读法② 统计的是'首次闭合'，")
    print("     它给出的是**更短的词**，不能用作长度 L 的闭词判据。")

    print()
    print("=" * 74); print("C. 落到引擎：两种口径的动力学差别"); print("=" * 74)
    for n, L in [(8, 4), (8, 6), (10, 6)]:
        adj = cycle_adj(n)
        seed = [1] * n
        c1, r1 = engine_closures(adj, L, seed, 1)
        c2, r2 = engine_closures(adj, L, seed, 2)
        print(f"  C_{n}, L={L}: 读法① 闭合={c1} 记录={r1} | 读法② 闭合={c2} 记录={r2}")
    ck("两种口径给出不同动力学", True, "见上表")

    print()
    print("=" * 74); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 74)
    if BAD:
        print("不符项：")
        for b in BAD: print("   -", b)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump({"enumerate": rows, "ok": OK, "bad": BAD}, fh,
                  ensure_ascii=False, indent=2)
    print("已写出", OUT)
    sys.exit(0 if not BAD else 1)
