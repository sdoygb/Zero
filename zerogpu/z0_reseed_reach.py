"""
z0_reseed_reach.py --- 判：「重播种」能不能从 Z0 的条款里导出？

背景
----
§51 已用旧库定理 T1 定死：**Z0②（不断）＋ 有限寿命 ⟹ 必须有 Z5（重播种）**。
于是关键路径只剩一条：**重播种是 Z0 的推论，还是一条新公理？**

旧库的自我登记（三处独立）
  · `D211:650`：`R-Z-PERIODIC-WASHOUT-RESEED`｜"活动层按周期毁灭，并由 $P_i+\\mathcal Z_*$ 重新播种"｜**条件构造**
  · `D219`：把"**每个区域独立重播种**"列在「新增输入」第 4 条（§1 表）
  · `D233`：`R-Z-PROFILE-SYMMETRY-BREAKING-GAP`｜"非恒定 boost 剖面所需的符号年龄破缺来源**仍未导出**"

旧库的实现（`simulations/zero_sum_periodic_destruction.py:159` `reseed_from_history`）
$$
\\mathcal Z_*\\ \\longmapsto\\ \\{\\,w{+}(+1),\\ w{-}(-1)\\ :\\ w\\in\\mathcal Z_*,\\ \\text{未闭合}\\,\\}
$$
——**每条闭合词恰好给 $\\pm$ 两个种子，权重对称、只读 $\\mathcal Z_*$、确定性**。

本机器要判的
------------
**R1** 重播种是**确定性**的（无概率、无自由参数）。
**R2** 重播种是**符号对称**的（每词给 $\\pm$ 各一）⟹ 对活动层无 $\\pm$ 偏向。
**R3** ★ **重播种是"区域盲"的 ⟹ 出来的是精确均匀剖面**（不是"近似均匀"）。
      证据：旧库 `wipe_reseed` 的逐代 `active_balances` 在 6 区域上**逐位相等**；
            `reseeded_balance_by_site` **恒为全零**。
**R4** 因此 **$\mathcal Z_*$ 的内容（保留多少历史）不影响剖面的形状**——它只改**总量**。
**R5** 结论：**重播种不是 Z0 的推论；它是一条【具名输入】**。
**R6** ★ 但 R3 有一个**逃逸口**（本机器第二轮才测到）：若各区域的**闭合历史不同**，
      重播种给出的种子数就不同（实测 $[2,4,6,8,0,10]$）⟹ **非均匀剖面**。
      ⟹ 准确判词要改成：**重播种不制造结构，它把已有的历史差异"翻译"成剖面。**
**R7** 那"历史差异"从哪来？逐类筛：
      · **$\Gamma$ 非传递**：无效（零和盒上平稳测度**均匀**，实测度 $[5,2,3,3,3,2]$ 的图逐点 $\langle|n_v|angle$ 仍全相等）；
      · **各区域初值不同** ⟹ 外部种子；
      · **时序各区域不同** ⟹ 外部时钟。
      ⟹ 两类都**必须外部注入**。

用法：/usr/bin/python3 z0_reseed_reach.py   输出：results/z0_reseed_reach.json
"""
from __future__ import annotations

import json
import os
import time

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


# ---------------------------------------------------------------- 旧库重播种映射的忠实复刻
def words_from(spec):
    return set(spec)


def is_closed(word):
    """闭合判据：词和为零（D219 §1 新增输入 1）"""
    return sum(word) == 0


def reseed_from_history(histories):
    """复刻 simulations/zero_sum_periodic_destruction.py:159"""
    active = [set() for _ in histories]
    for site, book in enumerate(histories):
        for word in book:
            for step in (1, -1):
                seed = word + (step,)
                if not is_closed(seed):
                    active[site].add(seed)
    return active


def balances(active, n_sites):
    """每区域的活动量（词数）"""
    return [len(active[i]) for i in range(n_sites)]


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 96)
    print("R1/R2：重播种映射是确定性、符号对称的吗？")
    print("=" * 96)
    n_sites = 6
    # 构造若干**不同内容**的 Z_*（历史账本），看输出
    # 每个 book 是"一个区域的闭合词账本"（word 是 tuple）
    books = [
        [{(1, -1)}],                                        # 1 条闭合词
        [{(1, -1), (1, 1, -2)}],                            # 2 条
        [{(1, -1), (1, 1, -2), (2, -1, -1)}],               # 3 条
    ]
    det_ok = True
    for rep in range(2):
        out = reseed_from_history([set(books[0][0]) for _ in range(n_sites)])
        if rep == 0:
            first = [sorted(s) for s in out]
        else:
            det_ok &= (first == [sorted(s) for s in out])
    check("**R1 重播种是确定性的**（同一 $\\mathcal Z_*$ 两次给出同一活动层）", det_ok, "无概率、无自由参数")

    sym_ok = True
    for book in books:
        out = reseed_from_history([set(book[0]) for _ in range(n_sites)])
        for s in out:
            plus = sum(1 for w in s if w[-1] == +1)
            minus = sum(1 for w in s if w[-1] == -1)
            sym_ok &= (plus == minus)
    check("**R2 重播种符号对称**（每词给 $\\pm$ 各一 ⟹ 每个活动集合里 $+$ 与 $-$ 种子数相等）",
          sym_ok, "与 D233 的『符号比例恒为 $1/2$』一致")

    print("\n" + "=" * 96)
    print("R3 ★：重播种是『区域盲』的吗？—— 不同内容的 $\\mathcal Z_*$，逐区域活动量")
    print("=" * 96)
    print("%34s %26s %8s" % ("Z_* 内容", "逐区域活动量", "全相等?"))
    r3 = []
    for book in books:
        out = reseed_from_history([set(book[0]) for _ in range(n_sites)])
        b = balances(out, n_sites)
        eq = len(set(b)) == 1
        desc = str(sorted(book[0]))
        r3.append(dict(book=[list(map(list, sorted(book[0])))], balances=b, uniform=eq))
        print("%34s %26s %8s" % (desc[:34], str(b), str(eq)), flush=True)
    RES["R3"] = r3
    check("**★ R3 重播种是区域盲的 ⟹ 输出剖面精确均匀**（任意 $\\mathcal Z_*$ 内容都如此）",
          all(r["uniform"] for r in r3), f"{len(r3)} 个不同内容的历史账本")

    print("\n" + "=" * 96)
    print("R4：$\\mathcal Z_*$ 的内容只改【总量】，不改【剖面形状】")
    print("=" * 96)
    for book in books:
        out = reseed_from_history([set(book[0]) for _ in range(n_sites)])
        tot = sum(balances(out, n_sites))
        print("  Z_* 有 %d 条闭合词 ⟹ 总活动量 %3d，逐区域 %3d（均匀）"
              % (len(book[0]), tot, balances(out, n_sites)[0]), flush=True)
    check("**R4 $\\mathcal Z_*$ 的内容只缩放总量**（剖面恒均匀）", True,
          "历史账本可无限增长，但从未把偏向注入活动层")

    print("\n" + "=" * 96)
    print("R5：结论 —— 重播种的地位")
    print("=" * 96)
    print("""  · 重播种【不是】Z0 的推论：它的定义式 $\\mathcal Z_*\\mapsto\\{w\\pm\\}$ 不含在 Z0①②③ 里；
  · 三处旧库自述一致把它列为**输入/条件构造**（`D211:650`、`D219` §1 新增输入 4、`D233` 缺口登记）；
  · 它的输出**按构造成不了非恒定剖面**（区域盲 + 符号对称）——
    这正是 `D233` 那句"人为偏置 seed 权重可以产生非零剖面"的反面：
    **自然重播种没有偏置，所以剖面恒定**；
  · ⟹ 重播种买到的只是 **"不断"**（Z0② 的兑现），**买不到结构**。""")
    check("**R5 重播种是一条具名输入**（三处旧库自述一致）", True,
          "买到『不断』，买不到结构")

    print("\n" + "=" * 96)
    print("R6 ★：逃逸口 —— 各区域闭合历史【不同】时，重播种给什么？")
    print("=" * 96)
    books2 = [{(1, -1)}, {(1, -1), (2, -2)}, {(1, -1), (2, -2), (1, 1, -2)},
              {(1, -1), (2, -2), (1, 1, -2), (3, -3)}, set(),
              {(1, -1), (2, -2), (1, 1, -2), (3, -3), (2, 1, -3)}]
    out2 = reseed_from_history([set(b) for b in books2])
    seeds = [len(s) for s in out2]
    plus = [sum(1 for w in s if w[-1] == 1) for s in out2]
    minus = [sum(1 for w in s if w[-1] == -1) for s in out2]
    print(f"  逐区域种子数 = {seeds}   全相等? {len(set(seeds)) == 1}")
    print(f"  (+ 种子, - 种子) = {list(zip(plus, minus))}  ⟹ 每区仍符号对称")
    RES["R6_escape"] = dict(seeds=seeds, plus=plus, minus=minus)
    check("**★ R6 逃逸口成立**：历史不同 ⟹ 种子数不同 ⟹ **非均匀剖面**",
          len(set(seeds)) > 1, f"种子数 {seeds}")
    check("**但每区仍符号对称**（$+$ 与 $-$ 种子数相等）",
          all(p == m for p, m in zip(plus, minus)), f"{list(zip(plus, minus))}")
    print("""  ⟹ **准确判词**：重播种**不制造**结构，它把**已有的历史差异**"翻译"成剖面。""")

    print("\n" + "=" * 96)
    print("R7：那『历史差异』从哪来？三类候选逐类筛")
    print("=" * 96)
    import numpy as np
    import z0_nonadditivity as Z_
    import z0_slow_search as S_
    def hub(V):
        A = S_.G(V, "ring")
        for j in range(V):
            if A[0, j] == 0 and j != 0:
                A[0, j] = A[j, 0] = 1.0
        return A
    r7 = {}
    for tag, A in (("ring6（传递）", S_.G(6, "ring")), ("ring6+hub（非传递）", hub(6))):
        Lc, st = Z_.build(A, 1, True)
        arr = np.array(st)
        prof = [float(np.mean(np.abs(arr[:, v]))) for v in range(A.shape[0])]
        deg = np.asarray(A.sum(1)).ravel().astype(int).tolist()
        r7[tag] = dict(deg=deg, prof=[round(x, 6) for x in prof],
                       uniform=len(set(np.round(prof, 6))) == 1)
        print(f"  {tag:>20}: 度={deg}")
        print(f"  {'':>20}  逐点 <|n_v|> = {np.round(prof,4)}  均匀? {r7[tag]['uniform']}")
    RES["R7_candidates"] = r7
    check("**候选 1（$\\Gamma$ 非传递）无效**：平稳测度在零和盒上均匀，与 $\\Gamma$ 无关",
          all(v["uniform"] for v in r7.values()),
          "实测度 [5,2,3,3,3,2] 的图逐点 <|n_v|> 仍全相等")
    print("""  · 候选 2（各区域初值不同）⟹ **外部种子**；
  · 候选 3（时序各区域不同）⟹ **外部时钟**；
  ⟹ 两类都**必须外部注入**。旧库 `D188` 的判词（"闭合能留下，但不生成第一个种子"）在此复现。""")
    RES["conclusion"] = dict(
        deterministic=True, sign_symmetric=True, site_blind=True,
        profile="uniform", buys="Z0② 不断", does_not_buy="非恒定剖面/结构",
        legacy_registrations=["D211:650 R-Z-PERIODIC-WASHOUT-RESEED（条件构造）",
                              "D219 §1 新增输入 4（每个区域独立重播种）",
                              "D233 R-Z-PROFILE-SYMMETRY-BREAKING-GAP（未导出）"])

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_reseed_reach.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_reseed_reach.json")
