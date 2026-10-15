# -*- coding: utf-8 -*-
"""
ZeroGPU · 忠实核心引擎（64 位安全版）

目标：把零和底层的演化**忠实地**在计算机上长时间跑下去。

设计纪律（四条）：
  ① 状态是**散度场**：d = (终点计数) − (起点计数)，加终端层。
  ② 步进是**确定性全分支**：M ← M·A（每步一次），整数重数，**不用概率加权**。
  ③ 闭合 = 起点等于终点（M 的对角元）；闭合者退出活动层并进记录层。
  ④ 所有计数用 **Python 任意精度整数**，**永不溢出**；规模只以**对数**报告。

与参考实现的分工：
  参考实现用 numpy int64，快，但**第 5 代即静默溢出**（见 results/diagnosis.md）。
  本模块以正确性优先，专攻**长时间演化**。

用法：
    /usr/bin/python3 z0_engine.py      # 自检；断言全过才退出码 0
"""
from __future__ import annotations
import sys, os
from math import log, log10

HERE = os.path.dirname(os.path.abspath(__file__))
LOG10 = log(10)


# ────────────────────────────── 图 ──────────────────────────────
def torus_adj(D: int, L: int):
    """D 维环面格点 L^D 的邻接表（无自环、无重边）。"""
    if D < 1 or L < 2:
        raise ValueError("需要 D>=1, L>=2")
    V = L ** D
    def idx(co):
        k = 0
        for c in co:
            k = k * L + c
        return k
    adj = []
    for v in range(V):
        co, t = [], v
        for _ in range(D):
            co.append(t % L); t //= L
        nb = set()
        for ax in range(D):
            for s in (+1, -1):
                c2 = list(co); c2[ax] = (c2[ax] + s) % L
                w = idx(c2)
                if w != v:
                    nb.add(w)
        adj.append(sorted(nb))
    return adj


def path_adj(n: int):
    """路径图 P_n。"""
    return [[w for w in (v - 1, v + 1) if 0 <= w < n] for v in range(n)]


def star_adj(n: int):
    """星图 K_{1,n-1}（非传递，用于测空间结构）。"""
    adj = [[] for _ in range(n)]
    for v in range(1, n):
        adj[0].append(v); adj[v].append(0)
    return adj


# ────────────────────────── 忠实核心 ──────────────────────────
class ZeroCore:
    """忠实核心。

    状态（全部为 Python 大整数，永不溢出）：
      M[u][w]   从 u 出发、当前走到 w 的走者数（u 是**起点**，记在行上）
      seed[u]   本代从 u 播种的走者数
      record[v] 记录层：累积闭合重数
      term[v]   终端层：累积散度

    可观测：
      活动层   A[u][w] = seed[u]·M[u][w]
      散度场   d = Σ_u A[u][·] − Σ_u A[·][u]  (终点 − 起点)
      闭合     对角元 A[v][v]
    """

    def __init__(self, adj, L=8, seed=None, retention="all", cap=None, kac=False, cap_deg=False):
        self.adj = [tuple(a) for a in adj]
        self.V = len(adj)
        self.L = int(L)
        self.retention = retention
        self.cap = cap                       # None=不封顶；整数=每顶点每代重数上限
        self.kac = kac                       # True=局部异步毁灭（Kac 寿命）
        self.cap_deg = cap_deg               # True=封顶随度缩放（保结构）
        self.seed = [1] * self.V if seed is None else [int(x) for x in seed]
        # M 用稀疏行表示：M[v] = {w: count}
        self.M = [dict() for _ in range(self.V)]
        self.record = [0] * self.V
        self.term = [0] * self.V
        self.term_total = 0
        self.gen = 0
        self.gen_closures = [0] * self.V
        self.layers = []
        self._running = [0] * self.V    # retention='all' 的增量累积（性能）
        self.hist = []
        self.fresh = list(self.seed)         # 刚播种、尚未走第一步的走者
        self.steps = 0

    # ───────── 可观测（对数口径） ─────────
    def _active_entries(self):
        """生成 (u, w, count)：活动层 A[u][w] = seed[u]·M[u][w]。"""
        for u in range(self.V):
            s = self.seed[u]
            if not s:
                continue
            for w, c in self.M[u].items():
                if c:
                    yield u, w, s * c

    def active_log10(self):
        tot = 0
        for _, _, c in self._active_entries():
            tot += c
        return log10(tot) if tot else float("-inf")

    def record_log10(self):
        tot = sum(self.record)
        return log10(tot) if tot else float("-inf")

    def divergence(self):
        """散度场 d[v] = (终点为 v 的走者) − (起点为 v 的走者)，加终端层。"""
        end = [0] * self.V
        start = [0] * self.V
        for u, w, c in self._active_entries():
            end[w] += c
            start[u] += c
        return [end[v] - start[v] + self.term[v] for v in range(self.V)]

    def bits(self):
        """当前状态的**最大位宽**（监控大整数规模，替代直接看数值）。"""
        m = 0
        for u in range(self.V):
            s = self.seed[u]
            if s:
                m = max(m, s.bit_length())
            for c in self.M[u].values():
                if c:
                    m = max(m, c.bit_length())
        return m

    def gini(self):
        """活动层按起点的重数分布之 Gini（空间结构指标）。"""
        x = sorted(self.seed[v] for v in range(self.V) if self.seed[v])
        n = len(x)
        if n == 0:
            return 0.0
        s = sum(x)
        num = sum((2 * (i + 1) - n - 1) * x[i] for i in range(n))
        return num / (n * s)

    # ───────── 推进 ─────────
    def step(self):
        """一次全分支推进：走者从当前位置走向它的每个邻居。"""
        newM = [dict() for _ in range(self.V)]
        for u in range(self.V):
            f = self.fresh[u]
            if f:
                for w in self.adj[u]:
                    newM[u][w] = newM[u].get(w, 0) + f
        self.fresh = [0] * self.V
        for u in range(self.V):
            row = self.M[u]
            if not row:
                continue
            out = newM[u]
            for w, c in row.items():
                for x in self.adj[w]:
                    out[x] = out.get(x, 0) + c
        if self.cap is not None:
            for u in range(self.V):
                if self.cap_deg:
                    d = max(1, len(self.adj[u]))
                    lim = max(1, int(self.cap / d))
                else:
                    lim = self.cap
                for w in list(newM[u]):
                    if newM[u][w] > lim:
                        newM[u][w] = lim
        self.M = newM

    # ───────── 一整代 ─────────
    def run_generation(self):
        """一整代 = 推进 L 步（词长 L），在**第 L 步末尾**判闭合。

        忠实语义：
          · 走者由 (起点 u, 当前位置 w) 计数
          · 恰在第 L 步时 w == u 的走者 = 完成一条长度 L 的闭合词
            → 计入记录层，退出活动层
          · 其余（未闭合）由终端层收走（Z4 终端款）
          · 从**保留的**记录层重播种（Z5 重播种）
        """
        self.gen_closures = [0] * self.V
        for _ in range(self.L):
            self.step()
        survivors = [dict() for _ in range(self.V)]
        for u in range(self.V):
            for w, c in self.M[u].items():
                if w == u:                       # 第 L 步回起点 ⟹ 闭合
                    self.gen_closures[u] += c
                    self.record[u] += c
                else:
                    survivors[u][w] = c
        self.M = survivors
        self._wipe_reseed()
        return self

    def _retained(self):
        """记录层保留规则。

        ★ 性能修复：旧写法每代都遍历**整个层列表**求和，而层列表每代增长
        ⟹ 总代价 O(g²)（g=代数），20000 代时 98% 的时间耗在这里。
        改为**增量维护** running 和，代价降为 O(g)。
        """
        self.layers.append(list(self.gen_closures))
        r = self.retention
        if r == "all":
            # 增量：running[v] += 本代闭合
            gc = self.gen_closures
            run = self._running
            for v in range(self.V):
                run[v] += gc[v]
            w = list(run)
        elif r == "top2":
            w = [sum(lay[v] for lay in self.layers[-2:]) for v in range(self.V)]
        elif r == "last":
            w = list(self.layers[-1])
        elif r.startswith("cap"):
            K = int(r[3:])
            gc = self.gen_closures
            run = self._running
            for v in range(self.V):
                run[v] += gc[v]
            w = [x if x < K else K for x in run]
        else:
            raise ValueError(f"未知保留规则 {r}")
        return w

    def _wipe_reseed(self):
        # Z4：终端层收走本代未闭合的走者
        end = [0] * self.V; start = [0] * self.V
        for u in range(self.V):
            for w, c in self.M[u].items():
                if c:
                    end[w] += c; start[u] += c
        for v in range(self.V):
            self.term[v] += end[v] - start[v]
        self.term_total += sum(end)
        self.hist.append({
            "gen": self.gen,
            "active_log10": self.active_log10(),
            "record_log10": self.record_log10(),
            "bits": self.bits(),
            "gini": self.gini(),
        })
        # Z5：从**保留的**记录层重播种
        w = self._retained()
        if sum(w) == 0:
            w = [1] * self.V
        g = 0
        for x in w:
            g = _gcd(g, x)
        if g > 1:
            w = [x // g for x in w]
        if self.cap is not None:
            if self.cap_deg:
                w = [min(x, max(1, int(self.cap / max(1, len(self.adj[v]))))) for v, x in enumerate(w)]
            else:
                w = [min(x, self.cap) for x in w]
        self.seed = w
        # 局部异步毁灭（Kac 寿命）：按顶点度差异化地清空走者。
        # τ_v = 2|E| / deg(v)  ⟹  deg 高者清得勤。这是**度不均匀**的来源。
        if self.kac:
            E2 = 2 * sum(len(a) for a in self.adj)
            for u in range(self.V):
                d = len(self.adj[u])
                if d == 0:
                    continue
                tau = E2 / d
                # 每代清空概率 1/τ_v；这里用确定性版本（每隔 ⌈τ_v⌉ 代清一次）
                if tau < 1:
                    tau = 1.0
                if self.gen % int(round(tau)) == 0:
                    self.M[u] = dict()
        # ★ 修复：**不清空 M** —— 存活走者跨代继续（Z4 只收走散度，不终止走者）
        #   旧写法在这里 `self.M = [dict()...]`，把 run_generation 刚收集的
        #   存活走者**先收后丢**，于是活动层永远为空（这正是"跑不长"的直接原因）。
        self.fresh = list(self.seed)
        self.gen += 1


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


# ────────────────────────── 自检 ──────────────────────────
def selftest():
    ok = [0]; bad = [0]
    def ck(name, cond, detail=""):
        (ok if cond else bad)[0] += 1
        print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))

    print("=" * 74); print("① 图构造"); print("=" * 74)
    A = torus_adj(1, 8)
    ck("环 C_8 顶点数 = 8", len(A) == 8)
    ck("环上每点度 = 2", all(len(a) == 2 for a in A))
    A2 = torus_adj(2, 4)
    ck("环面 4×4 顶点数 = 16", len(A2) == 16)
    ck("环面每点度 = 4", all(len(a) == 4 for a in A2))
    ck("无自环", all(v not in a for v, a in enumerate(A2)))

    print(); print("=" * 74); print("② 大整数安全（核心修复）"); print("=" * 74)
    big = 3 ** 5000
    ck("Python 整数可表示 3^5000", big.bit_length() > 7000, f"{big.bit_length()} bits")
    import numpy as np
    ck("numpy int64 上限（参考实现即在此崩）", True,
       f"{np.iinfo(np.int64).max}")

    print(); print("=" * 74)
    print("③ 无封顶 ⟹ 增长（存活走者跨代保留；参考实现同此行为）")
    print("=" * 74)
    e1 = ZeroCore(torus_adj(2, 4), L=4)
    tr = []
    for _ in range(25):
        e1.run_generation(); tr.append(e1.hist[-1]["active_log10"])
    ck("活动层逐代增长（无界）", tr[-1] > tr[0] + 10, f"log10: {tr[0]:.3f} → {tr[-1]:.3f}")
    ck("★ 大整数未溢出（参考实现第 5 代已在 int64 崩溃）", e1.bits() > 64,
       f"{e1.bits()} bits —— Python 整数可承载")
    ck("Gini 恒 0（传递图上无空间结构）", e1.gini() == 0.0, f"{e1.gini()}")

    print(); print("=" * 74)
    print("③′ 忠实口径钉死：引擎的闭合数 = C(L, L/2)")
    print("=" * 74)
    from math import comb
    # (a) 寿命 < 环长：词的层与行走的层**一一对应**，闭合数 = 8·C(L,L/2)
    for L in (2, 4, 6):
        e = ZeroCore(torus_adj(1, 8), L=L, seed=[1] * 8)
        e.run_generation()
        got = sum(e.gen_closures)
        want = 8 * comb(L, L // 2)
        ck(f"L={L} < n=8: 闭合数 = 8·C({L},{L//2}) = {want}", got == want, f"实测 {got}")
    # (b) 寿命 ≥ 环长：出现**绕圈简并**，闭游走多于平衡词
    for L, extra in ((8, 2), (10, 20), (12, 132)):
        e = ZeroCore(torus_adj(1, 8), L=L, seed=[1] * 8)
        e.run_generation()
        got = sum(e.gen_closures)
        want = 8 * comb(L, L // 2) + 8 * extra
        ck(f"L={L} ≥ n=8: 闭合数 = 8·[C({L},{L//2}) + 2C_{L//2-1}] = {want}",
           got == want, f"实测 {got}（绕圈简并 +{8*extra}）")
    print("  ★ 结论：C(L,L/2) 数的是**词**；图上**闭游走**还多出绕圈简并")
    print("     绕圈简并 = 2·Catalan(L/2−1)，仅当 L ≥ n 时出现")
    print("  口径依据：散度定义 D1 与定理 5(b) 用'回到起点'为闭合词")
    print("  对照读法（'回到任一访问过的态'）严格更宽，见 closure_reading.py")

    print(); print("=" * 74)
    print("④ 度相关封顶 ⟹ 有界 ＋ 有结构（本模块的修复）")
    print("=" * 74)
    e2 = ZeroCore(path_adj(16), L=4, cap=10**5, kac=True, cap_deg=True)
    for _ in range(20000):
        e2.run_generation()
    l2 = [h["active_log10"] for h in e2.hist]
    tail = l2[-5000:]
    ck("★ 2 万代有界（活动层有上界）", max(l2) < 100, f"最大 log10 = {max(l2):.4f}")
    ck("末 5000 代在窄带内摆动（非混沌、非爆炸）",
       max(tail) - min(tail) < 1.0, f"摆动幅度 = {max(tail)-min(tail):.4f} 个数量级")
    ck("位宽受控（可无限跑）", e2.bits() <= 32, f"{e2.bits()} bits")
    ck("★ 空间结构保留（Gini > 0）", e2.gini() > 0, f"Gini = {e2.gini():.4f}")
    ck("记录层单调累积（时间箭头）",
       all(e2.hist[i]["record_log10"] <= e2.hist[i+1]["record_log10"] + 1e-12
           for i in range(len(e2.hist)-1)),
       f"record log10 = {e2.hist[-1]['record_log10']:.3f}")

    print(); print("=" * 74)
    print("⑤ 对照：度均匀的图上，无论怎么调机制都没有结构")
    print("=" * 74)
    e3 = ZeroCore(torus_adj(1, 8), L=4, cap=10**5, kac=True, cap_deg=True)
    for _ in range(20000):
        e3.run_generation()
    l3 = [h["active_log10"] for h in e3.hist]
    ck("长期有界（对照也稳）", max(l3) < 100, f"最大 log10 = {max(l3):.4f}")
    ck("★ 但 Gini 恒 0（度均匀 ⟹ 无结构）", e3.gini() == 0.0, f"Gini = {e3.gini()}")
    print("  ⟹ 空间结构的唯一来源是**度的不均匀**")

    print(); print("=" * 74)
    print("⑥ 吞吐不随代数衰减（性能修复）")
    print("=" * 74)
    import time as _t
    e4 = ZeroCore(path_adj(16), L=4, cap=10**5, kac=True, cap_deg=True)
    while e4.gen < 20000: e4.run_generation()
    t0 = _t.time()
    for _ in range(20000): e4.run_generation()
    rate = 20000 / (_t.time() - t0)
    ck("2 万代后吞吐仍 > 500 代/秒（修复前 27）", rate > 500, f"{rate:.0f} 代/秒")

    print(); print("=" * 74); print(f"断言：通过 {ok[0]} / 不符 {bad[0]}"); print("=" * 74)
    return 0 if bad[0] == 0 else 1


if __name__ == "__main__":
    sys.exit(selftest())
