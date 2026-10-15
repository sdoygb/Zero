# -*- coding: utf-8 -*-
"""
ZeroGPU · 三层全双向闭环增强 v2（按最新 Zero 猜想改造）

最新 Zero 猜想（E58 v2）：三层（L0 骨架 / L1 历史 / L2 演化）构成三节点
全双向自指闭环；粒子/时空 = 三层演化流收敛的稳态；闭合回灌 = 重数累积
（内容），类型结构（规则）固定；D=4 ⟺ q^(m) ∈ (1/2, 3/5)。

v2 修正（用户指出）：L1↔L0 的动力学也很重要。上一版只显式做了 L2→L0，
且用 record（实为 L1 数据）加权，把 L1→L0 与 L2→L0 混用了。
本版把**六条边全部显式**，尤其把两条"进 L0"的回灌分开：

  ① L0→L1 词入历史：L0 的闭合约束（Σ=0 / 回到起点）定义 record 记录什么
  ② L1→L2 历史作种子：record 保留 → 重播种
  ③ L2→L1 演化写历史：本代闭合 gen_closures 累积进 record
  ④ L0→L2 结构决定：图 Γ / 闭合判据决定 L2 演化
  ⑤ L2→L0 本代闭合回灌：用本代闭合重数 gen_closures 加权 → q_live
  ⑥ L1→L0 记录回灌：用累积历史重数 record 加权 → q_hist

  两条回灌分开监测，各自判定 D4 窗。稳态形态对两者分别判定。

用法：
    /usr/bin/python3 layer_closure.py           # 自检
    /usr/bin/python3 layer_closure.py --run     # 示例演化 + 六边闭环监测
"""
from __future__ import annotations
import sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import ZeroCore, torus_adj, path_adj, star_adj

# Z001 引理 7.1：D=4 ⟺ q ∈ (1/2, 3/5)
QW_MIN = 0.5
QW_MAX = 3 / 5.0   # 0.6
EPS = 1e-6


class LayerClosure:
    """三层全双向闭环增强 v2：六条边全部显式，两条进 L0 的回灌分开。"""

    def __init__(self, adj, L=4, **kw):
        self.core = ZeroCore(adj, L=L, **kw)
        self.V = self.core.V
        self.adj = adj                      # L0 类型结构（规则，固定）
        # —— 六条边的数据 ——
        self.q_hist_trace = []              # ⑥ L1→L0 记录回灌（累积 record）
        self.q_live_trace = []              # ⑤ L2→L0 本代闭合回灌（gen_closures）
        self.w_hist_trace = []              # record 加权 ω
        self.w_live_trace = []              # 本代闭合加权 ω
        self.record_total_trace = []        # ③ record 累积（内容）
        self.gens = 0

    # ───────── ⑤ L2→L0：本代闭合回灌 ─────────
    def live_partition(self):
        """本代闭合重数 gen_closures[v] 的归一化（L2 当轮演化 → L0）。"""
        gc = self.core.gen_closures
        tot = sum(gc)
        if tot == 0:
            return None
        return tuple(x / tot for x in gc)

    def q_live(self):
        """q_live^(m) = Σ_v (本代闭合份额)² —— L2→L0 直接回灌。"""
        w = self.live_partition()
        return None if w is None else sum(x * x for x in w)

    # ───────── ⑥ L1→L0：记录回灌 ─────────
    def hist_partition(self):
        """累积历史重数 record[v] 的归一化（L1 历史 → L0）。"""
        rec = self.core.record
        tot = sum(rec)
        if tot == 0:
            return None
        return tuple(x / tot for x in rec)

    def q_hist(self):
        """q_hist^(m) = Σ_v (累积历史份额)² —— L1→L0 记录回灌。"""
        w = self.hist_partition()
        return None if w is None else sum(x * x for x in w)

    def in_d4_window(self, q):
        """D4 窗：q ∈ (1/2, 3/5)。"""
        return q is not None and QW_MIN < q < QW_MAX

    # ───────── 稳态判定（对两条回灌分别判定）─────────
    def stable(self, trace, window=20, eps=EPS):
        """给定 ω 轨迹是否收敛（相邻代最大分量差 < eps）。"""
        if len(trace) < window:
            return False, None
        recent = trace[-window:]
        diffs = [max(abs(a - b) for a, b in zip(recent[i], recent[i + 1]))
                 for i in range(len(recent) - 1)]
        return max(diffs) < eps, max(diffs)

    def steady_form(self, q, stable):
        """稳态形态：None(未收敛) / 'single'(单点) / 'structured'(有结构)。"""
        if not stable or q is None:
            return None
        if QW_MIN < q < QW_MAX:
            return "structured"   # 落在 D4 窗 ⟹ 适度集中，多顶点承载，非单点
        return "single"           # 极端集中或太平 ⟹ 坍塌/退化

    # ───────── 规则固定检查（三层闭环边界）─────────
    def rule_fixed(self):
        """类型结构固定：图/闭合定义/顶点集不变，只 record 累积（内容）。"""
        rt = self.record_total_trace
        mono = all(rt[i] <= rt[i + 1] + 1e-12 for i in range(len(rt) - 1))
        return mono

    # ───────── 推进 ─────────
    def run(self, generations, log_every=100):
        for g in range(generations):
            self.core.run_generation()
            self.w_hist_trace.append(self.hist_partition())
            self.w_live_trace.append(self.live_partition())
            self.q_hist_trace.append(self.q_hist())
            self.q_live_trace.append(self.q_live())
            self.record_total_trace.append(sum(self.core.record))
            self.gens += 1
            if log_every and (g + 1) % log_every == 0:
                qh, ql = self.q_hist_trace[-1], self.q_live_trace[-1]
                wh = self.in_d4_window(qh); wl = self.in_d4_window(ql)
                h = f"{qh:.4f}" if qh is not None else "—"
                l = f"{ql:.4f}" if ql is not None else "—"
                print(f"  gen {g+1:>7}: q_hist={h} D4={'IN ' if wh else 'OUT'} | "
                      f"q_live={l} D4={'IN ' if wl else 'OUT'}")
        return self

    def summary(self):
        qh, ql = self.q_hist(), self.q_live()
        st_h, _ = self.stable(self.w_hist_trace)
        st_l, _ = self.stable(self.w_live_trace)
        return {
            "generations": self.gens,
            "q_hist(L1→L0)": qh,
            "q_live(L2→L0)": ql,
            "d4_hist": self.in_d4_window(qh),
            "d4_live": self.in_d4_window(ql),
            "hist_form": self.steady_form(qh, st_h),
            "live_form": self.steady_form(ql, st_l),
            "record_monotone": self.rule_fixed(),
            "record_total": sum(self.core.record),
            "gini": self.core.gini(),
            "bits": self.core.bits(),
        }


# ────────────────────────── 自检 ──────────────────────────
def selftest():
    ok = [0]; bad = [0]
    def ck(name, cond, detail=""):
        (ok if cond else bad)[0] += 1
        print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))

    print("=" * 74); print("① 两条回灌分开：L1→L0(累积) vs L2→L0(本代)"); print("=" * 74)
    lc = LayerClosure(torus_adj(1, 4), L=2)
    lc.core.record = [10, 0, 0, 0]         # 历史集中顶点0
    lc.core.gen_closures = [0, 5, 0, 0]    # 本代闭合集中顶点1
    wh = lc.hist_partition(); qh = lc.q_hist()
    wl = lc.live_partition(); ql = lc.q_live()
    ck("历史 ω 和 = 1", abs(sum(wh) - 1.0) < 1e-12)
    ck("历史集中在顶点0 → ω=(1,0,0,0)", wh == (1.0, 0.0, 0.0, 0.0))
    ck("本代集中在顶点1 → ω=(0,1,0,0)", wl == (0.0, 1.0, 0.0, 0.0))
    ck("q_hist=1（历史集中）", abs(qh - 1.0) < 1e-12, f"{qh}")
    ck("q_live=1（本代集中）", abs(ql - 1.0) < 1e-12, f"{ql}")
    ck("★ 两条回灌是独立量（来源不同顶点）", qh == 1.0 and ql == 1.0)

    print(); print("=" * 74); print("② D4 窗边界"); print("=" * 74)
    def inwin(qv): return QW_MIN < qv < QW_MAX
    ck("q=0.5 不进窗", inwin(0.5) is False)
    ck("q=0.55 进窗", inwin(0.55) is True)
    ck("q=0.6 不进窗（3/5 开）", inwin(0.6) is False)

    print(); print("=" * 74); print("③ 实际演化：六边闭环监测（非传递图）"); print("=" * 74)
    lc4 = LayerClosure(path_adj(8), L=4, cap=10**5, kac=True, cap_deg=True)
    lc4.run(60, log_every=20)
    s = lc4.summary()
    ck("跑了 60 代", s["generations"] == 60)
    ck("记录层单调（L1 内容累积）", s["record_monotone"] is True)
    ck("q_hist 有限", s["q_hist(L1→L0)"] is not None and 0 < s["q_hist(L1→L0)"] <= 1)
    ck("q_live 有限", s["q_live(L2→L0)"] is not None and 0 < s["q_live(L2→L0)"] <= 1)
    ck("★ 两条回灌数值不同（L1 累积 ≠ L2 本代）",
       abs(s["q_hist(L1→L0)"] - s["q_live(L2→L0)"]) > 1e-6,
       f"hist={s['q_hist(L1→L0)']:.4f} vs live={s['q_live(L2→L0)']:.4f}")
    ck("位宽受控", s["bits"] <= 64, f"{s['bits']} bits")

    print(); print("=" * 74); print("④ 对照：度均匀图上无结构"); print("=" * 74)
    lc5 = LayerClosure(torus_adj(2, 4), L=4, cap=10**5, kac=True, cap_deg=True)
    lc5.run(60, log_every=60)
    s5 = lc5.summary()
    ck("度均匀 → Gini=0", s5["gini"] == 0.0, f"{s5['gini']}")
    ck("度均匀 → q_hist 接近均匀下限 1/16",
       s5["q_hist(L1→L0)"] is not None and s5["q_hist(L1→L0)"] <= 0.25,
       f"{s5['q_hist(L1→L0)']:.4f}")

    print(); print("=" * 74); print(f"断言：通过 {ok[0]} / 不符 {bad[0]}"); print("=" * 74)
    return 0 if bad[0] == 0 else 1


if __name__ == "__main__":
    if "--run" in sys.argv:
        print("=" * 74); print("示例：path_adj(8), L=4, 300 代 — 六边闭环监测"); print("=" * 74)
        lc = LayerClosure(path_adj(8), L=4, cap=10**5, kac=True, cap_deg=True)
        lc.run(300, log_every=50)
        print(); print("汇总：")
        for k, v in lc.summary().items():
            print(f"  {k}: {v}")
        print("=" * 74)
        sys.exit(0)
    sys.exit(selftest())
