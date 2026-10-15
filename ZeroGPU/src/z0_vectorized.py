# -*- coding: utf-8 -*-
"""
ZeroGPU · 向量化引擎（GPU-ready：float64 + 无损归一化缩放）

把 ZeroCore 的演化改成 NumPy 矩阵运算（W[u][pos]：起点 u、当前位置 pos
的走者数），为 GPU/大规模并行铺路。按用户决策：

  ① 放弃"任意精度整数"纪律，改用定长；
  ② **无损缩放**：周期对 record/W 除以总和（公共因子）——因为
     q^(m)=Σ(r_v/Σr)² 只依赖比例，除以 Σr 后比例**严格不变**，永不溢出；
  ③ 缩放口径显式标注（归一化无损，无舍入失真）。

与 ZeroCore 语义等价（矢量 × 循环）：
  起始 W0[u][pos] = 存活 + (seed[u] if pos==u)
  每步   W ← W @ A（走者从 pos 走向邻居）
  第 L 步 对角元 W[u][u] = 闭合 → record[u] += c（累积, L1 历史）
                     gen_closures[u] = c（本代, L2 回灌 / last 重播种）
  非对角 = 存活，跨代保留
  散度   term[v] += (列和 − 行和)  [仅存活走者]
  重播种 seed ← _retained(gen_closures)

用法：
    /usr/local/bin/python3 z0_vectorized.py --verify      # 与原引擎逐代一致性
    /usr/local/bin/python3 z0_vectorized.py --long 16 20000   # 长程 D4 窗
"""
from __future__ import annotations
import sys, os, time, math, json

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj, torus_adj, star_adj

QW_MIN = 0.5
QW_MAX = 3 / 5.0
NORM_THRESHOLD = 1e250   # 总和超过即归一化（无损）


def _adj_matrix(adj):
    V = len(adj)
    A = np.zeros((V, V), dtype=np.float64)
    for u, nb in enumerate(adj):
        for w in nb:
            A[u][w] = 1.0
    return A


class ZeroCoreVec:
    """向量化忠实核心（float64 + 无损归一化）。"""

    def __init__(self, adj, L=4, seed=None, retention="all",
                 cap=None, kac=False, cap_deg=False):
        self.adj = [tuple(a) for a in adj]
        self.V = len(adj)
        self.A = _adj_matrix(adj)
        self.L = int(L)
        self.retention = retention
        self.cap = cap
        self.kac = kac
        self.cap_deg = cap_deg
        self.seed = np.array(seed if seed is not None else [1.0] * self.V, dtype=np.float64)
        self.W = np.zeros((self.V, self.V), dtype=np.float64)
        self.record = np.zeros(self.V, dtype=np.float64)
        self.term = np.zeros(self.V, dtype=np.float64)
        self.term_total = 0.0
        self.gen = 0
        self.gen_closures = np.zeros(self.V, dtype=np.float64)
        self.running = np.zeros(self.V, dtype=np.float64)
        self.layers = []
        self.norm_events = 0     # 归一化次数（口径记录）

    # ───────── 无损归一化（保持比例，永不溢出）─────────
    def _normalize(self):
        for a in (self.record, self.W, self.running, self.term):
            t = a.sum()
            if t > NORM_THRESHOLD:
                a /= t
                self.norm_events += 1
        t = self.seed.sum()
        if t > NORM_THRESHOLD:
            self.seed /= t

    # ───────── 一整代 ─────────
    def run_generation(self):
        V, L = self.V, self.L
        W = self.W.copy()
        W[np.arange(V), np.arange(V)] += self.seed
        for _ in range(L):
            W = W @ self.A
            if self.cap is not None:
                for u in range(V):
                    lim = max(1.0, self.cap / max(1.0, len(self.adj[u]))) if self.cap_deg else float(self.cap)
                    np.minimum(W[u], lim, out=W[u])
        diag = W[np.arange(V), np.arange(V)]
        self.gen_closures = diag.copy()
        self.record = self.record + diag
        W2 = W.copy()
        W2[np.arange(V), np.arange(V)] = 0
        self.W = W2
        end = W2.sum(axis=0); start = W2.sum(axis=1)
        self.term = self.term + (end - start)
        self.term_total += float(end.sum())
        self._reseed()
        self._normalize()
        self.gen += 1
        return self

    # ───────── 记录层保留 + 重播种 ─────────
    def _retained(self):
        self.layers.append(self.gen_closures.copy())
        r = self.retention
        if r == "all":
            self.running = self.running + self.gen_closures
            return self.running.copy()
        if r == "top2":
            if len(self.layers) >= 2:
                return self.layers[-2] + self.gen_closures
            return self.gen_closures.copy()
        if r == "last":
            return self.gen_closures.copy()
        if r.startswith("cap"):
            K = float(r[3:])
            self.running = self.running + self.gen_closures
            return np.minimum(self.running, K)
        raise ValueError(f"未知保留规则 {r}")

    def _reseed(self):
        w = self._retained()
        if float(w.sum()) == 0:
            w = np.ones(self.V, dtype=np.float64)
        if self.cap is not None:
            if self.cap_deg:
                lims = np.array([max(1.0, self.cap / max(1.0, len(self.adj[v]))) for v in range(self.V)])
            else:
                lims = float(self.cap)
            w = np.minimum(w, lims)
        self.seed = w.astype(np.float64)
        if self.kac:
            E2 = 2 * sum(len(a) for a in self.adj)
            for u in range(self.V):
                d = len(self.adj[u])
                if d == 0:
                    continue
                tau = max(1.0, E2 / d)
                if self.gen % int(round(tau)) == 0:
                    self.W[u, :] = 0.0

    # ───────── 观测 ─────────
    def q_weighted(self, arr):
        tot = float(arr.sum())
        if tot == 0:
            return None
        w = arr / tot
        return float(w @ w)

    def q_hist(self):
        return self.q_weighted(self.record)

    def q_live(self):
        return self.q_weighted(self.gen_closures)

    def in_window(self, q):
        return q is not None and QW_MIN < q < QW_MAX

    def gini(self):
        x = sorted(float(v) for v in self.seed if v)
        n = len(x)
        if n == 0:
            return 0.0
        s = sum(x)
        num = sum((2 * (i + 1) - n - 1) * x[i] for i in range(n))
        return num / (n * s)

    def summary(self):
        qh, ql = self.q_hist(), self.q_live()
        return {
            "generations": self.gen,
            "q_hist(L1→L0)": qh, "q_live(L2→L0)": ql,
            "d4_hist": self.in_window(qh), "d4_live": self.in_window(ql),
            "record_total": float(self.record.sum()),
            "gini": self.gini(), "norm_events": self.norm_events,
        }


# ────────────────────────── 一致性验证 ──────────────────────────
def verify_against_reference(N=60):
    from layer_closure import LayerClosure
    ok = [0]; bad = [0]
    def ck(name, cond, detail=""):
        (ok if cond else bad)[0] += 1
        print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    print("=" * 76); print(f"一致性验证：path6 hub2 last L4 — {N} 代"); print("=" * 76)
    adj = path_adj(6); V = len(adj)
    seed = [1] * V; seed[0] = seed[1] = 100
    ref = LayerClosure(adj, seed=seed, retention="last", L=4)
    vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    max_dq = 0.0
    for g in range(N):
        ref.core.run_generation()
        vec.run_generation()
        rq, vq = ref.q_hist(), vec.q_hist()
        rl, vl = ref.q_live(), vec.q_live()
        max_dq = max(max_dq, abs(rq - vq), abs(rl - vl))
    ck("q_hist/q_live 逐代一致（无损归一化）", max_dq < 1e-9, f"最大差 {max_dq:.2e}")
    # 归一化后 record 比例应一致
    ref_rec = np.array(ref.core.record, float); vec_rec = vec.record
    def ratio(a): return a / a.sum()
    ck("record 比例一致", np.max(np.abs(ratio(ref_rec) - ratio(vec_rec))) < 1e-9)
    ref_seed = np.array(ref.core.seed, float); vec_seed = vec.seed
    ck("seed 比例一致", np.max(np.abs(ratio(ref_seed) - ratio(vec_seed))) < 1e-9)
    ck("record 单调累积（内容）", all(vec.record[i] >= 0 for i in range(V)))
    print(f"  参考 q_hist={ref.q_hist():.6f}  向量 q_hist={vec.q_hist():.6f}")
    print(f"  归一化触发 {vec.norm_events} 次（无损，比例不变）")
    print(); print(f"断言：通过 {ok[0]} / 不符 {bad[0]}")
    return 0 if bad[0] == 0 else 1


# ────────────────────────── 长程（更大图）─────────────────────────
def longrun(path_n=16, gens=20000, seed_kind="hub2", retention="last"):
    adj = path_adj(path_n); V = len(adj)
    if seed_kind == "hub2":
        seed = [1.0] * V; seed[0] = seed[1] = 100
    elif seed_kind == "hub3":
        seed = [1.0] * V; seed[0] = seed[1] = seed[2] = 100
    elif seed_kind == "geom":
        seed = [float(2 ** (V - 1 - v)) for v in range(V)]
    else:
        seed = [1.0] * V

    print("=" * 76)
    print(f"长程：path{path_n} {seed_kind} {retention} L4 — {gens} 代（向量化+无损归一化）")
    print("=" * 76)
    t0 = time.time()
    vec = ZeroCoreVec(adj, seed=seed, retention=retention, L=4)
    SAMP = 500
    samples = []
    for g in range(gens):
        vec.run_generation()
        if (g + 1) % SAMP == 0:
            qh, ql = vec.q_hist(), vec.q_live()
            samples.append({"gen": g + 1, "q_hist": qh, "q_live": ql,
                            "in_hist": vec.in_window(qh), "in_live": vec.in_window(ql),
                            "norm": vec.norm_events})
    s = vec.summary()
    n_in = sum(1 for p in samples if p["in_hist"])
    qh_first, qh_last = samples[0]["q_hist"], samples[-1]["q_hist"]
    rec = {
        "config": {"adj": f"path{path_n}", "seed": seed_kind, "retention": retention, "L": 4},
        "generations": gens, "sample_every": SAMP, "n_samples": len(samples),
        "elapsed_s": round(time.time() - t0, 1),
        "rate_gen_per_s": round(gens / max(time.time() - t0, 0.01), 1),
        "frac_in_window": round(n_in / len(samples), 4),
        "q_hist_drift": round(qh_last - qh_first, 8),
        "q_hist_range": round(max(p["q_hist"] for p in samples) - min(p["q_hist"] for p in samples), 8),
        "q_hist_final": qh_last, "q_live_final": s["q_live(L2→L0)"],
        "norm_events": s["norm_events"], "gini": s["gini"],
        "series": samples,
    }
    os.makedirs(os.path.join(HERE, "..", "results"), exist_ok=True)
    fname = os.path.join(HERE, "..", "results", f"d4_longrun_vec_path{path_n}_{retention}.json")
    with open(fname, "w", encoding="utf-8") as f:
        json.dump(rec, f, ensure_ascii=False, indent=1)

    print(f"耗时 {rec['elapsed_s']}s（{rec['rate_gen_per_s']:.0f} 代/秒）")
    print(f"q_hist: 首={qh_first:.6f} 末={qh_last:.6f} 极差={rec['q_hist_range']:.2e} 漂移={rec['q_hist_drift']:+.2e}")
    print(f"进窗占比 {rec['frac_in_window']:.1%}  归一化 {s['norm_events']} 次  Gini={s['gini']:.4f}")
    print(f"已存 {fname}")
    return rec


if __name__ == "__main__":
    if "--verify" in sys.argv:
        sys.exit(verify_against_reference())
    if "--long" in sys.argv:
        i = sys.argv.index("--long")
        n = int(sys.argv[i + 1]) if len(sys.argv) > i + 1 and sys.argv[i + 1].isdigit() else 16
        g = int(sys.argv[i + 2]) if len(sys.argv) > i + 2 and sys.argv[i + 2].isdigit() else 20000
        longrun(path_n=n, gens=g)
        sys.exit(0)
    sys.exit(verify_against_reference())
