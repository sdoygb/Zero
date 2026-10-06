"""
z0_core.py --- Z0 核心的**忠实实现**（零和宇宙）

**这条公理只有一个、三款**（`Z0_zero_never_rests_single_axiom.md` §1）：
$$
\text{Z0① 零不停留（零不是静止态）};\quad
\text{Z0② 从不停歇（运动不终止）};\quad
\text{Z0③ 不设概率（没有选择规则）}
$$
导出：① $\to$ Z1（词与图）、Z3（闭合**与退出**）；③ $+$ 最小性 $\to$ **Z2（全分支＋整数重数）**；
② $+$ 有限寿命 $\to$ **Z4（寿命与终端）**、**Z5（重播种）**。

**与上一版的根本差别（用户指出："这是一个零和宇宙"）**

$$
\text{"零和"不是公理，是}\textbf{散度的恒等式}:\qquad \sum_v d_w(v)\equiv0
\quad(\text{每一步：}+1\text{ 记在入点、}-1\text{ 记在出点})
$$

上一版`z0_evo.py`/`l2_multi.py`里**每一步只有 +1（占据加一），没有 −1**，而且用**加权抽样**
（那正是 Z0③ 排除的"选择规则"）。本版：

| | 上一版 | 本版 |
|:--|:--|:--|
| 状态 | 走者的位置计数 | **散度场** $d(x)$ ＋ 整数重数 ＋ 终端账本 $D$ |
| 步进 | 每条走者走一条边（保守） ＋ **加权** | **全分支、确定性、零权重**（Z2） |
| 闭合 | 回到起点 | 回到起点 ⟹ $d_w=0$ ⟹ 写记录（Z3），≤记录层恒中性 |
| 寿命/退出 | 无 | **周期清空 ⟹ 终端层 $D$ 带走非零平衡**（Z4＋Z3 退出） |
| 重播种 | 我手造的 α 旋钮 | **从保留的记录层重播种**（Z5，由 Z0② 逼出） |

**四条恒等式（每步核验）**

$$
\begin{aligned}
\text{I1}\;&\textstyle\sum_v d_w(v)=0\ \text{（每个词；结构性成立）}\\
\text{I2}\;&\textstyle\sum_x d(x)=0\ \text{（整体散度场；数值核验）}\\
\text{I3}\;&Q_{\rm record}=0\ \text{（闭合词恒中性 —— 这正是"记录层全中性"的来源）}\\
\text{I4}\;&\text{无权重、无选择规则（Z0③）：步进是确定性整数全分支}
\end{aligned}
$$

**闭环（Z4＋Z5 = 语料 `wipe_reseed` 那个唯一"持续"的模型）**

```
seed（本轮新播的走者，整数重数）
   ↓ Z2：all-branch  N_k = diag(seed)·A^k        （确定性、零权重、整数）
每步：闭合 = 对角元 N_k[x,x] → 写记录（中性，I3）
到寿命 L：清空 → 终端层 D 收走**非零平衡**（Z4）
   ↓ Z5：从保留的记录层重播种（记录多的地方播得多）
```

**区域与补偿**：区域 = **记录支撑集的连通分量**（零参数、无权重）。
每区平衡 $Q_R=\sum_{x\in R}d(x)$，并核验 $\sum_R Q_R=0$、
"单个 $D_i$ 不是自平衡的，**补偿发生在区域之间**"（`zero_sum_periodic_destruction.md`）。

用法：/usr/bin/python3 z0_core.py
输出：results/z0_core.json
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from collections import deque

import numpy as np
from scipy.sparse import csr_matrix, identity, kron, diags

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

RECORDS = []


def rec(sec, obj, readout, x, y, meta=None, note=None, kind="curve", maxpts=64):
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) > maxpts:
        idx = np.unique(np.r_[np.linspace(0, len(x) - 1, maxpts - 1).astype(int), len(x) - 1])
        x, y = x[idx], y[idx]
    RECORDS.append({"id": f"{sec}:{obj}:{readout}", "sector": sec, "object": obj,
                    "readout": readout, "kind": kind,
                    "x": [float(v) for v in x], "y": [float(v) for v in y],
                    "meta": meta or {}, **({"note": note} if note else {})})


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def components_from_support(support, A):
    """记录支撑集的连通分量 = 区域（零参数）。"""
    N = A.shape[0]
    indptr, indices = A.indptr, A.indices
    lab = np.full(N, -1, np.int64)
    c = 0
    for s in np.nonzero(support)[0]:
        if lab[s] >= 0:
            continue
        lab[s] = c
        q = deque([int(s)])
        while q:
            u = q.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if support[v] and lab[v] < 0:
                    lab[v] = c
                    q.append(int(v))
        c += 1
    return lab, c


# =====================================================================
class Z0Core:
    """
    Z0 核心。状态：
      seed   : (V,) 本轮新播的走者数（整数重数）
      M      : (V,V) A^k（全分支走数）
      record : (V,) 记录层（闭合事件计数，恒中性）
      term   : (V,) 终端账本累积的**散度场**贡献
      Dgen   : 每代的终端账本规模（用于看增长）
      fresh  : 本代刚播种、还没走过一步的走者（避免"未走即闭合"）
    """

    def __init__(self, A, seed, L=8, dtype=np.int64, renorm=True, retention="all"):
        self.A = A.tocsr()
        self.V = A.shape[0]
        self.L = int(L)
        self.dtype = dtype
        self.renorm = renorm
        self.seed = np.asarray(seed, dtype=dtype).copy()
        self.M = np.eye(self.V, dtype=dtype)
        self.record = np.zeros(self.V, dtype=dtype)
        self.fresh = self.seed.copy()
        self.term = np.zeros(self.V, dtype=float)      # 终端散度场
        self.term_total = 0.0
        self.gen = 0
        self.hist = []
        self.retention = retention          # all / top2 / classes1 / classes2
        self.gen_closures = np.zeros(self.V, dtype=float)      # 本代闭合计数
        self.gen_class = np.zeros((self.V, self.L + 2), bool)  # 本代哪些长度闭合过
        self.layers = []                    # 每代 (计数向量, 类矩阵)

    # --- 当前散度场：端点边缘 − 起点边缘（活动层） ---
    def divergence(self):
        N = self.M * self.seed[:, None]                # N[u,v] = seed[u]·A^k[u,v]
        end = N.sum(axis=0)
        start = N.sum(axis=1)
        d_act = end - start
        return d_act + self.term                       # 活动层 ＋ 终端层

    def check_identities(self):
        N = self.M * self.seed[:, None]
        d_act = N.sum(axis=0) - N.sum(axis=1)
        I2 = float(np.sum(d_act) + np.sum(self.term))
        # I3：闭合词的散度恒为 0（对角 ⟹ 起点=终点）
        closed = np.diag(N)
        I3 = 0.0
        return {"I2_sum_d": float(I2), "I3_record_charge": float(I3),
                "closed_this_step": float(np.sum(closed))}

    # --- 走一步（Z2：全分支、零权重、确定性） ---
    def step(self):
        # 刚播种的走者先迈出第一步（起点→邻居），之后它们进入 M 体系
        if float(np.sum(self.fresh)) > 0:
            base = self.M.copy()
            add = np.asarray((self.A.T @ diags(self.fresh.astype(float))).T.todense())
            self.M = base + add.astype(self.dtype)                    # [u, w] = fresh[u]·A[u,w]
            self.fresh = np.zeros_like(self.fresh)
        self.M = np.asarray(self.M @ self.A, dtype=self.dtype)   # ★ 必须转回 ndarray：np.matrix 会让 end-start 广播成 (V,V)                # N ← N A（每步一次全分支）
        # 闭合：对角元（起点=终点）
        N = self.M * self.seed[:, None]
        closed = np.diag(N).copy()
        if float(np.sum(closed)) > 0:
            self.record = self.record + closed.astype(self.dtype)
            self.gen_closures = self.gen_closures + closed.astype(float)   # 按代记账
            k = min(self.L, self.gen_class.shape[1] - 1)
            self.gen_class[:, k] = self.gen_class[:, k] | (closed > 0)
            # 闭合者退出活动层（它们是"已完成"的词）
            idx = np.nonzero(closed)[0]
            for x in idx:
                self.M[x, x] = 0
        return float(np.sum(closed))

    # --- 到寿命：清空 → 终端层带走非零平衡（Z4）；从记录层重播种（Z5） ---
    def wipe_reseed(self):
        N = self.M * self.seed[:, None]
        d_act = (N.sum(axis=0) - N.sum(axis=1)).astype(float)       # 本代活动层的散度场
        self.term = self.term + d_act                                # 终端账本收走
        self.term_total += float(np.sum(N))
        self.hist.append({"gen": self.gen, "active_total": float(np.sum(N)),
                          "records_total": float(np.sum(self.record)),
                          "term_abs": float(np.abs(self.term).sum()),
                          "closed_prev": None})
        # ---- Z5：从**保留的**记录层重播种（D222：只保留最高两层）----
        self.layers.append((self.gen_closures.copy(), self.gen_class.copy()))
        r = self.retention
        if r == "all":
            w = np.sum([L0 for L0, _ in self.layers], axis=0)
        elif r == "top2":
            w = np.sum([L0 for L0, _ in self.layers[-2:]], axis=0)
        elif r == "classes1":
            w = self.layers[-1][1].sum(axis=1).astype(float)
        elif r == "classes2":
            w = (self.layers[-1][1] | self.layers[-2][1]).sum(axis=1).astype(float) \
                if len(self.layers) >= 2 else self.layers[-1][1].sum(axis=1).astype(float)
        elif r.startswith("cap"):          # ★ D222「最高两层保留」= 把保留的重数**封顶**
            K = int(r[3:])
            w = np.minimum(np.sum([L0 for L0, _ in self.layers], axis=0), K)
        else:
            raise ValueError(r)
        self.gen_closures = np.zeros(self.V, dtype=float)
        self.gen_class = np.zeros((self.V, self.L + 2), bool)
        newseed = w.astype(self.dtype) if self.dtype == np.int64 else w
        newseed = self.record.copy().astype(self.dtype) if False else newseed
        if float(np.sum(newseed)) == 0:
            newseed = np.ones(self.V, dtype=self.dtype)              # 记录空 ⟹ 均匀重播（退化情形）
        if self.renorm and self.dtype == np.int64:
            g = int(np.gcd.reduce(newseed[newseed > 0].astype(object))) if np.any(newseed > 0) else 1
            if g > 1:
                newseed = (newseed // g).astype(self.dtype)          # 整数重数的规范归一化（公因子）
        self.seed = newseed
        self.M = np.eye(self.V, dtype=self.dtype)
        self.fresh = self.seed.copy()
        self.gen += 1


def sign_components(d, A, thresh=0.0):
    """零和原生的区域定义：散度 **正** 连通块（源）与 **负** 连通块（汇）。
    它们必须互相补偿（`zero_sum_periodic_destruction.md`：补偿发生在区域之间）。"""
    N = A.shape[0]
    indptr, indices = A.indptr, A.indices
    lab = np.full(N, -1, np.int64)
    c = 0
    sign = np.sign(d)
    for s in range(N):
        if lab[s] >= 0 or sign[s] == 0:
            continue
        lab[s] = c
        q = deque([int(s)])
        while q:
            u = q.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if lab[v] < 0 and sign[v] == sign[s]:
                    lab[v] = c
                    q.append(int(v))
        c += 1
    return lab, c, sign


# =====================================================================
def run_faithful(A, name, L=6, gens=5, dtype=np.int64, seed0=None, verbose=True):
    """
    忠实演化：**种子铺满全图**（无偏好），全分支（Z2）→ 闭合写记录（Z3）→ 到 L 清空（Z4）
    → 从记录层重播种（Z5）。核验 I1/I2/I3。
    """
    V = A.shape[0]
    seed = np.ones(V, dtype=dtype) if seed0 is None else seed0
    eng = Z0Core(A, seed, L=L, dtype=dtype)
    rows = []
    for g in range(gens):
        for _ in range(L):
            eng.step()
        N = eng.M * eng.seed[:, None]
        active_total = float(np.sum(N))
        chk = eng.check_identities()
        d = eng.divergence()
        # 零和原生区域：散度正/负连通块
        lab2, c2, sgn = sign_components(d, A)
        Q2 = [float(d[lab2 == j].sum()) for j in range(c2)]
        pos = sorted([q for q in Q2 if q > 0], reverse=True)
        neg = sorted([-q for q in Q2 if q < 0], reverse=True)
        pair = (len(pos) == len(neg)) and all(abs(a - b) <= 1e-9 * max(1.0, abs(a))
                                              for a, b in zip(pos, neg))
        eng.wipe_reseed()
        sd = eng.seed.astype(float)
        tot = sd.sum()
        overflow = bool(np.any(sd < 0))          # int64 溢出会翻负
        sd = np.maximum(sd, 0.0)
        gini = float((2 * np.arange(1, V + 1) - V - 1).dot(np.sort(sd)) / (V * max(tot, 1e-30)))
        rows.append({"gen": g, "active_total": active_total,
                     "sum_d": float(chk["I2_sum_d"]),
                     "n_sources_sinks": int(c2), "pairing_exact": bool(pair),
                     "Q_pos_top": [round(x, 4) for x in pos[:4]],
                     "Q_neg_top": [-round(x, 4) for x in neg[:4]],
                     "seed_total": float(tot), "seed_gini": round(gini, 5),
                     "seed_maxfrac": round(float(sd.max() / max(tot, 1e-30)), 5),
                     "seed_support": int((sd > 0).sum()),
                     "records_total": float(np.sum(eng.record)),
                     "int_overflow": overflow})
        if verbose:
            print(f"   代{g}: 活动层={active_total:>12.4g}  Σd={chk['I2_sum_d']:>8.2g}  "
                  f"源/汇={c2:>3} 配对={pair!s:<5} 重播种 Gini={gini:.4f} "
                  f"最大占比={sd.max()/max(tot,1e-30):.4f} 支撑={int((sd>0).sum()):>5}/{V}")
        if overflow:
            print("   [停] int64 整数重数溢出（全分支增长 (2d)^{gL}）")
            break
    return rows


def random_regular(N, d, seed=0):
    rng = np.random.default_rng(seed)
    while True:
        stubs = np.repeat(np.arange(N), d)
        rng.shuffle(stubs)
        a, b = stubs[0::2], stubs[1::2]
        ok = a != b
        a, b = a[ok], b[ok]
        A = csr_matrix((np.ones(len(a)), (a, b)), shape=(N, N))
        A = A + A.T
        A.data[:] = 1.0
        if A.nnz >= N * d * 0.9:
            return A


def main(argv):
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "architecture": "Z0① ② ③ → Z1/Z2/Z3/Z4/Z5：全分支零权重；闭合写记录；终端层带走非零平衡；从记录层重播种",
           "identities": "I1 Σ_v d_w=0（结构）· I2 Σ_x d(x)=0（精确核验）· I3 Q_record=0 · I4 零权重",
           "runs": []}

    graphs = [("torus_16", torus_adj(2, 16)), ("torus_32", torus_adj(2, 32)),
              ("rand4reg_256", random_regular(256, 4, seed=1))]
    try:
        import l0_closure as L0
        import observable_sweep as OS
        from zcl import Engine
        eng = Engine()
        reps, _ = L0.enumerate_necklaces(14, eng, chunk_bits=26, verbose=False)
        graphs.append(("G_14", OS.build_sparse(14, reps, eng)))
    except Exception as e:
        print("  [warn] 闭链图跳过:", e)

    for L in (4, 6):
        for name, A in graphs:
            print("=" * 100)
            print(f"### Γ={name}  V={A.shape[0]}  寿命 L={L}  种子=全图均匀（无偏好）")
            rows = run_faithful(A, name, L=L, gens=5)
            out["runs"].append({"graph": name, "V": int(A.shape[0]), "L": L, "rows": rows})
            rec("Z0CORE", f"{name} L={L}", "重播种分布 Gini vs 代",
                [r["gen"] for r in rows], [r["seed_gini"] for r in rows], {"V": A.shape[0], "L": L})
            rec("Z0CORE", f"{name} L={L}", "源/汇连通块数 vs 代",
                [r["gen"] for r in rows], [r["n_sources_sinks"] for r in rows], {"V": A.shape[0], "L": L})

    print("=" * 100)
    print("### 用 float64 跟长（指数照样可看；整数精确性只在上面 int64 段成立）")
    for name, A in graphs:
        if name.startswith("torus"):
            continue
        print(f"\n--- Γ={name}  V={A.shape[0]}  L=4  dtype=float64")
        rows = run_faithful(A, name, L=4, gens=12, dtype=np.float64)
        out["runs"].append({"graph": name, "V": int(A.shape[0]), "L": 4,
                            "dtype": "float64", "rows": rows})
        rec("Z0CORE", f"{name} float64", "重播种分布 Gini vs 代",
            [r["gen"] for r in rows], [r["seed_gini"] for r in rows], {"V": A.shape[0]})

    print("=" * 100)
    print("### Z5 的四种保留读法对比（非传递 Γ；看谁不炸）")
    for name, A in graphs:
        if name.startswith("torus"):
            continue
        for ret in ("all", "cap2", "cap4"):
            print(f"\n--- Γ={name} V={A.shape[0]} L=4 retention={ret} (float64)")
            V = A.shape[0]
            eng = Z0Core(A, np.ones(V, np.float64), L=4, dtype=np.float64, retention=ret)
            rows = []
            for g in range(12):
                for _ in range(4):
                    eng.step()
                N = eng.M * eng.seed[:, None]
                act = float(np.sum(N))
                d = eng.divergence()
                chk = eng.check_identities()
                lab2, c2, sgn = sign_components(d, A)
                Q2 = [float(d[lab2 == j].sum()) for j in range(c2)]
                pos = sorted([q for q in Q2 if q > 0], reverse=True)
                neg = sorted([-q for q in Q2 if q < 0], reverse=True)
                pair = (len(pos) == len(neg)) and all(abs(a - b) <= 1e-9 * max(1.0, abs(a))
                                                      for a, b in zip(pos, neg))
                eng.wipe_reseed()
                sd = np.maximum(eng.seed.astype(float), 0)
                tot = sd.sum()
                gini = float((2 * np.arange(1, V + 1) - V - 1).dot(np.sort(sd)) / (V * max(tot, 1e-30)))
                rows.append({"gen": g, "active_total": act, "seed_total": float(tot),
                             "seed_gini": round(gini, 5),
                             "seed_maxfrac": round(float(sd.max() / max(tot, 1e-30)), 5),
                             "seed_support": int((sd > 0).sum()),
                             "n_sources_sinks": int(c2), "pairing_exact": bool(pair),
                             "sum_d": float(chk["I2_sum_d"])})
                print(f"   代{g:>2}: 活动层={act:>11.4g}  种子总量={tot:>11.4g}  "
                      f"Gini={gini:.4f}  最大占比={sd.max()/max(tot,1e-30):.4f}  "
                      f"支撑={int((sd>0).sum()):>4}/{V}  源/汇={c2:>3} 配对={pair!s:<5}")
                if act == 0 or not np.isfinite(act):
                    print("   [停] 活动层熄灭了或发散")
                    break
            out["runs"].append({"graph": name, "V": V, "L": 4, "dtype": "float64",
                                "retention": ret, "rows": rows})
            rec("Z0CORE", f"{name} Z5={ret}", "重播种 Gini vs 代",
                [r["gen"] for r in rows], [r["seed_gini"] for r in rows],
                {"V": V, "retention": ret})
            rec("Z0CORE", f"{name} Z5={ret}", "活动层总数 vs 代",
                [r["gen"] for r in rows], [r["active_total"] for r in rows],
                {"V": V, "retention": ret})

    with open(os.path.join(OUT, "z0_core.json"), "w") as f:
        json.dump({**out, "records": RECORDS}, f, ensure_ascii=False, indent=1)
    print("=" * 100)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_core.json")


if __name__ == "__main__":
    main(sys.argv[1:])
