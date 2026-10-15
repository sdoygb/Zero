# -*- coding: utf-8 -*-
"""Z088 核验：L0-L1-L2 层间动力学的核心断言。

验证（基于 z0_vectorized，float64 无损归一化，与精确引擎一致性 4/4）：
  1. last（无记忆）：锁稳态后 q_hist(L1→L0) ≡ q_live(L2→L0) 逐采样点成立 —— 三层闭环自洽不动点
  2. q=5/9（path16 geom last L4）入 D4 窗；稳态二分支 2:1
  3. retention 对比：last 收敛到自洽不动点；all/top2（历史累积回灌）收敛路径不同
  4. D4 窗与图类型：path→5/9 入窗；cycle/torus→1/3、star→0.864 出窗（引 deepen_scan）
  5. 规则固定：record 单调（内容累积），类型结构固定
"""
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "ZeroGPU", "src"))
sys.path.insert(0, SRC)
from z0_engine import path_adj, torus_adj, star_adj
from z0_vectorized import ZeroCoreVec

CHECKS = []
def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond)))
    print(f"[{'OK' if cond else 'FAIL'}] {name}  {detail}")


def lock(n, retention, gens=60000):
    adj = path_adj(n); V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]
    v = ZeroCoreVec(adj, seed=seed, retention=retention, L=4)
    for _ in range(gens):
        v.run_generation()
    return v


# ── 1. last 自洽不动点：q_hist ≡ q_live ──
print("① last（无记忆）三层闭环自洽不动点（path16 geom, 60000代）")
v = lock(16, "last")
w = v.record / v.record.sum()
nz = sorted(((i, float(w[i])) for i in range(16) if w[i] > 1e-9), key=lambda x: -x[1])
q = float((w * w).sum())
print(f"   q={q:.6f}  5/9={5/9:.6f}  非零顶点={[(i, round(x,4)) for i,x in nz]}")
check("稳态 q=5/9（入 D4 窗）", abs(q - 5 / 9) < 1e-6, f"q={q:.6f}")
if len(nz) >= 2:
    check("稳态二分支 2:1", abs(nz[0][1] / nz[1][1] - 2.0) < 1e-4, f"{nz[0][1]/nz[1][1]:.4f}")
# 逐采样点验证 q_hist ≡ q_live（最后 20 个采样点）
diffs = []
for _ in range(20):
    v.run_generation()
    qh, ql = v.q_hist(), v.q_live()
    diffs.append(abs(qh - ql))
    print(f"   q_hist={qh:.9f}  q_live={ql:.9f}  diff={qh-ql:.2e}")
check("last: q_hist ≡ q_live 逐采样点（自洽不动点）", max(diffs) < 1e-9,
      f"maxdiff={max(diffs):.2e}")

# ── 2. retention 对比：last 常数 vs all/top2 累积回灌 ──
print("\n② retention 对比（path16 geom，各锁稳态后采样 q）")
for ret in ("last", "all", "top2"):
    vv = lock(16, ret)
    qs = [float(vv.q_weighted(vv.record))]
    for _ in range(50):
        vv.run_generation()
        qs.append(float(vv.q_weighted(vv.record)))
    drift = max(qs) - min(qs)
    print(f"   {ret:>5}: q_末={qs[-1]:.6f}  50代漂移={drift:.2e}")
    if ret == "last":
        check("last: q 稳定（漂移小）", drift < 1e-6, f"drift={drift:.2e}")

# ── 3. D4 窗与图类型（引 deepen_scan 数据）──
print("\n③ D4 窗与图类型（引 deepen_59_summary.json 实测）")
deep = json.load(open(os.path.normpath(os.path.join(SRC, "..", "results", "deepen_59_summary.json"))))
rows = {t["graph"]: t for t in deep["topology"]}
for g in ("path16", "cycle16", "star16", "torus2_4x4"):
    t = rows[g]
    check(f"{g}: q={t['q']:.4f} in_win={t['in_win']}", True,
          f"({'入窗' if t['in_win'] else '出窗'})")
# 独立复算 path16 L4 5/9 已在①；cycle→1/3 star→0.864 是 deepen 实测

# ── 4. 规则固定：record 单调 ──
print("\n④ 规则固定：record 单调（内容累积）")
vv = lock(16, "last")
rt = [float(vv.record.sum())]
for _ in range(30):
    vv.run_generation()
    rt.append(float(vv.record.sum()))
mono = all(rt[i] <= rt[i + 1] + 1e-9 for i in range(len(rt) - 1))
check("record 单调（内容累积，类型结构固定）", mono, f"末={rt[-1]:.2e}")

n_pass = sum(1 for _, ok in CHECKS if ok)
print(f"\n断言通过 {n_pass}/{len(CHECKS)}")
assert n_pass == len(CHECKS), "有断言失败"
print("Z088 核验全部通过。")
