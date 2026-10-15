# -*- coding: utf-8 -*-
"""
ZeroGPU · 稳态谱计算（能否产生粒子，先不管是什么粒子）

实现"稳态的谱"操作定义：
  ① 锁稳态：演化到 ω(t)=record/Σrecord 收敛（q=5/9 不动点）；
  ② 数值雅可比：在稳态，扰动重数 ω_j，演化 m 代，测响应 → J(V×V)；
  ③ 本征分解 J v_k = λ_k v_k，按 |λ|、arg λ、空间局域性分类；
  ④ 粒子候选 = 边缘稳定（|λ|≈1）+ 局域（v_k 集中少数顶点）的模式。

数据保留到 results/。

用法：
    /usr/local/bin/python3 spectrum.py            # path16 geom last L4
    /usr/local/bin/python3 spectrum.py path64     # path64 版本（更慢）
"""
from __future__ import annotations
import sys, os, time, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

RES = os.path.normpath(os.path.join(HERE, "..", "results"))


class SpectrumAnalyzer:
    """在锁定的稳态上计算谱（数值雅可比 + 本征分解）。"""

    def __init__(self, adj, seed, retention="last", L=4, n_steady=60000):
        self.vec = ZeroCoreVec(adj, seed=seed, retention=retention, L=L)
        self.V = len(adj)
        self.n_steady = n_steady

    def lock_steady(self):
        for _ in range(self.n_steady):
            self.vec.run_generation()
        self.omega_star = self.vec.q_weighted(self.vec.record)
        self.record_star = self.vec.record.copy()
        self.W_star = self.vec.W.copy()
        self.seed_star = self.vec.seed.copy()
        return self

    def q_w(self):
        return self.vec.q_weighted(self.vec.record)

    def reset_to_steady(self):
        self.vec.record = self.record_star.copy()
        self.vec.W = self.W_star.copy()
        self.vec.seed = self.seed_star.copy()

    def jacobian(self, eps=1e-4, m=500):
        """数值雅可比：扰动 ω_j，演化 m 代，测响应（双演化差分离本底）。"""
        V = self.V
        J = np.zeros((V, V))
        # 基准（无扰动，分离本底累积）
        self.reset_to_steady()
        for _ in range(m):
            self.vec.run_generation()
        base = self.q_w()
        for j in range(V):
            self.reset_to_steady()
            self.vec.record[j] *= (1.0 + eps)
            for _ in range(m):
                self.vec.run_generation()
            pert = self.q_w()
            J[:, j] = (pert - base) / eps
        return J

    def spectrum(self, eps=1e-4, m=500):
        """返回 (eigvals, eigvecs, J)。"""
        J = self.jacobian(eps, m)
        eigvals, eigvecs = np.linalg.eig(J)
        return eigvals, eigvecs, J

    @staticmethod
    def classify(eigvals, eigvecs, stable_tol=0.05):
        """按 |λ|、arg λ、局域性分类。"""
        V = len(eigvals)
        rows = []
        for k in range(V):
            lam = eigvals[k]
            v = np.abs(eigvecs[:, k])
            v = v / v.sum()                       # 归一化（局域性）
            gini = SpectrumAnalyzer.gini_of(v)
            mx = v.max()                          # 最大分量占比
            absl = abs(lam)
            argl = np.angle(lam)
            if abs(absl - 1.0) < stable_tol:
                kind = "stable"                   # 边缘稳定（|λ|≈1）
            elif absl < 1.0:
                kind = "decay"                    # 衰减
            else:
                kind = "unstable"                 # 不稳定
            osc = abs(argl) > 0.1                 # 振荡（相位转）
            local = mx > 0.5 or gini > 0.5        # 局域（集中少数顶点）
            rows.append(dict(k=k, lam=f"{lam.real:.4f}{lam.imag:+.4f}i",
                             abs=round(absl, 4), arg=round(argl, 4),
                             kind=kind, osc=bool(osc), gini=round(gini, 3),
                             mx=round(mx, 3), local=bool(local)))
        return rows

    @staticmethod
    def gini_of(x):
        s = sorted(float(v) for v in x)
        n = len(s)
        if n == 0 or sum(s) == 0:
            return 0.0
        ssum = sum(s)
        num = sum((2 * (i + 1) - n - 1) * s[i] for i in range(n))
        return num / (n * ssum)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "path16"
    n = int(which.replace("path", ""))
    n_steady = 60000 if n <= 16 else 80000
    adj = path_adj(n)
    V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]

    print("=" * 84)
    print(f"稳态谱计算：path{n} geom last L4（锁稳态 {n_steady} 代）")
    print("=" * 84)
    t0 = time.time()
    sa = SpectrumAnalyzer(adj, seed, retention="last", L=4, n_steady=n_steady)
    sa.lock_steady()
    q = sa.omega_star
    print(f"稳态锁定：q=5/9?  q={q:.6f}  差|q-5/9|={abs(q-5/9):.2e}  ({time.time()-t0:.1f}s)")
    t0 = time.time()
    eigvals, eigvecs, J = sa.spectrum(eps=1e-4, m=500)
    print(f"雅可比+本征分解完成 ({time.time()-t0:.1f}s, V={V})")
    rows = sa.classify(eigvals, eigvecs)

    # 分类统计
    stable = [r for r in rows if r["kind"] == "stable"]
    stable_local = [r for r in stable if r["local"]]
    decay = [r for r in rows if r["kind"] == "decay"]
    unstable = [r for r in rows if r["kind"] == "unstable"]
    osc = [r for r in rows if r["osc"]]
    print("\n模式分类：")
    print(f"  边缘稳定 |λ|≈1 : {len(stable)} 个  其中局域(粒子候选): {len(stable_local)} 个")
    print(f"  衰减 |λ|<1     : {len(decay)} 个")
    print(f"  不稳定 |λ|>1   : {len(unstable)} 个")
    print(f"  振荡 argλ≠0    : {len(osc)} 个")

    print("\n全部模式（k, |λ|, arg, 类型, 局域性 Gini/max）：")
    for r in sorted(rows, key=lambda x: abs(x["abs"] - 1.0)):
        tag = "★粒子候选" if r["kind"] == "stable" and r["local"] else r["kind"]
        osc_mark = " 振荡" if r["osc"] else ""
        print(f"  k={r['k']:>2}  |λ|={r['abs']:.4f}  arg={r['arg']:+.3f}  "
              f"{tag}{osc_mark}   Gini={r['gini']:.3f} mx={r['mx']:.3f}")

    out = {"graph": f"path{n}", "n_steady": n_steady, "q_steady": q,
           "mode_counts": {"stable": len(stable), "stable_local": len(stable_local),
                           "decay": len(decay), "unstable": len(unstable), "osc": len(osc)},
           "modes": rows}
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, f"spectrum_path{n}.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {RES}/spectrum_path{n}.json")


if __name__ == "__main__":
    main()
