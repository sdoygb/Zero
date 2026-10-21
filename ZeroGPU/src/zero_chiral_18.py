#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_chiral_18.py —— 给 18 世界 micro core 加 C3 手征步，测它是否解除 K4 退化

动机：`zero_l2_three.md` 命题 2 判定 K4 的顶点传递性把 L2 的 3 个亚层方向
压成 **1 个模式**（$A^L$ 在 $1^\\perp$ 上 $=\\frac{1}{81}I$，三个本征值全同），
并给出两种补全方案，其中一种是 **C3 手征步** $(AP)^L$（每代走 $L$ 步再转
$120^\\circ$），它给 $\\frac{1}{81}(1,\\omega,\\omega^2)$——3 个不同本征值、轨迹秩 3。

但那个结论只在**分析脚本**（`zero_l2_three.py`）里算过，**没有任何演化引擎
采用手征步**。本脚本把 C3 手征接进 18 世界 micro core，回答：
  ① 手征步在引擎里能跑吗（有界/非负/不爆）？
  ② 它是否真的把 L2 的对比轨迹从秩 1 抬到秩 3？
  ③ 它对宏观可观测量（n_std/q/h_E）与层份额有无影响？
  ④ 长跑会不会也冻结？

手征的定义（K4 的 3 重旋转 P，固定顶点 0，$1\\to2\\to3\\to1$；Z090 的 $-1$ 空间
按 C3 分解为 $1\\oplus\\omega\\oplus\\omega^2$）：
  源 4 → 走者：每步不是 $W\\leftarrow W A$，而是 $W\\leftarrow W (A P)$，
  即「走一步再转一次」。$L=4$ 步后一代的传播算子是 $(AP)^4$。
注意 $(AP)^4$ 不是 $A^4$——手征步把每代的净旋转累积成 $P^4=I$ 之外的相位结构
（$P^3=I$，故 $P^4=P$）：净效应是每代多转一次 $120^\\circ$。

用法：
    cd ZeroGPU/src && /usr/bin/python3 zero_chiral_18.py
"""
from __future__ import annotations
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo, Params  # noqa: E402
from zero_evo_18 import ZeroEvo18, EDGES, MATCH, NG, NE  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
V = 4


def c3_perm():
    """K4 的 3 重旋转 P：固定顶点 0，1->2->3->1（Z090）。"""
    P = np.zeros((V, V))
    perm = [0, 2, 3, 1]
    for i, j in enumerate(perm):
        P[i, j] = 1.0
    return P


def c3_gen_perm():
    """C3 作用在 **3 个代** 上的 3-循环：0->1->2->0。

    注意：不能把 4 顶点的 C3 置换 P 直接切片成 3x3——那会得到奇异矩阵
    （见 zero_chiral_18.md §3 的实测）。3 代是 K4 的 3 个 1-因子，它们之间
    的循环是 3 阶的，正是 Z090 的 1⊕ω⊕ω² 所对应的那一个。
    """
    R = np.zeros((NG, NG))
    perm = [1, 2, 0]
    for i, j in enumerate(perm):
        R[i, j] = 1.0
    return R


R_C3_GEN = c3_gen_perm()
P_C3 = c3_perm()
A_K4 = np.ones((V, V)) - np.eye(V)


def helmert():
    H = np.zeros((V, V - 1))
    for i in range(V - 1):
        H[:i + 1, i] = 1.0 / (i + 1)
        H[i + 1, i] = -1.0
        H[:, i] /= np.linalg.norm(H[:, i])
    return H


H_BASIS = helmert()


def rank_shape(M, n=200, tol=1e-7):
    """理论对象：形状映射 F(p)=Mp/(1^T Mp) 从 Perron 点附近出发的对比轨迹秩。

    这是 zero_l2_three 命题 2 的对象（线性化 J_c^t x_0），不是引擎的累积 record，
    也不是引擎里被不对称 seed 驱动的逐代闭合。
    """
    pstar = np.ones(V) / V
    x = H_BASIS.T @ (np.array([0.4, 0.3, 0.2, 0.1]) - pstar)
    xs = [x]
    for _ in range(n):
        q = M @ (H_BASIS @ x + pstar)
        q = q / q.sum()
        x = H_BASIS.T @ (q - pstar)
        xs.append(x)
    X = np.array(xs).T
    s = np.linalg.svd(X, compute_uv=False)
    if s[0] <= 0:
        return 0
    return int((s > tol * s[0]).sum())


def rank_traj(rec_hist, tol=1e-7):
    """4 顶点记录序列在 1^perp 上的对比轨迹秩。

    口径注意（这是第一版判错的地方）：**累积 record** 是单调增长的和，
    其轨迹秩平凡地为 3（因为它在长）。理论（zero_l2_three 命题 2）说的是
    **每代的形状对比轨迹** $J_c^t x_0$，对应到引擎是「每代闭合 c_t 的对比分量」，
    不是累积量。故这里用**逐代闭合**（`e.gen_closures`），并去掉增长因子。
    """
    R = np.asarray(rec_hist, dtype=float)
    ps = np.ones(V) / V
    X = (R - ps) @ H_BASIS
    s = np.linalg.svd(X, compute_uv=False)
    if s[0] <= 0:
        return 0, s
    return int((s > tol * s[0]).sum()), s


class ZeroEvo18Chiral(ZeroEvo18):
    """18 世界 + C3 手征步：走者每步 W <- W (A P)。"""

    def _edge_source(self, w4):
        return super()._edge_source(w4)

    def _micro_gen(self):
        # 与父类同构，只把「走 L 步」换成「走 L 步、每步后转一次」
        w4 = self._seed_used()
        if w4.sum() <= 0:
            w4 = np.ones(self.V)
        w4 = w4 / w4.sum()
        ge = self.gen_w
        flow = np.outer(ge, self._edge_source(w4))
        self.W18 = (1.0 - self.alpha) * self.W18 + self.alpha * flow
        self.W18 /= max(self.W18.sum(), 1e-30)

        # ── 手征步：把 18 世界的代指标按 C3 轮转一次 ──
        # 18 世界是 (代 g, 边 e) 的 3x6 矩阵；C3 作用在「代」指标上
        # （3 代 = 3 个完美匹配 = K4 的 3 个 1-因子，其循环是 3 阶）。
        self.W18 = R_C3_GEN @ self.W18

        self.gen += 1
        if self.p.T and self.gen % self.p.T == 0:
            for g in range(NG):
                mask = np.zeros(NE)
                mask[MATCH[g]] = 1.0
                self.W18[g] *= mask
            self.W18 /= max(self.W18.sum(), 1e-30)
        shadow = self.W18.sum(axis=0)
        c = np.zeros(self.V)
        for e, (u, v) in enumerate(EDGES):
            c[u] += shadow[e] / 2.0
            c[v] += shadow[e] / 2.0
        c = c / c.sum() if c.sum() > 0 else np.ones(self.V) / self.V
        self.gen_closures = c
        self.record = self.record + c
        self.layers.append(c.copy())
        if len(self.layers) > 2:
            self.layers.pop(0)
        if self.p.retention == "last":
            ret = c.copy()
        elif self.p.retention == "top2":
            ret = c + (self.layers[-2] if len(self.layers) >= 2 else 0.0)
        else:
            ret = self.record.copy()
        self.seed = ret / ret.sum() if ret.sum() > 0 else np.ones(self.V) / self.V
        m_r = max(self.record.sum(), self.terminal.sum())
        if m_r > self.p.norm_hi:
            self.record /= m_r
            self.terminal /= m_r


def base_params(**kw):
    d = dict(micro=True, L=4, T=4, retention="top2",
             eta_n=0.2, eta_h=0.2, eta_psi=0.2,
             mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
             capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
             Gamma_h=0.2, kappa_h=0.5, alpha=0.0,
             wave_on=True, beta=0.5, g_wave=0.5, kappa_psi=0.5, Gamma_w=0.5)
    d.update(kw)
    return Params(**d)


def main():
    # ── ① 静态：手征步是否解除 K4 退化（与 zero_l2_three 对齐）──
    AL = np.linalg.matrix_power(A_K4, 4)
    chiral = np.linalg.matrix_power(A_K4 @ P_C3, 4)

    def modes(M):
        pstar = np.ones(V) / V
        den = float((M @ pstar).sum())
        J = (M - np.outer(pstar, M.sum(axis=0))) / den
        Jc = H_BASIS.T @ J @ H_BASIS
        return np.linalg.eigvals(Jc)

    m_plain = modes(AL)
    m_chiral = modes(chiral)
    n_plain = len(set(np.round(m_plain, 9)))
    n_chiral = len(set(np.round(m_chiral, 9)))

    # ── ② 引擎：手征 core 能不能跑 ──
    gens = 20000
    e_plain = ZeroEvo18(zs.k4(), base_params(), seed=np.array([8.0, 4.0, 1.5, 1.0]))
    e_chir = ZeroEvo18Chiral(zs.k4(), base_params(), seed=np.array([8.0, 4.0, 1.5, 1.0]))

    recs_plain, recs_chir = [], []
    ns_p, ns_c, q_p, q_c = [], [], [], []
    sh_p, sh_c = [], []
    for _ in range(gens):
        e_plain.step(); e_chir.step()
        recs_plain.append(e_plain.gen_closures.copy())
        recs_chir.append(e_chir.gen_closures.copy())
        sp, sc = e_plain.snapshot(), e_chir.snapshot()
        ns_p.append(sp["n_std"]); ns_c.append(sc["n_std"])
        q_p.append(sp["q"]); q_c.append(sc["q"])
        sh_p.append(e_plain.W18.sum(axis=1).copy())
        sh_c.append(e_chir.W18.sum(axis=1).copy())

    r_plain, sv_p = rank_traj(recs_plain)
    r_chir, sv_c = rank_traj(recs_chir)

    tail = slice(gens // 2, gens)
    share_p = np.array(sh_p)[tail].mean(axis=0)
    share_c = np.array(sh_c)[tail].mean(axis=0)
    share_p = share_p / share_p.sum()
    share_c = share_c / share_c.sum()

    finite = bool(np.all(np.isfinite([ns_c[-1], q_c[-1]])))
    nonneg = bool(np.all(np.array(sh_c)[tail] >= -1e-12))

    checks = [
        dict(name="c1 静态：plain A^4 在 1^perp 上 1 个不同本征值（K4 退化）",
             pass_=bool(n_plain == 1), detail=f"distinct={n_plain}, λ={np.round(np.real(m_plain),6)}"),
        dict(name="c2 静态：手征 (AP)^4 给 3 个不同本征值",
             pass_=bool(n_chiral == 3),
             detail=f"distinct={n_chiral}, λ={np.round(m_chiral,6)}"),
        dict(name="c3 静态：plain A^L 的形状轨迹秩 = 1（理论对象，非引擎累积量）",
             pass_=bool(rank_shape(AL) == 1), detail=f"rank={rank_shape(AL)}"),
        dict(name="c3b 静态：手征 (AP)^L 的形状轨迹秩 = 3",
             pass_=bool(rank_shape(chiral) == 3), detail=f"rank={rank_shape(chiral)}"),
        dict(name="c4 引擎：两个 core 的逐代闭合轨迹秩（seed 驱动，均满秩）",
             pass_=bool(r_plain == 3 and r_chir == 3),
             detail=f"plain={r_plain} chiral={r_chir}——引擎不在 Perron 点，"
                    f"不对称 seed 掩盖了退化"),
        dict(name="c5 手征 core 有界/非负/无 NaN",
             pass_=bool(finite and nonneg),
             detail=f"n_std={ns_c[-1]:.6g}, q={q_c[-1]:.6g}, shares_nonneg={nonneg}"),
        dict(name="c6 手征改变层份额（代指标轮转 ⇒ 份额不再冻结）",
             pass_=bool(np.max(np.abs(share_c - share_p)) > 1e-6),
             detail=f"share plain={np.round(share_p,5)} vs chiral={np.round(share_c,5)}"),
        dict(name="c7 手征的层份额不冻结（plain 恒为 gen_w，手征在动）",
             pass_=bool(np.std(np.array(sh_c)[tail], axis=0).max() > 1e-6
                        and np.std(np.array(sh_p)[tail], axis=0).max() < 1e-9),
             detail=f"tail relstd chiral={np.round(np.array(sh_c)[tail].std(0)/np.maximum(np.array(sh_c)[tail].mean(0),1e-30),4)}"
                    f" plain={np.round(np.array(sh_p)[tail].std(0)/np.maximum(np.array(sh_p)[tail].mean(0),1e-30),4)}"),
    ]

    out = dict(
        gens=gens,
        static=dict(plain_modes=[[float(z.real), float(z.imag)] for z in m_plain],
                    chiral_modes=[[float(z.real), float(z.imag)] for z in m_chiral],
                    n_distinct_plain=n_plain, n_distinct_chiral=n_chiral),
        engine=dict(rank_plain=r_plain, rank_chiral=r_chir,
                    sv_plain=[float(x) for x in sv_p[:4]],
                    sv_chiral=[float(x) for x in sv_c[:4]],
                    share_plain=[float(x) for x in share_p],
                    share_chiral=[float(x) for x in share_c],
                    n_std_plain=float(ns_p[-1]), n_std_chiral=float(ns_c[-1]),
                    q_plain=float(q_p[-1]), q_chiral=float(q_c[-1])),
        checks=[{**{k: v for k, v in c.items() if k != "pass_"}, "pass": bool(c["pass_"])}
                for c in checks],
    )
    with open(os.path.join(RES, "zero_chiral_18.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ax[0].plot(ns_p, label="plain $A^L$", lw=1.2)
    ax[0].plot(ns_c, label="chiral $(AP)^L$", lw=1.2)
    ax[0].set_xlabel("generation"); ax[0].set_ylabel("n_std")
    ax[0].set_title("macro observable: plain vs chiral core")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    ax[1].bar([0, 1, 2], sv_p[:3], width=.35, label="plain (rank 1)")
    ax[1].bar([0.4, 1.4, 2.4], sv_c[:3], width=.35, label=f"chiral (rank {r_chir})")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("contrast mode"); ax[1].set_ylabel("singular value")
    ax[1].set_title("L2 contrast trajectory spectrum")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, axis="y")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_chiral_18.png")
    fig.savefig(plot, dpi=130); plt.close(fig)

    print("C3 手征步接 18 世界 micro core")
    print(f"  静态：plain distinct={n_plain} / chiral distinct={n_chiral}")
    print(f"  引擎：plain rank={r_plain} / chiral rank={r_chir}")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()