#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_l2_slots.py —— L2 三亚层的完整动力学：链式交换、3 个模式与相互影响
==========================================================================
上一轮（`zero_l2_three.py`）判定：K4 上 4 顶点占据的 3 个自由度在**线性行走**
下本征值全同（$1/81$），只有 1 个模式——因为 K4 顶点传递，$A^L$ 在 $1^\\perp$ 上
是标量。

本轮把 L2 的 3 个亚层当成**独立的 3 维动力学对象**，给出**完整**（闭合、3 模、
显式相互影响）的规律：

$$
\\frac{d f}{dt}=-\\kappa\\,BB^{\\mathsf T}(f-f_{\\rm seed}),
\\qquad f=B\\log p,
\\qquad
BB^{\\mathsf T}=
\\begin{pmatrix}2&-1&0\\\\-1&2&-1\\\\0&-1&2\\end{pmatrix}.
$$

$B$ 是路径 $0-1-2-3$ 的边界算子（相邻差）；$f_i=\\log(p_i/p_{i+1})$ 是第 $i$ 个
亚层的（带符号）对数比。三点：

  1. **闭合**：$f$ 是 3 维的，方程只含 $f$、$f_{\\rm seed}$、$\\kappa$。3 个亚层
     通过 $BB^{\\mathsf T}$ 的**非对角 $-1$** 直接交换（0↔1、1↔2），0 与 2 不直接
     耦合；但静态 Green 函数 $G=(BB^{\\mathsf T})^{-1}$ 是**满支撑**的
     （$G_{02}=1/4\\ne0$）——0 与 2 通过 1 间接影响。
  2. **3 个模式**：$BB^{\\mathsf T}$ 的本征值 $2-\\sqrt2,\\;2,\\;2+\\sqrt2$
     （$0.586,2,3.414$），三个模式是
     对称慢模 $(1,\\sqrt2,1)$、反对称模 $(1,0,-1)$、交替快模 $(1,-\\sqrt2,1)$。
     弛豫率 $\\kappa(2-\\sqrt2),\\;2\\kappa,\\;\\kappa(2+\\sqrt2)$ **全不同**。
  3. **输入只有一条**：3 个亚层的有效耦合图是**链**（相邻亚层共享一个顶点，
     形状相关的耦合只能沿链走）。$B$ 与 $BB^{\\mathsf T}$ 由链**唯一确定**，
     没有可调矩阵；唯一自由参数是速率 $\\kappa$（和 L1 种子给出的目标
     $f_{\\rm seed}$）。

对照：若把 4 顶点当作完整的 K4（所有对都耦合），槽投影是
$B L_{K4}B^{\\mathsf T}(BB^{\\mathsf T})^{-1}=4I$——**又是标量**，1 个模式。
所以"链"不是随手挑的：它是 K4 投影里**唯一非标量的形状部分**。

输出：results/zero_l2_slots.json、zero_l2_slots.png

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_l2_slots.py
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

RES = os.path.join(HERE, "..", "results")
SEED = np.array([8.0, 4.0, 1.5, 1.0])


def boundary(n: int = 4) -> np.ndarray:
    """路径 0-1-...-(n-1) 的边界算子 B（(n-1) x n）。"""
    B = np.zeros((n - 1, n))
    for i in range(n - 1):
        B[i, i] = 1.0
        B[i, i + 1] = -1.0
    return B


B = boundary(4)
M = B @ B.T
LAM, VEC = np.linalg.eigh(M)
G = np.linalg.inv(M)                    # 静态 Green 函数（/kappa）


def f_of(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    return B @ np.log(p)


def response(t: float, kappa: float = 1.0) -> np.ndarray:
    """R(t)=exp(-kappa M t)。"""
    return VEC @ np.diag(np.exp(-kappa * LAM * t)) @ VEC.T


def relax(f0, fstar, kappa=1.0, t=20.0):
    return fstar + response(t, kappa) @ (f0 - fstar)


def impulse_peaks() -> np.ndarray:
    ts = np.linspace(0.0, 8.0, 801)
    R = np.array([response(t) for t in ts])
    return np.abs(R).max(axis=0)


def t_half() -> np.ndarray:
    """每个 (i,j) 对的 50% 积分响应时间：∫_0^T R_ij dt 达到 G_ij/2 的 T。"""
    ts = np.linspace(0.0, 20.0, 2001)
    R = np.array([response(t) for t in ts])
    out = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            I = np.concatenate([[0.0], np.cumsum(0.5 * (R[1:, i, j] + R[:-1, i, j])
                                                * np.diff(ts))])
            k = int(np.argmax(I >= 0.5 * G[i, j]))
            out[i, j] = ts[k]
    return out


def multiplicative(f0, fstar, kappa=0.1, steps=400):
    """乘性（保正）版本：log p <- log p - kappa B^T (f-f*)，再归一化。"""
    # 从 f0 还原一个 p（在 log p 的零和规范下取最小范数解）
    p = np.exp(np.linalg.lstsq(B, f0, rcond=None)[0])
    p /= p.sum()
    mins = []
    for _ in range(steps):
        lp = np.log(p)
        f = B @ lp
        lp = lp - kappa * (B.T @ (f - fstar))
        lp -= lp.mean()
        p = np.exp(lp)
        p /= p.sum()
        mins.append(float(p.min()))
    return p, min(mins)


def main() -> None:
    seed_p = SEED / SEED.sum()
    f_seed = f_of(seed_p)

    # ── 模式 ──
    modes = [[float(x) for x in VEC[:, k]] for k in range(3)]
    # 统一符号：第一分量 > 0
    for k in range(3):
        if modes[k][0] < 0:
            modes[k] = [-x for x in modes[k]]

    # ── 任意初值 → 三个模式都激发 ──
    rng = np.random.default_rng(3)
    f0 = f_seed + 0.3 * rng.normal(size=3)
    amps = [float(VEC[:, k] @ (f0 - f_seed)) for k in range(3)]
    ts = np.linspace(0.0, 8.0, 161)
    traj = np.array([f_seed + response(t) @ (f0 - f_seed) for t in ts])
    rank = int((np.linalg.svd(traj - f_seed, compute_uv=False) > 1e-9).sum())

    # ── 相互影响 ──
    P = impulse_peaks()
    TH = t_half()
    # 脉冲：slot0 被扰动后，三个槽的响应曲线（用于画图）
    tgrid = np.linspace(0.0, 8.0, 401)
    curves = np.array([[response(t)[i, 0] for t in tgrid] for i in range(3)])

    # ── 乘性保正 ──
    p_end, p_min = multiplicative(f0, f_seed)

    # ── 驱动局部性：顶点扰动只改相邻槽 ──
    drive_support = {int(u): [int(i) for i in range(3) if abs(B[i, u]) > 0]
                     for u in range(4)}

    # ── K4 对照：投影 4I ──
    Lk = np.array([[3.0, -1, -1, -1], [-1, 3, -1, -1],
                   [-1, -1, 3, -1], [-1, -1, -1, 3]])
    K4proj = B @ Lk @ B.T @ np.linalg.inv(M)
    k4_eig = np.linalg.eigvals(K4proj)

    checks = [
        dict(name="BB^T 三个本征值不同（2-√2, 2, 2+√2）",
             pass_=bool(abs(LAM[0] - (2 - math.sqrt(2))) < 1e-12
                        and abs(LAM[1] - 2) < 1e-12
                        and abs(LAM[2] - (2 + math.sqrt(2))) < 1e-12),
             detail=f"{np.round(LAM, 6).tolist()}"),
        dict(name="三个模式都被任意初值激发（振幅非零）",
             pass_=all(abs(a) > 1e-3 for a in amps),
             detail=f"amps={np.round(amps, 4).tolist()}"),
        dict(name="对比轨迹秩 = 3", pass_=rank == 3, detail=f"rank={rank}"),
        dict(name="弛豫到 seed 目标（t=60，误差 <1e-10）",
             pass_=float(np.abs(relax(f0, f_seed, 1.0, 60.0) - f_seed).max()) < 1e-10,
             detail=f"err={float(np.abs(relax(f0, f_seed, 1.0, 60.0) - f_seed).max()):.2e}"),
        dict(name="静态 Green 函数满支撑（G_02 ≠ 0）",
             pass_=bool(abs(G[0, 2]) > 1e-9),
             detail=f"G=\n{np.round(G, 4).tolist()}"),
        dict(name="脉冲横向响应非零（0→1, 0→2）",
             pass_=bool(P[0, 1] > 1e-3 and P[0, 2] > 1e-3),
             detail=f"peak(0→1)={P[0, 1]:.4f}, peak(0→2)={P[0, 2]:.4f}"),
        dict(name="乘性版本保正（min p > 0）", pass_=p_min > 0,
             detail=f"min p={p_min:.3e}"),
        dict(name="驱动是局部的（内部顶点只驱动相邻槽）",
             pass_=all(len(v) <= 2 for v in drive_support.values()),
             detail=f"support={drive_support}"),
        dict(name="K4 全图投影仍是标量（1 模式，故链是必要的）",
             pass_=bool(np.allclose(np.eye(3) * 4.0, K4proj)),
             detail=f"eig={np.round(k4_eig, 6).tolist()}"),
    ]

    res = dict(
        operator=dict(B=B.tolist(), BBt=M.tolist(), eigen=[float(x) for x in LAM],
                      G=G.tolist()),
        f_seed=[float(x) for x in f_seed],
        modes=modes,
        amplitudes=amps,
        rank=rank,
        peak_matrix=P.tolist(),
        t_half=TH.tolist(),
        curves={"t": tgrid.tolist(), "from_slot0": curves.tolist()},
        checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                 "pass": bool(c["pass_"])} for c in checks],
        K4_projection=dict(matrix=K4proj.tolist(), eigen=[[float(z.real), float(z.imag)]
                                                          for z in k4_eig]),
        multiplicative=dict(p_end=[float(x) for x in p_end], min_p=p_min),
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    # (a) 三个模式
    xs = np.arange(3)
    for k in range(3):
        axes[0].plot(xs, modes[k], "o-", lw=1.6, label=f"mode {k+1} (λ={LAM[k]:.3f})")
    axes[0].axhline(0, color="k", lw=0.6)
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels(["slot 0", "slot 1", "slot 2"])
    axes[0].set_title("3 modes of the slot chain")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    # (b) 脉冲响应（slot0 → 三个槽）
    for i in range(3):
        axes[1].plot(tgrid, curves[i], lw=1.5, label=f"Δf{ i }")
    axes[1].axhline(0, color="k", lw=0.6)
    axes[1].set_xlabel("t (1/κ)")
    axes[1].set_title("impulse at slot 0")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    # (c) Green 函数
    im = axes[2].imshow(G, cmap="viridis")
    for i in range(3):
        for j in range(3):
            axes[2].text(j, i, f"{G[i, j]:.2f}", ha="center", va="center",
                         color="w" if G[i, j] < 0.6 else "k")
    axes[2].set_xticks(range(3))
    axes[2].set_yticks(range(3))
    axes[2].set_xticklabels(["slot 0", "slot 1", "slot 2"])
    axes[2].set_yticklabels(["slot 0", "slot 1", "slot 2"])
    axes[2].set_title("static influence G = (BBᵀ)⁻¹")
    fig.colorbar(im, ax=axes[2], fraction=0.046)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_l2_slots.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    def _j(o):
        if isinstance(o, complex):
            return [float(o.real), float(o.imag)]
        return str(o)

    with open(os.path.join(RES, "zero_l2_slots.json"), "w") as f:
        json.dump(res, f, indent=1, default=_j)

    print("L2 三亚层完整动力学（链式交换）")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  模式本征值 {np.round(LAM, 4).tolist()}")
    print(f"  f_seed {np.round(f_seed, 4).tolist()}")
    print(f"  峰响应矩阵\n{np.round(P, 4)}")
    print(f"  t50\n{np.round(TH, 2)}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
