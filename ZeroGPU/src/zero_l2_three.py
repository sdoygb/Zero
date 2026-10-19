#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_l2_three.py —— 去掉亚亚层后，L2 只剩 3 个亚层：动力学完整性审计
=====================================================================
问题：把 L2 的亚亚层（16/8/4）全部去掉，只留 **3 个亚层**（4 顶点占据
$p\\in\\Delta^3$ 的 3 个独立自由度），这套 3 维动力学是否**完整**？

"完整"拆成四项：
  A. **闭合**：3 维状态自足（$p_{t+1}=F(p_t)$，不需要亚亚层/纤维隐变量）；
  B. **模式数**：形状线性化在固定点处有 3 个不同的非平凡本征值；
  C. **非平凡**：无外驱动时存在非平凡内动力学（方向不恒定），而非只有平凡吸引子；
  D. **可完成**：若 B/C 不满足，最少加什么（槽层算符/手征步）才能补成 3 模。

结论（先说）：
  1. **闭合是完整的**：正则舞台上闭合序列严格满足 $c_{t+1}=A^L c_t$
     （引擎实测 $2.8\\times10^{-16}$），3 维状态无隐变量——**去掉亚亚层合法**。
  2. **模式数不完整**：K4 是顶点传递图，$A^L=\\alpha I+\\beta J$，形状线性化
     在 $1^\\perp$ 上 $=\\alpha/\\langle d\\rangle I$——三个方向**本征值全同**，
     对比轨迹秩 1。"3 亚层"仍是运动学计数。
  3. **非平凡性不完整**：无宏观回写时只有平凡吸引子（fresh→均匀，
     persistent→seed 冻结）；对比方向不变、无内振荡。
  4. **可完成**：给槽层/手征结构一个新的具名输入即可补成 3 模——
     C3 手征步 $(AP)^L$ 的形状本征值 $1,\\omega,\\omega^2$（秩 3）；
     槽拉普拉斯 $BB^{\\mathsf T}$ 的本征值 $2-\\sqrt2,2,2+\\sqrt2$（秩 3）。
     两者都是**输入**，不是从现有 L0/L1/引擎导出的。

输出：results/zero_l2_three.json、zero_l2_three.png

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_l2_three.py
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
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
L = 4
V = 4


# ───────────────────────── 基本工具 ─────────────────────────
def contrast_basis(n: int = V) -> np.ndarray:
    """1^⊥ 的正交基（Helmert 型，n x (n-1)）。"""
    H = np.zeros((n, n - 1))
    for j in range(n - 1):
        H[: j + 1, j] = 1.0
        H[j + 1, j] = -(j + 1)
    Q, _ = np.linalg.qr(H)
    return Q


Q = contrast_basis(V)


def adj(neighbors) -> np.ndarray:
    return zs.adj_matrix(neighbors).astype(float)


def perron(M: np.ndarray) -> np.ndarray:
    """M 的 Perron 向量（取最大实部本征向量），归一化 sum=1。"""
    ev, vec = np.linalg.eig(M)
    k = int(np.argmax(ev.real))
    v = np.real(vec[:, k])
    v = np.abs(v)
    return v / v.sum()


def shape_map(AL: np.ndarray, p: np.ndarray) -> np.ndarray:
    q = AL @ p
    return q / q.sum()


def shape_jacobian(AL: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """形状映射在 Perron 固定点的 Jacobian（4x4，含固定点方向）。

    F(p)=AL p / (1^T AL p)，J=(AL - p* 1^T AL)/(1^T AL p*)。
    J 有一个 0 本征值（固定点方向），其余是形状模。
    """
    pstar = perron(AL)
    den = float((AL @ pstar).sum())
    J = (AL - np.outer(pstar, AL.sum(axis=0))) / den
    return J, pstar


def shape_modes(J: np.ndarray, tol: float = 1e-9) -> list[complex]:
    ev = np.linalg.eigvals(J)
    return [complex(z) for z in ev if abs(z) > tol]


def n_distinct(ev, tol: float = 1e-6) -> int:
    cl: list[complex] = []
    for z in ev:
        if not any(abs(z - c) < tol for c in cl):
            cl.append(z)
    return len(cl)


def rank_traj(X: np.ndarray, tol: float = 1e-7) -> int:
    """列向量轨迹的秩（相对容差）。"""
    s = np.linalg.svd(X, compute_uv=False)
    if s[0] <= 0:
        return 0
    return int((s > tol * s[0]).sum())


def slots(p: np.ndarray) -> np.ndarray:
    return np.array([abs(math.log(p[i] / p[i + 1])) for i in range(3)])


def cos_step(X: np.ndarray, n: int = 16) -> float:
    """相邻两步对比方向夹角的余弦（只取前 n 步，避免衰减到 0 的数值噪声）。"""
    a, b = X[:, :n], X[:, 1:n + 1]
    num = (a * b).sum(axis=0)
    den = np.linalg.norm(a, axis=0) * np.linalg.norm(b, axis=0)
    ok = den > 1e-12 * max(1.0, float(np.abs(X).max()))
    return float(np.median(num[ok] / den[ok]))


def stage_case(name: str, AL: np.ndarray, gens: int = 400) -> dict:
    """给定一代映射 AL，返回模式谱与对比轨迹秩。"""
    J, pstar = shape_jacobian(AL)
    Jc = Q.T @ J @ Q
    modes = shape_modes(Jc)
    rng = np.random.default_rng(11)
    p0 = rng.dirichlet(np.ones(V))
    x = Q.T @ (p0 - pstar)
    xs = [x]
    for _ in range(gens):
        x = Jc @ x
        xs.append(x)
    X = np.array(xs).T                # 线性化对比轨迹（精确到 K4 的线性情形）
    Xn = X / np.maximum(np.linalg.norm(X, axis=0), 1e-300)
    ps = [p0]
    for _ in range(gens):
        ps.append(shape_map(AL, ps[-1]))
    ps = np.array(ps)
    F = np.array([slots(p) for p in ps])
    return dict(name=name,
                contrast_norm=Xn[:, :40].tolist(),
                pstar=[float(x) for x in pstar],
                modes=[[float(z.real), float(z.imag)] for z in modes],
                n_modes=n_distinct(modes),
                rank_contrast=rank_traj(X),
                cos_step=cos_step(X),
                p_final=[float(x) for x in ps[-1]],
                slots_trace=F.tolist())


# ───────────────────────── 引擎实测 ─────────────────────────
def engine_checks() -> dict:
    K4 = zs.k4()
    A = adj(K4)
    AL = np.linalg.matrix_power(A, L)
    out: dict = {}

    # 1) fresh（retention=last, T=0）：c_{t+1} = A^L c_t
    p = Params(micro=True, micro_semantics="fresh", L=L, T=0,
               retention="last", wave_on=False, poisson_on=False)
    e = ZeroEvo(K4, p, seed=np.array([8.0, 4.0, 1.5, 1.0]))
    cs = []
    for _ in range(20):
        e.step()
        cs.append(e.gen_closures.copy())
    cs = np.array(cs)
    dev = 0.0
    for t in range(len(cs) - 1):
        pred = AL @ cs[t]
        pred *= cs[t + 1].sum() / pred.sum()
        dev = max(dev, float(np.abs(pred - cs[t + 1]).max()
                             / max(1e-300, np.abs(cs[t + 1]).max())))
    out["closure_dev"] = dev

    # 2) fresh record 轨迹秩（T=0 与 T=4）
    J, pstar = shape_jacobian(AL)
    for T in (0, 4):
        p = Params(micro=True, micro_semantics="fresh", L=L, T=T,
                   retention="top2", wave_on=False, poisson_on=False)
        e = ZeroEvo(K4, p, seed=np.array([8.0, 4.0, 1.5, 1.0]))
        shapes = []
        for _ in range(400 if T else 200):
            e.step()
            shapes.append(e.record / e.record.sum())
        X = Q.T @ (np.array(shapes) - pstar).T
        out[f"fresh_T{T}_rank"] = rank_traj(X)
        out[f"fresh_T{T}_final"] = [float(x) for x in shapes[-1]]

    # 3) persistent：record=seed 冻结
    p = Params(micro=True, micro_semantics="persistent", L=L, T=0,
               retention="top2", wave_on=False, poisson_on=False)
    e = ZeroEvo(K4, p, seed=np.array([8.0, 4.0, 1.5, 1.0]))
    shapes = []
    for _ in range(100):
        e.step()
        shapes.append(e.record / e.record.sum())
    X = Q.T @ (np.array(shapes) - pstar).T
    out["persistent_rank"] = rank_traj(X)
    out["persistent_final"] = [float(x) for x in shapes[-1]]
    out["persistent_seed"] = [float(x) for x in (np.array([8.0, 4.0, 1.5, 1.0]) / 14.5)]
    return out


# ───────────────────────── 槽层算符与两种补全 ─────────────────────────
def slot_operator() -> dict:
    B = np.zeros((3, 4))
    for i in range(3):
        B[i, i] = 1.0
        B[i, i + 1] = -1.0
    Ls = B @ B.T
    ev = np.linalg.eigvalsh(Ls)
    kappa = 0.1
    rng = np.random.default_rng(5)
    f = rng.dirichlet(np.ones(3))
    traj = [f.copy()]
    for _ in range(200):
        f = f - kappa * (Ls @ f)
        traj.append(f.copy())
    traj = np.array(traj)
    return dict(B=B.tolist(), L_slot=Ls.tolist(),
                eigen=[float(x) for x in ev],
                rate=[float(1.0 - kappa * x) for x in ev],
                rank_traj=rank_traj(traj.T),
                min_entry=float(traj.min()),
                final=[float(x) for x in traj[-1]])


def chiral_case(gens: int = 400) -> dict:
    """C3 手征步：p <- P A^L p（先走 L 步，再 1->2->3->1 旋转）。"""
    A = adj(zs.k4())
    AL = np.linalg.matrix_power(A, L)
    P = np.zeros((V, V))
    P[0, 0] = 1.0
    P[2, 1] = P[3, 2] = P[1, 3] = 1.0
    ALc = P @ AL
    J, pstar = shape_jacobian(ALc)
    Jc = Q.T @ J @ Q
    modes = shape_modes(Jc)
    rng = np.random.default_rng(11)
    p0 = rng.dirichlet(np.ones(V))
    x = Q.T @ (p0 - pstar)
    xs = [x]
    for _ in range(gens):
        x = Jc @ x
        xs.append(x)
    X = np.array(xs).T
    Xn = X / np.maximum(np.linalg.norm(X, axis=0), 1e-300)
    ps = [p0]
    for _ in range(gens):
        ps.append(shape_map(ALc, ps[-1]))
    ps = np.array(ps)
    F = np.array([slots(p) for p in ps])
    return dict(pstar=[float(x) for x in pstar],
                contrast_norm=Xn[:, :40].tolist(),
                modes=[[float(z.real), float(z.imag)] for z in modes],
                n_modes=n_distinct(modes),
                rank_contrast=rank_traj(X),
                cos_step=cos_step(X),
                min_entry=float(ps.min()),
                slots_trace=F.tolist(),
                p_final=[float(x) for x in ps[-1]])


# ───────────────────────── 图与主程序 ─────────────────────────
def make_plot(res: dict) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    k4 = np.array(res["stages"]["K4"]["contrast_norm"])
    ch = np.array(res["chiral"]["contrast_norm"])
    for ax, M, title in ((axes[0], k4, "K4 undirected: contrast direction fixed"),
                         (axes[1], ch, "K4 + C3 chiral: contrast rotates")):
        for i in range(3):
            ax.plot(M[i], lw=1.6, label=f"x{i+1}")
        ax.set_xlabel("generation")
        ax.set_ylabel("normalized contrast component")
        ax.set_title(title)
        ax.legend(fontsize=8, ncol=3)
        ax.set_ylim(-1.05, 1.05)
        ax.grid(alpha=0.3)
    ax = axes[2]
    op = res["slot_operator"]
    ax.bar(range(3), op["eigen"], color=["#888", "#3a7", "#37a"])
    ax.set_xticks(range(3))
    ax.set_xticklabels(["2-√2", "2", "2+√2"])
    ax.set_title("slot Laplacian BBᵀ eigenvalues")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    path = os.path.join(RES, "zero_l2_three.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main() -> None:
    A = adj(zs.k4())
    AL = np.linalg.matrix_power(A, L)
    c4 = np.linalg.matrix_power(adj(zs.cycle(4)), L)
    p4 = np.linalg.matrix_power(adj(zs.path(4)), L)

    stages = {
        "K4": stage_case("K4", AL),
        "C4": stage_case("C4", c4),
        "P4": stage_case("P4", p4),
    }
    eng = engine_checks()
    op = slot_operator()
    ch = chiral_case()

    checks = [
        dict(name="闭合律 c_{t+1}=A^L c_t（引擎实测）", pass_=eng["closure_dev"] < 1e-12,
             detail=f"max rel dev = {eng['closure_dev']:.3e}"),
        dict(name="K4 形状模数 = 1（三方向本征值全同）", pass_=stages["K4"]["n_modes"] == 1,
             detail=f"modes={np.round(stages['K4']['modes'], 6).tolist()}"),
        dict(name="K4 对比轨迹秩 = 1（方向不变）", pass_=(stages["K4"]["rank_contrast"] == 1
                                                          and stages["K4"]["cos_step"] > 0.999),
             detail=f"rank={stages['K4']['rank_contrast']}, cos={stages['K4']['cos_step']:.6f}"),
        dict(name="对照 C4：也退化（对称舞台通病）", pass_=stages["C4"]["n_modes"] <= 2,
             detail=f"modes={np.round(stages['C4']['modes'], 4).tolist()}"),
        dict(name="对照 P4：n_modes=2（非 3）", pass_=stages["P4"]["n_modes"] == 2,
             detail=f"modes={np.round(stages['P4']['modes'], 4).tolist()}"),
        dict(name="persistent 记录冻结到 seed（秩 1）", pass_=eng["persistent_rank"] == 1,
             detail=f"rank={eng['persistent_rank']}"),
        dict(name="槽层算符 BB^T 有 3 个不同本征值", pass_=op["eigen"][0] < op["eigen"][1] < op["eigen"][2],
             detail=f"{np.round(op['eigen'], 4).tolist()}"),
        dict(name="槽扩散轨迹秩 = 3 且非负", pass_=(op["rank_traj"] == 3 and op["min_entry"] > 0),
             detail=f"rank={op['rank_traj']}, min={op['min_entry']:.3e}"),
        dict(name="C3 手征：3 个不同形状模（1,ω,ω²）", pass_=ch["n_modes"] == 3,
             detail=f"{np.round(ch['modes'], 6).tolist()}"),
        dict(name="C3 手征：对比轨迹秩 = 3、方向在转且非负",
             pass_=(ch["rank_contrast"] == 3 and ch["cos_step"] < 0.9 and ch["min_entry"] > 0),
             detail=f"rank={ch['rank_contrast']}, cos={ch['cos_step']:.4f}, min={ch['min_entry']:.2e}"),
    ]

    res = dict(L=L, checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                             "pass": bool(c["pass_"])} for c in checks],
               engine=eng, stages=stages, slot_operator=op, chiral=ch)

    plot = make_plot(res)

    def _j(o):
        if isinstance(o, complex):
            return [float(o.real), float(o.imag)]
        return str(o)

    with open(os.path.join(RES, "zero_l2_three.json"), "w") as f:
        json.dump(res, f, indent=1, default=_j)

    print(f"L2 3 亚层完整性审计  L={L}")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  引擎闭合偏差 {eng['closure_dev']:.2e}")
    for k in ("K4", "C4", "P4"):
        s = stages[k]
        print(f"  {k}: 模 {np.round(s['modes'], 4).tolist()} 秩 {s['rank_contrast']}")
    print(f"  C3: 模 {np.round(ch['modes'], 4).tolist()} 秩 {ch['rank_contrast']}")
    print(f"  槽算符 {np.round(op['eigen'], 4).tolist()} 秩 {op['rank_traj']}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)}  → {plot}")


if __name__ == "__main__":
    main()
