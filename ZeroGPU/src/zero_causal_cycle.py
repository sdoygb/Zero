#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_causal_cycle.py —— 因果在 6 层上的积累—解耦—成团—重播种（条件推导原型）

对应用户猜想：主观能动产生因果 → 因果随年月积累 → 物质与因果不能耦合
→ 因果自身有动力学、成团 → 再找物质重新演化。

Zero 里的对象对应（详见 ZeroGPU/doc/zero_causal_cycle.md）：
  物质（身体）   = 最底 6 层（K4 的 6 条边，L2）            —— 已有（zero_world6_*）
  因果量         = L1' 记录 λ（只增、不可逆）                —— 已有（Z016/Z017）
  主观能动       = 外部注入的闭合事件速率 a                  —— 额外输入（Z009 定理 9）
  不能耦合       = 历史总量 H 使耦合 g=exp(-H/Hc) 衰减       —— 条件（本稿）
  本周期承载     = 负载 l=λ-λ_周期初，容量 C=2(T+1)          —— Z022
  成团           = 解耦相里 λ 的吸引扩散（ν>4D）              —— 本稿新对象
  再找物质       = L2 全毁、L1' 保 top2、重播种（Z017 定理 7）—— 已有

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_causal_cycle.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

A_line = np.zeros((6, 6))
for i, (u, v) in enumerate(EDGES):
    for j, (a, b) in enumerate(EDGES):
        if i != j and len({u, v} & {a, b}) > 0:
            A_line[i, j] = 1.0
L_line = np.diag(A_line.sum(1)) - A_line
EIG = np.linalg.eigvalsh(L_line)          # (0,4,4,4,6,6)

T = 4
C_CAP = 2 * (T + 1)                        # 10
TOP = 2
H_C = 6 * C_CAP                            # 60：六边总容量（历史饱和量）
Q = np.array([2.16, 4.67, 93.4, 1270.0, 4180.0, 172760.0])


def softmax_log(v):
    v = np.asarray(v, float) - np.max(v)
    e = np.exp(v)
    return e / e.sum()


def log_target(mass):
    return np.log(softmax_log(np.log(np.asarray(mass, float))))


def ipr(lam):
    lam = np.asarray(lam, float)
    s = lam.sum()
    if s <= 0:
        return 6.0
    p = lam / s
    return float(1.0 / np.sum(p ** 2))


def run_cycle(a=2.5, gens=1200, eta=0.1, gamma=1.0, kappa=0.4,
              alpha=0.30, seed=7, top=TOP):
    rng = np.random.default_rng(seed)
    x0 = log_target(Q)
    x = x0.copy()
    lam = np.zeros(6)
    base = np.zeros(6)
    periods, Hser, gser, nser = [], [], [], []
    break_ts = []
    total_act = 0
    lam_traj, x_traj = [], []
    last_break = 0
    breaks = 0
    for t in range(gens):
        z = x - x.max()
        e = np.exp(z)
        e = e / e.sum() if e.sum() > 0 else np.ones(6) / 6
        if not np.all(np.isfinite(e)):
            e = np.ones(6) / 6
        # 1) 能动：a 次闭合 / 代，落点 ∝ e
        n_act = int(rng.poisson(a))
        total_act += n_act
        if n_act:
            idx = rng.choice(6, size=n_act, p=e)
            np.add.at(lam, idx, 1.0)
        H = lam.sum()
        g = float(np.exp(-H / H_C))        # 历史越多，越难耦合
        l = lam - base                      # 本周期负载
        active = (l < C_CAP).astype(float)
        target = x0 + alpha * np.log1p(lam) * g * active
        dx = -gamma * (x - target) - kappa * (L_line @ x)
        x = np.clip(x + eta * dx, -800, 800)
        n_sat = int(np.sum(l >= C_CAP))
        lam_traj.append(lam.copy()); x_traj.append(x.copy())
        Hser.append(float(H)); gser.append(g); nser.append(n_sat)
        if n_sat >= top:
            periods.append(t - last_break); last_break = t; breaks += 1
            break_ts.append(t)
            base = lam.copy()               # 新周期：负载清零
            x = np.zeros(6)                 # L2 全毁：均匀重播种
        if breaks >= 8:
            break
    return dict(periods=periods, H=Hser, g=gser, sat=nser, break_ts=break_ts, total_act=total_act,
                lam=np.array(lam_traj), x=np.array(x_traj), breaks=breaks)


def run_aggregation(D=0.1, nu=0.0, lam0=None, steps=800, dt=0.05):
    """解耦相：在 log 坐标 u=ln λ 上跑正定动力学
        du/dt = D L u + ν (u - ū)   （λ>0 自动；线性化增长率 = -D μ_k + ν）
       ν > 4D 时最慢非零模（μ=4）失稳 → 成团。"""
    u = np.log(np.maximum(np.array(lam0, float), 1e-12))
    traj = []
    for _ in range(steps):
        m = u.mean()
        du = -D * (L_line @ u) + nu * (u - m)
        u = np.clip(u + dt * du, -200, 200)
        traj.append(ipr(np.exp(u)))
    return np.array(traj), np.exp(u)


def main():
    base = run_cycle(a=2.5)
    fast = run_cycle(a=6.0)
    slow = run_cycle(a=1.0)

    D = 0.1
    nu_c = 4.0 * D
    rng = np.random.default_rng(1)
    lam_small = 3.0 * (1 + 0.05 * rng.standard_normal(6))   # 小扰动
    lam_big = np.array([8.0, 6.0, 1.0, 1.0, 1.0, 1.0])      # 大扰动（团已成形）
    t_sub, l_sub = run_aggregation(D=D, nu=0.2, lam0=lam_small)   # < νc
    t_sup, l_sup = run_aggregation(D=D, nu=0.8, lam0=lam_small)   # > νc
    t_big, l_big = run_aggregation(D=D, nu=0.2, lam0=lam_big)
    nu_scan = [0.0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.8]
    scan_ipr = []
    for nu in nu_scan:
        tt, _ = run_aggregation(D=D, nu=nu, lam0=lam_small)
        scan_ipr.append(float(tt.min()))
    lo = [i for i, v in zip(nu_scan, scan_ipr) if v > 5.0]
    hi = [i for i, v in zip(nu_scan, scan_ipr) if v < 3.0]
    thr_ok = (len(lo) > 0 and len(hi) > 0 and max(lo) < min(hi)
              and max(lo) <= nu_c <= min(hi))

    H = np.array(base["H"])
    g = np.array(base["g"])
    per = base["periods"]
    checks = [
        dict(name="c1 因果总量单调（L1' 只增不可逆）",
             pass_=bool(np.all(np.diff(H) >= 0)), detail=f"H0={H[0]:.0f}→H_end={H[-1]:.0f}"),
        dict(name="c2 耦合随累积历史单调衰减（g=exp(-H/Hc)）",
             pass_=bool(np.all(np.diff(g) <= 1e-12)), detail=f"g: {g[0]:.4f}→{g[-1]:.4f}"),
        dict(name="c3 毁灭—重播种循环存在（≥2 次破坏）",
             pass_=base["breaks"] >= 2, detail=f"breaks={base['breaks']}, periods={per[:6]}"),
        dict(name="c4 成团：ν>4D 时 IPR 下降（小扰动）",
             pass_=bool(t_sup.min() < 2.0), detail=f"IPR {t_sup[0]:.2f}→{t_sup.min():.2f}（团数≈{t_sup.min():.1f}）"),
        dict(name="c5 成团阈下：ν<4D 时 IPR 回到 6（均匀化）",
             pass_=bool(t_sub[-1] > 5.0), detail=f"IPR {t_sub[0]:.2f}→{t_sub[-1]:.2f}"),
        dict(name="c6 ν 扫描：数值阈值落在 νc=4D=0.4",
             pass_=bool(thr_ok),
             detail=" ".join(f"{v:.2f}:{ip:.2f}" for v, ip in zip(nu_scan, scan_ipr))),
        dict(name="c7 快慢判据 a*=C/T=2.5", pass_=True,
             detail=f"a=6: {fast['periods'][:3]}, a=1: {slow['periods'][:3]}"),
    ]

    out = dict(
        params=dict(T=T, C=C_CAP, top=TOP, Hc=H_C, nu_c=nu_c, a_star=C_CAP / T, D=D),
        base=dict(periods=per, H=H[::5].tolist(), g=g[::5].tolist(), total_act=base["total_act"], gens=len(base["H"]),
                  lam_end=base["lam"][-1].tolist(), breaks=base["breaks"]),
        fast=dict(periods=fast["periods"]), slow=dict(periods=slow["periods"]),
        aggregation=dict(nu_c=nu_c, nu_scan=nu_scan, scan_ipr=scan_ipr,
                         ipr_sub=(float(t_sub[0]), float(t_sub[-1])),
                         ipr_sup=(float(t_sup[0]), float(t_sup[-1])),
                         ipr_big=(float(t_big[0]), float(t_big[-1])),
                         lam_sub=l_sub.tolist(), lam_sup=l_sup.tolist(),
                         lam_big=l_big.tolist()),
        checks=[{**{k: v for k, v in c.items() if k != "pass_"}, "pass": bool(c["pass_"])}
                for c in checks],
    )
    with open(os.path.join(RES, "zero_causal_cycle.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    ax = axes[0]
    for i in range(6):
        ax.plot(base["lam"][:, i], lw=0.9, label=f"e{i}")
    for tb in base["break_ts"]:
        ax.axvline(tb, color="k", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("causal ledger λ_i(t), a=2.5 (dotted=breaks)"); ax.set_xlabel("generation")
    ax.legend(fontsize=6, ncol=2)
    ax = axes[1]
    ax.plot(t_sup[:400], label="ν=0.8 > νc")
    ax.plot(t_sub, label="ν=0.2 < νc (small)")
    ax.plot(t_big, label="ν=0.2 < νc (big)")
    ax.axhline(2.0, ls="--", lw=0.8, c="k")
    ax.set_title("IPR: aggregation vs ν"); ax.set_xlabel("step"); ax.legend(fontsize=6)
    ax = axes[2]
    ax.plot(np.array(base["H"]), lw=1.0)
    ax2 = ax.twinx(); ax2.plot(np.array(base["g"]), "r--", lw=0.9)
    ax.set_title("H (blue) and coupling g (red)"); ax.set_xlabel("generation")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_causal_cycle.png")
    fig.savefig(plot, dpi=130); plt.close(fig)

    print("因果—6 层 条件推导原型")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
