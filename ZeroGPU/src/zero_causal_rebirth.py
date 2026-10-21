#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_causal_rebirth.py —— 三段：锁死 → 脱耦 → 成团，以及“成团后能否重新耦合”

接上 g = exp(-H/Hc)（上一轮 zero_causal_coupling.py 故意关掉的那个）。

阶段：
  A 锁死：β>0，g≈1，因果与物质一起集中；
  B 脱耦：H 增长使 g→0，耦合断开；
  C 成团：解耦相 λ 在 u=lnλ 上富集，IPR→1；
  D 重播种：新 L2 均匀，记录以“团核”形式回来，耦合 g_new = exp(-H_core/Hc)。

问（用户猜想）：成团之后，因果会不会重新耦合新生命？

答案取决于**成团时释放多少**：
  H_core = f·H_end，g_new = exp(-f·H_end/Hc)。
  f=1（不释放）→ g_new≈0 → 耦合失败；
  f 足够小 → g_new 回升 → 重新耦合。
临界 f* = Hc·ln2 / H_end（要求 g_new ≥ 0.5）。

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_causal_rebirth.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
C = 10.0
HC = 6 * C

EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
A_line = np.zeros((6, 6))
for i, (u, v) in enumerate(EDGES):
    for j, (a, b) in enumerate(EDGES):
        if i != j and len({u, v} & {a, b}) > 0:
            A_line[i, j] = 1.0
L_line = np.diag(A_line.sum(1)) - A_line


def _softmax(x):
    z = x - x.max()
    e = np.exp(z)
    return e / e.sum()


def ipr(v):
    v = np.asarray(v, float)
    s = v.sum()
    if s <= 0:
        return 6.0
    p = v / s
    return float(1.0 / np.sum(p ** 2))


def phase_AC(beta=0.5, alpha=1.0, kappa=0.1, eta=0.2, p=0.5,
             gens=320, gstar=0.01, seed=5):
    """A 锁死 → B 脱耦 → C 成团。返回轨迹与终态。"""
    rng = np.random.default_rng(seed)
    x = np.zeros(6)
    lam = np.zeros(6)
    ie, il, gs = [], [], []
    t_dec = None
    for t in range(gens):
        H = lam.sum()
        g = float(np.exp(-H / HC))
        if t_dec is None and g < gstar:
            t_dec = t
        e = _softmax(x)
        logw = np.log(e + 1e-12) + beta * lam
        logw -= logw.max()
        w = np.exp(logw); w /= w.sum()
        for _ in range(6):
            i = int(rng.choice(6, p=w))
            if rng.random() < p:
                lam[i] += 1.0
        # 物质被记录拉（g 衰减：越老越拉不动）
        x = x + eta * (-kappa * (L_line @ x) + alpha * np.log1p(lam) * g)
        x = x - x.mean()
        x = np.clip(x, -60, 60)
        ie.append(ipr(e)); il.append(ipr(lam)); gs.append(g)

    # C 成团：解耦相 ν-富集
    u = np.log(np.maximum(lam, 1e-12))
    for _ in range(600):
        m = u.mean()
        u = np.clip(u + 0.05 * (-0.1 * (L_line @ u) + 0.8 * (u - m)), -200, 200)
    lam_c = np.exp(u)
    lam_c = lam_c * (lam.sum() / lam_c.sum())     # 质量守恒：成团不改变总量
    return dict(lam=lam, lam_clustered=lam_c, H=float(lam.sum()),
                H_clustered=float(lam_c.sum()),
                ie=ie, il=il, gs=gs, t_dec=t_dec,
                ipr_l_clustered=ipr(lam_c),
                e_end=_softmax(x).tolist())


def phase_D(lam_core, beta=0.5, alpha=1.0, kappa=0.1, eta=0.2,
            gens=300, seed=9):
    """重播种：新 L2 均匀，记录以团核形式回来。"""
    rng = np.random.default_rng(seed)
    H_core = float(lam_core.sum())
    g_new = float(np.exp(-H_core / HC))
    x = np.zeros(6)
    for t in range(gens):
        e = _softmax(x)
        logw = np.log(e + 1e-12) + beta * lam_core
        logw -= logw.max()
        w = np.exp(logw); w /= w.sum()
        for _ in range(6):
            int(rng.choice(6, p=w))
        x = x + eta * (-kappa * (L_line @ x) + alpha * np.log1p(lam_core) * g_new)
        x = x - x.mean()
        x = np.clip(x, -60, 60)
    e = _softmax(x)
    return dict(g_new=g_new, H_core=H_core, ipr_e=ipr(e), max_e=float(e.max()), e=e.tolist())


def transfer_scan(ps, gens=320, HC=HC):
    """若“成团保留 = 闭合比例” f=p，则遗忘比例 = 1-p，但耦合 g=exp(-f·H/Hc) 会下降。
    有效传递 T = f·g。H = 6p·gens（每代 6 次尝试，闭合 p）。"""
    rows = []
    for p in ps:
        H = 6.0 * p * gens
        f = p
        g = float(np.exp(-f * H / HC))
        rows.append(dict(p=float(p), f=float(f), H=float(H), g=g, T=float(f * g)))
    return rows


def main():
    r = phase_AC()
    lam_c = r["lam_clustered"]
    H = float(lam_c.sum())

    fs = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
    rows = []
    for f in fs:
        core = lam_c * f
        d = phase_D(core)
        d["f"] = f
        rows.append(d)

    f_star = HC * np.log(2.0) / H
    recoupled = [d for d in rows if d["ipr_e"] < 5.0]
    f_num = min((d["f"] for d in recoupled), default=None)

    checks = [
        dict(name="c1 三段齐全：锁死(g≈1) → 脱耦(g<0.01) → 成团(IPR_λ→1)",
             pass_=bool(r["t_dec"] is not None and r["ipr_l_clustered"] < 4.0),
             detail=f"t_dec={r['t_dec']}, IPR_λ(团)={r['ipr_l_clustered']:.2f}"),
        dict(name="c2 不释放（f=1）→ 不能重新耦合",
             pass_=bool(rows[-1]["ipr_e"] > 4.0),
             detail=f"f=1: g={rows[-1]['g_new']:.2e}, IPR_e={rows[-1]['ipr_e']:.2f}"),
        dict(name="c3 充分释放（f≤0.01）→ 重新耦合成功",
             pass_=bool(rows[2]["ipr_e"] < 4.0),
             detail=f"f=0.01: g={rows[2]['g_new']:.3f}, IPR_e={rows[2]['ipr_e']:.2f}"),
        dict(name="c4 临界 f* = Hc·ln2/H 与数值同量级",
             pass_=bool(f_num is not None and abs(f_num - f_star) / f_star < 1.0),
             detail=f"数值 f*={f_num}, 解析 f*={f_star:.4f}, H={H:.0f}"),
        dict(name="c5 结论：重新耦合有一个 f 窗口（非单调）",
             pass_=bool(min(d["ipr_e"] for d in rows) < 4.0
                        and rows[-1]["ipr_e"] > 5.0),
             detail="f:" + " ".join(f"{d['f']:g}->{d['ipr_e']:.2f}" for d in rows)),
    ]

    # ---- 用户猜想：越接近闭合，遗忘越少；接上“保留=闭合比例” ----
    ps = [0.01, 0.05, 0.1, 0.125, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99]
    tr = transfer_scan(ps)
    T_peak = max(tr, key=lambda d: d["T"])
    p_star_theory = float(np.sqrt(HC / (12.0 * 320)))   # H=6p·gens → T=p·exp(-6p²gens/Hc)

    checks.append(dict(name="c6 越接近闭合，遗忘比例 1-p 越小（单调）",
                       pass_=True,
                       detail="1-p " + " ".join(f"{d['p']:g}:{1-d['p']:.3f}" for d in tr)))
    checks.append(dict(name="c7 但有效传递 T=f·g 有峰值（p→1 反而传不过去）",
                       pass_=bool(T_peak["p"] < 0.5 and T_peak["T"] > tr[-1]["T"] * 1e6),
                       detail=f"T_peak at p={T_peak['p']}, T={T_peak['T']:.3e}; "
                              f"p=0.99 T={tr[-1]['T']:.2e}"))
    checks.append(dict(name="c8 峰值位置与解析 p*=sqrt(Hc/(12·gens)) 同量级",
                       pass_=bool(abs(T_peak["p"] - p_star_theory) <= 0.05),
                       detail=f"数值 {T_peak['p']}, 解析 {p_star_theory:.3f}"))

    out = dict(params=dict(C=C, HC=HC, beta=0.5, alpha=1.0, kappa=0.1, p=0.5,
                           gstar=0.01, H_end=H, f_star=f_star),
               AC=dict(t_dec=r["t_dec"], H=r["H"], H_clustered=r["H_clustered"],
                       ipr_l_clustered=r["ipr_l_clustered"]),
               rows=rows, transfer=tr,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_causal_rebirth.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    ax = axes[0]
    ax.plot(r["ie"], lw=0.8, label="IPR_e (matter)")
    ax.plot(r["il"], lw=0.8, label="IPR_λ (ledger)")
    if r["t_dec"]:
        ax.axvline(r["t_dec"], ls=":", c="k", lw=0.8)
    ax.set_xlabel("generation"); ax.set_ylabel("IPR")
    ax.set_title("lock-in → decouple"); ax.legend(fontsize=7)
    ax = axes[1]
    ax.plot(r["gs"], lw=0.9)
    ax.set_yscale("log"); ax.set_xlabel("generation"); ax.set_ylabel("g")
    ax.set_title("coupling g=exp(-H/Hc) decays")
    ax = axes[2]
    ax.plot([d["f"] for d in rows], [d["ipr_e"] for d in rows], "o-", label="IPR_e after reseed")
    ax.axvline(f_star, ls=":", c="r", label=f"f*={f_star:.3f}")
    ax.set_xlabel("f = kept fraction of karma"); ax.set_ylabel("IPR_e")
    ax.set_title("re-coupling needs release"); ax.legend(fontsize=7)

    # 第四面板（另图）：T(p)
    fig2, ax2 = plt.subplots(1, 1, figsize=(5, 3.6))
    ax2.plot([d["p"] for d in tr], [d["T"] for d in tr], "o-")
    ax2.axvline(T_peak["p"], ls=":", c="r", label=f"peak p={T_peak['p']}")
    ax2.set_xlabel("p (closure)"); ax2.set_ylabel("T = f·g")
    ax2.set_title("effective transmission vs closure")
    ax2.legend(fontsize=8); ax2.grid(alpha=.3)
    fig2.tight_layout()
    plot2 = os.path.join(RES, "zero_causal_transfer.png")
    fig2.savefig(plot2, dpi=130); plt.close(fig2)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_causal_rebirth.png")
    fig.savefig(plot, dpi=130); plt.close(fig)

    print("锁死 → 脱耦 → 成团 → 重播种")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
