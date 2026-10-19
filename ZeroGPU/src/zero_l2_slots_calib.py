#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_l2_slots_calib.py —— 链式模型对引擎的标定：κ、三种激发模式
=====================================================================
问题（用户 (a)）：把链路模型的 κ 与引擎的衰减标定，看它能否复现 Z161 的
几何/相位/容量三种激发模式。

做法：
  1. **κ 标定**：引擎 fresh 无宏观时，形状对比按 $1/81$ 每代衰减（实测算得
     $-\\ln81=-4.394$）。链模型的慢模速率是 $\\kappa(2-\\sqrt2)$，两者相等给
     $\\kappa=7.50$/代。核对链模型的衰减曲线与引擎的逐代槽值。
  2. **三种模式**：按 Z161/zero_layers 的原参数与原初值（rng(7) 逐例推进）
     复跑 geom/phase/capacity，取记录槽 $f_{\\rm obs}$ 与引擎的种子调制
     $d=\\eta_n x+\\eta_h h+\\eta_\\psi|\\psi|^2$；比较
       · 链模型的稳态 $f_{\\rm seed}+G\\,Bd$（源项口径）与 $f_{\\rm seed}+Bd$（目标口径）；
       · **局部**响应 $\\alpha\\,Bd$（相邻差，不经 $G$）。
  3. **判定**：链模型能否复现三种模式。

关键数字（先说）：
  · κ=7.50/代：链的慢模速率 $4.394$ 与引擎的 $1/81$ **精确对上**；
  · 但**瞬态形状不对**：第 1 代引擎槽 $(0.0135,0.0085,0.0017)$，链给
    $(0.0077,0.0109,0.0077)$（相对误差 43%/28%/350%）——引擎的 $f_1$ 是
    "种子被 $A^L$ 均匀化"的形状，链给的是慢模 $(1,\\sqrt2,1)$ 的形状；
  · **三种模式**：geom 与 phase 的 $f_{\\rm obs}$ 几乎正比于 $Bd$
    （$\\cos=0.9935,1.0000$，$\\alpha=0.353,0.271$）；capacity 不（$\\cos=-0.147$，
    它由 503 次硬容量事件驱动，不是种子调制）；
  · 链的 $G$ 口径 $\cos=0.740,0.826,0.084$，稳态预测误差 $0.74$–$1.23$——
    **比局部 $Bd$ 差**。

判定：**不能**。κ 标定只在"速率"上成功；链模型的**全局交换**（$G$）没有被引擎
采用——引擎的 L2 记录是**局部**的（$f\\approx\\alpha Bd$），三个亚层的相互影响
在引擎里只通过归一化/幸存者反馈出现，不是链上的直接交换。

输出：results/zero_l2_slots_calib.json、zero_l2_slots_calib.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_l2_slots_calib.py
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
SEED = (8.0, 4.0, 1.5, 1.0)
L = 4

B = np.zeros((3, 4))
for _i in range(3):
    B[_i, _i] = 1.0
    B[_i, _i + 1] = -1.0
M = B @ B.T
LAM, VEC = np.linalg.eigh(M)
G = np.linalg.inv(M)


def f_of(p):
    p = np.asarray(p, float)
    return B @ np.log(p)


def contrast(p):
    p = np.asarray(p, float)
    return p - p.mean()


def chain_decay(f0, kappa, gens):
    return np.array([VEC @ np.diag(np.exp(-LAM * kappa * t)) @ VEC.T @ f0
                     for t in range(gens)])


def cos(a, b):
    return float(np.dot(a, b) / max(1e-30, np.linalg.norm(a) * np.linalg.norm(b)))


def engine_decay(gens=8):
    """fresh 无宏观：逐代闭合形状的槽与对比。"""
    p = Params(micro=True, micro_semantics="fresh", L=L, T=0,
               retention="last", wave_on=False, poisson_on=False)
    e = ZeroEvo(zs.k4(), p, seed=np.array(SEED))
    fs, xs, ps = [], [], []
    for _ in range(gens):
        e.step()
        rec = e.gen_closures / e.gen_closures.sum()
        ps.append(rec)
        fs.append(f_of(rec))
        xs.append(contrast(rec))
    return np.array(fs), np.array(xs), np.array(ps)


def macro_cases(gens=10000):
    """按 zero_layers 的原参数/原 rng 序列复跑三种宏观。"""
    base = dict(micro=True, micro_semantics="fresh", L=L, T=4,
                retention="top2", eta_n=0.2, eta_h=0.2, eta_psi=0.2,
                mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
                capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
                Gamma_h=0.2, kappa_h=0.5, alpha=0.0)
    cases = {
        "no_macro": dict(eta_n=0.0, eta_h=0.0, eta_psi=0.0, g_src=0.0, kappa_h=0.0),
        "geom": dict(),
        "phase": dict(beta=0.5, g_wave=0.5, kappa_psi=0.5, Gamma_w=0.5, wave_on=True),
        "capacity": dict(capacity="hard", K_hard=1.0, reset_amp=0.05),
    }
    rng = np.random.default_rng(7)
    out = {}
    for name, kw in cases.items():
        p = Params(**{**base, **kw})
        x0 = rng.normal(size=4)
        x0 -= x0.mean()
        x0 *= 0.05 / max(1e-12, np.abs(x0).max())
        e = ZeroEvo(zs.k4(), p, seed=np.array(SEED), n0=1.0 + x0)
        e.run(gens, every=500)
        rec = e.record / e.record.sum()
        n = e.n / max(1e-30, e.n.mean()) - 1.0
        h = e.h - e.h.mean()
        a = np.abs(e.psi) ** 2
        a = a / max(1e-30, a.mean()) - 1.0
        d = p.eta_n * n + p.eta_h * h + p.eta_psi * a
        out[name] = dict(record=[float(x) for x in rec], f=f_of(rec),
                         d=d, Bd=B @ d, GBd=G @ (B @ d),
                         events=len(e.events))
    return out


def main() -> None:
    fts, xs, ps = engine_decay()
    # 对比衰减速率
    tt = np.arange(1, 6)
    slope = float(np.polyfit(tt, np.log(np.linalg.norm(xs, axis=1)[1:6]), 1)[0])
    rho = math.exp(slope)
    kappa = -slope / (2.0 - math.sqrt(2.0))

    f0 = fts[0]
    chain = chain_decay(f0, kappa, 8)
    rel1 = np.abs(chain[1] - fts[1]) / np.maximum(np.abs(fts[1]), 1e-12)

    # 三种宏观
    mac = macro_cases()
    fseed = f_of(np.array(SEED) / sum(SEED))

    rows = {}
    for name in ("geom", "phase", "capacity"):
        m = mac[name]
        f = m["f"]
        Bd = m["Bd"]
        GBd = m["GBd"]
        alpha = float(np.dot(f, Bd) / np.dot(Bd, Bd))
        rows[name] = dict(
            f=[float(x) for x in f],
            d=[float(x) for x in m["d"]],
            Bd=[float(x) for x in Bd],
            GBd=[float(x) for x in GBd],
            cos_Bd=cos(f, Bd),
            cos_GBd=cos(f, GBd),
            alpha=alpha,
            resid_local=float(np.linalg.norm(f - alpha * Bd)),
            pred_target=[float(x) for x in (fseed + Bd)],
            pred_source=[float(x) for x in (fseed + GBd)],
            err_target=float(np.linalg.norm(fseed + Bd - f)),
            err_source=float(np.linalg.norm(fseed + GBd - f)),
            events=m["events"],
        )

    # 判据（审计事实）
    checks = [
        dict(name="引擎对比衰减 = 1/81（斜率 -ln81）",
             pass_=abs(slope + math.log(81)) < 1e-6,
             detail=f"slope={slope:.6f}, rho={rho:.6f}"),
        dict(name="κ 标定：慢模速率 = 引擎速率",
             pass_=abs(kappa * (2 - math.sqrt(2)) + slope) < 1e-9,
             detail=f"κ={kappa:.4f}/gen, 慢模 {kappa*(2-math.sqrt(2)):.4f}"),
        dict(name="链只对上速率、瞬态形状不对（t=1 相对误差 >20%）",
             pass_=bool(np.max(rel1) > 0.2),
             detail=f"engine f1={np.round(fts[1],4).tolist()}, chain={np.round(chain[1],4).tolist()}, rel={np.round(rel1,3).tolist()}"),
        dict(name="geom 模式局部（cos(f,Bd) > 0.99）",
             pass_=rows["geom"]["cos_Bd"] > 0.99,
             detail=f"cos={rows['geom']['cos_Bd']:.4f}, α={rows['geom']['alpha']:.3f}"),
        dict(name="phase 模式局部（cos(f,Bd) > 0.99）",
             pass_=rows["phase"]["cos_Bd"] > 0.99,
             detail=f"cos={rows['phase']['cos_Bd']:.4f}, α={rows['phase']['alpha']:.3f}"),
        dict(name="capacity 不是种子调制模式（cos < 0）",
             pass_=rows["capacity"]["cos_Bd"] < 0.0,
             detail=f"cos={rows['capacity']['cos_Bd']:.4f}, events={rows['capacity']['events']}"),
        dict(name="链的 G 口径稳态预测误差 > 0.2（不能复现）",
             pass_=all(rows[k]["err_source"] > 0.2 for k in rows),
             detail=f"err_source={ {k: round(rows[k]['err_source'],3) for k in rows} }"),
    ]

    res = dict(
        decay=dict(slope=slope, rho=rho, kappa=kappa,
                   f_engine=[[float(x) for x in row] for row in fts[:4]],
                   f_chain=[[float(x) for x in row] for row in chain[:4]],
                   rel_t1=[float(x) for x in rel1]),
        f_seed=[float(x) for x in fseed],
        rows=rows,
        checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                 "pass": bool(c["pass_"])} for c in checks],
    )

    # ── 图 ──
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    # (a) 衰减：引擎 vs 链
    for i in range(3):
        axes[0].plot(fts[:5, i], "o-", lw=1.2, label=f"engine slot {i}")
        axes[0].plot(chain[:5, i], "--", lw=1.2, label=f"chain slot {i}")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("generation")
    axes[0].set_ylabel("|f_i|")
    axes[0].set_title(f"decay: engine vs chain (κ={kappa:.2f})")
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].grid(alpha=0.3)
    # (b) 三种模式：f_obs（柱） vs α·Bd（叉）
    xs3 = np.arange(3)
    w = 0.25
    for j in range(3):
        axes[1].bar(xs3 + (j - 1) * w,
                    [rows[n]["f"][j] for n in ("geom", "phase", "capacity")],
                    width=w, label=f"f_obs slot{j}")
    for j, mk in zip(range(3), ("x", "+", "*")):
        axes[1].plot(xs3,
                     [rows[n]["alpha"] * rows[n]["Bd"][j]
                      for n in ("geom", "phase", "capacity")],
                     mk, ms=9, mew=1.6, color="k",
                     label="α·Bd (local)" if j == 0 else None)
    axes[1].set_xticks(xs3)
    axes[1].set_xticklabels(["geom", "phase", "capacity"])
    axes[1].set_title("observed (bars) vs local α·Bd (marks)")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3, axis="y")
    # (c) cos 对比
    names = ("geom", "phase", "capacity")
    axes[2].plot([rows[n]["cos_Bd"] for n in names], "o-", label="cos(f, Bd)")
    axes[2].plot([rows[n]["cos_GBd"] for n in names], "s--", label="cos(f, GBd)")
    axes[2].axhline(0, color="k", lw=0.6)
    axes[2].set_xticks(range(3))
    axes[2].set_xticklabels(names)
    axes[2].set_title("local (Bd) vs chain (GBd)")
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_l2_slots_calib.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    def _j(o):
        if isinstance(o, complex):
            return [float(o.real), float(o.imag)]
        return str(o)

    with open(os.path.join(RES, "zero_l2_slots_calib.json"), "w") as f:
        json.dump(res, f, indent=1, default=_j)

    print("链式模型 × 引擎标定")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  κ={kappa:.4f}/gen, 引擎 f1={np.round(fts[1],4).tolist()}, 链 f1={np.round(chain[1],4).tolist()}")
    for n in ("geom", "phase", "capacity"):
        r = rows[n]
        print(f"  {n:9s} f={np.round(r['f'],3).tolist()} cos_Bd={r['cos_Bd']:.3f} "
              f"cos_GBd={r['cos_GBd']:.3f} α={r['alpha']:.3f} err_source={r['err_source']:.3f}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
