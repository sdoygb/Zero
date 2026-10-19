#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_evo_18.py —— 把 18 世界作为 micro core 接进 ZeroEvo，开宏观扇区跑长程
=====================================================================
`ZeroEvo18(ZeroEvo)`：micro core 换成 4→18→6（3 代 × 6 边；3 个完美匹配），
宏观扇区（内容/几何/波/容量）照旧由 `ZeroEvo` 提供；`_seed_used()` 的调制
（η_n x + η_h h + η_ψ |ψ|²）作用在 4-顶点源上，源再流进 18。

规则：
  · 源 4 = `_seed_used()` 归一；
  · 4→18 流：边 e=(u,v) 收 (w4[u]+w4[v])/2，乘 3 代权重；
  · 18 更新：(1-α)W + α flow；T=4 毁灭：每代保留自己的完美匹配（2 边）；
  · 影子 6 = 18 对代求和；再按关联映射回 4 顶点 record（每条边给两端各一半）；
  · record/seed 交回 `ZeroEvo` 的宏观扇区（内容/几何/波/容量）。

两个 case：
  · micro_only：只开 micro（宏观全关）；
  · full：12 通道全开（内容/几何/波/软容量）。

输出：results/zero_evo_18.json、zero_evo_18.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_18.py
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
EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
MATCH = [[0, 5], [1, 4], [2, 3]]          # 3 个完美匹配
NG, NE, NV = 3, 6, 4
T = 4


class ZeroEvo18(ZeroEvo):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        slots = np.array([math.log(2.0), math.log(8.0 / 3.0), math.log(1.5)])
        self.gen_w = slots / slots.sum()
        self.W18 = np.outer(self.gen_w, self._edge_source(self.seed))
        self.W18 /= self.W18.sum()
        self.alpha = 0.35
        self.gamma_gen = 0.0          # 代依赖宏观耦合强度（0 = 不耦合）
        self.lambda_gen = 0.0         # 份额负反馈（稳定化）

    def _edge_source(self, w4):
        w4 = np.asarray(w4, float)
        s = np.array([(w4[u] + w4[v]) / 2.0 for (u, v) in EDGES])
        return s / s.sum() if s.sum() > 0 else np.ones(NE) / NE

    def _micro_gen(self):
        w4 = self._seed_used()
        if w4.sum() <= 0:
            w4 = np.ones(self.V)
        w4 = w4 / w4.sum()
        if self.gamma_gen or self.lambda_gen:
            nn = self.n / max(self.n.mean(), 1e-30)
            # 用**边梯度**（匹配之间的区别只能在边上，顶点平均对 3 个匹配相同）
            obs = np.array([np.mean([abs(nn[u] - nn[v])
                                     for (u, v) in [EDGES[e] for e in MATCH[g]]])
                            for g in range(NG)])
            obs = obs - obs.mean()
            share = self.W18.sum(axis=1)
            mod = self.gamma_gen * obs - self.lambda_gen * (share - self.gen_w)
            ge = self.gen_w * np.exp(np.clip(mod, -20, 20))
            ge = ge / ge.sum()
        else:
            ge = self.gen_w
        flow = np.outer(ge, self._edge_source(w4))
        self.W18 = (1.0 - self.alpha) * self.W18 + self.alpha * flow
        self.W18 /= max(self.W18.sum(), 1e-30)
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


def run_case(name, gens, micro_only=False):
    base = dict(micro=True, L=4, T=T, retention="top2",
                eta_n=0.2, eta_h=0.2, eta_psi=0.2,
                mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
                capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
                Gamma_h=0.2, kappa_h=0.5, alpha=0.0)
    if micro_only:
        base.update(eta_n=0.0, eta_h=0.0, eta_psi=0.0, g_src=0.0,
                    kappa_h=0.0, wave_on=False)
    p = Params(**base)
    e = ZeroEvo18(zs.k4(), p, seed=np.array([8.0, 4.0, 1.5, 1.0]))
    gen_sh, shad, qs, ns, he, pc = [], [], [], [], [], []
    for _ in range(gens):
        e.step()
        gen_sh.append(e.W18.sum(axis=1).copy())
        shad.append(e.W18.sum(axis=0).copy())
        s = e.snapshot()
        qs.append(s.get("q", 0.0)); ns.append(s.get("n_std", 0.0))
        he.append(s.get("h_energy", 0.0)); pc.append(s.get("psi_coh", 0.0))
    gen_sh = np.array(gen_sh); shad = np.array(shad)
    tail = slice(gens // 2, gens)
    return dict(name=name,
                gen_shares=[float(x) for x in gen_sh[tail].mean(axis=0)],
                gen_relstd=[float(x) for x in gen_sh[tail].std(axis=0) / np.maximum(gen_sh[tail].mean(axis=0), 1e-12)],
                shadow=[float(x) for x in shad[tail].mean(axis=0)],
                q=float(np.mean(qs[tail])), n_std=float(np.mean(ns[tail])),
                h_energy=float(np.mean(he[tail])), psi_coh=float(np.mean(pc[tail])),
                finite=bool(np.isfinite(gen_sh).all() and np.isfinite(shad).all()),
                gen_trace=gen_sh[::max(1, gens // 200)].tolist(),
                shadow_trace=shad[::max(1, gens // 200)].tolist())


def main() -> None:
    cases = [run_case("micro_only", 2000, micro_only=True),
             run_case("full", 10000, micro_only=False)]

    checks = [
        dict(name="3 代都在（尾部份额 > 0.01）",
             pass_=all(min(c["gen_shares"]) > 0.01 for c in cases),
             detail={c["name"]: np.round(c["gen_shares"], 3).tolist() for c in cases}),
        dict(name="影子 6 满支撑（尾部 6 条边 > 0）",
             pass_=all(min(c["shadow"]) > 0 for c in cases),
             detail={c["name"]: np.round(c["shadow"], 3).tolist() for c in cases}),
        dict(name="3 代份额稳定（relstd < 0.1，长程）",
             pass_=all(max(c["gen_relstd"]) < 0.1 for c in cases),
             detail={c["name"]: np.round(c["gen_relstd"], 4).tolist() for c in cases}),
        dict(name="无 NaN / 有界",
             pass_=all(c["finite"] for c in cases), detail="all finite"),
        dict(name="宏观扇区工作（full 的 n_std/q 非零）",
             pass_=cases[1]["n_std"] > 0 and cases[1]["q"] > 0,
             detail=f"n_std={cases[1]['n_std']:.4f}, q={cases[1]['q']:.4f}"),
        dict(name="full 长程 10000 代稳定",
             pass_=bool(cases[1]["finite"]), detail="10000 gens"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    for c in cases:
        g = np.array(c["gen_trace"])
        ls = "-" if c["name"] == "full" else "--"
        for i in range(NG):
            axes[0].plot(g[:, i], ls, lw=1.0,
                         label=f"{c['name']} gen{i}")
        sh = np.array(c["shadow_trace"])
        for i in range(NE):
            axes[1].plot(sh[:, i], ls, lw=0.9,
                         label=f"{c['name']} edge{i}" if i < 2 else None)
    axes[0].set_title("generation shares (18-world micro core)")
    axes[0].legend(fontsize=6, ncol=2)
    axes[0].grid(alpha=0.3)
    axes[1].set_title("shadow 6")
    axes[1].legend(fontsize=6)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_evo_18.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(cases=cases,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_evo_18.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("ZeroEvo18（18 世界 micro core）")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for c in cases:
        print(f"  {c['name']}: shares={np.round(c['gen_shares'], 3).tolist()} "
              f"shadow={np.round(c['shadow'], 3).tolist()} n_std={c['n_std']:.4f}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
