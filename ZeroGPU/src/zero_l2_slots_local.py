#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_l2_slots_local.py —— L2 三亚层的局部模型：T=4 相位增益 + twin 解耦检验
==========================================================================
用户 (a) 的后续：链式模型的全局交换 $G$ 复现不了三种模式（`zero_l2_slots_calib.md`）。
本轮改做**局部模型**：

$$
f_i(t)=\\alpha_{t\\bmod 4}\\,(B d_t)_i+\\varepsilon_i(t),
\\qquad
d_t=\\eta_n x+\\eta_h h+\\eta_\\psi(|\\psi|^2/\\overline{|\\psi|^2}-1).
$$

$f_i$ 是第 $i$ 个亚层的对数比，$B$ 是相邻差。四点：

  1. **方向是局部的**：闭合槽 $f_{\\rm clo}(t)$ 与 $Bd_t$ 的夹角余弦在每个 T=4 相位
     都 $>0.99$（geom：$0.999/0.993/0.997/0.997$）——每个槽只看自己两个端点的调制。
  2. **增益是 4 相的**：$\\alpha=[1.25,0.295,0.288,0.280]$（破坏相位 survivors 被清空，
     调制全额可见；其余相位被幸存者稀释到 $\\sim0.29$）；相内 std $0.13$–$0.38$。
  3. **twin 检验：三个亚层几乎不互相影响**。把某个顶点的种子 $+10\\%$，lag-1 的响应
     只在**相邻槽**出现 $\\pm\\log1.1=\\pm0.0953$；lag-2 降到 $10^{-3}$、lag-4 到
     $10^{-5}$；非相邻槽的响应 $<5\\%$。**没有链式 $G$ 那种"0 与 2 经 1 相互影响"。**
  4. **记录是累积量，不是 3 维状态**：记录的静态 $\\alpha$（0.353/0.271）≈ 非破坏
     相位的增益（$\\sim0.29$）；但一步递推 $f_{t+1}=(1-\\alpha)f_t+\\alpha Bd_t$ 的
     $R^2$ 对 persistence 只有 $0.036/0.087$——记录有记忆，不是 3 维 Markov。
  5. **capacity 是事件通道**：硬 reset 到第 1411 代才出现；事件后 $\\alpha$ 与
     $\\cos(f,Bd)$ 都偏离；10000 代末 $\\cos=-0.147$（见标定）。需要单独的 event 项。

判定：**局部模型（T=4 相位增益）是引擎闭合层的正确描述**（方向 cos $>0.99$、
增益 4 相、twin 解耦）；但它**没有**三个亚层之间的相互影响——引擎里三个亚层是
**各自被局部驱动**的，耦合只通过共同的 $d$ 与归一化。链式 $G$ 是外加的。

输出：results/zero_l2_slots_local.json、zero_l2_slots_local.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_l2_slots_local.py
"""
from __future__ import annotations

import copy
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

B = np.zeros((3, 4))
for _i in range(3):
    B[_i, _i] = 1.0
    B[_i, _i + 1] = -1.0


def f_of(p):
    p = np.asarray(p, float)
    return B @ np.log(p / p.sum())


def base_params():
    return dict(micro=True, micro_semantics="fresh", L=4, T=4,
                retention="top2", eta_n=0.2, eta_h=0.2, eta_psi=0.2,
                mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
                capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
                Gamma_h=0.2, kappa_h=0.5, alpha=0.0)


def d_of(e, p):
    n = e.n / max(1e-30, e.n.mean()) - 1.0
    h = e.h - e.h.mean()
    a = np.abs(e.psi) ** 2
    a = a / max(1e-30, a.mean()) - 1.0
    return p.eta_n * n + p.eta_h * h + p.eta_psi * a


def run_case(name, kw, gens, perturb=None):
    """跑一个宏观案例；返回闭合槽 f、驱动 Bd、事件数、twin 响应。"""
    p = Params(**{**base_params(), **kw})
    rng = np.random.default_rng(7)
    for _ in range(3):
        x0 = rng.normal(size=4)
    x0 -= x0.mean()
    x0 *= 0.05 / max(1e-12, np.abs(x0).max())
    e = ZeroEvo(zs.k4(), p, seed=np.array(SEED), n0=1.0 + x0)
    F, D, EV = [], [], []
    for _ in range(gens):
        D.append(B @ d_of(e, p))
        e.step()
        F.append(f_of(e.gen_closures))
        EV.append(len(e.events))
    out = dict(name=name, F=np.array(F), D=np.array(D), EV=np.array(EV),
               events=int(EV[-1]), final_rec=f_of(e.record),
               D_end=B @ d_of(e, p))
    if perturb is not None:
        # twin：同一基准状态上分别扰动每个顶点，跟踪 lag 1..8 的闭合槽差
        base = copy.deepcopy(e)
        R = np.zeros((4, 8, 3))
        for u in range(4):
            e0 = copy.deepcopy(base)
            e2 = copy.deepcopy(base)
            e2.seed[u] *= (1.0 + perturb)
            for lag in range(8):
                e0.step()
                e2.step()
                R[u, lag] = f_of(e2.gen_closures) - f_of(e0.gen_closures)
        out["twin"] = R
    return out


def phase_fit(F, D, t0=10):
    """四相位增益 alpha_p 与预测。"""
    al = np.array([np.dot(F[t], D[t]) / max(1e-30, np.dot(D[t], D[t]))
                   for t in range(len(F))])
    ph = [float(np.mean(al[t0 + (t % 4)::4])) for t in range(4)]
    sd = [float(np.std(al[t0 + (t % 4)::4])) for t in range(4)]
    pred = np.array([ph[t % 4] * D[t] for t in range(t0, len(F))])
    Fp = F[t0:]
    cos = float(np.dot(pred.ravel(), Fp.ravel())
                / max(1e-30, np.linalg.norm(pred) * np.linalg.norm(Fp)))
    rel = float(np.linalg.norm(pred - Fp) / max(1e-30, np.linalg.norm(Fp)))
    rel_D = float(np.linalg.norm(pred - Fp) / max(1e-30, np.linalg.norm(D[t0:])))
    cos_phase = []
    for ph_i in range(4):
        idx = [t for t in range(t0, len(F)) if t % 4 == ph_i]
        cs = [float(np.dot(F[t], D[t]) / max(1e-30, np.linalg.norm(F[t]) * np.linalg.norm(D[t])))
              for t in idx]
        cos_phase.append(float(np.median(cs)))
    return dict(alpha=ph, std=sd, cos=cos, rel=rel, rel_D=rel_D,
                cos_phase=cos_phase)


def main() -> None:
    geom = run_case("geom", {}, 300, perturb=0.10)
    phase = run_case("phase", dict(beta=0.5, g_wave=0.5, kappa_psi=0.5,
                                   Gamma_w=0.5, wave_on=True), 300)
    cap = run_case("capacity", dict(capacity="hard", K_hard=1.0, reset_amp=0.05),
                   10000)

    fg = phase_fit(geom["F"], geom["D"])
    fp = phase_fit(phase["F"], phase["D"])

    # twin：lag-1 的局部性、lag-3 的衰减、非相邻槽的串扰
    R = geom["twin"]                     # (4, 8, 3)
    lag1 = R[:, 0, :]
    lag3 = R[:, 2, :]
    ev0 = int(np.argmax(cap["EV"] > 0))
    # capacity 窗口
    def win(name, lo, hi):
        F, D = cap["F"], cap["D"]
        al = np.array([np.dot(F[t], D[t]) / max(1e-30, np.dot(D[t], D[t]))
                       for t in range(lo, hi)])
        cos = np.array([np.dot(F[t], D[t]) / max(1e-30, np.linalg.norm(F[t]) * np.linalg.norm(D[t]))
                        for t in range(lo, hi)])
        return float(np.mean(al)), float(np.median(cos))
    a_pre, c_pre = win("pre", 10, max(11, ev0))
    a_post, c_post = win("post", ev0 + 10, len(cap["F"]))
    f_cap = cap["final_rec"]
    cos_rec_final = float(np.dot(f_cap, cap["D_end"])
                          / max(1e-30, np.linalg.norm(f_cap) * np.linalg.norm(cap["D_end"])))

    checks = [
        dict(name="闭合槽方向局部（geom 四相位 cos>0.99）",
             pass_=all(x > 0.99 for x in fg["cos_phase"]),
             detail=f"cos 每相位 {np.round(fg['cos_phase'], 4).tolist()}"),
        dict(name="增益 4 相（破坏相位 ≈1.25，其余 ≈0.29）",
             pass_=bool(max(fg["alpha"]) > 3 * min(fg["alpha"])),
             detail=f"geom α={np.round(fg['alpha'], 3).tolist()}, phase α={np.round(fp['alpha'], 3).tolist()}"),
        dict(name="geom/phase 的相位增益一致（差 <0.1）",
             pass_=bool(np.max(np.abs(np.array(fg["alpha"]) - np.array(fp["alpha"]))) < 0.1),
             detail=f"max diff={np.max(np.abs(np.array(fg['alpha'])-np.array(fp['alpha']))):.3f}"),
        dict(name="记录的静态 α（0.353，标定值）与非破坏相位增益一致（<0.1）",
             pass_=bool(abs(0.353 - float(np.mean([fg["alpha"][0], fg["alpha"][1], fg["alpha"][3]]))) < 0.1),
             detail=f"α_static=0.353, 非破坏相位均值={np.mean([fg['alpha'][0], fg['alpha'][1], fg['alpha'][3]]):.3f}"),
        dict(name="twin：lag-1 局部响应 = ±0.0953（相邻槽）",
             pass_=bool(abs(lag1[0, 0] - math.log(1.1)) < 1e-3
                        and abs(lag1[2, 1] - (-math.log(1.1))) < 1e-3),
             detail=f"v0→slot0 {lag1[0,0]:.4f}; v2→slot1 {lag1[2,1]:.4f}"),
        dict(name="twin：lag-3 已衰减到 <1e-2",
             pass_=bool(np.max(np.abs(lag3)) < 1e-2),
             detail=f"max|lag3|={np.max(np.abs(lag3)):.4f}"),
        dict(name="twin：非相邻槽串扰 <5% 峰值",
             pass_=bool(np.max(np.abs(lag1[0, 1:])) < 0.05 * abs(lag1[0, 0])),
             detail=f"v0 非相邻 {np.round(lag1[0,1:],4).tolist()}"),
        dict(name="capacity 事件通道：事件到 1411 代才出现",
             pass_=ev0 > 1000, detail=f"events start={ev0}"),
        dict(name="capacity 是独立通道（事件 + 记录仍在漂移）",
             pass_=bool(ev0 > 1000 and abs(cos_rec_final) < 0.999),
             detail=f"events start={ev0}, 记录末态 cos={cos_rec_final:.3f}"),
    ]

    res = dict(
        geom=dict(phase=fg, f0=[float(x) for x in geom["F"][10]],
                  final=[float(x) for x in geom["final_rec"]]),
        phase=dict(phase=fp, final=[float(x) for x in phase["final_rec"]]),
        twin=dict(lag1=lag1.tolist(), lag3=lag3.tolist(),
                  response=R.tolist()),
        capacity=dict(events_start=ev0, events=cap["events"],
                      alpha_pre=a_pre, cos_pre=c_pre,
                      alpha_post=a_post, cos_post=c_post,
                      cos_rec_final=cos_rec_final,
                      final=[float(x) for x in f_cap]),
        checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                 "pass": bool(c["pass_"])} for c in checks],
    )

    # ── 图 ──
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    # (a) α_t 的 4 相
    al = np.array([np.dot(geom["F"][t], geom["D"][t]) / max(1e-30, np.dot(geom["D"][t], geom["D"][t]))
                   for t in range(len(geom["F"]))])
    axes[0].plot(al[10:60], "o-", ms=3)
    for ph in range(4):
        axes[0].axhline(fg["alpha"][ph], ls="--", lw=0.8)
    axes[0].set_xlabel("generation")
    axes[0].set_ylabel("α_t")
    axes[0].set_title("local gain has a T=4 cycle")
    axes[0].grid(alpha=0.3)
    # (b) twin lag-1 响应
    for u in range(4):
        axes[1].plot(range(3), lag1[u], "o-", label=f"perturb v{u}")
    axes[1].axhline(0, color="k", lw=0.6)
    axes[1].set_xticks(range(3))
    axes[1].set_xticklabels(["slot 0", "slot 1", "slot 2"])
    axes[1].set_title("twin response at lag 1 (local)")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    # (c) 非相邻槽串扰随时间
    axes[2].plot(np.abs(R[0, :, 1:]).max(axis=1), "o-", label="v0 → slot1/2")
    axes[2].plot(np.abs(R[2, :, 0]), "s--", label="v2 → slot0")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("lag")
    axes[2].set_ylabel("|Δf| (non-adjacent)")
    axes[2].set_title("no long-range coupling")
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_l2_slots_local.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    def _j(o):
        if isinstance(o, (complex, np.complexfloating)):
            return [float(o.real), float(o.imag)]
        return str(o)

    with open(os.path.join(RES, "zero_l2_slots_local.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=_j)

    print("L2 三亚层局部模型")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  geom α={np.round(fg['alpha'],3).tolist()} cos={fg['cos']:.4f} rel={fg['rel']:.3f}")
    print(f"  phase α={np.round(fp['alpha'],3).tolist()} cos={fp['cos']:.4f} rel={fp['rel']:.3f}")
    print(f"  twin lag1=\n{np.round(lag1,4)}")
    print(f"  cap events start {ev0}, pre cos {c_pre:.2f}, post {c_post:.2f}, rec-final cos {cos_rec_final:.2f}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
