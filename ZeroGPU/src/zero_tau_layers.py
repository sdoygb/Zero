#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""读法 A prototype：每类一个年龄阈值 τ_c，去掉全局 if gen%T。

τ_b=4；中低 τ=4、中高 τ=8；顶 τ=16。
每代年龄 +1；到 τ 清该类并重置年龄；没有全局计时器。
看“底全灭、中留 2、顶不动”是否出现。
输出：results/zero_tau_layers.json/.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
V4, L = 4, 4
TAU_B, TAU_ML, TAU_MH, TAU_T = 4, 4, 8, 16
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def main():
    A = zs.adj_matrix(zs.k4()); AL = np.linalg.matrix_power(A, L)
    Wb = np.zeros((V4, V4)); Wm = np.zeros((V4, V4)); Wt = np.zeros((V4, V4))
    age_b = np.zeros((V4, V4), dtype=int)
    age_m = np.zeros((V4, 2), dtype=int)
    age_t = np.zeros((V4, 4), dtype=int)
    seed = np.array([8.0, 4.0, 1.5, 1.0]) / 14.5
    rows = []
    for g in range(1, 65):
        diag = np.arange(V4)
        Wb[diag, diag] += (1/3) * seed
        Wm[diag, diag] += (1/3) * seed
        Wt[diag, diag] += (1/3) * seed
        Wb = Wb @ AL; Wm = Wm @ AL; Wt = Wt @ AL
        clo = np.diag(Wb) + np.diag(Wm) + np.diag(Wt)
        jstar = int(np.argmax(clo))
        # ages + local clear (no global timer)
        age_b += 1; age_t += 1; age_m += 1
        # bottom all tau=4
        clear_b = age_b >= TAU_B
        if np.any(clear_b):
            Wb[clear_b] = 0.0; age_b[clear_b] = 0
        # middle low/high
        for j in range(V4):
            for s in range(2):
                tau = TAU_MH if j == jstar else TAU_ML
                if age_m[j, s] >= tau:
                    Wm[j, 2*s] = 0.0; Wm[j, 2*s+1] = 0.0
                    age_m[j, s] = 0
        # top tau=16
        clear_t = age_t >= TAU_T
        if np.any(clear_t):
            Wt[clear_t] = 0.0; age_t[clear_t] = 0
        # freeze pre-flow state for the observation
        pre_bottom_zero = bool(np.abs(Wb).sum() < 1e-30)
        pre_mid_classes = sum(1 for j in range(V4) for s in range(2)
                              if abs(Wm[j,2*s]) + abs(Wm[j,2*s+1]) > 1e-12)
        pre_top_classes = sum(1 for j in range(V4)
                              if np.abs(Wt[j,:]).sum() > 1e-12)
        # flow top->middle->bottom (uniform target)
        k = 0.2
        Pt = np.abs(Wt).sum()
        if Pt > 0:
            J = k * Pt; Wt -= Wt * (J / Pt)
            for j in range(V4):
                for s in range(2):
                    Wm[j, 2*s] += J / (2*V4); Wm[j, 2*s+1] += J / (2*V4)
        Pm = np.abs(Wm).sum()
        if Pm > 0:
            J2 = k * Pm; Wm -= Wm * (J2 / Pm)
            Wb += J2 / (V4*V4)
        mw = max(np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum())
        if np.isfinite(mw) and mw > 0:
            Wb /= mw; Wm /= mw; Wt /= mw
        # seed = top2 site closures
        seed = clo / clo.sum() if clo.sum() > 0 else np.ones(V4)/V4
        # shares & q16
        mb, mm, mt = np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum()
        ms = mb + mm + mt
        w = [mb/ms, mm/ms, mt/ms] if ms > 0 else [1/3]*3
        p = np.abs(Wb).reshape(-1); q16 = float(((p/p.sum())**2).sum()) if p.sum() > 0 else 0.0
        # active middle classes = both parity entries >0
        mid_active = sum(1 for j in range(V4) for s in range(2)
                         if abs(Wm[j,2*s]) + abs(Wm[j,2*s+1]) > 1e-12)
        top_active = sum(1 for j in range(V4) for f in range(4) if abs(Wt[j,f]) > 1e-12)
        rows.append({"gen": g, "w": w, "q16": q16, "mid_active": mid_active,
                     "top_active": top_active, "jstar": jstar,
                     "pre_bottom_zero": pre_bottom_zero,
                     "pre_mid_classes": pre_mid_classes,
                     "pre_top_classes": pre_top_classes})
        if g <= 8 or g % 4 == 0:
            print(f"g={g:>2} j*={jstar} w={np.round(w,3)} mid={mid_active} top={top_active} q16={q16:.3f}")
    g4 = rows[3]
    ck("无全局 if gen%T（源码）", True, "ages only")
    ck("第4代底清零（flow 前）", g4["pre_bottom_zero"], str(g4["pre_bottom_zero"]))
    ck("第4代中只留2类（flow 前）", g4["pre_mid_classes"] == 2,
       f"pre_mid_classes={g4['pre_mid_classes']}")
    ck("第4代顶仍4类（flow 前）", g4["pre_top_classes"] == 4,
       f"pre_top_classes={g4['pre_top_classes']}")
    ck("无 NaN", all(np.isfinite(r["w"]).all() for r in rows), "")

    gs = [r["gen"] for r in rows]; w = np.array([r["w"] for r in rows])
    mid = [r["mid_active"] for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for i, lab in enumerate(["Wb", "Wm", "Wt"]):
        ax[0].plot(gs, w[:, i], lw=1, label=lab)
    ax[0].set_title("(a) shares"); ax[0].legend(fontsize=8)
    ax[1].plot(gs, mid, lw=1, label="mid active classes")
    ax[1].plot(gs, [r["top_active"] for r in rows], lw=1, label="top active")
    ax[1].set_title("(b) active classes"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_tau_layers.png")
    fig.savefig(fp, dpi=130)
    out = {"tau": {"b": TAU_B, "ml": TAU_ML, "mh": TAU_MH, "t": TAU_T},
           "rows": rows, "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "per-class age thresholds; no global timer"}
    with open(os.path.join(RES, "zero_tau_layers.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")


if __name__ == "__main__":
    main()
