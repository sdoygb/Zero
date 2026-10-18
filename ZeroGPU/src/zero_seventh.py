# -*- coding: utf-8 -*-
"""第七条边（L0 内容 → L2 seed）的受控判定。

定义：在 seed 上加一个无量纲偏置
    seed_v <- seed_v * (1 + λ (m_v/mean(m) - 1))，
m 取边⑥的记录（default）或边⑤的本代闭合；归一化保持 seed 总量不变。
λ=0 回到 Z161/Z162 的六边基线；λ≠0 是新的“第七条边”，不是 Z088 的六边之一。

输出：results/zero_seventh.json
"""
from __future__ import annotations
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo, Params  # noqa: E402
from z0_engine import path_adj  # noqa: E402
import zero_sim as zs  # noqa: E402

OUT = os.path.join(HERE, "..", "results")
CK: list[dict] = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def _seed_for(V):
    if V == 4:
        return np.array([8.0, 4.0, 1.5, 1.0])
    return np.array([2.0 ** (V - 1 - i) for i in range(V)])


def run_case(adj, sem, ret, T, norm, lam, gens, mode="hist"):
    V = len(adj)
    e = ZeroEvo(adj, Params(micro=True, micro_semantics=sem, L=4, T=T,
                            retention=ret, norm_mode=norm, l0_back=lam,
                            l0_back_mode=mode, wave_on=False,
                            poisson_on=False), seed=_seed_for(V))
    e.run(gens, every=max(1, gens // 10))
    s = e.snapshot()
    rec = e.record / e.record.sum()
    slots = []
    if V == 4:
        slots = sorted(abs(math.log(max(rec[i], 1e-300) / max(rec[i + 1], 1e-300)))
                       for i in range(3))
    return {"q_live": s["q_live"], "q_hist": s["q_hist"],
            "support": int((rec > 1e-4).sum()), "max": float(rec.max()),
            "record": [float(x) for x in rec], "slots": [float(x) for x in slots]}


def run_stability(adj, sem, ret, T, norm, lam, gens=20000, mode="hist"):
    V = len(adj)
    e = ZeroEvo(adj, Params(micro=True, micro_semantics=sem, L=4, T=T,
                            retention=ret, norm_mode=norm, l0_back=lam,
                            l0_back_mode=mode, wave_on=False,
                            poisson_on=False), seed=_seed_for(V))
    qs = []
    for g in range(gens):
        e.step()
        if g >= gens - 1000 and g % 10 == 0:
            qs.append(e.snapshot()["q_live"])
    return {"q_last": float(qs[-1]), "q_range": float(max(qs) - min(qs)),
            "support": int((e.record / e.record.sum() > 1e-4).sum())}


def main():
    scans = {}
    specs = [
        ("K4_persistent", zs.k4(), "persistent", "top2", 4, "common",
         [0.0, 0.001, 0.1, 1.0], 5000),
        ("path16_persistent", path_adj(16), "persistent", "last", 0, "separate",
         [0.0, 0.001, 0.1], 20000),
        ("K4_fresh", zs.k4(), "fresh", "top2", 4, "common",
         [0.0, 1.0], 5000),
        ("path16_fresh", path_adj(16), "fresh", "top2", 4, "common",
         [0.0, 0.01, 0.1, 0.5, 1.0, 5.0], 10000),
    ]
    for name, adj, sem, ret, T, norm, lams, gens in specs:
        rows = []
        print(f"{name}（{sem}, {ret}, T={T}, {norm}）")
        for lam in lams:
            r = run_case(adj, sem, ret, T, norm, lam, gens)
            r["lambda"] = lam
            rows.append(r)
            print(f"  λ={lam:<6} q_live={r['q_live']:.6f} q_hist={r['q_hist']:.6f} "
                  f"support={r['support']} max={r['max']:.4f}")
        scans[name] = rows

    stab = []
    print("path16 fresh 稳定性（20k 代，末 1k 代 q 范围）")
    for lam in (0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0):
        r = run_stability(path_adj(16), "fresh", "top2", 4, "common", lam)
        r["lambda"] = lam
        stab.append(r)
        print(f"  λ={lam:<6} q_last={r['q_last']:.6f} range={r['q_range']:.2e} "
              f"support={r['support']}")

    print("=" * 74)
    print("第七条边判据")
    print("=" * 74)
    k4p = {r["lambda"]: r for r in scans["K4_persistent"]}
    p16 = {r["lambda"]: r for r in scans["path16_persistent"]}
    k4f = {r["lambda"]: r for r in scans["K4_fresh"]}
    f16 = {r["lambda"]: r for r in scans["path16_fresh"]}
    st = {r["lambda"]: r for r in stab}

    ck("λ=0 回到 K4 persistent 基线（q=0.395957）",
       abs(k4p[0.0]["q_live"] - 0.395957) < 1e-5, f"{k4p[0.0]['q_live']:.6f}")
    ck("K4 persistent 任意 λ>0 都 winner-take-all（q=1, 支撑 1）",
       k4p[0.1]["q_live"] > 0.999 and k4p[0.1]["support"] == 1,
       f"λ=0.1 q={k4p[0.1]['q_live']:.6f}")
    ck("path16 persistent λ>0 也 winner-take-all",
       p16[0.1]["q_live"] > 0.999 and p16[0.1]["support"] == 1,
       f"λ=0.1 q={p16[0.1]['q_live']:.6f}")
    ck("K4 fresh 对 λ 不敏感（m 均匀）",
       abs(k4f[0.0]["q_live"] - k4f[1.0]["q_live"]) < 1e-9
       and abs(k4f[1.0]["q_live"] - 0.25) < 1e-9,
       f"{k4f[0.0]['q_live']:.6f}")
    ck("path16 fresh λ=0.01 发生选择（q 升、支撑 8）",
       f16[0.01]["q_live"] > 0.15 and f16[0.01]["support"] == 8,
       f"q={f16[0.01]['q_live']:.6f} support={f16[0.01]['support']}")
    ck("path16 fresh λ≥0.5 进入振荡（末 1k 代 q 范围 >0.01）",
       st[0.5]["q_range"] > 0.01 and st[1.0]["q_range"] > 0.01,
       f"λ=1 range={st[1.0]['q_range']:.3f}")

    out = {"scans": scans, "stability": stab, "checks": CK,
           "n_fail": sum(1 for c in CK if not c["pass"]),
           "definition": "seed_v *= 1 + λ (m_v/mean(m) - 1); m=record (hist) or gen_closures (live)",
           "anchor": "λ dimensionless; no new dimensionful anchor"}
    with open(os.path.join(OUT, "zero_seventh.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(OUT, "zero_seventh.json"))


if __name__ == "__main__":
    main()
