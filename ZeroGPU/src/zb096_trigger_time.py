#!/usr/bin/env python3
"""P7′：记录越阈值触发事件的首次穿越时刻 t* 与延迟闭式

模型：R_{t+1}=R_t+g(m_t)，g(m)=c·m，m_t=r（等幅脉冲）；
     seed=K·R（对角核 K_eff 或一跳核 ΣK）；
     触发 = seed 首次 > τ。
检验闭式 t* = ceil(τ/(c·K_eff))。
"""
import json, math
from pathlib import Path
import numpy as np

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb096_trigger_time.json")

def simulate(tau, c=1.0, K_eff=1.0, n=500, g="linear", sat=None):
    R, t_star = 0.0, None
    for t in range(1, n + 1):
        m = c
        if g == "linear":
            R += m
        else:
            R += m / (1.0 + R / (sat or 100.0))
        seed = K_eff * R
        if seed > tau:
            t_star = t
            break
    return t_star, R

def main():
    rows = []
    for tau in (10, 50, 100):
        for c in (0.5, 1.0, 2.0):
            for K_eff in (0.5, 1.0, 2.0):
                t_sim, _ = simulate(tau, c=c, K_eff=K_eff)
                t_pred = math.floor(tau / (c * K_eff)) + 1
                rows.append({"tau": tau, "c": c, "K_eff": K_eff, "t_sim": t_sim, "t_pred": t_pred,
                             "match": t_sim == t_pred})
    ok = sum(r["match"] for r in rows)
    print(f"闭式 t*=ceil(τ/(c·K_eff)) 与模拟一致：{ok}/{len(rows)}")
    for r in rows[:6]:
        print(f"  τ={r['tau']:>3} c={r['c']:<4} K={r['K_eff']:<4} 模拟 {r['t_sim']:>3}  闭式 {r['t_pred']:>3}  {'✓' if r['match'] else '✗'}")

    # 非线性 g：闭式失效但 t* 仍存在且唯一
    t_non, _ = simulate(50, c=1.0, K_eff=1.0, g="sat", sat=20.0)
    print(f"\n非线性 g（饱和 20）：τ=50 的首次穿越时刻 {t_non}（存在）")
    out = {"rows": rows, "match_count": ok, "total": len(rows), "nonlinear_t_star": t_non}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)

if __name__ == "__main__":
    raise SystemExit(main())
