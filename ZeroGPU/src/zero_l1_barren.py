#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L1 词层的原生梯度方差：面积相位 + 双覆盖读出。

原生对象（全部来自库内，不引入外部线路）
  · L1 词层：长度 L、和为 0 的 ±1 词 W_L，|W_L| = C(L,L/2)（Z161）。
  · 高度 h(w)=max_k |前缀和|；层 b_h（Z161）。
  · 面积 a(w)=Σ_k 前缀和_k；a(-w)=-a(w)（Z048）。
  · 双覆盖 R|w>=|-w>；h(-w)=h(w)（Z047）。

原生参数化族与读出
  |psi_theta> = N^{-1/2} Σ_w e^{i theta a(w)} |w>
  C(theta) = <psi_theta| R |psi_theta> = N^{-1} Σ_w cos(2 a(w) theta)
  梯度 dC/dtheta = -2 N^{-1} Σ_w a(w) sin(2 a(w) theta)

定理（本脚本核验）
  Var_theta[dC/dtheta] = 4 N^{-2} Σ_{a≠0} a^2 n_a^2，n_a = #{w: a(w)=a}。
  逐层：同一公式限制在 h(w)=h，N=b_h。

判据
  N1 精确公式 vs θ 网格（L=4/8/12/16）。
  N2 全局方差随 L 增长（不是 2^{-L}）。
  N3 平均 |梯度| 保持 O(1)（≥1）。
  N4 顶层 b_H=2，顶层方差最大。
  N5 log Var vs log N 斜率 >0（没有 1/N 抑制）。
  N6 底层（b_1=2L）方差也随 L 增长。
  N7 L=16 的逐层方差全正，顶层最大。

输出：results/zero_l1_barren.json
"""
from __future__ import annotations

import json
import math
import pathlib
import sys
from collections import Counter

import numpy as np

HERE = pathlib.Path(__file__).resolve()
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
import zero_sim as zs  # noqa: E402

OUT = ROOT / "results" / "zero_l1_barren.json"
OK, BAD, ROWS = [], [], []


def ck(name: str, cond: bool, detail: str = "") -> bool:
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def exact_var(a: np.ndarray, N: int) -> float:
    """4 N^{-2} Σ_{a≠0} a^2 n_a^2."""
    c = Counter(a.tolist())
    s = sum(aa * aa * n * n for aa, n in c.items() if aa != 0)
    return 4.0 * s / N ** 2


def grid_var(a: np.ndarray, ntheta: int = 4001) -> tuple[float, float]:
    """θ 网格上的 Var[dC/dθ] 与 mean|dC/dθ|。"""
    th = np.linspace(0.0, 2.0 * np.pi, ntheta)[:-1]
    # C(theta) = mean cos(2 a theta); 逐 θ 向量化
    C = np.array([np.mean(np.cos(2.0 * t * a)) for t in th])
    g = np.gradient(C, th)
    return float(np.var(g)), float(np.mean(np.abs(g)))


def main():
    print("=" * 72)
    print("zero_l1_barren.py — L1 词层的原生梯度方差（面积相位 + 双覆盖读出）")
    print("=" * 72)
    table = []
    for L in range(4, 21, 2):
        ws = zs.balanced_words(L)
        N = len(ws)
        a = np.array([zs.word_area(w) for w in ws], dtype=np.int64)
        h = np.array([zs.layer_of(w) for w in ws], dtype=np.int64)
        var = exact_var(a, N)
        mg = None
        if L <= 16:
            vg, mg = grid_var(a)
            rel = abs(vg - var) / var
            ck(f"N1 L={L} exact vs grid", rel < 0.01, f"rel={rel:.2e}")
        c = Counter(a.tolist())
        n_a = max(c.values()) if c else 0
        top = h == L // 2
        bot = h == 1
        vtop = exact_var(a[top], int(top.sum()))
        vbot = exact_var(a[bot], int(bot.sum()))
        table.append({"L": L, "N": N, "var": var, "mean_abs_grad": mg,
                      "var_top": vtop, "b_top": int(top.sum()),
                      "var_bottom": vbot, "b_bottom": int(bot.sum()),
                      "max_area": int(np.max(np.abs(a))) if a.size else 0,
                      "n_area_values": len(c)})
        print(f"    L={L:2d} N={N:7d}  Var={var:8.3f}  mean|g|={('%.3f' % mg) if mg is not None else '—'}"
              f"  Var_top={vtop:9.1f} (b={int(top.sum())})  Var_bot={vbot:.2f} (b={int(bot.sum())})")

    ck("N2 global variance grows with L", table[-1]["var"] > 2.0 * table[0]["var"],
       f"Var(20)/Var(4)={table[-1]['var']/table[0]['var']:.2f}")
    mgs = [r["mean_abs_grad"] for r in table if r["mean_abs_grad"] is not None]
    ck("N3 mean |gradient| stays O(1)", min(mgs) > 1.0, f"min={min(mgs):.3f}")
    ck("N4 top layer b_H=2 and largest", table[0]["b_top"] == 2 and
       all(r["var_top"] > r["var_bottom"] for r in table),
       "b_H=2; Var_top > Var_bottom for all L")
    ln = np.log([r["N"] for r in table]); lv = np.log([r["var"] for r in table])
    slope_n = float(np.polyfit(ln, lv, 1)[0])
    ck("N5 no 1/N suppression", slope_n > 0.0, f"d log Var / d log N = {slope_n:.3f}")
    ck("N6 bottom-layer variance grows", table[-1]["var_bottom"] > table[0]["var_bottom"],
       f"{table[0]['var_bottom']:.2f} -> {table[-1]['var_bottom']:.2f}")

    # L=16 逐层
    L = 16
    ws = zs.balanced_words(L)
    a = np.array([zs.word_area(w) for w in ws], dtype=np.int64)
    h = np.array([zs.layer_of(w) for w in ws], dtype=np.int64)
    layers = []
    N16 = len(ws)
    for hh in range(1, L // 2 + 1):
        m = h == hh
        bh = int(m.sum())
        vh = exact_var(a[m], bh)
        w = bh / N16
        layers.append({"h": hh, "b_h": bh, "var": vh, "weight": w,
                       "eff_var": (w ** 2) * vh,
                       "max_area": int(np.max(np.abs(a[m]))) if m.any() else 0})
    ck("N7 L=16 layer variances positive, top largest",
       all(x["var"] > 0 for x in layers) and layers[-1]["var"] == max(x["var"] for x in layers),
       f"top var={layers[-1]['var']:.1f}")
    peak = max(layers, key=lambda x: x["eff_var"])
    ck("N8 effective layer variance peaks inside (not at top/bottom)",
       1 < peak["h"] < L // 2, f"peak h={peak['h']}, eff={peak['eff_var']:.3f}")
    ck("N9 total layer-phase variance is O(1)+",
       sum(x["eff_var"] for x in layers) > 1.0,
       f"sum eff={sum(x['eff_var'] for x in layers):.3f}")
    print("    L=16 layers (h,b_h,Var,weight,eff):",
          [(x["h"], x["b_h"], round(x["var"], 2), round(x["weight"], 4), round(x["eff_var"], 4))
           for x in layers])

    out = {"script": "zero_l1_barren.py", "checks": ROWS,
           "n_pass": len(OK), "n_fail": len(BAD),
           "table": table, "layers_L16": layers, "slope_logvar_logn": slope_n,
           "formula": "Var = 4 N^-2 sum_{a!=0} a^2 n_a^2"}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
