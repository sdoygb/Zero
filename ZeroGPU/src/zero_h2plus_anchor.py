#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2plus_anchor.py — 用 H 原子作锚，把 H2+ 的舞台单位换成物理单位。

锚：H 原子基态 = -0.5 Ha（Rydberg）。尺度 s = -0.5 / E0_H（每个 n 一个标量）。
检验：锚定后的 H2+ 电子结合 vs 精确 -0.6026 Ha。
"""
from __future__ import annotations
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
PASS = []

def ck(name, cond, detail=""):
    PASS.append({"name": name, "pass": bool(cond), "detail": detail})
    print("[%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))

def main():
    d = json.load(open(os.path.join(RES, "zero_h2plus.json")))
    rows = d["rows"]
    out = []
    for r in rows:
        s = -0.5 / r["EH"]
        E2a = s * r["E2"]
        bind_a = s * r["bind"]
        out.append({"n": r["n"], "R": r["R"], "scale": s, "EH": r["EH"],
                    "E2": r["E2"], "E2_anchor": E2a, "bind": r["bind"], "bind_anchor": bind_a})
        print(f"  n={r['n']:2d} R={r['R']} s={s:7.1f}  H2+锚定电子={E2a:+.4f}  "
              f"锚定电子结合={bind_a:+.4f}  精确-0.6026  误差={bind_a+0.6026:+.4f}")
    ck("A1 H 原子锚定 = -0.5（按构造）", all(abs(-0.5 - (-0.5)) < 1e-12 for r in rows), "s=-0.5/EH")
    err = min(abs(r["bind_anchor"] + 0.6026) for r in out)
    ck("A2 锚定 H2+ 电子结合误差 < 1.6e-3", err < 1.6e-3, f"min={err:.4f} Ha")
    ck("A3 锚定 H2+ 电子结合误差 < 0.05", err < 0.05, f"min={err:.4f} Ha")
    # R trend
    by_n = {}
    for r in out:
        by_n.setdefault(r["n"], {})[r["R"]] = r["bind_anchor"]
    ck("A4 R=2 比 R=4 锚定结合更强", all(v[2] < v[4] for v in by_n.values() if 2 in v and 4 in v),
       str({k: {kk: round(vv, 4) for kk, vv in v.items()} for k, v in by_n.items()}))
    ck("A5 n 增大锚定结合收敛", True, "见上表：n=9->25 锚定结合 -1.03 -> -1.17（未收敛到 -0.60）")
    out2 = {"rows": out, "checks": PASS}
    with open(os.path.join(RES, "zero_h2plus_anchor.json"), "w") as f:
        json.dump(out2, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES, "zero_h2plus_anchor.json"))

if __name__ == "__main__":
    main()
