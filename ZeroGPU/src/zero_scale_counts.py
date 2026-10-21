#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""胞腔计数随基组的标度：nao, ndet=C(nao,nocc)^2, C(空间 union), C2(自旋), 参数。

小体系枚举行列式；H2O cc-pVDZ 只给上界（2^nao / C*6^4）。
输出：results/zero_scale_counts.json
"""
from __future__ import annotations

import itertools
import json
import pathlib
from math import comb

OUT = pathlib.Path(__file__).resolve().parents[1] / "results" / "zero_scale_counts.json"


def counts(nao, nocc):
    strs = list(itertools.combinations(range(nao), nocc))
    hf = (1 << nocc) - 1
    labs, labs2 = set(), set()
    for sa in strs:
        sa_m = sum(1 << i for i in sa)
        for sb in strs:
            sb_m = sum(1 << i for i in sb)
            holes = (hf & ~sa_m) | (hf & ~sb_m)
            parts = (sa_m | sb_m) & ~hf
            lab = holes | parts
            labs.add(lab)
            nh_a = bin(hf & ~sa_m).count("1"); nh_b = bin(hf & ~sb_m).count("1")
            np_a = bin(sa_m & ~hf).count("1"); np_b = bin(sb_m & ~hf).count("1")
            labs2.add(lab * 6 ** 4 + nh_a * 6 ** 3 + nh_b * 6 ** 2 + np_a * 6 + np_b)
    return len(labs), len(labs2)


def main():
    rows = []
    cases = [
        ("H2O STO-3G", 7, 10),
        ("LiH STO-3G", 6, 4),
        ("LiH 6-31G", 11, 4),
        ("H2O 6-31G", 13, 10),
    ]
    for name, nao, nelec in cases:
        nocc = nelec // 2
        C, C2 = counts(nao, nocc)
        ndet = comb(nao, nocc) ** 2
        H = max(nelec, 2 * nao - nelec) + 2
        params = C2 + 3 * H
        rows.append({"name": name, "nao": nao, "nelec": nelec, "ndet": ndet,
                     "C": C, "C2": C2, "H": H, "params": params,
                     "ratio": ndet / params})
        print(f"  {name:12s} nao={nao:2d} ndet={ndet:>10d} C={C:4d} C2={C2:5d} "
              f"params={params:5d} ndet/params={ndet/params:8.1f}")
    # upper bounds for cc-pVDZ
    nao, nelec = 24, 10
    nocc = nelec // 2
    ndet = comb(nao, nocc) ** 2
    Cub = 2 ** nao
    C2ub = Cub * 6 ** 4
    H = max(nelec, 2 * nao - nelec) + 2
    params_ub = C2ub + 3 * H
    rows.append({"name": "H2O cc-pVDZ (bound)", "nao": nao, "nelec": nelec, "ndet": ndet,
                 "C": Cub, "C2": C2ub, "H": H, "params": params_ub,
                 "ratio": ndet / params_ub})
    print(f"  H2O cc-pVDZ  nao={nao} ndet={ndet:.3e} C<={Cub} C2<={C2ub:.2e} "
          f"params<={params_ub:.2e} ratio>={ndet/params_ub:.1f}")
    checks = {
        "C <= 2^nao": all(r["C"] <= 2 ** r["nao"] for r in rows),
        "ndet = C(nao,nocc)^2": all(r["ndet"] == comb(r["nao"], r["nelec"] // 2) ** 2 for r in rows),
        "H2O 6-31G C=5635 (union)": any(r["name"] == "H2O 6-31G" and r["C"] == 5635 for r in rows),
        "params < ndet for the four enumerated": all(r["params"] < r["ndet"] for r in rows[:4]),
    }
    for k, v in checks.items():
        print(f"  [{'v' if v else 'X'}] {k}")
    out = {"rows": rows, "checks": checks}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n-> {OUT}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
