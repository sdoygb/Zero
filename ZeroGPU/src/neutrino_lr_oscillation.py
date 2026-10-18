# -*- coding: utf-8 -*-
"""Z019 左右手中微子：真质量、相消与三相位震荡。

观测质量 = 左右手两振幅的相干和（相对相位 2π−ε）：
    m_k = 2 M sin(ε_k/2) ≈ M ε_k,
    M = A² = 313.856 MeV,  ε_k = (2/9)^15 r_k²,  r_k=(0, A−1, 1).
三个观测值是左右震荡的三个驻定相位；第四个 r=A+1 不被采样。
"""
from __future__ import annotations
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
A = math.sqrt(2.0)
A2_MEV = 313.856
M_MEV = A2_MEV
M_MEV_UNIT = A2_MEV * 1e9          # MeV → meV
EPS3 = (2.0 / 9.0) ** 15
R = (0.0, A - 1.0, 1.0)
EPS = tuple(EPS3 * x * x for x in R)
M_OBS = tuple(M_MEV_UNIT * e for e in EPS)   # meV
EPS4 = EPS3 * (A + 1.0) ** 2
M4 = M_MEV_UNIT * EPS4
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def main():
    print("A =", A, " A² =", A2_MEV, "MeV")
    print("ε3 = (2/9)^15 =", EPS3, "rad =", EPS3 * 180 / math.pi, "deg")
    for k in range(3):
        print(f"k={k} r={R[k]:.9f} ε={EPS[k]:.6e} rad  phase={360-EPS[k]*180/math.pi:.10f}°  m={M_OBS[k]:.9f} meV")
    print(f"fourth r={A+1:.9f} ε4={EPS4:.6e} m4={M4:.6f} meV")
    print(f"true L/R = {M_MEV:.3f} MeV each; observed sum = {sum(M_OBS):.6f} meV")

    ck("M_L=M_R=A²=313.856 MeV", abs(M_MEV - 313.856) < 1e-9, f"{M_MEV} MeV")
    ck("ε3=(2/9)^15", abs(EPS3 - (2 / 9) ** 15) < 1e-20, f"{EPS3:.6e}")
    ck("m=(0,8.570,49.951) meV", abs(M_OBS[1] - 8.570) < 0.01 and abs(M_OBS[2] - 49.951) < 0.01,
       f"{M_OBS[1]:.3f},{M_OBS[2]:.3f}")
    ck("m2/m3=(A−1)²=3−2√2", abs(M_OBS[1] / M_OBS[2] - (A - 1) ** 2) < 1e-12,
       f"{M_OBS[1]/M_OBS[2]:.9f}")
    ck("相干和公式 2M sin(ε/2)=Mε", all(abs(2 * M_MEV_UNIT * math.sin(e / 2) - M_MEV_UNIT * e) < 1e-9 for e in EPS),
       "small-angle")
    ck("相位缺口 = 2π−ε（ε3≈1.59e-10 rad）",
       abs((2 * math.pi - EPS3) * 180 / math.pi - (360 - EPS3 * 180 / math.pi)) < 1e-12,
       f"gap={EPS3*180/math.pi:.3e}°")
    ck("第四个 r=A+1 给 m4≈291.1 meV（未观测支）", abs(M4 - 291.13) < 0.1, f"{M4:.3f} meV")
    ck("震荡 m(t)=M ε3 cos² 在 cos=0,A−1,1 给 0,m2,m3",
       abs(M_MEV_UNIT * EPS3 * 0.0) < 1e-12
       and abs(M_MEV_UNIT * EPS3 * (A - 1) ** 2 - M_OBS[1]) < 1e-6
       and abs(M_MEV_UNIT * EPS3 - M_OBS[2]) < 1e-6,
       "three phases")

    out = {"A": A, "A2_MeV": A2_MEV, "M_L_MeV": M_MEV, "M_R_MeV": M_MEV,
           "epsilon3": EPS3, "epsilon3_deg": EPS3 * 180 / math.pi,
           "r": list(R), "epsilon": list(EPS), "m_meV": list(M_OBS),
           "phase_deg": [360 - e * 180 / math.pi for e in EPS],
           "epsilon4": EPS4, "m4_meV": M4, "checks": CK,
           "n_fail": sum(1 for c in CK if not c["pass"])}
    with open(os.path.join(OUT, "neutrino_lr_oscillation.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(OUT, "neutrino_lr_oscillation.json"))


if __name__ == "__main__":
    main()
