#!/usr/bin/env python3
"""用 AME2020 实测核质量表检验 Zb042 半经验式的**形状**（系数由拟合定，见 Zb042 §5）

步骤：
 1) 从 YAML 读 (Z, N, A, mass_excess_keV)；
 2) 由质量过剩算结合能 B = Z·Δ_H + N·Δ_n − Δ(A,Z)；
 3) 只用偶-偶核（避开配对项 δ），拟合 Zb042 (4.1) 的四个系数；
 4) 看残差随 A、Z、(N−Z)/A 的结构。
"""
import numpy as np, re, json
from pathlib import Path

YAML = Path("/tmp/ame.yaml")
OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb042_ame_residuals.json")

# AME2020 的氢原子与中子质量过剩（keV）
D_H, D_n = 7288.971064, 8071.3181

def load():
    txt = YAML.read_text()
    rows = []
    for blk in txt.split("\n- nuclide_id:"):
        if "Z:" not in blk: continue
        def g(k, cast=float):
            m = re.search(rf"^\s+{k}:\s*(\S+)\s*$", blk, re.M)
            return None if m is None else cast(m.group(1))
        Z = g("Z", int); N = g("N", int); A = g("A", int); D = g("mass_excess_keV")
        if None in (Z, N, A, D): continue
        rows.append((Z, N, A, D))
    return np.array(rows, float)

def main():
    R = load()
    Z, N, A, D = R.T
    B = Z*D_H + N*D_n - D          # keV
    print(f"读到 {len(R)} 个核（AME2020 实测行）")
    ee = (Z.astype(int)%2==0) & (N.astype(int)%2==0) & (A>=40)   # 液滴模型只在中重核可靠
    print(f"其中偶-偶核 {ee.sum()} 个 —— 用它们避开配对项")
    Zc, Nc, Ac, Bc = Z[ee], N[ee], A[ee], B[ee]

    # 设计矩阵：[A, A^{2/3}, Z(Z-1)/A^{1/3}, (A-2Z)^2/A]
    X = np.column_stack([Ac, Ac**(2/3), Zc*(Zc-1)/Ac**(1/3), (Ac-2*Zc)**2/Ac])
    coef, *_ = np.linalg.lstsq(X, Bc, rcond=None)
    pred = X @ coef
    res = Bc - pred
    names = ["a_V", "a_S", "a_C", "a_sym"]
    print("\n拟合系数（Zb042 §5 说明：四个系数是外部拟合输入，库里只给形状）")
    for n,v in zip(names, coef): print(f"  {n:6s} = {v:9.3f} keV")
    print(f"\n残差：RMS = {np.sqrt(np.mean(res**2)):.2f} keV，中位 |res| = {np.median(np.abs(res)):.2f} keV，最大 {np.abs(res).max():.1f} keV")

    # 残差的结构：按 A 分箱、按 |N-Z| 分箱
    print("\n按 A 分箱（残差均值 ± 标准差）：")
    bins = [(40,60),(60,80),(80,100),(100,120),(120,150),(150,180),(180,210),(210,260)]
    for lo,hi in bins:
        m = (Ac>=lo)&(Ac<hi)
        if m.sum()>3:
            print(f"  A∈[{lo:3d},{hi:3d}): n={m.sum():4d}  均值 {res[m].mean():+7.2f}  σ {res[m].std():6.2f} keV")
    asym = (Nc-Zc)/Ac
    print("\n按 |(N−Z)/A| 分箱：")
    for lo,hi in [(0.0,0.05),(0.05,0.1),(0.1,0.15),(0.15,0.2),(0.2,0.3)]:
        m = (np.abs(asym)>=lo)&(np.abs(asym)<hi)
        if m.sum()>3:
            print(f"  |asym|∈[{lo:.2f},{hi:.2f}): n={m.sum():4d}  均值 {res[m].mean():+7.2f}  σ {res[m].std():6.2f} keV")
    # 幻数附近
    print("\n幻数附近（N 或 Z 落在 8,20,28,50,82,126）的残差：")
    magic = {8,20,28,50,82,126}
    mm = np.array([(int(n) in magic) or (int(z) in magic) for z,n in zip(Zc,Nc)])
    print(f"  幻数行 n={mm.sum()}  均值 {res[mm].mean():+7.2f}  σ {res[mm].std():6.2f} keV")
    print(f"  其它行 n={(~mm).sum()}  均值 {res[~mm].mean():+7.2f}  σ {res[~mm].std():6.2f} keV")
    # 拟合优度
    ss = 1 - np.sum(res**2)/np.sum((Bc-Bc.mean())**2)
    print(f"\nR² = {ss:.6f}")
    json.dump(dict(n_total=int(len(R)), n_even_even=int(ee.sum()),
                   coef={n: float(v) for n,v in zip(names, coef)},
                   rms=float(np.sqrt(np.mean(res**2))), r2=float(ss),
                   res_A={f"{lo}-{hi}": float(res[(Ac>=lo)&(Ac<hi)].mean()) for lo,hi in bins if ((Ac>=lo)&(Ac<hi)).sum()>3},
                   res_magic_mean=float(res[mm].mean()), res_other_mean=float(res[~mm].mean())),
              open(OUT,"w"), indent=2)
    print("落盘", OUT)

if __name__=="__main__":
    main()
