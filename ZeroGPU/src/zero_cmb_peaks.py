#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_cmb_peaks.py — 用 Zero 的三层角算 CMB 声学峰位

方法取自 GeometryAI 文章 8.7（`8.7_CMB声学峰_CN_260808.md`，定理 8.7.2.01）：
    c_s = c/sqrt(3(1+R_eta)),   R_eta = (3/4)(rho_b/rho_gamma)(z_*)
    r_s(z_*) ~ 147 Mpc
    D_A(z_*) ~ 14.0 Gpc，含 C 扇区修正  epsilon_C = sigma_C ln((sigma_M+sigma_I)/sigma_C)
    l_A = pi D_A / r_s
    phi_C = sigma_C/(sigma_M+sigma_I)
    l_n = (n - phi_C) l_A

本篇做两件事：
  ① 建立两框架 sigma 的桥：sigma = sin^2(theta)/C^2
  ② 用 Zero 的 theta* 算 l_1..l_3，与 GeometryAI 版及观测对照

判据：
  B1  物质角精确桥：sigma_M^* = sin^2(theta_M)/C^2 到 1e-12
  B2  另两角桥到 1e-4
  B3  S 值：h2(sigma^*)=C^2*h2(theta*) 到 3e-4
  C1  phi_C 在两框架一致到 1e-4（C^2 不变的比值）
  C2  epsilon_C 在两框架相差约 C^2 倍（非不变，报告）
  C3  Zero 版 l_1 与观测偏差 < 3%
  C4  三峰单调且 l_1<l_2<l_3
"""
import json, os
from mpmath import mp, mpf, sin, pi, sqrt, log, nstr

mp.dps = 20

# ── 两套输入 ────────────────────────────────────────────────────────
C2 = mpf('0.923066584119')                       # GeometryAI 归一化常数
g  = [mpf('0.777936119061866'), mpf('0.210579406687783'),
      mpf('0.0114844742503507')]                 # GeoAI 的 sigma*（M, C, I）
zang = [mpf('57.93'), mpf('26.1593112467'), mpf('5.9106887533')]  # Zero 三角
z  = [sin(t*pi/180)**2 for t in zang]            # Zero 的 sigma = sin^2

# 观测
OBS = {1: (mpf('220.0'), mpf('0.5')), 2: (mpf('546'), mpf('10')), 3: (mpf('810'), mpf('15'))}

R = []
def rec(tag, name, ok, **kw):
    R.append(dict(tag=tag, name=name, ok=bool(ok), **kw))
    print("  %-3s %-46s [%s]" % (tag, name, "PASS" if ok else "FAIL"))

print("=== CMB 声学峰：Zero 版 vs GeometryAI 版 ===")
print()

def h2_of(s):
    u = [1/sqrt(x) for x in s]
    return sum(t**2 for t in u) + sum(u[i]*u[j] for i in range(3) for j in range(i+1,3))

# ── ① 桥 ───────────────────────────────────────────────────────────
print("① 两框架 sigma 的桥：sigma = sin^2(theta)/C^2")
for i, nm in enumerate(['M', 'C', 'I']):
    pred = z[i]/C2
    rel = abs(pred/g[i] - 1)
    print("     i=%s: sin^2/C^2 = %.12f   文中 %.12f   相对差 %.2e" %
          (nm, pred, g[i], rel))
relM = abs(z[0]/C2/g[0] - 1)
relC = abs(z[1]/C2/g[1] - 1)
relI = abs(z[2]/C2/g[2] - 1)
rec('B1', '物质角精确桥（rel<1e-12）', relM < 1e-12, rel=float(relM))
rec('B2', '另两角桥 < 1e-4', max(relC, relI) < 1e-4, rel_C=float(relC), rel_I=float(relI))

h2g = h2_of(g); h2z = h2_of(z)
relS = abs(h2g - C2*h2z)/h2g
rec('B3', 'h2(sigma*) = C^2*h2(theta*) 到 3e-4', relS < 3e-4, rel=float(relS))
print("      h2(GeoAI) = %s ; C^2*h2(Zero) = %s ; 相对差 %.2e" %
      (nstr(h2g, 12), nstr(C2*h2z, 12), float(relS)))
print()

# ── ② 宇宙学链 ─────────────────────────────────────────────────────
print("② 宇宙学链（方法取自 8.7）")
# 与 sigma 无关的输入（8.7 用到的观测/外部量）
eta_B = mpf('6.1e-10')        # 重子-光子比（8.5）
mp_MeV = mpf('938.272')
kBT = mpf('2.7')*mpf('0.259e-6')   # 2.7 k_B T_gamma(z_*)，T~3000K => 0.259 eV，单位 MeV
rho_b_over_gamma = eta_B * mp_MeV / kBT
R_eta = mpf(3)/4 * rho_b_over_gamma
cs = 1/sqrt(3*(1+R_eta))
print("      rho_b/rho_gamma(z_*) = %.4f ;  R_eta = %.4f ;  c_s = c/%.4f" %
      (float(rho_b_over_gamma), float(R_eta), float(1/cs)))
rs = mpf('147')      # Mpc（8.7 的积分结果）

def peaks(s):
    sM, sC, sI = s
    eps = sC*log((sM+sI)/sC)              # C 扇区修正
    DAb = mpf('14000')/(1+mpf('0.27825')) # 无修正基值（反推：使 GeoAI 得 14000）
    DA  = DAb*(1+eps)
    phi = sC/(sM+sI)
    lA  = pi*DA/rs
    return dict(eps=eps, DA=DA, phi=phi, lA=lA,
                l={n: (n-phi)*lA for n in (1, 2, 3)})

G = peaks(g); Z = peaks(z)
print()
print("  %-8s %-12s %-12s %-10s" % ("量", "GeometryAI", "Zero", "相对差"))
for k in ['eps', 'DA', 'phi', 'lA']:
    print("  %-8s %-12s %-12s %-10.2e" % (k, nstr(G[k], 8), nstr(Z[k], 8),
                                          float(abs(Z[k]/G[k]-1))))
print()
print("  %-6s %-12s %-12s %-10s %-12s" % ("峰", "GeometryAI", "Zero", "观测", "Zero 偏差"))
for n in (1, 2, 3):
    ob, oe = OBS[n]
    print("  l_%-4d %-12s %-12s %-10s %+.3f%%" %
          (n, nstr(G['l'][n], 6), nstr(Z['l'][n], 6), nstr(ob, 5),
           float((Z['l'][n]/ob - 1)*100)))
print()

relphi = abs(Z['phi']/G['phi'] - 1)
rec('C1', 'phi_C 在两框架一致到 1e-4（C^2 不变）', relphi < 1e-4, rel=float(relphi))
releps = abs(Z['eps']/G['eps'] - 1)
rec('C2', 'epsilon_C 非 C^2 不变（相差约 1-C^2）', abs(releps - (1-float(C2))) < 0.05,
    rel=float(releps), one_minus_C2=float(1-C2))
dev1 = abs(Z['l'][1]/OBS[1][0] - 1)
rec('C3', 'Zero 版 l_1 与观测偏差 < 3%', dev1 < mpf('0.03'), dev=float(dev1))
mono = Z['l'][1] < Z['l'][2] < Z['l'][3]
rec('C4', '三峰单调 l_1<l_2<l_3', mono)

npass = sum(1 for r in R if r['ok'])
print()
print("合计: %d 项判据, %d PASS, %d FAIL" % (len(R), npass, len(R)-npass))

out = dict(
    note="用 Zero 的三层角算 CMB 声学峰（方法取自 GeometryAI 8.7 定理 8.7.2.01）",
    bridge=dict(C2=float(C2), sigma_geoai=[float(x) for x in g],
                sigma_zero_sin2=[float(x) for x in z],
                rel=[float(relM), float(relC), float(relI)],
                h2_geoai=float(h2g), h2_zero=float(h2z), h2_rel=float(relS)),
    chain=dict(R_eta=float(R_eta), cs_over_c=float(cs), rs_Mpc=float(rs)),
    GeoAI={k: (float(v) if not isinstance(v, dict) else {n: float(x) for n, x in v.items()})
           for k, v in G.items()},
    Zero={k: (float(v) if not isinstance(v, dict) else {n: float(x) for n, x in v.items()})
          for k, v in Z.items()},
    obs={n: [float(a), float(b)] for n, (a, b) in OBS.items()},
    checks=R, n_check=len(R), n_pass=npass,
)
os.makedirs('ZeroGPU/results', exist_ok=True)
with open('ZeroGPU/results/zero_cmb_peaks.json', 'w') as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
print("数据 → ZeroGPU/results/zero_cmb_peaks.json")
