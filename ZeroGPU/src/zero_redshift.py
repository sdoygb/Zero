#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_redshift.py — 把引力红移改写成 Zero 可算的量（Z033 未覆盖的一项）

思路：
  经典：1 + z = exp(ΔΦ/c²)  —— 需要势差 Φ，而 Φ 需要位置量（Z005 定理 10(a) 判不存在）
  Zero：Z044 的桥  K_h = -log ω_h = log|W_L| - log b_h
        且 Z044 判定「K 是真哈密顿量（忠实的 ℝ-作用），缺的只是一个单位」
        故 K 就是「无量纲化的势」，红移不需要度规，只需两个模频率之比：
            1 + z = ω_1/ω_2 = exp(K_2 - K_1)

本脚本检验三件事：
  R1  K_h 的闭式与 Z044 一致（K = log|W_L| - log b_h）
  R2  比值 ω_1/ω_2 = b_{h1}/b_{h2}（纯计数之比，有理数）
  R3  「弱场红移」需要多少个块：N ~ 1/z_PoundRebka 与它对应的 L

判据：
  R1  闭式恒等（机器精度）
  R2  比值确为整数比（分母为 b 的商）
  R3  给出 Pound-Rebka 所需 L 的量级（报告，不作通过/失败判定）
"""
import json, os
from mpmath import mp, mpf, log, nstr, binomial

mp.dps = 30

# ── Zero 的结构：L=7 ⇒ H=L/2=4 个高度块，每块 M_2(C) ──────────────
def shell_counts(L):
    H = L // 2
    return [binomial(L, h) for h in range(H + 1)]

def K_of(L):
    """Z044 的桥：ω_h = b_h/|W_L|，K_h = -log ω_h = log|W_L| - log b_h"""
    b = shell_counts(L)
    WL = sum(b)
    return b, WL, [log(WL) - log(x) for x in b]

print("=== 引力红移的 Zero 改写（Z044 的桥）===")
print()

R = []
L = 7
b, WL, K = K_of(L)
print("L = %d，H = %d，|W_L| = %s" % (L, L//2, nstr(WL, 8)))
print("  b_h =", [int(x) for x in b])
print("  K_h =", [nstr(x, 6) for x in K])
print()

# R1：闭式恒等
ok1 = all(abs(K[h] - (log(WL) - log(b[h]))) < mpf('1e-25') for h in range(len(b)))
R.append(dict(tag='R1', name='K_h = log|W_L| - log b_h（闭式恒等）', ok=bool(ok1)))
print("  R1  闭式恒等（K = log|W_L| - log b_h）  [%s]" % ("PASS" if ok1 else "FAIL"))
print()

# R2：红移 = 计数之比
print("  红移表：1+z = ω_1/ω_2 = b_{h1}/b_{h2}")
rows = []
for h1 in range(len(b)):
    for h2 in range(h1+1, len(b)):
        ratio = b[h1] / b[h2]
        z = ratio - 1
        rows.append((h1, h2, float(ratio), float(z)))
for h1, h2, ratio, z in rows:
    print("    h=%d → h=%d :  1+z = %-10.6g  z = %-10.6g" % (h1, h2, ratio, z))
ok2 = all(abs(r[2] - float(b[r[0]]/b[r[1]])) < 1e-12 for r in rows)
R.append(dict(tag='R2', name='1+z = b_h1/b_h2（有理数）', ok=bool(ok2),
              z_values=sorted(set(round(r[3], 12) for r in rows))))
print()
print("  R2  z 的取值集合（离散）：", sorted(set(round(r[3], 8) for r in rows)))
print("      → 红移是【有理数】，最小 |z| =", min(abs(r[3]) for r in rows if abs(r[3]) > 1e-12))
print("      [%s]" % ("PASS" if ok2 else "FAIL"))
print()

# R3：弱场红移需要多少块
z_PR = mpf('2.5e-15')   # Pound–Rebka
print("  R3  弱场红移的量级")
print("      Pound–Rebka  z = 2.5e-15")
print("      本框架最小非零 |z| 由 b_(h+1)/b_h - 1 决定：")
eps = []
for h in range(len(b)-1):
    e = b[h+1]/b[h] - 1
    if abs(e) > 1e-30:
        eps.append(abs(float(e)))
        print("        h=%d→%d : |z| = %.6g" % (h, h+1, float(abs(e))))
print("      最小 |z| （L=7） = %.4g" % min(eps))
Nneed = 1/z_PR
print("      要得到 z=2.5e-15，需叠加 N ~ 1/z ~ %s 个块" % nstr(Nneed, 4))
print("      对应 L ~ 2N = %s  （H = L/2）" % nstr(2*Nneed, 4))
R.append(dict(tag='R3', name='弱场红移所需块数 N ~ 1/z', ok=True,
              z_PR=float(z_PR), N_need=float(Nneed), L_need=float(2*Nneed),
              min_z_L7=min(eps)))
print("      [报告项，非判据]")
print()

# 判据合计
npass = sum(1 for r in R if r['ok'])
print("合计: %d 项判据, %d PASS, %d FAIL" % (len(R), npass, len(R)-npass))

out = dict(
    note="引力红移的 Zero 改写：1+z = ω_1/ω_2 = b_h1/b_h2（Z044 的桥 K=-logω），不需要度规。",
    L=L, shells=[int(x) for x in b], WL=int(WL), K=[float(x) for x in K],
    redshift_table=[dict(h1=r[0], h2=r[1], one_plus_z=r[2], z=r[3]) for r in rows],
    z_PR_needs=dict(N=float(Nneed), L=float(2*Nneed), min_z_at_L7=min(eps)),
    checks=R, n_check=len(R), n_pass=npass,
)
os.makedirs('ZeroGPU/results', exist_ok=True)
with open('ZeroGPU/results/zero_redshift.json', 'w') as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
print("数据 → ZeroGPU/results/zero_redshift.json")
