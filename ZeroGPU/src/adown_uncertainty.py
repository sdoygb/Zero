#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adown_uncertainty.py —— A_dn 到底被钉得多紧？两个效应必须分开报
================================================================================
背景（Z126）：
  质量律  m = u_s A^2 M^2，  M = 1 + A cos(theta)，三代相位相隔 2π/3。
  下型    A_dn = 2.682850，u_dn = 0.98798，最轻代相位 114.7387°。
  只给比值：整体标度全进 u_s，**A 只由质量比值定**。

两个**不同性质**的效应，本篇分开报：
  【一】分支退化（离散、精确）：同一组质量比值有多个精确解 (A, theta)，u_s 各自补偿。
  【二】输入敏感度（连续）：固定一条分支，把 m_d,m_s,m_b 在 PDG 误差内扰动，
        看 A 连续移动多少。

【检验】
  【V1】复现 Z126 的 A_dn、u_dn、theta
  【V2】比值由 A、theta 定，与 u_s 无关
  【V3】分支退化：列出全部精确分支
  【V4】各分支生成同一个质量比值
  【V5】各分支各配一个 u_s 后生成同一组质量
  【V6】哪些分支的 A <= 8/3（即不"超出上界"）
  【V7】连续敏感度：固定本支，PDG 误差内 A 的范围
  【V8】判据分辨力：0.032%／0.61%／1.20% 与敏感度的信噪比

用法：/usr/local/bin/python3 adown_uncertainty.py
输出：ZeroGPU/results/adown_uncertainty.json
"""
from __future__ import annotations
import json, math, os

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "adown_uncertainty.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


A2_MEV = 313.856
A_DN = 2.682850
U_DN = 0.98798
TH_DN = math.radians(114.7387)
MASS = {"d": 4.67, "s": 93.4, "b": 4180.0}

# PDG 量级的**非对称**相对误差（下型最不确定；m_b 最准）
ERR = {"d": (-0.036, +0.103), "s": (-0.036, +0.092), "b": (-0.005, +0.007)}


def M_of(theta, A):
    return 1.0 + A * math.cos(theta)


def masses(A, th, u):
    return [u * A2_MEV * M_of(th + 2 * math.pi * i / 3, A) ** 2 for i in range(3)]


def ratios(A, th):
    m = masses(A, th, 1.0)
    return np.array(m) / m[0]


def res(p, tr):
    A, th = p
    r = ratios(A, th)
    if np.any(r <= 0) or not np.all(np.isfinite(r)):
        return np.array([1e6, 1e6])
    return np.log(r) - np.log(tr)


def all_branches(tr, A_lo=1.0, A_hi=6.0, step=0.2):
    """全精确分支（残差 < 1e-14，按正 A 去重）。"""
    sols = []
    for A0 in np.arange(A_lo, A_hi + 1e-9, step):
        for th0 in np.arange(0.0, 2 * math.pi, step):
            try:
                r = least_squares(res, [A0, th0], args=(tr,),
                                  xtol=1e-15, ftol=1e-15, gtol=1e-15)
            except Exception:
                continue
            if float(np.sum(r.fun ** 2)) < 1e-14:
                A, th = r.x[0], r.x[1] % (2 * math.pi)
                if A < 0:                       # A→-A, θ→θ+π 等价
                    A, th = -A, (th + math.pi) % (2 * math.pi)
                ang = lambda x, y: abs((x - y + math.pi) % (2 * math.pi) - math.pi)
                if all(abs(A - a) > 1e-3 or ang(th, t) > 1e-3
                       for _, a, t in sols):
                    sols.append((float(np.sum(r.fun ** 2)), A, th))
    return sorted(sols)


print("=" * 94)
print("【V1】复现 Z126 的下型参数")
print("=" * 94)
m = masses(A_DN, TH_DN, U_DN)
print(f"  A_dn={A_DN}, u_dn={U_DN}, theta={math.degrees(TH_DN):.4f}°")
print(f"  ⟹ 质量 d,s,b = {m[0]:.3f}, {m[1]:.2f}, {m[2]:.1f} MeV")
chk(abs(m[0] - MASS["d"]) / MASS["d"] < 0.02, "复现 m_d", f"{m[0]:.3f} vs {MASS['d']}")
chk(abs(m[1] - MASS["s"]) / MASS["s"] < 0.03, "复现 m_s", f"{m[1]:.2f} vs {MASS['s']}")
chk(abs(m[2] - MASS["b"]) / MASS["b"] < 0.02, "复现 m_b", f"{m[2]:.1f} vs {MASS['b']}")

tr = np.array(m) / m[0]
print()
print("=" * 94)
print("【V2】比值只由 (A, theta) 定，与 u_s 无关")
print("=" * 94)
for u in (0.5, 1.0, 7.3):
    r = masses(A_DN, TH_DN, u)
    r = np.array(r) / r[0]
    print(f"  u={u:>4}: 比值 = [{r[0]:.6f}, {r[1]:.4f}, {r[2]:.3f}]")
r1 = ratios(A_DN, TH_DN)
chk(np.allclose(r1, tr), "比值与 u_s 无关", f"{r1[1]:.4f} vs {tr[1]:.4f}")

print()
print("=" * 94)
print("【V3】分支退化：同一组质量比值的全部精确解")
print("=" * 94)
brs = all_branches(tr)
print(f"  精确分支数（正 A，去重）：{len(brs)}")
print(f"  {'#':>2} {'A':>12} {'theta':>11} {'< 8/3 ?':>9}")
for i, (e, A, th) in enumerate(brs, 1):
    print(f"  {i:>2} {A:>12.6f} {math.degrees(th):>10.4f}° {'是' if A < 8/3 else '否':>9}")
chk(len(brs) >= 2, "存在多于一条精确分支", f"{len(brs)} 条")

print()
print("=" * 94)
print("【V4】【V5】各分支生成同一比值；各配 u_s 后生成同一质量")
print("=" * 94)
print(f"  {'A':>12} {'比值(相对第 3 代归一的差)':>26} {'u_s（使 m_b 吻合）':>20} {'m_d 误差':>10}")
for e, A, th in brs:
    rr = ratios(A, th)
    dev = float(np.max(np.abs(rr / tr - 1)))
    u = MASS["b"] / (A2_MEV * M_of(th + 4 * math.pi / 3, A) ** 2)
    mm = masses(A, th, u)
    print(f"  {A:>12.6f} {dev:>26.2e} {u:>20.6f} "
          f"{100*(mm[0]-MASS['d'])/MASS['d']:>+9.3f}%")
chk(all(float(np.max(np.abs(ratios(A, th) / tr - 1))) < 1e-6 for _, A, th in brs),
    "全部分支给出同一个比值")

print()
print("=" * 94)
print("【V6】哪些分支不超出 8/3")
print("=" * 94)
below = [A for _, A, _ in brs if A <= 8 / 3]
print(f"  A <= 8/3 的分支：{len(below)}/{len(brs)} 个  ⟹ {['%.6f' % x for x in below]}")
chk(len(below) >= 1, "至少有一条分支的 A 不超过 8/3（故'超出上界'是分支相关）",
    f"{len(below)} 条")

print()
print("=" * 94)
print("【V7】连续敏感度：固定本支（A=2.682850），PDG 误差内的 A 范围")
print("=" * 94)
print(f"  非对称误差：m_d {ERR['d']}, m_s {ERR['s']}, m_b {ERR['b']}")
lo, hi = [], []
for sd in (ERR["d"][0], 0.0, ERR["d"][1]):
    for ss in (ERR["s"][0], 0.0, ERR["s"][1]):
        for sb in (ERR["b"][0], 0.0, ERR["b"][1]):
            tgt = np.array([MASS["d"] * (1 + sd), MASS["s"] * (1 + ss),
                            MASS["b"] * (1 + sb)])
            tgt = tgt / tgt[0]
            r = least_squares(res, [A_DN, TH_DN], args=(tgt,),
                              xtol=1e-15, ftol=1e-15, gtol=1e-15)
            A = r.x[0]
            if abs(A - A_DN) < 0.5 * A_DN:        # 只留不换支的
                lo.append(A)
print(f"  留在本支的角点：{len(lo)}/27  （其余换支，已排除）")
lo = np.array(lo)
print(f"  A ∈ [{lo.min():.6f}, {lo.max():.6f}]")
rel = 100 * (lo.max() - lo.min()) / 2 / A_DN
print(f"  相对半宽 = ±{rel:.3f}%")
chk(rel > 1.0, "本支的连续敏感度 > 1%", f"±{rel:.3f}%")

print()
print("=" * 94)
print("【V8】判据分辨力")
print("=" * 94)
gaps = [("36/5 的 +0.032%", 0.032), ("8/3 的 -0.603%", 0.603), ("A² 尺的 1.20%", 1.20)]
print(f"  {'判据给出的偏差':>20} {'信噪比 = 偏差/±%.2f%%' % rel:>22} {'可分辨？':>10}")
for nm, g in gaps:
    snr = g / rel
    print(f"  {nm:>20} {snr:>22.3f} {'是' if snr > 3 else '否':>10}")
chk(all(g / rel < 3 for _, g in gaps), "三条判据的信噪比都 < 3")

out = {
    "A_dn": A_DN, "u_dn": U_DN, "theta_deg": math.degrees(TH_DN),
    "masses": m,
    "branches": [{"A": A, "theta_deg": math.degrees(th),
                  "below_8_3": bool(A <= 8 / 3)} for _, A, th in brs],
    "branch_count": len(brs),
    "sens_rel_half_pct": rel,
    "sens_range": [float(lo.min()), float(lo.max())],
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for x in bad:
    print("  ✗", x)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)