#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_uud.py —— 重新奠基: 质子=uud 三色荷, 偏心由电荷定(不等份), 非线性合成
======================================================================================
用户: 三个夸克不是三等份, 偏心不等份, 由 uud 性质(电荷 2/3,2/3,-1/3)定.
      用 Zero 的量(电荷/格点/闭合), 不用经典拟合. 非线性合成.
机制候选: 偏心角 θ = 3π/4(零点) + 偏移(由电荷×格点构建).
合成: 三色荷 M 的非线性组合. 找 S=m_p/A²=2.9895 (若共享标度 A).
"""
from __future__ import annotations
import os, sys, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

me, m_p = 0.51099895, 938.272088
theta_e = 2*np.pi/3 + 2/9
A = np.sqrt(me)/(1+np.sqrt(2)*np.cos(theta_e)); A2 = A**2
S0 = m_p/A2   # 2.989468
P_("=" * 96)
P_("质子=uud 偏心不等份+非线性 · k4_proton_uud.py")
P_("=" * 96)
P_("目标合成 S0=%.6f (m_p=A²·S0, 共享标度 A²=%.3f)" % (S0, A2))
def M(t): return 1 + np.sqrt(2)*np.cos(t)

Z = 3*np.pi/4                 # 零点(光子)
qu, qd = 2/3, -1/3             # 电荷
# 格点/内生角
pi=np.pi; d0=np.pi/12-2/9; two9=2/9
grid = [0, d0, 2*d0, two9, pi/12, pi/8, pi/6, pi/4, pi/3, pi/2,
        pi/12+d0, two9+d0, pi/6+d0, 2*two9, pi/8+d0, pi/4+d0]

# ── ① 偏心偏移候选: 由电荷 × 格点构建 ──
# u 电荷 2/3 相关偏移, d 电荷 1/3 相关偏移 (用 |q|)
def off_cands(qfrac):
    """由 |q| 与格点构建的偏移候选 (相对零点)."""
    a = abs(qfrac)  # u:2/3, d:1/3
    c = []
    for g in grid:
        c += [a*g, -a*g, a*g/two9*0, 0]
    c += [a*np.pi, (1-a)*np.pi, a*2*np.pi/3, (1-a)*2*np.pi/3]
    return [x for x in c if abs(x)<np.pi]

off_u = sorted(set(round(x,9) for x in off_cands(2/3)))
off_d = sorted(set(round(x,9) for x in off_cands(1/3)))
P_("  u 偏移候选数=%d, d 偏移候选数=%d" % (len(off_u), len(off_d)))

# ── ② 搜索: 三色荷 θ=(Z+ou, Z+ou, Z+od), 非线性合成 → S0 ──
best = []
synth_names = {
    "2M_u²+M_d²": lambda Mu,Md: 2*Mu**2+Md**2,
    "M_u²·M_d":   lambda Mu,Md: Mu**2*Md,
    "M_u·M_d":     lambda Mu,Md: Mu*Md,
    "(2M_u+M_d)²/某": None,
}
# 显式合成函数
def synths(Mu, Md):
    return {
        "2M_u²+M_d²": 2*Mu**2+Md**2,
        "M_u²·M_d":   Mu**2*Md,
        "M_u·Md":     Mu*Md,
        "ΣM²/2":      (2*Mu**2+Md**2)/2,
        "(ΣM)²/3":    (2*Mu+Md)**2/3,
        "2M_u+Md":    2*Mu+Md,
        "M_u·Mu·Md":  Mu*Mu*Md,
        "M_u²+Md":    Mu**2+Md,
    }
for ou in off_u:
    for od in off_d:
        Mu = M(Z+ou); Md = M(Z+od)
        for sname, fn in synths(Mu,Md).items():
            rel = abs(fn-S0)/S0*100
            best.append((rel, sname, ou, od, Mu, Md, fn))
best.sort(key=lambda x:x[0])
P_("\n【Top12 最接近 S0=%.4f (m_p=A²·S, 偏差)】" % S0)
for rel, sname, ou, od, Mu, Md, S in best[:12]:
    P_("  合成=%-10s S=%.5f (偏差 %.3f%%)  ou=%.4f° od=%.4f°  Mu=%.4f Md=%.4f" %
       (sname, S, rel, np.degrees(ou), np.degrees(od), Mu, Md))

# ── ③ 精确命中检查 (rel < 0.1%) ──
P_("\n【精确命中 (偏差<0.1%)】")
exact = [b for b in best if b[0]<0.1]
if exact:
    for b in exact[:10]:
        P_("  ★ %-10s S=%.6f  ou=%.5f° od=%.5f°" % (b[1], b[6], np.degrees(b[2]), np.degrees(b[3])))
else:
    P_("  无偏差<0.1% 的精确命中")

out={"s0":S0,"best":[b[0:7] for b in best[:12]]}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES,"k4_proton_uud.json"),"w",encoding="utf-8") as f_:
    json.dump(out,f_,ensure_ascii=False,indent=1)
with open(os.path.join(RES,"k4_proton_uud.txt"),"w",encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_uud.json / .txt")
