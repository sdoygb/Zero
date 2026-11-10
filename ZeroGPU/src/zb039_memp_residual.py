#!/usr/bin/env python3
"""Zb039 推论 6.5 的 m_e/m_p 预言与实测的残差，以及三条修法的代价。"""
import math, json
th = math.pi/12 - 2/9
Me = 1 + math.sqrt(2)*math.cos(2*math.pi/3 + 2/9)
F0 = 1 - math.sqrt(5)*th*th
meas = 1/1836.15267343
pred = Me*Me/(3*F0)
r = pred/meas - 1
c = (F0*r)/th**4
out = dict(theta_kappa=th, M_e=Me, F0=F0, pred=pred, meas=meas, rel_residual=r,
           c_required=c, minus_1p5sqrt5=-1.5*math.sqrt(5),
           mp_residual_before=2.0e-6, mp_residual_after=2.0e-6 + (F0*r)/F0)
print(json.dumps(out, indent=2))
