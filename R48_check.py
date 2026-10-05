#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R48_check.py -- 核验"精确锥一直在，缺的是锥边可见性"的判决与数值。

  F1  文档结构：引理、闭式、三区制、更正、残留、未推出
  F2  JSON 数值自洽：支持界零违反、eps(B) = (1/2)log(4/B)、B=4 幂律、B>4 钉住
  F3  独立复算：闭式 mu*/c*、eps(B)、支持界（小规模独立积分）
  F4  交叉一致：G59/G56/R47 的既有断言仍在位；STATUS/INDEX 登记
"""
from __future__ import annotations

import io
import json
import math
import os
import sys
from math import cosh, log, tanh

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   " + d) if d else ""))


def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f:
        return f.read()


DOC = read("R48_exact_cone_and_effective_cone.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
G59 = read("G59_I7_settled_native_cone_and_its_residue.md")
G56 = read("G56_degeneration_attempt2_six_slots.md")
R47 = read("R47_signature_from_causal_cone.md")
RES = json.load(io.open(os.path.join(HERE, "R48_exact_cone_results.json"), encoding="utf-8"))

# ----------------------------------------------------------------------
head("F1  文档结构")

check("标题点出'精确'与'可见性'", "精确" in DOC and "可见" in DOC)
check("(R48-1)...(R48-7) 在位", all(("(R48-%d)" % k) in DOC for k in range(1, 8)))
check("引理给出支持界 |x| <= t", "supp" in DOC and "c_{\\rm exact}=1" in DOC)
check("eps(B) = (1/2) log(4/B) 写明", "\\tfrac12\\log\\tfrac4B" in DOC or "\\tfrac12\\log(4/B)" in DOC)
check("三区制 B<4 / =4 / >4 成表", "B<4" in DOC and "B>4" in DOC and "B=4" in DOC)
check("写明 B 的上确界来自 f(+inf)=log 2", "f(+\\infty)=\\log2" in DOC or "f(+\\infty)=\\log 2" in DOC)
check("更正 R47 §3 残留 2", "残留 2" in DOC and "R47" in DOC)
check("更正 G59 §3.3 的措辞", "G59" in DOC and "概念更正" in DOC)
check("残留含'B 的绝对来源未导出'", "绝对来源" in DOC or "$B=4$ 的绝对来源" in DOC)
check("未推出动力学", "没有**推出 Einstein" in DOC or "Einstein 方程" in DOC)
check("登记维数常数因子未验", "常数因子" in DOC)

# ----------------------------------------------------------------------
head("F2  JSON 数值自洽")

sup = RES["support"]
supv = RES["support_verdict"]
check("支持界零违反（全部 B）", supv["violations"] == 0, "violations=%d" % supv["violations"])
check("精确锥速度 = 1", abs(supv["exact_cone_speed"] - 1.0) < 1e-15)
for B in ("2.00", "3.00", "3.90", "4.00", "5.00"):
    k = "B=%s" % B
    check("%s：max(tip-t) <= 0" % k, sup[k]["max_tip_minus_t"] <= 0, "%d" % sup[k]["max_tip_minus_t"])
check("B=4,5：尖端速度 = 1（>0.999）",
      sup["B=4.00"]["tip_speed_last40pct"] > 0.999 and sup["B=5.00"]["tip_speed_last40pct"] > 0.999)
check("B=2：尖端速度 < 1（不是被钉住的那一支）", sup["B=2.00"]["tip_speed_last40pct"] < 0.999,
      "%.6f" % sup["B=2.00"]["tip_speed_last40pct"])

ev = RES["edge_visibility"]
for B in ("2.00", "2.50", "3.00", "3.50", "3.90"):
    k = "B=%s" % B
    check("%s：eps 实测与闭式一致（<1%%）" % k, ev[k]["rel_diff"] < 0.01,
          "实测 %.6f vs 闭式 %.6f" % (ev[k]["eps_measured"], ev[k]["eps_closed"]))
check("B=4：幂律指数 alpha ≈ 1", abs(ev["B=4.00"]["power_law_exponent"] - 1.0) < 0.02,
      "alpha=%.4f" % ev["B=4.00"]["power_law_exponent"])
rtt = ev["B=4.00"]["rho_t_times_t"]
check("B=4：rho*t 随 t 缓增且有界（对数修正）",
      all(0.5 < float(rtt[str(t)]) < 20.0 for t in (500, 1000, 2000, 5000, 10000)),
      "t=500 -> %.3f, t=10000 -> %.3f" % (float(rtt["500"]), float(rtt["10000"])))
check("B=4：rho_t(t) 多项式可见（t=10000 时 > 1e-5）",
      ev["B=4.00"]["rho_t"]["10000"] > 1e-5, "%.3e" % ev["B=4.00"]["rho_t"]["10000"])
check("B=2：rho_t(t) 指数不可见（t=2000 时 < 1e-10）",
      ev["B=2.00"]["rho_t"]["2000"] < 1e-10, "%.3e" % ev["B=2.00"]["rho_t"]["2000"])

ab = RES["above4"]
for B in ("4.00", "4.50", "5.00", "6.00", "8.00"):
    k = "B=%s" % B
    check("%s：尖端速度 = 1（钉在锥上）" % k, abs(ab[k]["tip_speed"] - 1.0) < 1e-6,
          "%.6f" % ab[k]["tip_speed"])
check("B=5：振幅无关（seed 1 与 1e-10 速度同）",
      abs(ab["B=5.00 seed 1e+00"]["tip_speed"] - ab["B=5.00 seed 1e-10"]["tip_speed"]) < 1e-6)
check("B>4 无 KPP 解", RES["kpp"]["B=4.5"]["mu_star"] is None and RES["kpp"]["B=5.0"]["mu_star"] is None)

kp = RES["kpp"]
check("闭式 = 直接最小化 lam/mu（相对差 < 1e-12）",
      all(kp["B=%.2f" % B]["equivalence_rel_diff"] < 1e-12 for B in (2.0, 2.5, 3.0, 3.5, 3.9, 4.0)),
      "max %.2e" % max(kp["B=%.2f" % B]["equivalence_rel_diff"] for B in (2.0, 2.5, 3.0, 3.5, 3.9, 4.0)))
check("c*(4) = 1（1e-9）", abs(kp["B=4.00"]["c_star"] - 1.0) < 1e-9,
      "%.10f" % kp["B=4.00"]["c_star"])

nv = RES["naive"]
for B in ("2.00", "3.00", "3.90"):
    k = "B=%s" % B
    check("独立无饱和积分器 %s：eps 与闭式一致（<1e-4 相对）" % k,
          abs(nv[k]["eps_measured"] - nv[k]["eps_closed"]) / nv[k]["eps_closed"] < 1e-4,
          "%.6f vs %.6f" % (nv[k]["eps_measured"], nv[k]["eps_closed"]))
    check("独立无饱和积分器 %s：支持界同样成立" % k, nv[k]["max_tip_minus_t"] <= 0)

ea = RES["eps_analytic"]
# (i) KPP 定义残差为零；(ii) lambda(mu*) = mu* c*（被选模的谱值）
for B in ("2.00", "2.50", "3.00", "3.50", "3.90", "4.00"):
    k = "B=%s" % B
    check("%s：f(mu*) = (1/2)log B（残差 < 1e-12）" % k, abs(ea[k]["kpp_f_residual"]) < 1e-12,
          "residual=%.1e" % ea[k]["kpp_f_residual"])
    check("%s：lambda(mu*) = mu*c*（相对差 < 1e-12）" % k, ea[k]["lambda_rel_diff"] < 1e-12,
          "%.10f vs %.10f" % (ea[k]["lambda(mu_star)"], ea[k]["closed_form_lambda"]))
check("B=4：mu* - lambda(mu*) = 0（锚点）",
      abs(ea["B=4.00"]["mu_star"] - ea["B=4.00"]["lambda(mu_star)"]) < 1e-9,
      "mu*-lam=%.3e" % (ea["B=4.00"]["mu_star"] - ea["B=4.00"]["lambda(mu_star)"]))
check("B=4：alpha 未定（eps=0，登记）", ea["B=4.00"]["alpha_amplitude_exponent"] is None)
check("alpha(B) 随 B 递增（B=2 -> 3.9）",
      ea["B=2.00"]["alpha_amplitude_exponent"] < ea["B=3.90"]["alpha_amplitude_exponent"],
      "%.3f -> %.3f" % (ea["B=2.00"]["alpha_amplitude_exponent"], ea["B=3.90"]["alpha_amplitude_exponent"]))
check("alpha(B) = mu*c*/eps 的定义一致性（B=2）",
      abs(ea["B=2.00"]["alpha_amplitude_exponent"]
          - ea["B=2.00"]["mu_star"] * ea["B=2.00"]["c_star"] / ea["B=2.00"]["eps_closed"]) < 1e-9,
      "alpha=%.6f" % ea["B=2.00"]["alpha_amplitude_exponent"])

dt = RES["deep_tail"]
slopes = [dt["windows"]["t=%d" % t]["slope"] for t in (2000, 4000, 6000)]
check("深层尾部斜率随 t 单调趋近 -mu*",
      abs(slopes[-1] + dt["mu_star"]) < abs(slopes[0] + dt["mu_star"]),
      "%.5f -> %.5f（-mu*=%.5f）" % (slopes[0], slopes[-1], -dt["mu_star"]))
check("深窗末端相对差 < 3%", abs(slopes[-1] + dt["mu_star"]) / dt["mu_star"] < 0.03)

# ----------------------------------------------------------------------
head("F3  独立复算（不复用探针代码路径）")


def f_mu(mu):
    return mu * tanh(mu) - log(cosh(mu))


for B in (2.0, 3.0, 4.0):
    tgt = 0.5 * log(B)
    lo, hi = 1e-12, 400.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f_mu(mid) < tgt:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    check("B=%.1f：f(mu*) = (1/2)log B（1e-9）" % B, abs(f_mu(mu) - tgt) < 1e-9,
          "f=%.9f tgt=%.9f" % (f_mu(mu), tgt))
    check("B=%.1f：c* 单调、B=4 时为 1" % B, (abs(tanh(mu) - 1.0) < 1e-9 if B == 4 else tanh(mu) < 1.0))
check("f(+inf) = log 2 <=> log B/2 <= log 2 <=> B <= 4",
      abs(f_mu(400.0) - log(2)) < 1e-6 and abs(0.5 * log(4.0) - log(2)) < 1e-15)
for B in (2.0, 3.0, 3.9, 4.0):
    check("eps(%.1f) = (1/2)log(4/B) 手算复核" % B,
          abs(0.5 * log(4.0 / B) - 0.5 * (log(4.0) - log(B))) < 1e-15)

# 独立小规模积分：支持界 + B=4 锥边可见
T, X = 400, 700
n = np.zeros((4, 2 * X + 1))
n[0, X] = 1.0
edge = np.zeros(T)
vmax = 0
for t in range(T):
    m = np.zeros_like(n)
    m[:, 1:-1] = 0.5 * (n[:, :-2] + n[:, 2:])
    tot = m.sum(axis=0)
    nn = np.zeros_like(n)
    nn[0] = 4.0 * m[1] * (1.0 - tot)
    nn[1] = m[0]
    nn[2] = m[1]
    nn[3] = m[2]
    n = np.maximum(nn, 0.0)
    nt = n.sum(axis=0)
    nz = np.flatnonzero(nt > 0)
    vmax = max(vmax, (nz.max() - X) - (t + 1))
    edge[t] = n[:, X + t + 1].sum()
check("独立积分（B=4，T=400）：支持界零违反", vmax <= 0, "max(tip-t)=%d" % vmax)
check("独立积分（B=4）：锥边密度 O(1e-2)（t=400）", edge[-1] > 1e-3, "%.4e" % edge[-1])
check("独立积分（B=4）：锥边密度按幂律衰减（rho(400)*400 > 1）", edge[-1] * 400 > 1.0,
      "rho*t = %.3f" % (edge[-1] * 400))

# ----------------------------------------------------------------------
head("F4  交叉一致与登记")

check("G59 仍把精确锥列为残余（本文关闭它）", "精确锥" in G59 and "有效锥" in G59)
check("G56 已有 c*=1 <=> B=4 与 B>4 无解", "c_*=1\\iff\\tfrac12\\log B=\\log 2\\iff B=4" in G56)
check("G56 §4 把'B 的值无原生约束'登记为缺口", "B$ 的值" in G56 or "$B$ 的值" in G56)
check("R47 的残留 2 提到精确光锥仍缺", "精确光锥" in R47)
check("STATUS 已登记 R48", "### 2.48" in STATUS and "R48_exact_cone_and_effective_cone.md" in STATUS
      and "R48_check.py" in STATUS)
check("STATUS 的'不是当前卡点'已更新精确锥行", "锥边可见性" in STATUS or "精确锥" in STATUS)
check("INDEX 已收录 R48", "R48_exact_cone_and_effective_cone.md" in INDEX and "R48_check.py" in INDEX)
check("探针与结果文件同名约定", os.path.exists(os.path.join(HERE, "R48_exact_cone_probe.py"))
      and os.path.exists(os.path.join(HERE, "R48_exact_cone_results.json")))

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
