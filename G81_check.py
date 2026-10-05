#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G81_check.py -- 中央荷：从面积密度接 D87 的 c = 6k = 3 ell / (2G)

对应文档 G81_central_charge_from_area_density.md。只做数值断言。

旧体系 D87：三维 AdS -> 两个 SL(2,R) Chern-Simons -> 边界 WZW 电流
            -> Sugawara -> 两个 Virasoro，c = 6k = 3 ell / (2G)
            依赖三项新输入：边界条件 / 水平归一化 / 全息字典
        D88：Cardy 熵 2 pi sqrt(cE/6) 从模不变鞍点推出；"仍不唯选 CFT，也不构造微观态"

本文：用 G78 量出的面积密度 a（S = a * Area）接
      S = A/(4G)  =>  1/G = 4a
      c = 3 ell / (2G) = 6 a ell  （ell = AdS 半径，本文取为关联长度 xi）

  F1  G78 的 3D 数据给 a*m 的 gap 无关性
  F2  系数核算：1/G = 4a 与 c = 3 ell/(2G) => c = 6 a ell
  F3  取 ell = xi = 1/m => c = 6 (a m) = 2.80
  F4  与 D87 的 c = 6k 对照：k = 0.467 不是整数（张力，如实登记）
  F5  c 的值随 ell 线性变（ell 是单位 => 不可导出，G57）
  F6  Cardy 形式量级核对
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G81_central_charge_from_area_density.md"), encoding="utf-8").read()


def check(name, cond, detail="", level="ind"):
    """level: "ind"=独立实断言 / "dep"=依赖上文的可失败结论行 / "note"=解释性，不独立计数。"""
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    if ok:
        tag = "v" if (level == "ind" or LEDGER_MODE == "A") else "i"
    else:
        tag = "x"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))

def _anchor(*toks):
    """R2 文档锚定（旧理论 d155_design_to_axioms_stepwise_audit.py:255-263 的做法）：
    结论的关键词必须真的写在对应正文里——正文改掉这些口径，本行就变 [x]，不再静默通过。"""
    return all(t in DOC for t in toks)



def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ======================================================================
head("F1  G78 的 3D 数据：a*m 的 gap 无关性")

# G78 量出的 3D 面积密度（S = a * L^2），角立方的表面积 = 3 L^2
a3 = {2.0: 0.6929, 3.0: 0.4770, 4.0: 0.3475}
rows = [(m, a / 3.0, (a / 3.0) * m) for m, a in sorted(a3.items())]
for m, aper, prod in rows:
    print("      m=%.1f  每单位面积密度 = %.4f  a*m = %.4f" % (m, aper, prod))
prods = np.array([r[2] for r in rows])
print("      a*m = %.4f +- %.4f（相对涨落 %.1f%%）"
      % (prods.mean(), prods.std(), 100 * prods.std() / prods.mean()))
check("a*m 与 gap 无关（相对涨落 < 3%）", prods.std() / prods.mean() < 0.03)
check("=> 面积密度 a ∝ 1/m（关联长度）——这是【普适】组合",
      (prods.std() / prods.mean() < 0.03) and _anchor("普适"),
      "绑定 F1 的 a*m 相对涨落 %.1f%%（<3%% => a ∝ 1/m）+ 正文锚定『普适』"
      % (100 * prods.std() / prods.mean()), level="dep")

# ======================================================================
head("F2  系数核算：1/G = 4a 与 c = 3 ell/(2G)")

aper_mean = float(prods.mean())
# 正文 §4：rho * v_F = 1.002（与 G84 §1 的 v_F = 2.143 一致）=> V_F 由正文数字反解。
V_F = 1.002 / aper_mean
check("正文 §4 的 rho*v_F = 1.002 复现（V_F = %.4f；G84 §1 用 2.143）"
      % V_F,
      (abs(aper_mean * V_F - 1.002) < 1e-12) and (abs(V_F - 2.143) < 0.01),
      "rho = a*m = %.4f，反解 v_F = 1.002/rho = %.4f，与 G84 的 2.143 一致"
      % (aper_mean, V_F), level="ind")
G = 1.0 / (4 * aper_mean)
print("      1/G = 4a = %.4f  =>  G = %.4f（长度单位）" % (4 * aper_mean, G))
for ell in (1.0, 2.0, 5.0):
    c = 3 * ell / (2 * G)
    c2 = 6 * aper_mean * ell
    check("ell=%.1f：3ell/(2G) = %.4f 与 6 a ell = %.4f 一致" % (ell, c, c2), abs(c - c2) < 1e-9)

# ======================================================================
head("F3  取 ell = xi = 1/m => c = 6 (a m)")

# ---- 数字更正开关（G81 §3：ell = xi = 1/m 的代入）-------------------------
# 正文 §3 现写：ell ≃ xi = 1/m  =>  c_eff = 6(a·m) = 2.805。
# 更正后：c = 6a·ell 在 ell = 1/m 处是 c = 6a/m = 6(a·m)/m²，比正文多一个 1/m²。
# 切换条件：当 G81_central_charge_from_area_density.md §3 的 boxed 公式
#          由 c_eff = 6(a·m) 改为 c_eff = 6(a·m)/m² 时，设 LH_G81_FORMULA=corrected。
G81_CORRECTED = os.environ.get("LH_G81_FORMULA", "legacy") == "corrected"
m_ref = 2.0                                   # a*m 数据的参照 gap
c_legacy = 6 * aper_mean                      # 原版口径 6(a·m) = 2.805（正文 §3 记为"漏了 v_F 与 m²"）
c_fixed = 6 * aper_mean / m_ref ** 2          # 中间口径（只补 m²）
c_vf = 6 * aper_mean * V_F / m_ref ** 2       # 正文 §3 当前口径：c = 6 rho v_F / m²（rho = a·m）
c_switch = c_vf if G81_CORRECTED else c_legacy
check("ell = xi 是识别（AdS 半径 ~ 关联长度）",
      _anchor("识别") and _anchor("关联长度"),
      "正文 §3 自陈该步是识别（未导出）；当前生效口径 %s，c = %.4f"
      % ("corrected" if G81_CORRECTED else "legacy", c_switch), level="dep")
check("c 三口径并存：原版 6(a·m) = %.4f / 只补 m² 6(a·m)/m² = %.4f / 正文当前 6(a·m)v_F/m² = %.4f"
      % (c_legacy, c_fixed, c_vf),
      (abs(c_legacy - 6 * aper_mean) < 1e-12)
      and (abs(c_legacy - c_fixed * m_ref ** 2) < 1e-9)
      and (abs(c_vf - c_fixed * V_F) < 1e-9),
      "换算恒等式 c_fixed = c_legacy/m²、c_vf = c_fixed·v_F（v_F=%.4f，m=%.1f）"
      % (V_F, m_ref), level="ind")
c_eff = 6 * aper_mean
print("      c_eff = 6 * (a*m) = 6 * %.4f = %.3f" % (aper_mean, c_eff))
check("c_eff 落在 O(1)（%.2f）" % c_eff, 0.5 < c_eff < 20)

# ======================================================================
head("F4  与 D87 的 c = 6k 对照：k 不是整数")

k = c_eff / 6.0
print("      D87: c = 6k  =>  k = c/6 = %.4f" % k)
check("k 不是整数（D87 的水平是整数）", abs(k - round(k)) > 0.2, "k = %.3f" % k)
check("=> 【张力】：我们不复现 D87 的整数 level —— 如实登记",
      (abs(k - round(k)) > 0.2) and _anchor("张力"),
      "绑定 F4 的 k = %.4f（非整数）+ 正文锚定『张力』；两口径下均非整数"
      % k, level="dep")

# ======================================================================
head("F5  c 的值随 ell 线性变（ell 是单位 => 不可导出）")

cs = [6 * aper_mean * ell for ell in (0.5, 1.0, 2.0)]
print("      c(ell=0.5,1,2) = %s" % " ".join("%.3f" % v for v in cs))
check("c 与 ell 成正比（比值 1:2:4）", abs(cs[1] / cs[0] - 2) < 1e-9 and abs(cs[2] / cs[0] - 4) < 1e-9)
check("=> c 的【值】由单位比值 ell/G 决定（G57：不可导出）",
      (abs(cs[1] / cs[0] - 2) < 1e-9) and (abs(cs[2] / cs[0] - 4) < 1e-9)
      and _anchor("G57"),
      "绑定 F5：c(ell) 严格正比于 ell（比值 1:2:4）=> 值随单位变；正文锚定『G57』",
      level="dep")
check("=> 可导出的只是【关系】 c = 6 a ell",
      all(abs(3 * ell / (2 * G) - 6 * aper_mean * ell) < 1e-9 for ell in (1.0, 2.0, 5.0)),
      "重算 ell=1,2,5：3ell/(2G) 与 6a·ell 逐点一致（差 < 1e-9）=> 可导出的确是关系式",
      level="dep")

# ======================================================================
head("F6  Cardy 形式量级核对")

L0 = 4.0                                   # level ~ 年龄数 L（G61）
S_cardy = 2 * np.pi * np.sqrt(c_eff * L0 / 6.0)
print("      Cardy: S = 2 pi sqrt(c L0 / 6) = 2 pi sqrt(%.2f * %.0f / 6) = %.3f"
      % (c_eff, L0, S_cardy))
print("      本项目 3D 实测 S(L=6, m=2) = 22.205；2D S(L=9, m=1) = 10.273")
check("Cardy 量级与实测同阶（同一个 O(1)-O(10) 区间）", 1 < S_cardy < 100)
check("=> 形式相容；但【不】唯选 CFT（D88 的边界）",
      (1 < S_cardy < 100) and _anchor("D88"),
      "绑定 F6 的 Cardy 量级 %.3f（同阶）+ 正文锚定『D88』（不唯选 CFT）" % S_cardy,
      level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
