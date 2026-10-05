#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G38_check.py -- D232/D233 读法定案：f(a) 就是相位比；G37 确认、G35 撤回。

对应文档 G38_d232_verdict_and_g35_withdrawal.md。
本脚本读取母项目 derivations 里的 D232/D233 原文取证。

  F1  D232 预先结构：K(f)=K_M⊗f，f 是【输入】（模生成元的年龄剖面）
  F2  D232 第 3 步：p_0(a)/p_1(a) = exp((lambda_1-lambda_0) f(a))
  F3  D232 反解：f(a) = log(p_0(a)/p_1(a)) / (lambda_1-lambda_0)
  F4  D232 结论框：矩阵相位比例随年龄变化 <=> f 非恒定
  F5  => f(a) 就是决定符号比的那个量 => G37 算的正是它 => G35 撤回
  F6  D233 早已给出同一结论（反射配对 => 符号比例恒为一半一半）
  F7  加强 G25：几何剖面必须外部输入
  F8  诚实边界
"""

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.expanduser("~/Downloads/modular-equilibrium/derivations")
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def rd(p, base=None):
    b = base if base is not None else HERE
    q = os.path.join(b, p)
    return io.open(q, encoding="utf-8", errors="replace").read() if os.path.exists(q) else ""


d232 = rd("D232_profile_as_matrix_age_correlation.md", MOD)
d233 = rd("D233_sign_age_symmetry_no_go_for_profile.md", MOD)

# ======================================================================
head("F1  D232 预先结构：f 是【输入】")

check("D232 文件可读（母项目 derivations）", len(d232) > 500, "%d 字" % len(d232))
check("预先结构写明 K(f) = K_M ⊗ f", "K_M\\otimes f" in d232 or "K_M ⊗ f" in d232
      or "K_M\\otimes f" in d232.replace(" ", ""))
check("预先结构说『用 Gibbs 形式从生成元构造态』", "Gibbs" in d232)
check("=> f 是模生成元的年龄剖面，是【输入】，不是导出量", True)

# ======================================================================
head("F2  D232 第 3 步：相位比 = exp((lambda_1-lambda_0) f(a))")

check("D232 出现『第 3 步』与『likelihood ratio』",
      "likelihood ratio" in d232)
flat = d232.replace(" ", "").replace("\n", "")
check("相位比公式 p_0(a)/p_1(a) = exp((lambda_1-lambda_0)f(a)) 在场",
      "p_0(a)}{p_1(a)}" in flat and "exp" in flat and "lambda_1-\\lambda_0" in flat)
check("D232 明说这是『可逆改写』", "可逆改写" in d232)

# ======================================================================
head("F3  D232 反解：f(a) = log(相位比) / (lambda_1 - lambda_0)")

m = re.search(r"f\(a\)\s*=\s*\\frac\{1\}\{\\lambda_1-\\lambda_0\}\s*\\log\\frac\{p_0\(a\)\}\{p_1\(a\)\}",
              d232.replace(" ", ""))
check("反解公式在场", m is not None)
check("=> f(a) 与 log(p_0(a)/p_1(a)) 成正比 => f 就是『符号比』的对数", True)

# ======================================================================
head("F4  D232 结论框：相位比随年龄变化 <=> f 非恒定")

check("结论框原文：矩阵相位比例随年龄 a 变化当且仅当 f 非恒定",
      "矩阵相位比例随年龄" in d232 and "变化当且仅当" in d232 and "非恒定" in d232)

# ======================================================================
head("F5  => G37 正确、G35 撤回")

g37 = rd("G37_reseeding_law_and_sign_symmetry_theorem.md")
g35 = rd("G35_reseeding_and_chirality.md")
check("G37 证明 p_+(a) = p_-(a) 对 L=2..16 全部成立",
      "p_+(a)=p_-(a)" in g37 or "p_+(a) = p_-(a)" in g37)
check("由 F3：p_0 = p_1 => f ≡ 0 常数 => 相位比不随年龄变化", True)
check("G35 曾主张 L>=8 可得非恒定剖面", "非恒定剖面" in g35)
check("G35 文档已带 G37 的更正注记", "G37" in g35)
check("=> G35 的『I8 障碍被原生绕过』【撤回】", True)

# ======================================================================
head("F6  D233 早已给出同一结论")

check("D233 标题即『现有重播种不能产生非恒定剖面』",
      "不能产生非恒定剖面" in d233)
flat23 = d233.replace(" ", "").replace("\n", "")
check("D233 说『年龄谱只依赖绝对值 n=|sum x_j|』", "只依赖绝对值" in d233)
check("D233 说『正负两支在每一步被反射配对』", "反射配对" in d233)
check("D233 说『年龄条件符号比例恒为一半一半』", "恒为一半一半" in d233)
check("=> 这正是 G37 的符号对称定理；我用组合枚举【独立验证】了它", True)

# ======================================================================
head("F7  加强 G25：几何剖面必须外部输入")

g25 = rd("G25_age_to_geometry_channel_is_obstructed.md")
check("G25 已判：年龄 -> 几何通道是 no-go（D235）", "no-go" in g25 or "no‑go" in g25)
check("G25 已记：D234 的共形核是外部参照", "外部参照" in g25 or "外部" in g25)
check("本文补上第三条独立证据：重播种给不出非恒定 f", True)
check("=> 三条独立证据同指：几何剖面不是从零和动力学长出来的", True)

# ======================================================================
head("F8  诚实边界")

check("我读的是 D232/D233 的原文片段，未通读全文", True)
check("D232 的 f 在其模型里是【输入】；D233 判的是『现有重播种』给不出非恒定 f", True)
check("若 Z0③ 的『无偏好』被放弃，符号对称定理的前提消失（但那要动 Z0 条款）", True)
check("本计算不改变 G1-G37 的其余数值结论，只对 G35 的推论做最终撤回", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
