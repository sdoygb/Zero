#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R8_check.py -- Jacobson 2016 补前提尝试的独立核验。

对应文档 R8_jacobson_entanglement_equilibrium_completion.md。

  F1  外部论文与首选目标在位
  F2  条件定理 C1-C4 与失败树 L1-L4 在位
  F3  任意忠实态精确熵差恒等式 S(sigma)-S(rho)=Tr(Delta K)-D(sigma||rho)
  F4  一阶第一定律的 O(eps^2) 余项
  F5  Einstein 系数 eta=1/(4G) 的代数消去
  F6  当前 gap/面积律硬障碍被明确登记
  F7  文档没有越过“条件恢复”
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
NCHECK = 0


def check(name, cond, detail=""):
    global FAIL, NCHECK
    NCHECK += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    with io.open(os.path.join(HERE, name), encoding="utf-8") as f:
        return f.read()


DOC = read("R8_jacobson_entanglement_equilibrium_completion.md")

# ======================================================================
head("F1  外部论文与首选目标在位")

for token in [
    "Entanglement Equilibrium and the Einstein Equation",
    "10.1103/PhysRevLett.116.201101",
    "Thermodynamics of Spacetime",
    "10.1103/PhysRevLett.75.1260",
    "Jacobson 2016",
    "Jacobson 1995",
    "首选目标",
    "次选目标",
    "Cao-Carroll 2018",
    "历史控制",
]:
    check("文档锚点：%s" % token, token in DOC)

check("明确拒绝把 Oh-Park-Sin 或 Gorard 当第一目标",
      "Oh-Park-Sin" in DOC and "Gorard" in DOC and "第一目标" in DOC)

# ======================================================================
head("F2  条件定理 C1-C4 与失败树 L1-L4 在位")

for token in ["C1｜连续局域完成", "C2｜强图／预解意义的几何模流极限",
              "C3｜普适面积密度", "C4｜固定体积平衡",
              "J1｜BW/几何 boost 极限", "J2｜普适面积密度",
              "J3｜四维 Lorentzian 与固定体积局域化", "J4｜固定体积面积变分",
              "J5｜低维相关算符污染"]:
    check("条件/缺口锚点：%s" % token, token in DOC)

check("写明 C1-C4 下的 EFE 形式", "R8-EFE" in DOC and "8\\pi G" in DOC)
check("写明 Lambda 不是本步导出", "Lambda g_{ab}" in DOC and "不能由本引理给出数值" in DOC)

# ======================================================================
head("F3  任意忠实态精确熵差恒等式")

rng = np.random.default_rng(802)


def rand_density(dim, rng, floor=0.03):
    z = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    a = z @ z.conj().T + floor * np.eye(dim)
    return a / np.trace(a)


def entropy(rho):
    w = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    w = np.clip(w.real, 1e-14, None)
    return float(-np.sum(w * np.log(w)))


def rel_entropy(sigma, rho):
    ws, vs = np.linalg.eigh((sigma + sigma.conj().T) / 2)
    wr, vr = np.linalg.eigh((rho + rho.conj().T) / 2)
    ws = np.clip(ws.real, 1e-14, None)
    wr = np.clip(wr.real, 1e-14, None)
    logsig = vs @ np.diag(np.log(ws)) @ vs.conj().T
    logrho = vr @ np.diag(np.log(wr)) @ vr.conj().T
    return float(np.real(np.trace(sigma @ (logsig - logrho))))


worst = 0.0
for dim in (3, 5, 8):
    rho = rand_density(dim, rng)
    sigma = rand_density(dim, rng)
    wr, vr = np.linalg.eigh((rho + rho.conj().T) / 2)
    K = vr @ np.diag(-np.log(np.clip(wr.real, 1e-14, None))) @ vr.conj().T
    delta = sigma - rho
    lhs = entropy(sigma) - entropy(rho)
    rhs = float(np.real(np.trace(delta @ K))) - rel_entropy(sigma, rho)
    worst = max(worst, abs(lhs - rhs))
    check("dim=%d 精确熵差恒等式" % dim, abs(lhs - rhs) < 1e-10,
          "|left-right|=%.3e" % abs(lhs - rhs))

check("全随机样本最大误差 < 1e-10", worst < 1e-10, "worst=%.3e" % worst)

# ======================================================================
head("F4  一阶第一定律的 O(eps^2) 余项")

dim = 8
rho = rand_density(dim, rng)
rho_evals, rho_vecs = np.linalg.eigh((rho + rho.conj().T) / 2)
rho_evals = np.clip(rho_evals.real, 1e-6, None)
rho = rho_vecs @ np.diag(rho_evals) @ rho_vecs.conj().T
K = rho_vecs @ np.diag(-np.log(rho_evals)) @ rho_vecs.conj().T

tau = rand_density(dim, rng)

rows = []
for eps in (0.02, 0.01, 0.005):
    sigma = (1 - eps) * rho + eps * tau
    dS = entropy(sigma) - entropy(rho)
    first = float(np.real(np.trace((sigma - rho) @ K)))
    rem = abs(dS - first)
    rows.append((eps, rem))
    print("      eps=%.3f  |remainder|=%.6e" % (eps, rem))

check("rem0/rem1 > 3.5", rows[0][1] / rows[1][1] > 3.5,
      "ratio=%.3f" % (rows[0][1] / rows[1][1]))
check("rem1/rem2 > 3.5", rows[1][1] / rows[2][1] > 3.5,
      "ratio=%.3f" % (rows[1][1] / rows[2][1]))
check("小扰动余项远小于一阶量尺度", rows[-1][1] < 1e-3)

# ======================================================================
head("F5  Einstein 系数 eta=1/(4G)")

for G in (0.2, 1.0, 7.5):
    eta = 1 / (4 * G)
    coeff = 2 * np.pi / eta
    check("G=%.3f 时 2 pi / eta = 8 pi G" % G,
          abs(coeff - 8 * np.pi * G) < 1e-12,
          "coeff=%.12f" % coeff)

check("文档写明 eta=1/(4G) 是单位约定", "eta=\\frac1{4G}" in DOC and "单位约定" in DOC)

# ======================================================================
head("F6  gap/面积律硬障碍被登记")

for token in ["m_{\\rm stag}=0.23534171", "\\xi\\simeq\\frac{v_F}{m}",
              "2.157\\,L", "8.498", "\\eta_{\\rm face}m\\simeq0.467",
              "\\xi\\ll L"]:
    check("障碍锚点：%s" % token, token in DOC)

check("明确写出旧 1.078L 漏掉费米速度，修正为约 2.157L",
      "$1.078L$ 漏掉了费米速度" in DOC and "2.157\\,L" in DOC)
check("明确写出 8.498 只是有限窗口算术差距，不是不可能证明",
      "算术差距，不是结构性不可能证明" in DOC)
check("明确写出 eta_face 依赖 gap，稳定的是 eta_face*m",
      "稳定的是 $\\eta_{\\rm face}m" in DOC and "依赖 gap" in DOC)
check("登记 C3 当前不被支持及寿命细化条件下的结构障碍",
      "C3 当前不被现有数据支持" in DOC
      and "若寿命与区域尺度同步细化，则有结构性障碍" in DOC)
check("独立 L2 审计与核验文件在位",
      "R8_L2_area_density_audit.md" in DOC and "R8_L2_check.py" in DOC)
check("L1 独立审计裁决在位：现有 Zero 不能补 C2，只能条件桥",
      "R8_L1_refutation_attempt.md" in DOC
      and "R8_L1_check.py" in DOC
      and "现有 Zero 基础不能补上 C2" in DOC
      and "BGL／Hislop–Longo" in DOC)
check("C2 不再以无界交换子的字面算子范数作为目标",
      "强图／预解意义的几何模流极限" in DOC
      and "字面写法一般不成立" in DOC)
check("R9 低维相关算符批评已登记为 L5",
      "Casini–Galante–Myers" in DOC
      and "R^{2\\Delta}" in DOC
      and "先控制 L5" in DOC)
check("R10 Cao-Carroll 次选条件桥已链接",
      "R10_cao_carroll_bulk_entanglement_completion.md" in DOC
      and "只推出弱场 Einstein 方程" in DOC)

# ======================================================================
head("F7  文档没有越过“条件恢复”")

for token in ["不把条件恢复写成无条件导出", "条件恢复", "没有从 Zero 无条件导出四维 GR",
              "J1 未证时", "不能说本项目已经补上 Jacobson"]:
    check("诚实边界锚点：%s" % token, token in DOC)

forbidden = ["本项目已经补上 Jacobson 了", "无条件导出四维 GR 已完成"]
check("没有出现无条件的成功宣告", not any(x in DOC for x in forbidden))

# ======================================================================
head("汇总")
print("  断言 %d 项，不符项：%d" % (NCHECK, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
