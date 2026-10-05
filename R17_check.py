#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R17_check.py —— 【L1 必要性审计与 L5 前置门】的独立核验
========================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：L1、Z14-Z16 非必经、R17-STOP、L5 优先、L1 开放
  F2  临界路径分类：选择型 / 分析型 / 常数型
  F3  停止规则：不再下钻表示层
  F4  宇称 no-go：中央宇称清除单费米算符，但不清除二费米标量算符
  F5  L5 先于 L1 的顺序与 R12 边界
  F6  路线切换：L5-CERT 主攻，Z-CORE+Z-TAIL 支线
  F7  跨文档状态一致：R16/STATUS/R0/INDEX
"""
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = "v" if ok else "x"
    if level == "sup" and ok:
        tag = "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    return io.open(os.path.join(HERE, name), encoding="utf-8").read()


DOC = read("R17_L1_critical_path_and_L5_gate.md")
R0 = read("R0_publication_theorem.md")
R8 = read("R8_jacobson_entanglement_equilibrium_completion.md")
R9 = read("R9_external_GR_derivations_landscape.md")
R12 = read("R12_zero_native_gap_filling.md")
R13 = read("R13_L1_strong_resolvent_attempt.md")
R16 = read("R16_direction_audit_reduction_tree.md")
Z14 = read("Z14_closure_cyclic_order_base_theorem.md")
Z15 = read("Z15_zcar_no_go_and_jordan_wigner_readout.md")
Z16 = read("Z16_zunif_balanced_regular_module.md")
STATUS = read("STATUS.md")


# ---------------------------------------------------------------- F1
head("F1  文档锚点")

check("文档锁定 L1 为唯一锚点",
      "唯一目标" in DOC and "K_B" in DOC and "2\\pi B_B" in DOC)
check("Z14-Z16 非必经结论在位",
      "命题 R17.1" in DOC and "not\\subset" in DOC
      and "CriticalPath(L1)" in DOC)
check("R17-STOP 停止规则在位",
      "R17-STOP" in DOC and "不得继续增加新的" in DOC)
check("L5 优先于 L1 写明",
      "L5 是否应先于 L1" in DOC and "**是**" in DOC)
check("J1 未关闭写明",
      "L1 当前状态" in DOC and "**仍未关闭**" in DOC)
check("四个不做的事在位",
      "不再下钻" in DOC and "不在没有 L5 结论时硬攻" in DOC
      and "不再增加新的表示层" in DOC)


# ---------------------------------------------------------------- F2
head("F2  L1 临界路径分类")

SELECTION = ["Z-CRIT-DER", "Z-SCALE", "Z-HILB", "Z-CONF"]
ANALYTIC = ["Z-CORE", "Z-TAIL"]
CONSTANT = ["Z-STRESS"]
ALL = SELECTION + ANALYTIC + CONSTANT

check("七个 R13 缺口全部出现在 R17",
      all(x in DOC for x in ALL), "计数=%d" % sum(x in DOC for x in ALL))
check("选择型四类完整",
      all(x in DOC for x in SELECTION) and "选择型" in DOC)
check("分析型两类完整",
      all(x in DOC for x in ANALYTIC) and "分析型" in DOC)
check("常数型一类完整",
      all(x in DOC for x in CONSTANT) and "常数型" in DOC)
check("分类无重复：4+2+1=7",
      len(set(ALL)) == 7 and len(SELECTION) + len(ANALYTIC) + len(CONSTANT) == 7)
check("分析型与常数型可直接改变 L1 可证性",
      "可以直接改变 L1 的可证性" in DOC)


# ---------------------------------------------------------------- F3
head("F3  停止规则")

check("新子问题必须关闭 Z-CORE/Z-TAIL 或 L5",
      "关闭 `Z-CORE`／`Z-TAIL`" in DOC or "关闭 `Z-CORE`" in DOC)
check("继续 Z-UNIF 物理必然性被标为停",
      "继续解释 `Z-UNIF` 的物理必然性" in DOC and "**停**" in DOC)
check("自旋结构与模式相位保留为 Z-READ 条件",
      "作为 `Z-READ` 条件保留" in DOC)
check("Z-TAIL uniform 界被标为可做",
      "证明 `Z-TAIL` 的 uniform 界" in DOC and "**可做**" in DOC)
check("Z-CAR 分支被冻结为条件输入",
      "`Z-CAR` 分支冻结为条件输入" in DOC)


# ---------------------------------------------------------------- F4
head("F4  中央宇称 no-go")


def kron_all(mats):
    out = np.array([[1.0]], dtype=complex)
    for mat in mats:
        out = np.kron(out, mat)
    return out


NM = 3
ident = np.eye(2, dtype=complex)
sig_z = np.diag([1.0, -1.0]).astype(complex)
sig_minus = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
sig_plus = sig_minus.conj().T

parity = kron_all([sig_z] * NM)
c = []
cdag = []
for j in range(NM):
    prefix = [sig_z] * j
    suffix = [ident] * (NM - j - 1)
    c.append(kron_all(prefix + [sig_minus] + suffix))
    cdag.append(kron_all(prefix + [sig_plus] + suffix))

single_odd = all(
    np.allclose(parity @ c[j] @ parity, -c[j])
    and np.allclose(parity @ cdag[j] @ parity, -cdag[j])
    for j in range(NM)
)
check("中央宇称使单费米算符为奇", single_odd)

bilinear_even = True
pairing_even = True
for i in range(NM):
    for j in range(NM):
        o1 = cdag[i] @ c[j]
        o2 = c[i] @ c[j] + cdag[i] @ cdag[j]
        if not np.allclose(parity @ o1 @ parity, o1):
            bilinear_even = False
        if not np.allclose(parity @ o2 @ parity, o2):
            pairing_even = False

check("二费米数守恒双线性保持偶", bilinear_even)
check("二费米配对双线性保持偶", pairing_even)
check("文档明确不证明危险算符一定存在",
      "并不证明危险算符一定存在" in DOC)
check("文档明确只证明中央宇称不够",
      "不足以构成 L5 的对称保护证书" in DOC)


# ---------------------------------------------------------------- F5
head("F5  L5 先于 L1")

check("R8 明写先控制 L5 再攻 L1",
      "先控制 L5" in R8 and "再攻 L1" in R8)
check("R9 把低维相关算符列为必须解释项",
      "低维相关算符" in R9 and "固定体积" in R9)
check("R12 登记 L5 有限维 no-go",
      "L5-NG" in R12 and "D\\ne o(R^d)" in R12)
check("R17 把 L5 写成上游门槛",
      "上游门槛" in DOC and "不是 L1 旁边的技术注脚" in DOC)
check("R17 给出 L5-CERT 的 go/no-go",
      "L5-CERT" in DOC and "go/no-go" in DOC)


# ---------------------------------------------------------------- F6
head("F6  支线与路线排序")

check("1+1D 支线限定 Z-CORE+Z-TAIL",
      "1+1D 的 `Z-CORE + Z-TAIL`" in DOC)
check("支线不得冒充四维 L1",
      "禁止把这条 1+1D 条件定理写成四维 L1" in DOC)
check("Z-CONF 仍单列开放",
      "`Z-CONF` 仍必须单独开放" in DOC)
check("路线排序 L5 第一、Z-CORE+Z-TAIL 第二、Cao-Carroll 第三、Z-CAR 第四",
      all(x in DOC for x in ["**L5-CERT 的 go/no-go**", "**1+1D `Z-CORE+Z-TAIL`**",
                            "**Cao-Carroll 对冲**", "继续 `Z-CAR` 表示层细分"]))
check("L1 条件恢复判定不被抬高",
      "Conditional GR 判定不变" in DOC)


# ---------------------------------------------------------------- F7
head("F7  跨文档状态一致")

check("R16 已增加 R17 后记",
      "R17" in R16 and "L5" in R16)
check("STATUS 已登记 R17",
      "R17" in STATUS and "L5-CERT" in STATUS)
check("STATUS 仍保持 L1 开放",
      "L1" in STATUS and "仍未关闭" in STATUS)
check("R0 已把外部与方向审计范围扩到 R18 及以后",
      any(("R8`–`R%d" % n) in R0 or ("R8-R%d" % n) in R0
          for n in range(18, 40)))
check("Z14/Z15/Z16 都没有被写成关闭 L1",
      "J1 仍未关闭" in Z14 and "L1" in Z15 and "L1" in Z16)


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
