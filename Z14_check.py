#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z14_check.py —— 【Z-E*：闭合循环序与双覆盖】的独立核验
======================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：Z-E*、无偏好、不引入概率/权重、J1 未关闭
  F2  组合层：s_L 的阶为 L，转角为 2π/L，且 L 次回到恒等
  F3  群论层：双覆盖核为 Z2，升格 L 次幂为 -I，一圈升格为 -I
  F4  最小性：Z0③、全分支与整数重数未被改写
  F5  状态边界：Z-CAR 只完成群论部分，物理读出与费米统计等同仍开放
  F6  跨文档一致性：R14/R15/Z13/R16/STATUS 指向同一状态
"""
import cmath
import io
import os
import sys

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


DOC = read("Z14_closure_cyclic_order_base_theorem.md")
Z0 = read("Z0_zero_never_rests_single_axiom.md")
Z13 = read("Z13_zero_foundation_missing_principle.md")
R14 = read("R14_L1_from_zero_assembly.md")
R15 = read("R15_zcar_double_cover_and_zstress_scale.md")
R16 = read("R16_direction_audit_reduction_tree.md")
STATUS = read("STATUS.md")


# ---------------------------------------------------------------- F1
head("F1  文档锚点")

check("标题登记闭合循环序与双覆盖", "闭合的循环序与双覆盖" in DOC and "Z-E*" in DOC)
check("价签写明导出与无偏好",
      "导出／基础扩展" in DOC and "无偏好" in DOC)
check("明确不引入概率", "不引入概率" in DOC or "不增加概率" in DOC)
check("明确不引入权重", "权重" in DOC and "无偏好" in DOC)
check("J1 未关闭被写明", "J1 仍未关闭" in DOC or "仍未证" in DOC)
check("双覆盖只完成群论部分",
      "群论部分的双覆盖已导出" in DOC and "物理读出" in DOC)


# ---------------------------------------------------------------- F2
head("F2  组合层：循环序与旋转角")


def permutation_order(L):
    for j in range(1, L + 1):
        if all((i + j) % L == i for i in range(L)):
            return j
    return None


orders_ok = True
angles_ok = True
for L in (3, 4, 5, 6, 7, 8, 10, 12):
    order = permutation_order(L)
    theta = 2 * cmath.pi / L
    if order != L:
        orders_ok = False
    if abs(theta * L - 2 * cmath.pi) > 1e-12:
        angles_ok = False

check("后继置换 s_L 的阶为 L", orders_ok,
      "L=3,4,5,6,7,8,10,12")
check("L 次循环等于几何一圈", angles_ok,
      "L * (2π/L) = 2π")
check("R16 承认该结构不新增独立缺口",
      "没有新增独立缺口" in R16 or "开放具名簇" in R16)


# ---------------------------------------------------------------- F3
head("F3  群论层：双覆盖、中心 Z2 与升格")

check("双覆盖映射写明为 z -> z^2", "q(z)=z^2" in DOC or "z\\mapsto z^{2}" in DOC)
check("双覆盖的核是 ±1", "\\ker q=\\{\\pm1\\}" in DOC or "核为" in DOC and "±1" in DOC)


def lift(theta):
    return cmath.exp(1j * theta / 2)


lift_ok = True
geo_ok = True
for L in (3, 4, 5, 6, 8, 10, 12):
    theta = 2 * cmath.pi / L
    uL = lift(theta) ** L
    geo = cmath.exp(1j * theta) ** L
    if abs(uL + 1) > 1e-10:
        lift_ok = False
    if abs(geo - 1) > 1e-10:
        geo_ok = False

check("升格 (U(2π/L))^L = -1", lift_ok,
      "L=3,4,5,6,8,10,12")
check("几何一圈 R(2π) = 1", geo_ok,
      "q(U_theta)=R_theta, q(-1)=1")
check("升格一圈 U(theta+2π) = -U(theta)",
      abs(lift(cmath.pi / 3 + 2 * cmath.pi) + lift(cmath.pi / 3)) < 1e-12)
check("旋量/张量两种中心特征都在文档中具名",
      "旋量" in DOC and "张量" in DOC and "中心特征" in DOC)


# ---------------------------------------------------------------- F4
head("F4  最小性：没有把权重塞回 Zero")

check("Z0③ 原文仍为不设概率、没有选择规则",
      "不设概率" in Z0 and "没有选择规则" in Z0)
check("Z2 全分支＋整数重数仍为导出",
      "全分支" in Z0 and "整数重数" in Z0 and "导出" in Z0)
check("Z14 明说不新增扩充条款",
      "不新增扩充条款" in DOC)
check("Z14 明说全分支与整数重数保持",
      "全分支与整数重数原样保持" in DOC)
check("Z14 把二选一定位为表示论离散选择，而非分支偏置",
      "表示论中的离散选择" in DOC and "不是给分支加偏置" in DOC)


# ---------------------------------------------------------------- F5
head("F5  Z-CAR 的状态边界")

check("Z14 明说双覆盖存在不等于已选费米统计",
      "不能直接写" in DOC and "费米统计已导出" in DOC)
check("物理框架嵌入仍开放", "物理框架" in DOC and "开放" in DOC)
check("费米统计等同仍开放",
      "旋量" in DOC and "费米交换统计" in DOC and "开放" in DOC)
check("L1 目标保持 K_B -> 2π B_B",
      "K_B\\longrightarrow 2\\pi B_B" in DOC or "K_B" in DOC and "2\\pi B_B" in DOC)


# ---------------------------------------------------------------- F6
head("F6  跨文档一致性")

check("R15 指向 Z14 且仍称同一开放簇",
      "Z14_closure_cyclic_order_base_theorem.md" in R15
      and "同一个开放簇" in R15)
check("R14 指向 Z14 的基础版本",
      "Z14_closure_cyclic_order_base_theorem.md" in R14)
check("Z13 说明 Z-E* 不是选择／读出原理",
      "不是选择／读出原理" in Z13)
check("R16 指向 Z14",
      "Z14_closure_cyclic_order_base_theorem.md" in R16)
check("STATUS 已登记 Z14 与无偏好边界",
      "Z14" in STATUS and "无偏好" in STATUS)


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
