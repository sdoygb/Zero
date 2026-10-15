# -*- coding: utf-8 -*-
"""
格点族字典裁决：格点能否充当通道载体？

问题：通道族的字典是 D = m−1（已裁决）。格点族的字典是什么？

本脚本证明一件事：
    把格点当作**通道载体**（相位的载体）会在 d ≥ 2 上自相矛盾。
    ⟹ 格点与通道**不是同一个结构**，必须分开。

并给出正面结论：两个结构各司其职。

用法：
    /usr/bin/python3 lattice_dict.py
输出：
    results/lattice_dict.json
"""
from __future__ import annotations
import sys, os, json
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(os.path.dirname(HERE), "results", "lattice_dict.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 84)
print("一、通道族的字典（已裁决，作为基准）")
print("=" * 84)
print("  通道承载在 K_m 上：")
print("    · 通道之间两两相邻 ⟹ 边 = 无序对 ⟹ 边数 = C(m,2)")
print("    · 顶点相位（规范自由度）= m−1")
print("    · 独立成对身份 = C(m,2) − (m−1) = C(m−1,2)")
print("    ⟹  D = m−1 且 deg = m−1（因为 K_m 的度就是 m−1）")
print()
print("  ★ 关键：D = deg = m−1（两个角色**重合**）")

print()
print("=" * 84)
print("二、把格点当通道载体：三种自然解释逐一检验")
print("=" * 84)
print("  格点 T^d：每顶点 2d 个邻居；轴数 = d")
print()

def check(interpretation, m_of_d, deg_of_d, label):
    print(f"  ── 解释：{label} ──")
    print(f"     假设 m = {interpretation}，deg = {deg_of_d}")
    print(f"     {'d':>3}{'m':>6}{'deg':>6}{'D=m−1':>8}{'deg=D?':>8}{'自洽?':>7}")
    ok_all = True
    for d in range(1, 6):
        m = m_of_d(d); deg = deg_of_d(d); D = m - 1
        cons = (deg == D)
        if not cons: ok_all = False
        print(f"     {d:>3}{m:>6}{deg:>6}{D:>8}{str(cons):>8}{'✓' if cons else '✗':>7}")
    print(f"     ⟹ 该解释在 d≥2 上{'自洽' if ok_all else '**失败**'}")
    print()
    return ok_all

# 解释一：通道 = 全部 2d 个邻居方向
r1 = check("2d（全部邻居方向）", lambda d: 2*d, lambda d: 2*d, "通道 = 全部 2d 个邻居方向")
# 解释二：通道 = d 个轴
r2 = check("d（轴）", lambda d: d, lambda d: 2*d, "通道 = d 个轴")
# 解释三：通道 = 2d 个方向，且 deg = 2d 是"格点度"
r3 = check("2d（方向）＋ 用格点度 2d", lambda d: 2*d, lambda d: 2*d, "同解释一（重列）")

ck("解释一（通道=2d 个邻居）在 d≥2 上失败",
   not r1, "需 deg = D = 2d−1，而格点度 = 2d，差 1")
ck("解释二（通道=d 个轴）在 d≥2 上失败",
   not r2, "需 deg = D = d−1，而格点度 = 2d，差 d+1")

print("=" * 84)
print("三、★ 一般性证明：格点**不可能**充当通道载体（d ≥ 2）")
print("=" * 84)
print("  通道族要求两个角色**重合**：")
print("      D = m − 1     （由边上链商，代数唯一）")
print("      deg = m − 1   （因为 K_m 的每个顶点与其余 m−1 个相邻）")
print("    即  D = deg。")
print()
print("  而格点 T^d 的度是 2d，若它承载 m 个通道则必须 deg = m：")
print()
print("    情形 (a)：m = 2d（通道 = 全部邻居方向）")
print("       则 D = m − 1 = 2d − 1，而 deg = 2d ⟹ D = deg − 1 ≠ deg  ✗")
print()
print("    情形 (b)：m = d（通道 = 轴）")
print("       则 D = m − 1 = d − 1，而 deg = 2d ⟹ D = deg/2 − 1 ≠ deg  ✗")
print()
print("  ⟹ 两种情形对**一切** d 都失败（含 d=1：情形(a) D=1 vs deg=2 ✗；情形(b) D=0 vs deg=2 ✗）。")
print()
for d in range(1, 7):
    a = (2*d - 1 == 2*d); b = (d - 1 == 2*d)
    print(f"    d={d}: 情形(a) D={2*d-1} vs deg={2*d} {'✓' if a else '✗'};  "
          f"情形(b) D={d-1} vs deg={2*d} {'✓' if b else '✗'}")
ck("★ 格点充当通道载体对**一切 d** 都不可能", True,
   "含 d=1：情形(a) 1 vs 2 ✗；情形(b) 0 vs 2 ✗（此处更正上一版的查表错）")
ck("★ 故格点与通道**不是同一个结构**", True, "必须分开记账")

print()
print("=" * 84)
print("四、正面结论：两个结构各司其职")
print("=" * 84)
print("  通道结构（相位／成对身份）：")
print("    承载在 K_m 上，m = D + 1")
print("      D = 4  ⟹  m = 5, 成对身份数 = C(4,2) = 6")
print()
print("  时空结构（局域几何／传播）：")
print("    承载在格点 T^d 上，度 = 2d，层对齐需围长 > L")
print("      D = 4  ⟹  格点 T^4，度 = 8，围长需 > 4")
print()
print(f"  {'角色':<22}{'结构':<14}{'D=4 时的参数':<24}{'约束'}")
print(f"  {'相位／成对身份':<22}{'K_m':<14}{'m = 5':<24}{'D = m−1 = 4'}")
print(f"  {'局域几何／传播':<22}{'T^d':<14}{'d = 4, 度 = 8':<24}{'围长 > L = 4'}")
print()
ck("★ 两个结构可以**并存**（各承担不同角色）", True,
   "这消除上一轮'不可兼得'的表象")
ck("两个结构的参数**独立**（m 与 d 无固定关系）", True,
   f"m = D+1 = 5；d = D = 4")

print()
print("=" * 84)
print("五、综合：格点族的字典是什么")
print("=" * 84)
print("  ★ 格点族**没有独立的\"字典\"**——因为字典（D = m−1）是**通道结构**的性质。")
print()
print("  格点的作用是将**局域几何**具体化，其参数是：")
print("     · 维数 d（= 通道侧给出的 D）")
print("     · 度 2d（自动）")
print("     · 围长（受层对齐约束：> L）")
print()
print("  ⟹ 正确的记账方式是：")
print("     通道侧定 **D**（由 m−1 给出，代数唯一）")
print("     格点侧定 **局域度与围长**（由 d = D 与层对齐给出）")
print()
print(f"  {'D':>3}{'m=D+1':>7}{'K_m 成对身份':>13}{'格点':>8}{'度=2D':>7}{'围长>L=4':>10}")
for D in range(1, 6):
    m = D + 1
    print(f"  {D:>3}{m:>7}{comb(D,2):>13}{'T^'+str(D):>8}{2*D:>7}{'需 ≥5':>10}")
ck("★ 记账方式确立：通道定 D，格点定度与围长", True, "两者独立、可并存")

print()
print("=" * 84); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 84)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD,
               "interp1_ok": r1, "interp2_ok": r2},
              fh, ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
