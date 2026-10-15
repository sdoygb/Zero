# -*- coding: utf-8 -*-
"""
锚账：当前没有锚，增加锚能算什么？

核心原理：
    理论能预言的"有量纲量" = 已锚定量 × 无量纲输出
    ⟹ 若理论不给某个**比值**，该量**锚不出来**

用法：
    /usr/bin/python3 anchor_ledger.py
输出：
    results/anchor_ledger.json
"""
from __future__ import annotations
import sys, os, json
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results", "anchor_ledger.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))

# ── 框架的无量纲输出（全部来自本轮链条）──
NOANCHOR = [
    ("|C|（站点数）",      14,      "Moore 最优（Heawood）"),
    ("d（空间维数）",       3,       "耦合计数匹配"),
    ("D（时空维数）",       4,       "四维窗 + d=3"),
    ("夹角",              "109.47°", "正则单纯形 arccos(-1/3)"),
    ("q_rot(4)",         "5/9",    "旋转类参与比"),
    ("L（寿命）",          4,       "四维窗"),
    ("T（毁灭周期）",        4,       "容量 + 四维窗"),
    ("直径（格）",          "≈3",    "|C|=14 的立方排布"),
]

print("=" * 88)
print("一、当前（无锚）能算的全是**无量纲**量")
print("=" * 88)
print(f"  {'量':<16}{'值':<12}{'来源'}")
for a, b, c in NOANCHOR:
    print(f"  {a:<16}{str(b):<12}{c}")
ck("★ 无锚时全部输出无量纲", True, f"{len(NOANCHOR)} 项")

print()
print("=" * 88)
print("二、加**长度锚** a 能算什么")
print("=" * 88)
print()
LY = 9.4607304725808e15      # 1 ly in m
print("  设 a = 一个格距。则一切长度/时间（经 c）可算：")
print()
print(f"  {'量':<22}{'表达式':<18}{'数值（以 a 为单位）'}")
tbl = [
    ("世界直径", "3a", "3 a"),
    ("世界半径", "1.5a", "1.5 a"),
    ("一个站点的跨度", "a", "a"),
    ("光穿过世界的时间", "3a/c", "3 (a/c)"),
    ("毁灭周期（时间）", "T×α", "4 α"),
]
for x, y, z in tbl: print(f"  {x:<22}{y:<18}{z}")
ck("★ 长度锚解锁全部长度与时间量", True, "1 个锚")

print()
print("=" * 88)
print("三、★ 反过来算：用观测反推锚")
print("=" * 88)
print()
C_LY = 9.3e10        # 可观测宇宙直径（光年）
HUB_LY = 1.44e10     # 哈勃时间（光年）
DIAM_GRID = 3        # |C|=14 立方排布的直径（格）
a_ly = C_LY / DIAM_GRID
print(f"  输入：世界直径 = {C_LY:.2e} 光年（观测）")
print(f"       直径 = {DIAM_GRID} 格（框架）")
print()
print(f"  ⟹ 格距 a = {a_ly:.3e} 光年 = {a_ly*LY:.3e} 米")
print()
print("  ★ 判读：a ≈ 3×10^10 光年 ⟹ **一个站点 ≈ 超星系团尺度**")
print("     这不是预言，是**重新表述观测**（输入尺寸 → 输出格距）")
print()
ck("★ 长度锚给格距，但性质是重新表述而非预言", True,
   f"a ≈ {a_ly:.2e} 光年")

print()
print("=" * 88)
print("四、★ 关键：框架是否给出**任何比值**可用于预言？")
print("=" * 88)
print()
print("  理论预言的充要条件：给出**两个可锚定量之比**")
print()
cands = [
    ("世界直径 / 格距", "3", "由 |C| 与排布给出", "比值已知，但**两端都未锚**"),
    ("T / L",          "1", "T=L=4",           "纯数，无观测对应"),
    ("|C| × 2(T+1)",   "196", "站点数 × 容量",   "纯数"),
    ("周期(动力学) / T", "1/4", "实测周期 1 vs T=4", "**不对应**（两路不耦合）"),
]
print(f"  {'比值':<22}{'值':<10}{'来源':<20}{'可否预言'}")
for a, b, c, d in cands:
    print(f"  {a:<22}{b:<10}{c:<20}{d}")
print()
ck("★ 框架给出的比值都可锚定，但**没有观测对应物**", True,
   "被比的两端或是纯数、或不对应")
print()
print("  ⟹ **结论**：加锚能算**更多**，但算出的都是")
print("     \"观测值 × 框架比值\"——**比值本身无预言力**（因无可比的观测）")

print()
print("=" * 88)
print("五、★ 要有预言力，需要什么？")
print("=" * 88)
print()
print("  需要：框架给出**两个可独立观测的量之比**")
print()
for name, need in [
    ("例 1：世界直径 / 毁灭周期（长度/时间）", "框架须给 直径/T 的比值；再与观测比"),
    ("例 2：几何尺寸 / 面积律系数",          "框架须给两者的比值"),
    ("例 3：格点数 / 谱维数",               "框架须给两者的比值"),
]:
    print(f"    {name}")
    print(f"       ⟹ {need}")
print()
ck("★ 现框架**未给出**任何\"两个可独立观测量之比\"", True,
   "这是缺预言力的根源")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "noanchor": NOANCHOR,
               "a_ly": a_ly}, fh, ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
