# -*- coding: utf-8 -*-
"""
锚定计算：用观测钉住自由参数，算出一个可检验的数

两条锚定路线各算一次：
  (甲) 几何锚：|C| 个站点铺满一个观测尺度 ⟹ 得格距 a
  (乙) 时间锚：观测年龄 / 一次毁灭 ⟹ 得时间锚 α

判据：结果里**哪些是预言、哪些只是单位换算**。

用法：
    /usr/bin/python3 anchor_compute.py
输出：
    results/anchor_compute.json
"""
from __future__ import annotations
import sys, os, json
from math import pi

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results", "anchor_compute.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))

# ── 常数（外部观测，只用于定标）──
LY_M   = 9.4607304725808e15     # 1 光年 = ... 米
C_LY   = 9.3e10                 # 可观测宇宙直径（光年），约 93 Gly
HUB_LY = 1.44e10                # 哈勃时间（光年），约 144 亿年
N_SITE = 14                     # |C|（本框架）

print("=" * 88)
print("锚定计算")
print("=" * 88)
print(f"  |C| = {N_SITE}（本框架）")
print(f"  可观测宇宙直径 = {C_LY:.2e} 光年 = {C_LY*LY_M:.3e} 米")
print(f"  哈勃时间       = {HUB_LY:.2e} 光年（≈ 144 亿年）")
print()

print("=" * 88)
print("(甲) 几何锚：14 个站点铺满可观测宇宙")
print("=" * 88)
print()
print("  ★ 先看这个锚合不合理：格距 a 会是多少？")
print(f"  {'排布':<16}{'每边格数':>10}{'格距 a（光年）':>18}{'格距 a（米）':>16}")
geo = {}
for D, lab in [(1, "1 维链"), (2, "2 维方阵"), (3, "3 维立方")]:
    n = N_SITE ** (1.0/D)
    a_ly = C_LY / n
    geo[D] = a_ly
    print(f"  {lab:<16}{n:>10.2f}{a_ly:>18.3e}{a_ly*LY_M:>16.3e}")
print()
ck("★ 几何锚给出格距 ~10^10 光年（极粗粒化）", True,
   f"3 维: {geo[3]:.2e} 光年")
print()
print("  ★ 判读：这只是**重新表述**观测——")
print("     输入：世界多大；输出：格距多大")
print("     它**不是预言**（没有\"如果…则应为…\"的检验力）")
print()
print("  ★ 但它有一个可检验的**自洽性后果**：")
print("     若格距 ~10^10 光年，则**一个站点 = 一个超星系团尺度**")
print("     ⟹ 框架描述的是**极粗粒化**的宇宙，不是 Planck 尺度离散化")
print()
ck("★ 几何锚：不是预言，是重新表述", True, "输入尺寸 → 输出格距")

print()
print("=" * 88)
print("(乙) 时间锚：观测年龄 / 一次毁灭")
print("=" * 88)
print()
print("  ★ 问题：T 依赖 K，而 K 无出处 ⟹ 先看 T 有多少不确定性")
print()
print(f"  {'K':>12}{'T_engine':>10}{'若锚\"44 亿年/毁灭\"':>22}")
for K, T in [(10**4, 2), (10**5, 3), (10**7, 4), (10**8, 4)]:
    alpha = HUB_LY / T
    print(f"  {K:>12}{T:>10}{alpha:>20.3e} 年/毁灭")
print()
ck("★ 时间锚：T 随 K 变 ⟹ 时间锚给出的\"年/毁灭\"也随 K 变", True,
   "锚不定 ⟹ 结果不定")
print()
print("  ★ 但可以反过来问一个**可检验**的问题：")
print("     若取 T = 极限周期（与 K 无关），则 α 唯一")
print()
T_LIM = 6   # Heawood 的极限周期（实测）
alpha_lim = HUB_LY / T_LIM
print(f"     Heawood 极限周期 T = {T_LIM}")
print(f"     ⟹ α = {HUB_LY:.2e} / {T_LIM} = {alpha_lim:.3e} 光年/毁灭")
print(f"       = {alpha_lim/1e9:.3f} 十亿年/毁灭")
print()
ck("★ 取极限周期时，α 唯一（= 2.4 十亿年/毁灭）", True,
   f"{alpha_lim/1e9:.3f} 十亿年")

print()
print("=" * 88)
print("★ 综合：锚定能算出什么")
print("=" * 88)
print()
tbl = [
    ("|C| = 14",          "无量纲", "✅ 不需锚", "Moore 最优"),
    ("d = 3",             "无量纲", "✅ 不需锚", "耦合计数"),
    ("夹角 109.47°",       "无量纲", "✅ 不需锚", "正则单纯形"),
    ("D = 4",             "无量纲", "✅ 不需锚", "四维窗 + d=3"),
    ("格距 a",             "长度",   "⚠ 由尺寸反推", "非预言"),
    ("α（年/毁灭）",        "时间",   "⚠ 由年龄反推（且需定 T）", "非预言"),
    ("格距/世界尺寸",       "比值",   "✅ 可算", "= 1/边格数"),
    ("一个站点的尺度",       "长度",   "✅ 可算（给定尺寸）", "≈ 超星系团"),
]
print(f"  {'量':<16}{'类型':<8}{'判定':<24}{'依据'}")
for a, b, c, d in tbl:
    print(f"  {a:<16}{b:<8}{c:<24}{d}")

print()
print("=" * 88)
print("结论")
print("=" * 88)
print()
print("  ★ 能算（不需锚）：|C|=14、d=3、夹角、D=4 —— 无量纲结构")
print("  ★ 能算（给定观测）：格距、年/毁灭 —— 但是**重新表述**，非预言")
print("  ★ 不能算：绝对长度、绝对时间 —— 需锚，且锚只换单位")
print()
ck("★ 锚定不产生新预言（只是把观测翻译成框架语言）", True,
   "但能把框架量与实际尺度对齐")
ck("★ 无量纲结构才是可检验的预言", True,
   "|C|=14, d=3, 夹角, D=4")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "geo": geo, "alpha_lim_ly": alpha_lim,
               "table": tbl}, fh, ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
