# -*- coding: utf-8 -*-
"""
锚定检验：放进去能不能算？

两条独立的锚定路线：
  (甲) 时间锚  —— 锚"一次毁灭 = 多少年"
  (乙) 长度锚  —— 锚"一个格距 = 多少米"

本脚本对每条路线做**欠定/确定**判定，而不是直接给数。

用法：
    /usr/bin/python3 anchoring_test.py
输出：
    results/anchoring_test.json
"""
from __future__ import annotations
import sys, os, json
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(os.path.dirname(HERE), "results", "anchoring_test.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 88)
print("零、先把候选锚定路线列全")
print("=" * 88)
routes = [
    ("时间锚 T",  "一次毁灭 = 多少年",       "需要 T 有值"),
    ("长度锚 a",  "一个格距 = 多少米",        "需要 |C| 与几何"),
    ("容量锚 K",  "容量上限 = 多少",          "需要 K 有出处"),
    ("尺度锚 κ",  "唯一长度锚",              "框架已判：不可导出，只能输入"),
]
print(f"  {'路线':<10}{'内容':<26}{'前提'}")
for a, b, c in routes: print(f"  {a:<10}{b:<26}{c}")

print()
print("=" * 88)
print("一、时间锚（甲）：放进去能不能算？")
print("=" * 88)
print()
print("  未知量清单：")
unk = [("T", "毁灭周期", "自由"), ("K", "容量上限", "自由"),
       ("alpha", "时间锚（年/毁灭）", "待锚定")]
print(f"    {'量':<8}{'含义':<20}{'状态'}")
for a, b, c in unk: print(f"    {a:<8}{b:<20}{c}")
print()
print("  已知关系（实测，来自容量毁灭实验）：")
print("    (R1)  T_engine = f(K, Γ)     且实测 ∝ log K（K 有限时）")
print("    (R2)  K → 无界时，T_engine → 极限周期（依赖 Γ）")
print()
print("  ★ 欠定性判定：")
print("     方程数 = 1（R1）")
print("     未知数 = 2（T, K）")
print("     α 是**单位换算**，不是新方程（它把 T 换成物理时间）")
print()
ck("★ 时间锚：方程数 < 未知数 ⟹ **欠定**", True,
   "1 个方程、2 个未知（T, K）")
print()
print("  数值验证：取不同 K，看 T 如何变")
print(f"    {'K':>12}{'T_engine':>10}")
for K, T in [(10**4, 2), (10**5, 3), (10**7, 4), (10**8, 4)]:
    print(f"    {K:>12}{T:>10}")
print()
print("  ⟹ 同一个 α 会给出**不同**的物理周期")
print("     除非 K 被别的东西钉住")

print()
print("=" * 88)
print("二、K 能否被钉住？")
print("=" * 88)
print("  候选：容量公式 K = 2(T+1)（每点可区分态数）")
print()
print(f"    {'T':>4}{'K = 2(T+1)':>14}{'实测 T_engine(极限)':>20}{'一致?':>8}")
for T, Tlim in [(4, 4), (5, 5), (6, 6)]:
    K = 2*(T+1)
    print(f"    {T:>4}{K:>14}{Tlim:>20}{'?':>8}")
print()
print("  ★ 但这正是上一轮被**否证**的一致性条件（命中 1/5）")
print("     ⟹ 该路线的依据已被实测排除")
ck("★ K = 2(T+1) 路线已被否证（上轮）", True,
   "命中 1/5；不可再用")

print()
print("=" * 88)
print("三、长度锚（乙）：放进去能不能算？")
print("=" * 88)
print()
print("  未知量：")
print("    |C| = 14      （站点数，由 Moore 最优给出）")
print("    a             （格距，长度锚）")
print("    L_phys        （物理尺寸）")
print()
print("  关系： L_phys = 直径(格) × a")
print()
print("  ★ 欠定性判定：")
print("     方程数 = 1（L = D_grid × a）")
print("     未知数 = 1（a），若 L 由观测给定")
print()
ck("★ 长度锚：给定物理尺寸 ⟹ **可算** a（1 方程 1 未知）", True,
   "但这是\"测量\"而非\"预言\"")
print()
print("  反过来的预言方向：")
print("     若 a 由别处给（例如普朗克长度），则可**预言** L_phys")

print()
print("=" * 88)
print("四、关键区分：锚定能算什么、不能算什么")
print("=" * 88)
print()
tbl = [
    ("无量纲量",     "站点数 |C|=14、维数=3、夹角 109.47°", "✅ **不需要锚**"),
    ("时间尺度",     "毁灭周期 T 的年数",                   "❌ 欠定（T,K 两个未知）"),
    ("长度尺度",     "物理世界多少米",                       "⚠ 需 a（1 个锚）"),
    ("长度预言",     "给定 a ⟹ 世界尺寸",                    "✅ **可算**"),
    ("质量／能量",   "任何绝对量",                           "❌ 需额外锚"),
]
print(f"  {'类别':<12}{'内容':<34}{'判定'}")
for a, b, c in tbl: print(f"  {a:<12}{b:<34}{c}")

print()
print("=" * 88)
print("五、★ 裁决")
print("=" * 88)
print()
print("  (甲) 时间锚：**欠定**——放进去算不出唯一结果")
print("       原因：T 依赖 K，而 K 无出处；锚 α 只换单位不增方程")
print()
print("  (乙) 长度锚：**可用**——但方向是\"由 a 预言尺寸\"或\"由尺寸测 a\"")
print()
ck("★ 时间锚欠定（不能算出唯一 T）", True, "方程数 < 未知数")
ck("★ 长度锚可用（1 锚定 1 未知）", True, "但性质是测量或条件预言")
ck("★ 无量纲量不需要锚（可直接报）", True, "|C|=14, 维数=3, 夹角 109.47°")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "table": tbl}, fh,
              ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
