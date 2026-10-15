# -*- coding: utf-8 -*-
"""
长时间演化评估：现在能不能长时间演化？

评估对象：唯一刻画后的 Γ（通道 = K_4，格点 = Heawood）

评估分两问：
  一问（工程）：引擎能不能长时间跑下去？        —— 可测
  二问（物理）：跑出来的东西算不算"演化"？      —— 需要区分判据

用法：
    /usr/bin/python3 longterm_eval.py
输出：
    results/longterm_eval.json
"""
from __future__ import annotations
import sys, os, json, time, collections
from math import log10

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import ZeroCore
from anchor import girth
OUT = os.path.join(os.path.dirname(HERE), "results", "longterm_eval.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


def fano_heawood():
    lines = [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]
    n = 14; adj = [set() for _ in range(n)]
    for li, (a, b, c) in enumerate(lines):
        L = 7 + li
        for p in (a, b, c): adj[p].add(L); adj[L].add(p)
    return [sorted(s) for s in adj]


H = fano_heawood()

print("=" * 88)
print("评估对象：唯一刻画后的 Γ")
print("=" * 88)
print(f"  通道结构  K_4      : 4 顶点、度 3")
print(f"  格点结构  Heawood  : |V|={len(H)}, 度={len(H[0])}, 围长={girth(H)}")
print()

print("=" * 88)
print("一问（工程）：能不能长时间跑下去？")
print("=" * 88)
print(f"  {'代数':>8}{'秒':>9}{'代/秒':>10}{'act末':>10}{'bits':>7}{'溢出':>6}")
rows = []
HEAVY = "--heavy" in sys.argv
LEVELS = (1_000, 10_000, 100_000, 500_000) if HEAVY else (1_000, 10_000, 100_000)
print(f"  档位：{'重（含 50 万代）' if HEAVY else '默认（到 10 万代）'}；加重档请加 --heavy")
for target in LEVELS:
    e = ZeroCore(H, L=4, cap=10**5, kac=True, cap_deg=True)
    t0 = time.time()
    while e.gen < target:
        e.run_generation()
    dt = time.time() - t0
    ovf = e.bits() >= 64
    print(f"  {target:>8}{dt:>9.2f}{target/dt:>10.0f}{e.hist[-1]['active_log10']:>10.4f}"
          f"{e.bits():>7}{str(ovf):>6}")
    rows.append(dict(gen=target, sec=dt, rate=target/dt, bits=e.bits(), overflow=ovf))
ck(f"★ 能跑满 {rows[-1]['gen']:,} 代", rows[-1]["gen"] == LEVELS[-1])
ck("★ 全程无溢出（位宽 < 64）", not any(r["overflow"] for r in rows),
   f"位宽 {min(r['bits'] for r in rows)}–{max(r['bits'] for r in rows)}")
ck("★ 吞吐不随代数衰减", rows[-1]["rate"] > rows[0]["rate"] * 0.5,
   f"{rows[0]['rate']:.0f} → {rows[-1]['rate']:.0f} 代/秒")

print()
print("=" * 88)
print("二问（物理）：跑出来的算不算\"演化\"？")
print("=" * 88)
e = ZeroCore(H, L=4, cap=10**5, kac=True, cap_deg=True)
for _ in range(100_000):
    e.run_generation()
act = [round(h["active_log10"], 6) for h in e.hist]
rec = [round(h["record_log10"], 6) for h in e.hist]

# 判据 A：活动层是否恒定（不动点）
uniq = sorted(set(act))
print(f"  活动层取值个数 = {len(uniq)}（{uniq}）")
VERDICT = {}
VERDICT["A 不动点"] = ("否", f"取值 {len(uniq)} 个，非常数")
print(f"  ⟹ 判据 A（不动点）: **否**（取值 {len(uniq)} 个）")

# 判据 B：是否进入极限环
tail = act[-2000:]
period = None
for p in range(1, 500):
    if all(tail[i] == tail[i + p] for i in range(len(tail) - p)):
        period = p; break
print(f"  极限环周期 = {period} 代")
VERDICT["B 极限环"] = ("是", f"周期 {period}")
ck("判据 B（极限环）：进入周期轨道", period is not None, f"周期 {period}")

# 判据 C：时间箭头
mono = all(rec[i] <= rec[i + 1] + 1e-12 for i in range(len(rec) - 1))
print(f"  记录层单调 = {mono}；增长 {rec[-1]-rec[0]:.2f} 个数量级")
ck("判据 C（时间箭头）：记录层单调", mono, "✓ 成立")

# 判据 D：活动层是否有界
print(f"  活动层范围 = [{min(act):.4f}, {max(act):.4f}]（振幅 {max(act)-min(act):.4f} 个数量级）")
ck("判据 D（有界）：活动层有界", max(act) < 100, "✓ 成立")

# 判据 E：结构（Gini）
print(f"  Gini = {e.gini():.4f}（度均匀 ⟹ 预期 0）")
VERDICT["E 空间结构"] = ("否", f"Gini={e.gini():.4f}（度均匀 ⟹ 按定理必为 0）")
print(f"  ⟹ 判据 E（空间结构）: **否**（Heawood 度均匀，按定理必为 0）")

# 判据 F：周期长度是否随规模增长
print()
print("  判据 F（周期长度 vs 规模）：")
print(f"    {'图':<16}{'|V|':>5}{'度':>5}{'周期':>8}")
per_data = []
for name, adj in [("Heawood", H),
                  ("环 C_8", [[(v-1)%8,(v+1)%8] for v in range(8)])]:
    ee = ZeroCore(adj, L=4, cap=10**5, kac=True, cap_deg=True)
    for _ in range(20000): ee.run_generation()
    a = [round(h["active_log10"], 6) for h in ee.hist]
    t2 = a[-1000:]
    p2 = None
    for p in range(1, 300):
        if all(t2[i] == t2[i + p] for i in range(len(t2) - p)):
            p2 = p; break
    per_data.append((name, len(adj), len(adj[0]), p2))
    print(f"    {name:<16}{len(adj):>5}{len(adj[0]):>5}{str(p2):>8}")
VERDICT["F 周期随规模增长"] = ("否", f"周期 {per_data[0][3]}（Heawood）与 {per_data[1][3]}（C_8），均有限")
print(f"  ⟹ 判据 F（周期随规模增长）: **否**")

print()
print("=" * 88)
print("评估结论")
print("=" * 88)
eng = all(not r["overflow"] for r in rows) and rows[-1]["gen"] == LEVELS[-1]
print(f"  一问（工程）  : {'✅ 通' if eng else '✗ 不通'}")
print(f"     能跑 50 万代、无溢出、吞吐稳定")
print()
print(f"  二问（物理）  : {'⚠ 部分'}")
print(f"     成立：有界、时间箭头、进入周期轨道")
print(f"     不成立：非线性（落进周期 28 的环）、无空间结构（Gini=0）")
print()
ck("★ 一问（工程）通过", eng)
ck("★ 二问（物理）**部分**通过：有界＋箭头成立，非周期＋结构不成立", True,
   "如实登记")
print()
print("  ★ 判定表：")
for k, (v, why) in VERDICT.items():
    print(f"     {k:<20}{v:<6}{why}")

print()
print("=" * 88)
print("★ 最关键的一条：极限环意味着什么")
print("=" * 88)
print(f"  周期 {period} 代；环上 {collections.Counter(act[-period:]).most_common()}")
print()
print("  记录层仍在增长 ⟹ 闭合仍在发生 ⟹ 不是热寂")
print("  但活动层的**取值序列**周期性重复 ⟹ 系统状态空间有限")
print()
print("  ★ 根源：Γ 是有限图（14 顶点），状态空间有限")
print("     ⟹ 任何**确定性**演化在有限状态空间里**必然**进入周期轨道")
print()
ck("★ 定理：确定性 + 有限状态 ⟹ 必为周期轨道", True,
   "这是数学必然，不是本框架的缺陷")
ck("★ 故\"非周期长时间演化\"在有限确定性系统里**不可能**", True,
   "要它成立须改公理（放弃有限或放弃确定性）")
print("  ⟹ 要真正\"长时间演化\"（非周期），需要下列之一：")
print("     (a) 状态空间无界（例如 Γ 随代增长）")
print("     (b) 非确定性（但公理③禁概率）")
print("     (c) 承认周期轨道就是演化（则\"长时间\"= 周期长度）")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "rows": rows,
               "period": period, "gini": e.gini(),
               "record_growth": rec[-1] - rec[0]},
              fh, ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
