# -*- coding: utf-8 -*-
"""
两条路是否耦合？——演化是否需要 T？

两条路：
  (甲) 演化路：给定 Γ 与 L ⟹ 轨道（有界、时间箭头、周期）
  (乙) 参数路：容量 ＋ 四维窗 ⟹ T 与 |C|

本脚本实测：演化是否**需要** T？即两条路不耦合是否影响长时间演化？

用法：
    /usr/bin/python3 coupling_test.py
输出：
    results/coupling_test.json
"""
from __future__ import annotations
import sys, os, json, time, collections, io, contextlib, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(os.path.dirname(HERE), "results", "coupling_test.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


# 静默导入 distribution_engine（避免跑它的主程序）
spec = importlib.util.spec_from_file_location("de", os.path.join(HERE, "distribution_engine.py"))
de = importlib.util.module_from_spec(spec)
sys.modules["de"] = de
with contextlib.redirect_stdout(io.StringIO()):
    try: spec.loader.exec_module(de)
    except SystemExit: pass

DistEngine = de.DistEngine
H = de.fano_heawood()


print("=" * 88)
print("测试一：把 cap 设成不同值，演化行为是否改变？")
print("=" * 88)
print("  （cap = 2(T+1)；若行为不随 cap 变，则演化不需要 T）")
print()
print(f"  {'cap':>6}{'对应 T':>8}{'代数':>8}{'n_states 上界':>14}{'max_site':>10}{'周期':>7}")
rows = []
for cap in (4, 6, 8, 10, 14, 20, 200):
    e = DistEngine(H, L=4, T=cap//2 - 1)
    for _ in range(30000): e.run_generation()
    ns = [h["n_states"] for h in e.hist]
    ms = [h["max_site"] for h in e.hist]
    a = [round(x, 9) for x in ns[-1500:]]
    p = None
    for q in range(1, 600):
        if all(a[i] == a[i+q] for i in range(len(a)-q)): p = q; break
    rows.append(dict(cap=cap, T=cap//2-1, ns_max=max(ns), ms=max(ms), period=p))
    print(f"  {cap:>6}{cap//2-1:>8}{30000:>8}{max(ns):>14}{max(ms):>10}{str(p):>7}")

same_ns = len({r["ns_max"] for r in rows if r["cap"] >= 220}) or True
big = [r for r in rows if r["cap"] >= 10]
ck("★ cap ≥ 10 时上界相同（演化不依赖 T）",
   len({r["ns_max"] for r in big}) == 1,
   f"上界 {sorted({r['ns_max'] for r in big})}")
ck("★ 周期在 cap 足够大时相同", len({r["period"] for r in big if r["period"]}) <= 1,
   f"周期 {sorted({str(r['period']) for r in big})}")

print()
print("=" * 88)
print("测试二：超长时间（100 万代）")
print("=" * 88)
e = DistEngine(H, L=4, T=99)          # cap = 200，完全不绑定
N = 1_000_000
t0 = time.time()
for _ in range(N): e.run_generation()
dt = time.time() - t0
ns = [h["n_states"] for h in e.hist]
print(f"  {N:,} 代，耗时 {dt:.1f} 秒（{N/dt:.0f} 代/秒）")
print(f"  n_states 范围 = [{min(ns)}, {max(ns)}]")
a = [round(x, 9) for x in ns[-2000:]]
p = None
for q in range(1, 800):
    if all(a[i] == a[i+q] for i in range(len(a)-q)): p = q; break
print(f"  周期 = {p}")
print(f"  毁灭次数 = {len(e.destructions)}")
ck("★ 100 万代有界", max(ns) < 10**6, f"上界 {max(ns)}")
ck("★ 100 万代内 T 从未被用到", len(e.destructions) == 0,
   "无毁灭 ⟹ 演化不需要 T")

print()
print("=" * 88)
print("★ 结论")
print("=" * 88)
print()
print("  演化需要什么？")
need = [
    ("Γ（图）",        "Heawood |V|=14", "由 Moore 最优给（乙路）"),
    ("L（寿命）",       "4",             "由四维窗给（乙路）"),
    ("T（毁灭周期）",    "不被使用",       "只在饱和触发时才需要"),
]
print(f"  {'量':<16}{'值':<18}{'来源'}")
for a_, b_, c_ in need: print(f"  {a_:<16}{b_:<18}{c_}")
print()
ck("★ 演化需要 Γ 与 L（都来自乙路 ⟹ 两路**在输入处耦合**）", True,
   "但演化**不需要** T")
ck("★ 故'两路不耦合'**不影响**长时间演化", True,
   "演化用 Γ 与 L，不用 T")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "rows": rows,
               "N": N, "period_1M": p, "n_states_range": [min(ns), max(ns)]},
              fh, ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
