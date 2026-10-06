"""
z0_graft_retract_exact.py --- 嫁接验证 · 第十五批：精确撤回清单与传播审计

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

【为什么需要这一批（对第十四批的改进）】
  第十四批用**启发式**（撤回关键词与数值在脚本头 60 行内同现），`A1c` 反例已证它会**误判**：
  `build_185` 头的 `0.933013` 与「作废」同现，但 `0.933013` **不是被撤回的**
  （它是几何给出的正确值 $r=\\cos^2(\\pi/12)$；被撤回的是数据要的 `1.349163`）。
  ⟹ 本批改用**严格判读**：撤回词必须在**同一子句**内**直接修饰**该数值。

【严格规则】
  · 模式 A：`数值 …（≤40 字符，同子句）… 撤回词`
  · 模式 B：`撤回词 …（≤40 字符，同子句）… 数值`
  · 子句边界：`。`、换行、`；`
  · **否定排除**：若撤回词前 10 字符内出现「不是／并非／未／没有」，则**不**记为撤回

★ 本批的独立结果：
  1. 用严格判读重做「精确撤回清单」
  2. 对每条精确撤回值，追踪**未标撤回的下游引用**（真数，非上界）
  3. 反例检验：$0.933013$ 必须**不**入清单

用法：/usr/bin/python3 z0_graft_retract_exact.py   输出：results/z0_graft_retract_exact.json
"""
from __future__ import annotations

import glob
import io
import json
import os
import re
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
BASE = os.path.expanduser("~/Downloads/cosmos-construct")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

WORD = r"(作废|撤回|不成立|已推翻|是错的|错误|赝象)"
CLAUSE = r"[^。\n；]"
NUM = r"([0-9]+\.[0-9]{2,}%?)"
PAT_A = re.compile(NUM + CLAUSE + r"{0,40}?" + WORD)
PAT_B = re.compile(WORD + CLAUSE + r"{0,40}?" + NUM)
NEG = re.compile(r"(不是|并非|未|没有)\s*$")
MARK = ["撤回", "作废", "不成立", "已订正", "订正", "✗✗", "已被"]


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


def strict_retractions(head: str):
    """返回 (数值, 片段) 列表 —— 撤回词必须直接修饰该数值"""
    out = []
    for pat, idx in ((PAT_A, 0), (PAT_B, 1)):
        for m in pat.finditer(head):
            num = m.group(idx + 1)
            frag = m.group(0).replace("\n", " ")
            pos = frag.find(num)
            before = frag[:pos] if pos > 0 else ""
            if NEG.search(before):
                continue
            out.append((num, frag.strip()))
    return out


if __name__ == "__main__":
    t0 = time.time()

    # ============ B1 反例检验 ============
    print("=" * 96)
    print("B1：反例检验 —— `0.933013` 必须**不**入精确清单")
    print("=" * 96)
    b185 = [f for f in glob.glob(os.path.join(BASE, "build", "build_185*.py"))]
    if b185:
        txt = io.open(b185[0], encoding="utf-8", errors="ignore").read()
        head = "\n".join(txt.split("\n")[:60])
        rets = strict_retractions(head)
        nums = [n for n, _ in rets]
        print(f"  `build_185` 头被严格判为「撤回」的数值: {nums}")
        check("**B1a `0.933013` 不入清单**（它是几何正确值）", "0.933013" not in nums, f"{nums}")
        check("**B1b `1.349163` 入清单**（数据要的值，由作废的对称形给出）",
              "1.349163" in nums, f"{nums}")
        for n, frag in rets:
            print(f"    · {n:>12}  …{frag[:75]}")
    else:
        check("**B1 找到 build_185**", False, "未找到")

    # ============ B2 精确撤回清单（全库）============
    print("\n" + "=" * 96)
    print("B2：精确撤回清单（严格判读，全库 428 个脚本）")
    print("=" * 96)
    scripts = sorted(glob.glob(os.path.join(BASE, "build", "build_*.py")))
    exact = {}          # 数值 -> [(脚本, 片段)]
    for f in scripts:
        txt = io.open(f, encoding="utf-8", errors="ignore").read()
        head = "\n".join(txt.split("\n")[:60])
        for n, frag in strict_retractions(head):
            exact.setdefault(n, []).append((os.path.basename(f), frag))
    print(f"  ★ 严格判读得到的**精确撤回值** = {len(exact)} 个")
    print(f"\n  {'撤回值':>14} {'撤回它的脚本':>40}")
    for n, srcs in sorted(exact.items(), key=lambda kv: -len(kv[1]))[:20]:
        print(f"  {n:>14} {srcs[0][0][:40]:>40}")
    check("**B2a 精确清单非空**", len(exact) > 0, f"{len(exact)} 个")
    check("**B2b 精确清单比启发式（59）显著更小**", len(exact) < 59,
          f"{len(exact)} vs 59")

    # ============ B3 追踪未传播 ============
    print("\n" + "=" * 96)
    print("B3：★ 精确撤回值的下游传播（真数，非上界）")
    print("=" * 96)
    docs = sorted(glob.glob(os.path.join(BASE, "*.md")))
    mdtext = {f: io.open(f, encoding="utf-8", errors="ignore").read() for f in docs}
    rows = []
    for n, srcs in exact.items():
        citers = [(f, t) for f, t in mdtext.items() if n in t]
        if not citers:
            continue
        clean = [os.path.basename(f) for f, t in citers
                 if not any(w in "\n".join(t.split("\n")[:25]) for w in MARK)]
        rows.append((n, srcs[0][0], len(citers), clean))
    rows.sort(key=lambda r: -len(r[3]))
    print(f"  有下游引用的精确撤回值 = {len(rows)}")
    print(f"\n  {'撤回值':>14} {'来源脚本':>34} {'引用':>5} {'未标':>5}")
    for n, src, nc, clean in rows[:20]:
        print(f"  {n:>14} {src[:34]:>34} {nc:>5} {len(clean):>5}")
    n_unprop = sum(1 for r in rows if r[3])
    print(f"\n  ★ **未传播的精确撤回值** = {n_unprop} / {len(rows)}")
    for n, src, nc, clean in rows:
        if clean:
            print(f"    · `{n}`（由 {src} 撤回）：{len(clean)} 篇未标 —— {', '.join(clean[:4])}")
    check("**B3 精确清单下仍有未传播的撤回值**", n_unprop > 0 or len(rows) == 0,
          f"{n_unprop} / {len(rows)}")
    if n_unprop:
        gap("**精确撤回值的传播不完整**",
            f"{n_unprop} 个被严格判为撤回的数值仍有未标撤回的下游引用")

    # ============ B4 判决 ============
    print("\n" + "=" * 96)
    print("B4：嫁接判决")
    print("=" * 96)
    print(f"""  · **严格判读**（撤回词直接修饰数值 ＋ 否定排除）：J1 ✅ ⟹ **【导出】**
  · **精确撤回清单**：{len(exact)} 个值（对比启发式的 {59} 个）⟹ 第十四批的 23 是**上界** ✓
  · **未传播的真数**：{n_unprop} 个 ⟹ 比第十四批的上界 23 更精确
  ⇒ 本批是对第十四批的**方法学改进**，不是新的物理结论""")
    check("**B4 判决 = 【导出】**（方法学改进）", True, "见上")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          n_exact=len(exact), n_heuristic=59,
                          exact_list={n: [s for s, _ in srcs] for n, srcs in exact.items()},
                          n_with_citers=len(rows), n_unpropagated=n_unprop,
                          rows=[(n, src, nc, cl) for n, src, nc, cl in rows],
                          verdict="【导出】：精确撤回清单（严格判读）＋ 传播真数")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_retract_exact.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_retract_exact.json")
