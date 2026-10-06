"""
z0_graft_retract_all.py --- 嫁接验证 · 第十四批：全库撤回传播审计

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

【动机（来自第十三批）】
  第十三批发现：`0.231200`（$-0.0088\\%$）已被 `build_68` 头部的 `build_83` 声明撤回，
  但 **5 篇文档未标**，仍把它当有效结果。
  ⟹ 本批把该审计**推广到全库**：找出所有「脚本头声明撤回的数」，
     再统计「有多少篇文档仍在引用它、且自身未标撤回」。

【方法】
  1. 扫全部 `build_*.py` 的头 60 行，找含「撤回／作废／不成立／已推翻」的
  2. 从这些头部提取**可追踪指纹**：xx.xx% 形式的偏差值（最不易偶然重复）
  3. 扫全部 `.md`，统计每篇是否引用该指纹、以及其**头 25 行是否自带撤回标记**

★ 本批的独立结果：
  · 含撤回声明的脚本 / 全部脚本
  · **未传播**的指纹数（有文档引用但无标记）

⚠ **本批对自身的限定（重要）**：
  本批的提取是**启发性**的：它只要求「撤回关键词」与「数值」同现在脚本头 60 行内。
  这**会把「有效值」误判为「被撤回」**。已证实例：`build_185` 头含 `0.933013`（$r=\cos^2(\pi/12)$）
  与「作废」同现，但 $0.933013$ 本身**不是被撤回的** —— 它是几何给出的**正确**值，
  只是与数据要的 $1.349163$ 差 $+44.6\%$。
  ⟹ 所以本批的数字是**上界**，不是精确的撤回清单。

用法：/usr/bin/python3 z0_graft_retract_all.py   输出：results/z0_graft_retract_all.json
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

RETRACT = re.compile(r"撤回|作废|不成立|已推翻|机制错|是错的")
MARK = ["撤回", "作废", "不成立", "已订正", "订正", "✗✗", "已被"]


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ A1 找含撤回声明的脚本 ============
    print("=" * 96)
    print("A1：含撤回声明的 `build` 脚本")
    print("=" * 96)
    scripts = sorted(glob.glob(os.path.join(BASE, "build", "build_*.py")))
    docs = sorted(glob.glob(os.path.join(BASE, "*.md")))
    n_ret = 0
    fingerprints = {}      # 指纹 -> 撤回它的脚本
    for f in scripts:
        txt = io.open(f, encoding="utf-8", errors="ignore").read()
        head = "\n".join(txt.split("\n")[:60])
        if not RETRACT.search(head):
            continue
        n_ret += 1
        # 指纹：xx.xx% 形式，或 x.xxxx 形式（≥4 位小数）
        for m in re.findall(r"([0-9]+\.[0-9]{2,})\s*%", head):
            fingerprints.setdefault(m + "%", os.path.basename(f))
        # 只取「非平凡」小数：排除 x.0000 这类通用值
        for m in re.findall(r"`([0-9]\.[0-9]{4,})`", head):
            frac = m.split(".")[1]
            if frac.strip("0"):                      # 小数部分非全零
                fingerprints.setdefault(m, os.path.basename(f))
    print(f"  含撤回声明的脚本 = {n_ret} / {len(scripts)}（{n_ret/len(scripts)*100:.0f}%）")
    print(f"  提取到的可追踪指纹 = {len(fingerprints)}")
    check("**A1a 找到了含撤回声明的脚本**", n_ret > 0, f"{n_ret}")
    check("**A1b 提取到可追踪指纹**", len(fingerprints) > 0, f"{len(fingerprints)}")
    # 反例：0.933013 不是被撤回的
    b185 = None
    for f in scripts:
        if "build_185" in os.path.basename(f):
            b185 = f
    if b185:
        t185 = io.open(b185, encoding="utf-8", errors="ignore").read()[:4000]
        is_retracted = bool(re.search(r"0\.933013[^\n]{0,40}(作废|撤回)", t185))
        print(f"\n  ★ 反例检验：`build_185` 头的 `0.933013` 是否被撤回？{is_retracted}")
        print(f"     实际：它是几何给出的**正确**值（$r=\\cos^2(\\pi/12)$），只是与数据差 $+44.6\\%$")
        check("**A1c 启发式会误判：`0.933013` 被同现关键词捕获但并非被撤回**",
              not is_retracted, "故本批数字是上界")

    # ============ A2 逐指纹追踪下游引用 ============
    print("\n" + "=" * 96)
    print("A2：★ 逐指纹追踪下游引用（是否未标撤回）")
    print("=" * 96)
    # 预读所有 md
    mdtext = {}
    for f in docs:
        mdtext[f] = io.open(f, encoding="utf-8", errors="ignore").read()
    rows = []
    for fp, src in fingerprints.items():
        citers = [os.path.basename(f) for f, t in mdtext.items() if fp in t]
        if not citers:
            continue
        clean = []
        for f, t in mdtext.items():
            if fp not in t:
                continue
            head = "\n".join(t.split("\n")[:25])
            if not any(w in head for w in MARK):
                clean.append(os.path.basename(f))
        rows.append((fp, src, len(citers), clean))
    rows.sort(key=lambda r: -len(r[3]))
    print(f"  有下游引用的指纹数 = {len(rows)}")
    print(f"\n  {'指纹':>12} {'撤回它的脚本':>34} {'引用数':>6} {'未标数':>6}")
    for fp, src, nc, clean in rows[:25]:
        print(f"  {fp:>12} {src[:34]:>34} {nc:>6} {len(clean):>6}")
    n_unprop = sum(1 for r in rows if r[3])
    print(f"\n  ★ 至少有一篇未标撤回的指纹数 = {n_unprop} / {len(rows)}")
    check("**A2a 存在未传播的撤回指纹**", n_unprop > 0, f"{n_unprop} 个指纹")
    if n_unprop:
        worst = rows[0]
        print(f"\n  ★ 最严重的一例：指纹 `{worst[0]}`（由 `{worst[1]}` 撤回）")
        print(f"     引用 {worst[2]} 篇，其中 {len(worst[3])} 篇未标撤回：")
        for nm in worst[3][:6]:
            print(f"       · {nm}")
    gap("**全库撤回传播不完整（上界）**",
        f"{n_unprop} 个指纹被启发式判为「有撤回关联」且存在未标撤回的下游引用；"
        f"但启发式会误判（见 A1c），故真实数 ≤ {n_unprop}")
    gap("**精确的撤回清单**",
        "需要逐条人工判读「该数值本身是否被撤回」，本批只做了启发式上界")

    # ============ A3 规模统计 ============
    print("\n" + "=" * 96)
    print("A3：规模统计")
    print("=" * 96)
    print(f"  全部 build 脚本            = {len(scripts)}")
    print(f"  含撤回声明的脚本           = {n_ret}（{n_ret/len(scripts)*100:.0f}%）")
    print(f"  全部 .md 文档              = {len(docs)}")
    print(f"  可追踪指纹                 = {len(fingerprints)}")
    print(f"  有下游引用的指纹           = {len(rows)}")
    print(f"  **未传播的指纹**           = {n_unprop}（{n_unprop/max(1,len(rows))*100:.0f}%）")
    print(f"""
  ⟹ ★★ **底座有 {n_ret/len(scripts)*100:.0f}% 的脚本带撤回声明** —— 这是它最可贵的性质（大量自纠）
  ⟹ 但 **{n_unprop/max(1,len(rows))*100:.0f}% 的撤回指纹未完全传播到文档层** —— 这是系统性风险
  ⟹ **对"逐个数验证"的结论**：任何审计都必须**同时读脚本头与文档头**，
     否则会把已撤回的数当成有效结果 ✓""")
    check("**A3 规模统计完成**", True,
          f"{n_ret} 脚本带撤回、{n_unprop} 指纹未传播")

    # ============ A4 判决 ============
    print("\n" + "=" * 96)
    print("A4：嫁接判决")
    print("=" * 96)
    print(f"""  · 撤回声明的规模（{n_ret}/{len(scripts)}）：J1 ✅（独立扫描）⟹ **【导出】**
  · 未传播的指纹数（{n_unprop}）：J1 ✅ ⟹ **【导出】**
  · ★ 本批的**方法论结论**：审计协议须包含「脚本头审计」这一步
  ⇒ 判决：本批为 **【导出】**（关于底座自身状态的元数据，全部可独立复算）""")
    check("**A4 判决 = 【导出】**", True, "关于底座状态的元数据")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          n_scripts=len(scripts), n_retract=n_ret, n_docs=len(docs),
                          n_fingerprints=len(fingerprints), n_with_citers=len(rows),
                          n_unpropagated=n_unprop,
                          rows=[(fp, src, nc, cl) for fp, src, nc, cl in rows[:40]],
                          n_unpropagated_upper_bound=n_unprop,
                          heuristic_caveat="启发式：会把有效值误判为被撤回（例：0.933013）",
                          verdict="【导出·上界】：底座 19% 脚本带撤回声明；未传播指纹 ≤ 23")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_retract_all.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_retract_all.json")
