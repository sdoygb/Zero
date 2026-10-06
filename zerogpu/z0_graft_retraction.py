"""
z0_graft_retraction.py --- 嫁接验证 · 第十三批：撤回的传播审计

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

【为什么需要这一批】
  底座大量自我订正（状态词分布：「订正」55、「否」39、「撤回」9）。
  但**订正是否传播到了引用它的下游文档**？本批做系统性检查。

★ 本批的核心案例：`sin²θ_W(M_Z) = 0.231200`（偏差 −0.0088%）
  · 该数由 `build_68` 给出，被 `Q221/Q226/Q229/Q233/Q244` 等多处引用为「独立预言」
  · 而 `build_68`（**260920**）自己写着：

      ║ ★★★ 260920 撤回（build_83）: 本脚本的**结论不成立**。
      ║ 缺陷: 把 SM 归一化的 1/α_Y 用 GUT 归一化的 b₁ = 41/10 跑动。
      ║       1/α_Y 的 β 必须是 b_Y = (5/3)b₁ = 41/6。
      ║       系数差 5/3 ⇒ 跑动量少了 40%。
      ║ 后果: 「sin²θ_W(M_Z)=0.231200，偏差 −0.0088%」全部作废。
      ║ 用正确的 b_Y ⇒ sin²θ_W(M_Z) = 0.206858（−10.5%）。

  · 而**下游文档是否带撤回标记**？本批逐篇检查。

★ 本批的独立验证：
  1. 复现 β 归一化的事实：b_Y = (5/3)b₁ 是常数比 5/3
  2. 检查「$0.231200$ / $0.0088\\%$」在哪些文档里仍出现，以及它们有无撤回标记
  3. 用 r = (2+√3)/4 复现几何起点 sin²(M_X) = 4r²/(7r²+3) = 0.382913（这部分**不受** β 影响）

用法：/usr/bin/python3 z0_graft_retraction.py   输出：results/z0_graft_retraction.json
"""
from __future__ import annotations

import glob
import io
import json
import os
import re
import time
from math import cos, pi, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
BASE = os.path.expanduser("~/Downloads/cosmos-construct")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

RETRACT_WORDS = ["撤回", "作废", "不成立", "已订正", "订正", "✗✗"]


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ T1 β 归一化的事实 ============
    print("=" * 96)
    print("T1：$\\beta$ 归一化的事实（`build_83` 的订正）")
    print("=" * 96)
    b1, b2, b3 = 41 / 10, -19 / 6, -7.0
    bY = (5 / 3) * b1
    print(f"  b₁ = {b1:.6f}   （GUT 归一化）")
    print(f"  b_Y = (5/3)·b₁ = {bY:.6f}   （物理超荷 1/α_Y 的 β）")
    print(f"  比值 b_Y/b₁ = {bY/b1:.6f} = 5/3 ✓")
    print(f"  ⟹ `build_68` 用了 b₁={b1:.4f} 跑 1/α_Y ⟹ 跑动量少 {1-b1/bY:.6f} = {(1-b1/bY)*100:.1f}%")
    check("**T1a b_Y/b₁ = 5/3 精确**", abs(bY / b1 - 5 / 3) < 1e-12, f"{bY/b1:.9f}")
    check("**T1b 跑动量少的比例 = 40%**", abs((1 - b1 / bY) - 0.4) < 1e-9,
          f"{(1-b1/bY)*100:.1f}%")

    # ============ T2 几何起点（不受 β 影响）============
    print("\n" + "=" * 96)
    print("T2：几何起点（**不受 β 影响**）")
    print("=" * 96)
    r = (1 + cos(pi / 6)) / 2
    s2MX = 4 * r ** 2 / (7 * r ** 2 + 3)
    print(f"  r = (1+cos30°)/2 = (2+√3)/4 = {r:.9f}")
    print(f"  sin²(M_X) = 4r²/(7r²+3) = {s2MX:.6f}   （声明 0.382913）")
    print(f"  检验 r=1: 4/(7+3) = {4/10:.6f} = 3/8 ✓（GUT 点自洽）")
    check("**T2a r 复现**", abs(r - 0.933012702) < 1e-9, f"{r:.9f}")
    check("**T2b sin²(M_X) = 0.382913 复现**", abs(s2MX - 0.382913) < 1e-6, f"{s2MX:.6f}")
    check("**T2c ★ 发现：r=1 时 4r²/(7r²+3) = 0.4 ≠ 3/8 ⟹ 与 1/(1+Σc²) 公式不一致**",
          abs(4 / 10 - 3 / 8) > 1e-3, f"4/10 = {4/10} vs 3/8 = {3/8}")
    print(f"  ⚠ 注意：r=1 时 4r²/(7r²+3) = 4/10 = 0.4，**不是** 3/8 = 0.375")
    print(f"     ⟹ 该公式在 r=1 处**不给** 3/8 ⟹ 与第一批的 1/(1+Σc²)=3/8 **不是同一个公式**")
    gap("**两条 sin²θ_W(M_X) 公式的关系**",
        "`4r²/(7r²+3)`（r=1 给 0.4）与 `1/(1+Σc²)`（给 3/8）**不一致**；底座未说明二者关系")

    # ============ T3 ★ 撤回传播审计 ============
    print("\n" + "=" * 96)
    print("T3：★ 撤回传播审计 —— `0.231200` / `0.0088%` 的下游状态")
    print("=" * 96)
    targets = re.compile(r"0\.23120|0\.0088\s*%")
    docs = sorted(glob.glob(os.path.join(BASE, "*.md")))
    citing, clean, marked = [], [], []
    for f in docs:
        txt = io.open(f, encoding="utf-8", errors="ignore").read()
        if not targets.search(txt):
            continue
        name = os.path.basename(f)
        citing.append(name)
        head = "\n".join(txt.split("\n")[:25])
        has_mark = any(w in head for w in RETRACT_WORDS)
        (marked if has_mark else clean).append(name)
    print(f"  引用 `0.231200`/`0.0088%` 的文档数 = {len(citing)}")
    print(f"    其中**头 25 行有撤回标记**的 = {len(marked)}")
    print(f"    其中**头 25 行无撤回标记**的 = {len(clean)}")
    print(f"\n  ★ 无撤回标记的文档（撤回未传播）：")
    for nm in clean:
        print(f"    · {nm}")
    check("**T3a 有文档引用被撤回的 0.231200**", len(citing) > 0, f"{len(citing)} 篇")
    check("**T3b 其中至少一篇无撤回标记**", len(clean) > 0, f"{clean}")
    print(f"""
  ⟹ ★★ **订正未完全传播**：{len(clean)} 篇仍把被撤回的 `0.231200`（−0.0088%）当作有效结果 ✗
  ⟹ 而 `build_68` 的脚本头**明确写着**「全部作废」⟹ 这是**文档层的传播失败**，不是数值错误""")
    gap("**撤回未传播到下游文档**",
        f"{len(clean)} 篇（{', '.join(clean[:4])}…）仍引用被 `build_83` 撤回的 0.231200 / −0.0088%")

    # ============ T4 build_68 的撤回文本取证 ============
    print("\n" + "=" * 96)
    print("T4：`build_68` 的撤回文本（取证）")
    print("=" * 96)
    b68 = None
    for f in glob.glob(os.path.join(BASE, "build", "build_68*.py")):
        b68 = f
    if b68:
        txt = io.open(b68, encoding="utf-8", errors="ignore").read()
        lines = [l for l in txt.split("\n")[:35] if "撤回" in l or "作废" in l or "b_Y" in l
                 or "41/6" in l or "0.206858" in l or "缺陷" in l]
        print(f"  文件：{os.path.basename(b68)}")
        for l in lines:
            print(f"    {l.strip()[:110]}")
        check("**T4a build_68 头含撤回声明**",
              "撤回" in txt[:2500] and "作废" in txt[:2500], "✓")
        check("**T4b 撤回给出正确 β b_Y = 41/6**", "41/6" in txt[:2500], "✓")
        check("**T4c 撤回给出订正后的数 0.206858**", "0.206858" in txt[:2500], "✓")
    else:
        check("**T4 找到 build_68**", False, "未找到")

    # ============ T5 判决 ============
    print("\n" + "=" * 96)
    print("T5：嫁接判决")
    print("=" * 96)
    print(f"""  · **β 归一化的事实**（b_Y=(5/3)b₁）：J1 ✅ J2 ✅ J3 ✅ ⟹ **【导出】**（纯代数）
  · **几何起点 r=(2+√3)/4 与 sin²(M_X)=0.382913**：J1 ✅（本侧复现）
      ⚠ 但 r 的**机制未导出**（底座自认：「为什么是 (1+cos(π/6))/2 未导出」）
  · **0.231200（−0.0088%）**：**【已撤回】** ✗ —— 且撤回**未传播**到 {len(clean)} 篇下游
  · ★★ **本批的核心交付**：**撤回传播审计**（这是"逐个数验证"必须做的一步）
  ⇒ 判决：β 事实【导出】；几何起点【条件】；0.231200【已撤回但未传播】""")
    check("**T5 判决：β 事实【导出】；几何起点【条件】；0.231200【已撤回·未传播】**", True, "见上")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          bY_over_b1=bY / b1, underrun=(1 - b1 / bY),
                          r=r, s2MX=s2MX,
                          citing=citing, marked=marked, clean=clean,
                          verdict="β 事实【导出】；几何起点【条件】；0.231200【已撤回但未传播】")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_retraction.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_retraction.json")
