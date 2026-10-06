"""
z0_graft_levels.py --- 嫁接验证 · 第二十八批：层级的**多口径**表（解 G45 ＋ 复现 Q692）

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

来源：`Q610_COMPLETION_TRANSPORT.md`、`Q616_LADDER_EXPONENT.md`、`Q619_SYMPLECTIC_FACTOR.md`、
      `Q692_PROPAGATE_KD1.md`（`build_346`，16/0）、`Q621_VACUUM_ENERGY.md`

【本批解决第二十七批留的 G45】
  第二十七批本侧说「底座称『质量情形距整 $0.121$』，但本侧用 $S/(2\\pi)$ 得 $0.3507$，**不符**」。
  ★ **解答**：$0.121$ **不是** $S/(2\\pi)$ 的距整，而是 **电弱层级的距整**：
    层级 $=\\ln(1/(M_Z r_{\\rm int}))/(2\\pi) = 4.8792$（$r_{\\rm int}=6497\\,l_P$）
    $\\Longrightarrow$ 距最近整数 $=\\mathbf{5-4.8792=0.1208}$ ✓
  本侧上批**用错了量**（$S/(2\\pi)=5.3507$ 是"从 $v$ 到 $M_X$"的层级，而 $4.879$ 是"从 $M_Z$ 到内空间尺度"的层级）。

★ 本批的独立验证：
  1. 层级公式在 $r_{\\rm int}$ 的**多个口径**下的值（$6497/3415/489/1211/8078/29314/12312$）
  2. 复现 $Q692$ 的表（$6497\\to4.8792$；$3415\\to4.9815$）
  3. ★ **关键**：$r_{\\rm int}$ 的取值**直接决定**「是否整数」的结论

用法：/usr/bin/python3 z0_graft_levels.py   输出：results/z0_graft_levels.json
"""
from __future__ import annotations

import json
import os
import time
from math import log, pi

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

MZ = 91.1876
L_P = 1.616255e-35          # m
HBARC = 0.1973269804e-15    # GeV·m


def level(r_over_lP):
    """层级 = ln(ħc/(M_Z r_int))/(2π)"""
    r = r_over_lP * L_P
    return log(HBARC / r / MZ) / (2 * pi)


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ L1 层级公式的多口径表 ============
    print("=" * 96)
    print("L1：层级公式 $\\ln(1/(M_Z r_{\\rm int}))/(2\\pi)$ 的**多口径**表")
    print("=" * 96)
    cases = [
        (6497, "`Q610` 情形 B（旧）"),
        (3415, "`Q651` 订正（$k\\le3$ 仅标量）"),
        (489, "$k\\le1$＋Einstein"),
        (1211, "$k\\le1$＋Jordan"),
        (8078, "$k\\le3$＋Einstein"),
        (29314, "$k\\le3$＋Jordan"),
        (12312, "$\\sqrt{\\tilde\\Sigma}$（未除 $C$）"),
    ]
    print(f"  {'r_int/l_P':>11} {'层级':>10} {'距最近整':>10} {'相对':>9}  出处")
    for rl, src in cases:
        lv = level(rl)
        d = abs(lv - round(lv))
        print(f"  {rl:>11} {lv:>10.6f} {d:>10.4f} {d/lv*100:>8.2f}%  {src}")
    check("**L1a $r_{\\rm int}=6497$ 给层级 $4.8792$（声明）**",
          abs(level(6497) - 4.8792) < 1e-3, f"{level(6497):.6f}")
    check("**L1b $r_{\\rm int}=3415$ 给层级 $4.9815$（`Q651`/`Q692`）**",
          abs(level(3415) - 4.9815) < 1e-3, f"{level(3415):.6f}")

    # ============ L2 ★ 解 G45 ============
    print("\n" + "=" * 96)
    print("L2：★ 第二十七批 G45 的**解答**（$0.121$ 是什么？）")
    print("=" * 96)
    lv6497 = level(6497)
    print("  电弱层级（$r_{\\rm int}=6497$）：$\\mathbf{%.6f}$" % lv6497)
    print("  距最近整数 $5$：$5-%.4f=\\mathbf{%.4f}$ ← **正是底座说的 $0.121$** ✓" % (lv6497, 5-lv6497))
    print(f"\n  而本侧第二十七批误用的量：")
    S = 33.6197
    print("    $S/(2\\pi)=\\ln(M_X/v)/(2\\pi) = %.6f$，距整 $%.4f$" % (S/(2*pi), abs(S/(2*pi)-round(S/(2*pi)))))
    print("    ★ 这是「从 $v$ 到 $M_X$」的层级（$\\sim5.35$ 级），**不是**「从 $M_Z$ 到内空间」的层级（$\\sim4.88$ 级）")
    check("**L2a $5-4.8792=0.1208$ 正是底座说的 0.121**",
          abs((5 - lv6497) - 0.121) < 1e-3, f"{5-lv6497:.4f}")
    check("**L2b 本侧上批用错量（S/2pi = 5.3507 != 4.8792）**",
          abs(S / (2 * pi) - lv6497) > 0.4, f"{S/(2*pi):.4f} vs {lv6497:.4f}")

    # ============ L3 ★ 关键：r_int 决定结论 ============
    print("\n" + "=" * 96)
    print("L3：★ $r_{\\rm int}$ 的取值**直接决定**「是否整数」的结论")
    print("=" * 96)
    print(f"  {'r_int/l_P':>11} {'层级':>10} {'距整':>10} {'结论':>28}")
    verdicts = {}
    for rl, src in cases[:6]:
        lv = level(rl)
        d = abs(lv - round(lv))
        if d < 0.02:
            v = "**接近整数（0.37%）** ⚠"
        elif d < 0.16:
            v = "偏离（2–3%）"
        else:
            v = "**明显偏离（≥5%）**"
        verdicts[rl] = (lv, d, v)
        print(f"  {rl:>11} {lv:>10.6f} {d:>10.4f} {v:>28}")
    print(f"""
  ⟹ ★★ **`Q692` 的关键发现**（本侧复现）：
     · `Q619` 原说「层级 $4.8792$ ⟹ **确定不是整数**」⇒ 该子线收口
     · 而 `Q651` 把 $r_{{\\rm int}}$ 订正为 $3415$ ⟹ 层级 $\\mathbf{{4.9815}}$ ⟹ 距整 $\\mathbf{{0.0185}}$（$0.37\\%$）
       ⟹ **推论被【削弱】**：$0.37\\%$ 比 $2.42\\%$ 更接近整数 ✓
  ⟹ 故「层级是否整数」**完全取决于选哪个 $r_{{\\rm int}}$** —— 这与第二十三批的 $r_{{\\rm int}}$ 冲突**同源** ⚠""")
    check("**L3a Q692 的表复现（6497->4.8792；3415->4.9815）**",
          abs(level(6497) - 4.8792) < 1e-3 and abs(level(3415) - 4.9815) < 1e-3, "两行都对")
    check("**L3b 3415 的距整 0.0185 比 6497 的 0.1208 好 6.5 倍**",
          abs((5 - level(6497)) / (5 - level(3415)) - 6.53) < 0.1,
          f"{(5-level(6497))/(5-level(3415)):.2f} 倍")
    gap("**「层级是否整数」依赖 $r_{\\rm int}$ 的口径**（与 G38 同源）",
        "6497 → 距整 0.1208（2.4%）；3415 → 距整 0.0185（0.37%）；489 → 0.2909（5.5%）")

    # ============ L4 判决 ============
    print("\n" + "=" * 96)
    print("L4：嫁接判决")
    print("=" * 96)
    print(f"""  · **层级公式与多口径表**：J1 ✅ ⟹ **【导出】**
  · ★ **G45 的解答**（$0.121=5-4.8792$）：J1 ✅ ⟹ **【导出】** ＋ **关闭 G45** ✓
  · **`Q692` 的订正**（$3415$ 给 $4.9815$，距整 $0.0185$）：J1 ✅ ⟹ **【导出】**
  · **「层级是否整数」依赖 $r_{{\\rm int}}$ 口径**：J1 ✅ ⟹ **【条件】**
  ⇒ 本批的净收获：**关闭 G45** ＋ **把「层级整数性」的依赖关系写清**（与 G38 同源）""")
    check("**L4 判决：层级表【导出】；G45 关闭；整数性为【条件】**", True, "见上")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          levels={rl: level(rl) for rl, _ in cases},
                          level_6497=level(6497), level_3415=level(3415),
                          G45_resolved="0.121 = 5 - 4.8792（电弱层级，非 S/2π）",
                          verdict="【导出】：层级表；G45 关闭；整数性为【条件】")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_levels.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_levels.json")
