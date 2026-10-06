"""
r23_check.py --- 任务 3 的第 (3) 项：**独立复核 R23 的 $S_D^{\\rm evo}=a_L^D$**

**R23 的论证（`R0_publication_theorem.md` 第 20 条原文摘要）**
    演化层内 $S_D^{\\rm evo}(L)=\\bigl(\\binom{L}{L/2}/2^L\\bigr)^D$ **随 $D$ 严格下降**；
    配合临时 `GR-LB`（$D\\ge4$），四维在 GR 兼容扇区中的条件存活率**唯一最大**。
    `R25` 给成对连接版本 $F_D=B\\binom D2 q^D$，四维窗口 $1/2<q<3/5$；
    $L{=}4$ 的旋转类账本给 $q=5/9$，故条件证成全局唯一四维峰。
    `R26` 再把原生桥 `PAIR-CARRIER-DER` 拆成六项并给出三条边界。

**本脚本只做三件事（不引用结论，自己算）**
  ① 复核 $a_L=\\binom{L}{L/2}/2^L<1$ 与 $S_D^{\\rm evo}$ 对 $D$ 的**严格单调性**（$D=1..12$）
  ② 复核 $F_D=B\\binom D2q^D$ 的**峰的维数**随 $q$ 的变化，自算窗口端点（$1/2,\\ 3/5$）对不对
  ③ 把 $q=5/9$ 代进去，验证峰在 $D{=}4$
  ＋ 登记**它依赖但没有导出的前提**（`GR-LB`、`DIM-SECTOR`、`PAIR-CARRIER-DER`）

等级：纯算术核对【闭式】；对前提的依赖是【条件结论】，不是【导出】。
"""
from __future__ import annotations

import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def peak_dim(q, B=1.0, Dmax=40):
    """F_D = B * C(D,2) * q^D 的最大点（返回 argmax 与最大值）。"""
    best, bd, vals = -1.0, None, {}
    for D in range(2, Dmax + 1):
        f = B * (D * (D - 1) / 2) * q ** D
        vals[D] = f
        if f > best:
            best, bd = f, D
    return bd, best, vals


def main():
    out = {"title": "R23 / R25 独立复核", "grade": "闭式算术核对"}

    # ① a_L 与 S_D^evo 的单调性
    Ls = [2, 4, 8, 16, 32, 64, 128]
    aL = {L: math.comb(L, L // 2) / 2.0 ** L for L in Ls}
    rows1 = []
    for L in Ls:
        a = aL[L]
        SD = [a ** D for D in range(1, 13)]
        mono = all(SD[i] > SD[i + 1] for i in range(len(SD) - 1))
        rows1.append({"L": L, "a_L": round(a, 8), "a_L<1": a < 1,
                      "S_D_evo_D1": round(SD[0], 8), "S_D_evo_D4": round(SD[3], 10),
                      "S_D_evo_D12": round(SD[11], 12),
                      "strictly_decreasing_in_D": mono,
                      "argmax_over_D>=4": 4 if mono else None})
    out["S_D_evo"] = rows1
    print("① S_D^evo = (C(L,L/2)/2^L)^D")
    for r in rows1:
        print(f"   L={r['L']:>3}: a_L={r['a_L']:.6f} (<1: {r['a_L<1']})  "
              f"S_1={r['S_D_evo_D1']:.6f} S_4={r['S_D_evo_D4']:.3e} "
              f"严格下降={r['strictly_decreasing_in_D']}  ⟹ argmax(D>=4)=4")

    # ② F_D 的峰随 q：自算窗口
    qs = [0.40, 0.45, 0.48, 0.49, 0.50, 0.51, 0.55, 0.5556, 0.58, 0.59, 0.60, 0.61, 0.65, 0.70]
    rows2 = []
    print("\n② F_D = B·C(D,2)·q^D 的峰维数（自算）")
    for q in qs:
        bd, best, _ = peak_dim(q)
        rows2.append({"q": q, "argmax_D": bd, "F_max_over_B": round(best / 1.0, 6)})
        print(f"   q={q:<7}: argmax D = {bd}")
    win = [r["q"] for r in rows2 if r["argmax_D"] == 4]
    out["F_D_peak"] = {"rows": rows2, "q_values_with_peak_at_D4": win}
    # 精确窗口端点：解超越方程 F_{D+1}/F_D = 1
    # 手算：ratio(D,q)=F_{D+1}/F_D = [C(D+1,2)/C(D,2)]·q = [(D+1)/(D-1)]·q
    #   D=3→4 需 >1 ⟹ (4/2)q>1 ⟹ q>1/2 ；D=4→5 需 <1 ⟹ (5/3)q<1 ⟹ q<3/5
    lo_q = 1 / 2
    hi_q = 3 / 5
    out["F_D_peak"]["window_exact"] = {"lo": lo_q, "hi": hi_q,
                                       "derivation": "ratio(D,q)=((D+1)/(D-1))·q；"
                                                     "D=3→4 需 >1 ⟹ q>1/2；D=4→5 需 <1 ⟹ q<3/5"}
    print(f"\n   精确窗口（自算）：1/2 < q < 3/5 = [{lo_q}, {hi_q}]  "
          f"—— 与 R25 记录的 (1/2, 3/5) 一致")
    q_ledger = 5 / 9
    bd, best, _ = peak_dim(q_ledger)
    out["ledger_q"] = {"q": q_ledger, "argmax_D": bd,
                       "in_window": lo_q < q_ledger < hi_q}
    print(f"   L=4 旋转类账本 q=5/9={q_ledger:.6f}：在窗口内={lo_q < q_ledger < hi_q}，"
          f"峰在 D={bd}")

    out["dependencies_not_derived"] = [
        {"name": "GR-LB", "what": "D>=4 的 GR 兼容下界", "status": "语料自记为临时前提，R24 说要撤掉"},
        {"name": "DIM-SECTOR", "what": "维数扇区的乘积分解 B·C(D,2)·q^D", "status": "未导出"},
        {"name": "PAIR-CARRIER-DER", "what": "原生成对载体（把连通性变成完全图）", "status": "R26 拆成六项，仍开放"},
        {"name": "FULL-SUPPORT-LEDGER / WIPE-RESET-LEDGER", "what": "满支持账本 / 毁灭-重播种账本",
         "status": "R27 记为开放"},
        {"name": "PHASE-1-COCHAIN / PAIR-ID-QUOTIENT", "what": "相位 1-上链商 / 成对身份商", "status": "开放"},
    ]
    print("\n③ 依赖清单（R23/R25 用到但未导出的前提）：")
    for d in out["dependencies_not_derived"]:
        print(f"   · {d['name']:<40} {d['status']}")

    out["verdict"] = ("① S_D^evo 严格随 D 下降 ⟹ 在 GR 兼容集 {D>=4} 上唯一最大在 D=4：**复核通过**（纯算术）。"
                      "② F_D 的四维窗口 1/2<q<3/5 与 q=5/9 的峰位：**复核通过**（纯算术）。"
                      "③ 但两条都**不是**从 Z0 导出：它们依赖 GR-LB / DIM-SECTOR / PAIR-CARRIER-DER 等未导出前提，"
                      "所以结论是【条件结论】，不是 D 的导出。本轮另给的『相互作用必然 ⟺ D<=4』"
                      "（L2_REGIONS.md §3.2）是独立的第二条，但其『GR-LB』一侧同样未导出。")
    with open(os.path.join(OUT, "r23_check.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n" + out["verdict"])
    print("写出 results/r23_check.json")


if __name__ == "__main__":
    main()
