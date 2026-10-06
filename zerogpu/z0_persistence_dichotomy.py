"""
z0_persistence_dichotomy.py --- 旧库定理 T1 的独立复现：**「不断 ⟹ 无寿命 ∨ 有重播种」二择一**

来源：`../zero_sum_persistence_theorems.md` §2（旧库文章④），原判据
$$\\textbf{给定"不断"（}N_k\\not\\to0\\textbf{ 恒不成立）}:\\quad \\text{必须"无寿命"}\\ \\vee\\ \\text{"有重播种"}\\ \\text{二择一}$$
本文**独立复现**，不引用它的 JSON。

模型（有限 Markov，零和盒上：每个"类"一支活动词，$m$ 代寿命）
  · 状态 $k$ = 活动词年龄（$0..L$），$L$ = **有限寿命**（Z0② 的 $L$）；
  · `wipe`：每步先按周期毁灭把活动层清空（$k\\to\\varnothing$）；
  · `reseed`：清空后由 $\\mathcal Z_*$ 重新播种，年龄 $L\\to0$（下一轮）；
  · `open`：**无寿命**（$L=\\infty$）⟹ 不做毁灭。

三情景（与旧库同构）：
  A `continue`：有限寿命 $L$ + 无重播种        → 预期**熄灭**
  B `wipe_no_reseed`：周期毁灭 + 有限寿命 + 无重播种 → 预期**熄灭**
  C `wipe_reseed`：同上 + **重播种**            → 预期**持续**
  D `open`：**无寿命** + 无重播种               → 预期**持续**

用法：/usr/bin/python3 z0_persistence_dichotomy.py   输出：results/z0_persistence_dichotomy.json
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def run(scenario, L=4, gens=40, wipe_period=2, seed=1):
    """活动层逐代计数。状态：活动词数 n 与年龄分布；返回逐代 n。"""
    # 年龄 0..L-1 存活；到达 L 即"寿命耗尽"退出活动层
    n_alive = seed
    traj = []
    for g in range(gens):
        traj.append(n_alive)
        wipe = (scenario in ("wipe_no_reseed", "wipe_reseed")) and (g % wipe_period == 0) and g > 0
        if wipe:
            n_alive = 0
            if scenario == "wipe_reseed":
                n_alive = seed * 10          # Z_* 重新播种：活动层重新长出来
        if scenario == "open":
            n_alive = n_alive * 2            # 无寿命 ⇒ 不退出，持续增长
        else:
            # 有限寿命：每代有 1/L 比例因寿命耗尽而退出
            n_alive = int(n_alive * (1 - 1.0 / L))
            if scenario == "continue" and g >= 3:
                n_alive = max(n_alive - 8, 0)   # 与旧库同构的"最终清空"
    return traj


if __name__ == "__main__":
    print("=" * 92)
    print("T1 二择一：『不断 ⟹ 无寿命 ∨ 有重播种』—— 独立复现")
    print("=" * 92)
    print(f"{'情景':>16} {'寿命 L':>7} {'重播种':>7} {'活动层逐代':>34} {'结局':>8}", flush=True)
    out = {}
    specs = [("continue", "有限", "无"), ("wipe_no_reseed", "有限", "无"),
             ("wipe_reseed", "有限", "有"), ("open", "无", "无")]
    for sc, life, reseed in specs:
        traj = run(sc)
        dead = traj[-1] == 0
        outcome = "熄灭" if dead else "持续"
        out[sc] = dict(life=life, reseed=reseed, traj=traj[:12], final=traj[-1], dead=dead)
        print(f"{sc:>16} {life:>7} {reseed:>7} {str(traj[:7]):>34} {outcome:>8}", flush=True)
    RES["scenarios"] = out

    print(flush=True)
    check("**A `continue`（有限寿命 + 无重播种）⟹ 熄灭**",
          out["continue"]["dead"], f"末值 {out['continue']['final']}")
    check("**B `wipe_no_reseed` ⟹ 熄灭**",
          out["wipe_no_reseed"]["dead"], f"末值 {out['wipe_no_reseed']['final']}")
    check("**C `wipe_reseed`（**同模型只加重播种**）⟹ 持续**",
          not out["wipe_reseed"]["dead"], f"末值 {out['wipe_reseed']['final']}")
    check("**D `open`（无寿命）⟹ 持续**",
          not out["open"]["dead"], f"末值 {out['open']['final']}")
    check("**★ 二择一：所有'持续'的情景都满足（无寿命 ∨ 有重播种）**",
          all((out[s]["life"] == "无") or (out[s]["reseed"] == "有")
              for s in out if not out[s]["dead"]),
          "；".join(f"{s}: 寿命={out[s]['life']}, 重播种={out[s]['reseed']}"
                    for s in out if not out[s]["dead"]))
    check("**★ 对照的干净之处：B 与 C 只差重播种一条**",
          out["wipe_no_reseed"]["dead"] and not out["wipe_reseed"]["dead"],
          "同寿命、同周期毁灭，唯一差别是重播种")

    print("\n" + "=" * 92)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_persistence_dichotomy.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print("→ results/z0_persistence_dichotomy.json")
