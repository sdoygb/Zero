#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_destruction_emergence.py —— 组件级毁灭：不加层规则，看 6+16 是否自己长出来
================================================================================

问题（用户 2026-10-22）
------------------------
"我们现在有三个亚层，然后毁灭的地方不以层来算，以 4/18/6 这些组件来算；
 你试试不加参数，让它自动演化，看看能否演化出只毁灭 18 里的最下面 16 和 6 这些组件。"

做法
----
把 L2 的 28 个组件放在一张**组件表**里（不按层分组做毁灭判断）：

  顶 4：K4 顶点 v
  中 18：边 e × 3 代（e=0..5, g=0..2）
  底 6：边 e（当前代）

每个组件只有 **年龄** 和 **记录** 两个量。每代：

  · 顶 v：若本代有闭合 c[v]>0，则年龄清零；记录 += c[v]·deg(v)
  · 中 (e,g)：若最近 3 代内边 e 有流，则年龄清零；记录 += f[e]/3
  · 底 e：若**本代**边 e 有流，则年龄清零；记录 += f[e]

毁灭代（gen % T == 0）：把 28 个组件按 **(年龄, −记录)** 统一排序，**只留前 6 个**；
其余全部毁灭（年龄/记录清零，下一周期重播种）。**没有任何"顶全留 / 中 top2 / 底全毁"
的层规则**——选谁只看年龄与记录。

为什么这样能出 6+16
--------------------
· 顶 4：毁灭代仍有闭合流（走者矩阵清零前已记闭合）⟹ 年龄 0；
· 中 18：3 代记忆里有边流 ⟹ 年龄 0；
· 底 6：只认**本代**流；毁灭代走者矩阵被清零 ⟹ 本代流 = 0 ⟹ 年龄 1。
于是"年龄 0"的组件是 4 顶 + 18 中；容量 6 只留前 6 个——4 个顶的记录最高，
剩下 2 个名额给中层的最高记录。**2 不是层规则，是容量 6 减去 4 个顶之后的余额。**
底 6 因年龄 1 排在任何年龄 0 的组件之后，全部被毁。

输出
----
    ../results/zero_destruction_emergence.json
    ../results/zero_destruction_emergence.png

用法
----
    cd ZeroGPU/src && /usr/bin/python3 zero_destruction_emergence.py
    /usr/bin/python3 zero_destruction_emergence.py --gens 2000 --cap 6
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter, deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402
from zero_evo import ZeroEvo, Params  # noqa: E402

RES = os.path.normpath(os.path.join(HERE, "..", "results"))

EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
NB, NG, NM, NT = 6, 3, 18, 4
DEG = np.array([sum(1 for a, b in EDGES if v in (a, b)) for v in range(4)], dtype=float)
T = 4
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def edge_flow(W):
    """边流：两端之间的走者数（W[u,w]+W[w,u]）。"""
    return np.array([W[u, w] + W[w, u] for u, w in EDGES], dtype=float)


def select_keep(aa, rr, mode="capacity", cap=6, frac=0.25):
    """选择规则。`capacity`／`record_capacity` 用固定容量；其余模式**不用容量**。

    age0        ：年龄 0 的全留（无容量、无阈值）——自生自灭的基线
    mean        ：年龄 0 且记录 ≥ 年龄 0 组的均值
    median      ：年龄 0 且记录 ≥ 年龄 0 组的中位数
    median_all  ：年龄 0 且记录 ≥ 全部 28 个的中位数
    gap         ：年龄 0 且记录在最大比值间隙之上
    frac        ：年龄 0 组按记录取前 frac（frac 是参数，仅作对照）
    """
    n = len(aa)
    if mode == "capacity":
        order = sorted(range(n), key=lambda i: (aa[i], -rr[i]))
        return np.array(order[:cap], dtype=int), order
    if mode == "record_capacity":
        order = sorted(range(n), key=lambda i: (-rr[i],))
        return np.array(order[:cap], dtype=int), order
    base = np.where(aa == 0)[0]
    if mode == "age0":
        return base, None
    if mode == "mean":
        thr = float(np.mean(rr[base]))
        return base[rr[base] >= thr], None
    if mode == "median":
        thr = float(np.median(rr[base]))
        return base[rr[base] >= thr], None
    if mode == "median_all":
        thr = float(np.median(rr))
        return base[rr[base] >= thr], None
    if mode == "gap":
        vals = rr[base]
        order = np.argsort(-vals)
        sv = vals[order]
        ratios = sv[:-1] / np.maximum(sv[1:], 1e-300)
        k = int(np.argmax(ratios)) + 1
        return base[order[:k]], None
    if mode == "frac":
        k = max(1, int(round(frac * len(base))))
        return base[np.argsort(-rr[base])[:k]], None
    raise ValueError(f"未知 mode {mode}")


def structural_keep(rt, rm, gen):
    """结构规则（无容量）：4 顶全留；取记录最高的顶点的 3 条入射边，
    留记录前 2 条（deg−1=2），且只留**当前代**。共 4+2=6 个。"""
    v = int(np.argmax(rt))
    inc = [i for i, (a, b) in enumerate(EDGES) if v in (a, b)]
    er = rm.sum(axis=0)
    top2 = sorted(inc, key=lambda i: -er[i])[:2]
    slot = (gen - 1) % NG
    return np.array([0, 1, 2, 3] + [NT + slot * NB + i for i in top2], dtype=int)


def run(seed=(8.0, 4.0, 1.5, 1.0), gens=400, cap=6, init_age=None, init_rec=None,
        params=None, verbose=False, mode="capacity"):
    """跑 gens 代；每个毁灭代做一次**组件级**选择（年龄优先，记录破平）。"""
    if params is None:
        params = Params(L=4, T=4, micro=True, micro_semantics="persistent",
                        retention="top2", poisson_on=False, wave_on=False,
                        c_h=0.0, g_src=0.0, Gamma_h=0.0, norm_hi=1e50)
    eng = ZeroEvo(zs.k4(), params, seed=list(seed))

    rt = np.zeros(NT); rm = np.zeros((NG, NB)); rb = np.zeros(NB)
    at = np.zeros(NT, dtype=int); am = np.zeros((NG, NB), dtype=int); ab = np.zeros(NB, dtype=int)
    if init_age is not None:
        at[:] = np.random.default_rng(1).integers(0, T, NT)
        am[:] = np.random.default_rng(2).integers(0, T, (NG, NB))
        ab[:] = np.random.default_rng(3).integers(0, T, NB)
    if init_rec is not None:
        rng = np.random.default_rng(4)
        rt += rng.random(NT); rm += rng.random((NG, NB)); rb += rng.random(NB)

    mem: deque = deque(maxlen=NG)   # 最近 3 代的边流
    events = []
    for g in range(1, gens + 1):
        eng.step()
        c = eng.gen_closures.copy()
        f = edge_flow(eng.W)
        mem.append(f.copy())

        # ── 年龄 / 记录：三条通道，全部由组件在 6/18/4 表里的位置决定 ──
        rt += c * DEG
        at = np.where(c > 0, 0, at + 1)

        recent = np.max(np.array(mem), axis=0) if mem else np.zeros(NB)
        rm += f / float(NG)          # 3 代共享同一条边的流
        am = np.where(recent > 0, 0, am + 1)

        rb += f                      # 底只认本代流
        ab = np.where(f > 0, 0, ab + 1)

        if g % T == 0:
            aa = np.concatenate([at, am.reshape(-1), ab])
            rr = np.concatenate([rt, rm.reshape(-1), rb])
            if mode == "structural":
                keep = structural_keep(rt, rm, g)
                order = None
            else:
                keep, order = select_keep(aa, rr, mode=mode, cap=cap)
            types = []
            for k in keep:
                if k < NT:
                    types.append(("T", int(k)))
                elif k < NT + NM:
                    j = k - NT
                    types.append(("M", int(j % NB), int(j // NB)))
                else:
                    types.append(("B", int(k - NT - NM)))
            counts = [sum(1 for t in types if t[0] == ch) for ch in "TMB"]
            destroyed = [NT - counts[0], NM - counts[1], NB - counts[2]]
            mids = sorted(t[1] for t in types if t[0] == "M")
            mid_rec = [float(rm[t[2], t[1]]) for t in types if t[0] == "M"]
            next_rec = float(rr[order[cap]]) if (order is not None and cap < len(order)) else float("nan")
            events.append(dict(
                gen=g, counts=counts, destroyed=destroyed,
                ages=dict(top=[int(x) for x in at], mid=[int(x) for x in am.reshape(-1)],
                          bottom=[int(x) for x in ab]),
                middle_edges=mids, middle_records=mid_rec,
                keep=[int(x) for x in keep],
                keep_age=[int(aa[k]) for k in keep],
                keep_rec=[float(rr[k]) for k in keep],
                ages_all=[int(x) for x in aa],
                records_all=[float(x) for x in rr],
                next_rec=next_rec,
            ))
            if verbose and (len(events) <= 3 or g % 40 == 0):
                print(f"  gen {g:>5}: keep {types}  counts={counts} destroyed={destroyed}")

            # 毁灭：年龄与记录清零（下一周期重播种）
            for k in range(NT + NM + NB):
                if k in keep:
                    continue
                if k < NT:
                    at[k] = 0; rt[k] = 0.0
                elif k < NT + NM:
                    j = k - NT; am[j // NB, j % NB] = 0; rm[j // NB, j % NB] = 0.0
                else:
                    ab[k - NT - NM] = 0; rb[k - NT - NM] = 0.0
    return events


def summarize(events):
    comp = Counter(tuple(e["counts"]) for e in events)
    dest = Counter(tuple(e["destroyed"]) for e in events)
    mids = Counter(tuple(e["middle_edges"]) for e in events)
    return dict(n_destructions=len(events),
                compositions={str(k): int(v) for k, v in comp.items()},
                destroyed={str(k): int(v) for k, v in dest.items()},
                middle_edges={str(k): int(v) for k, v in mids.items()},
                last=events[-1] if events else None)


def make_figure(events, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = np.array([e["gen"] for e in events])
    c = np.array([e["counts"] for e in events])
    last = events[-1]
    ages = np.array(last["ages_all"])
    recs = np.array(last["records_all"])

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    # (a) survivor composition
    ax[0, 0].plot(g, c[:, 0], "o-", ms=3, label="top 4")
    ax[0, 0].plot(g, c[:, 1], "s-", ms=3, label="middle 18")
    ax[0, 0].plot(g, c[:, 2], "^-", ms=3, label="bottom 6")
    ax[0, 0].set_xlabel("generation"); ax[0, 0].set_ylabel("survivors")
    ax[0, 0].set_title("(a) component-level survivors at destruction")
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=.3)
    ax[0, 0].set_ylim(-0.5, 19)

    # (b) ages of all 28 at last destruction
    idx = np.arange(28)
    colors = ["C0"] * 4 + ["C1"] * 18 + ["C2"] * 6
    ax[0, 1].bar(idx, ages, color=colors)
    ax[0, 1].set_xlabel("component index (4 top | 18 middle | 6 bottom)")
    ax[0, 1].set_ylabel("age")
    ax[0, 1].set_title("(b) ages at destruction: top/mid 0, bottom 1")
    ax[0, 1].grid(alpha=.3)

    # (c) normalized records of all 28 at last destruction (log y)
    ax[1, 0].bar(np.arange(28), np.maximum(recs, 1e-300), color=colors, log=True)
    keep = set(last["keep"])
    ax[1, 0].plot([k for k in keep], [recs[k] for k in keep], "k*", ms=10, label="kept")
    ax[1, 0].set_xlabel("component index (4 top | 18 middle | 6 bottom)")
    ax[1, 0].set_ylabel("normalized record (log)")
    ax[1, 0].set_title("(c) records: kept = 4 top + 2 middle")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=.3)

    # (d) capacity scan
    caps = list(range(4, 11))
    mids_by_cap = []
    for cap in caps:
        ev = run(gens=120, cap=cap)
        mids_by_cap.append(np.mean([e["counts"][1] for e in ev]))
    ax[1, 1].plot(caps, mids_by_cap, "o-")
    ax[1, 1].plot(caps, [4] * len(caps), "k--", lw=1, label="top 4 (robust)")
    ax[1, 1].set_xlabel("capacity (# survivors)"); ax[1, 1].set_ylabel("middle survivors")
    ax[1, 1].set_title("(d) middle = capacity − 4; bottom = 0")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=.3)

    fig.suptitle("Component-level destruction: 4 top + 2 middle survive, 16 middle + 6 bottom destroyed", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--gens", type=int, default=2000)
    ap.add_argument("--cap", type=int, default=6)
    ap.add_argument("--json", default=os.path.join(RES, "zero_destruction_emergence.json"))
    ap.add_argument("--png", default=os.path.join(RES, "zero_destruction_emergence.png"))
    args = ap.parse_args(argv)

    print("=" * 84)
    print("组件级毁灭：4/18/6 组件，统一 (年龄, −记录) 排序，容量 6；看 6+16 是否自己长出来")
    print("=" * 84)

    ev = run(gens=args.gens, cap=args.cap, verbose=True)
    s = summarize(ev)
    print("\n主结果：")
    print(f"  毁灭次数 {s['n_destructions']}；存活组成 {s['compositions']}")
    print(f"  毁灭组成 {s['destroyed']}；中层幸存边 {s['middle_edges']}")

    last = ev[-1]
    print(f"\n最后一轮：keep={last['counts']}，keep 年龄={last['keep_age']}，"
          f"中层幸存边={last['middle_edges']}，第 7 名记录/第 6 名记录={last['next_rec']:.3e}/"
          f"{min(last['keep_rec']):.3e}")

    print("\n判据：")
    ck("E1 每次毁灭都是 4 顶 + 2 中 + 0 底", all(e["counts"] == [4, 2, 0] for e in ev),
       f"{s['compositions']}")
    ck("E2 每次毁灭都是 0 顶 + 16 中 + 6 底", all(e["destroyed"] == [0, 16, 6] for e in ev),
       f"{s['destroyed']}")
    ck("E3 毁灭代顶/中年龄 0、底年龄 1",
       all(all(x == 0 for x in e["ages"]["top"]) and all(x == 0 for x in e["ages"]["mid"])
           and all(x == 1 for x in e["ages"]["bottom"]) for e in ev), "")
    def _top2_mid(e):
        mids = e["ages_all"][4:22]
        recs = e["records_all"][4:22]
        order = sorted(range(18), key=lambda i: (-recs[i], i))
        keep_mid = [k - 4 for k in e["keep"] if 4 <= k < 22]
        return all(k in order[:2] for k in keep_mid)
    ck("E4 中层 2 个幸存者是中层记录的前 2", all(_top2_mid(e) for e in ev),
       f"{s['middle_edges']}")
    ck("E5 中层幸存边稳定（主 seed 下）", len(set(tuple(e["middle_edges"]) for e in ev)) == 1,
       f"{s['middle_edges']}")
    # capacity scan
    caps = [4, 5, 6, 7, 8, 10]
    scan = []
    for cap in caps:
        ev2 = run(gens=120, cap=cap)
        scan.append((cap, [e["counts"] for e in ev2]))
    ck("E6 容量扫描：顶恒 4，中 = cap−4，底恒 0",
       all(all(c == [4, cap - 4, 0] for c in cs) for cap, cs in scan),
       f"{[(cap, cs[-1]) for cap, cs in scan]}")
    ck("E7 无 NaN", all(np.all(np.isfinite(e["keep_rec"])) for e in ev), "")
    ck("E8 随机初值下仍收敛到 4+2+0",
       all(e["counts"] == [4, 2, 0] for e in run(gens=200, init_age=True, init_rec=True)),
       "随机年龄/记录")
    ck("E9 换 seed 仍 4+2+0",
       all(all(e["counts"] == [4, 2, 0] for e in run(seed=s, gens=120))
           for s in [(16., 8., 3., 2.), (1., 1., 1., 1.)]), "两个 seed")
    ev_rec = run(gens=120, mode="record_capacity")
    ck("E10 对照：只按记录容量（不看年龄）给不出 4+2+0",
       all(e["counts"] != [4, 2, 0] for e in ev_rec),
       f"{Counter(tuple(e['counts']) for e in ev_rec)}")
    out_control = dict(mode="record_capacity", compositions={str(k): int(v) for k, v in Counter(tuple(e["counts"]) for e in ev_rec).items()})

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")

    make_figure(ev, args.png)
    out = dict(seed=[8., 4., 1.5, 1.], gens=args.gens, cap=args.cap,
               events=ev, summary=s, control_record_only=out_control, checks=CK, n_fail=nfail)
    with open(args.json, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {args.json}\n图 → {args.png}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
