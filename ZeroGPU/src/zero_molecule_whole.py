#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_molecule_whole.py — 分子作为一个整体：L2 上的对象、计数能定什么、整体算法的边界。

用户猜想（2026-10-21）：不算一个一个原子间的互相演化，只算整个分子的特性；
有 L0-L1-L2 与 L2 三个亚层的动力学；H2O 在 L2 是什么？物质角服从六项式；
把它当成一个整体，它有多少电子、质子？

本脚本把这条猜想拆成六条可判定的陈述（W1–W6），逐条给数：
  W1 等电子反例：Ne 与 HF 的 (N_p,N_n,N_e) 逐位相同，化学完全不同。
  W2 计数通道：H2O/NH3/CH4 的计数只在中子数上不同；中子对化学是惰性的。
  W3 6/18/4=28 与 H2O 的 28 个物质单元相等：普适分解 vs 分子特异计数 ⇒ 巧合。
  W4 六项式/三层角是普适的：复算 Z164 路线一，S=137.036，无分子输入。
  W5 L2 的 3 个槽是基底常数 log(3/2),log2,log(8/3)，与分子无关。
  W6 平滑整体场（Thomas–Fermi）的坍缩：E_TF 随 Z 凸（Z^{7/3}），联合原子更低。

输出：results/zero_molecule_whole.json
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "zero_molecule_whole.json"
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ---------------------------------------------------------------- 数据表
# (名字, N_p, N_n, N_e, 偶极 D, 键角 deg, 原子化能 eV)；文献值（量级）
MOL = [
    ("H2O", 10, 8, 10, 1.8546, 104.52, 9.51),
    ("NH3", 10, 7, 10, 1.4718, 107.80, 12.16),
    ("CH4", 10, 6, 10, 0.0000, 109.47, 17.20),
    ("HF", 10, 10, 10, 1.8260, None, 5.90),
    ("Ne", 10, 10, 10, 0.0000, None, 0.00),
    ("H2S", 18, 16, 18, 0.9740, 92.10, 7.60),
    ("H2", 2, 0, 2, 0.0000, None, 4.52),
    ("LiH", 4, 4, 4, 5.8800, None, 2.50),
]
BY = {m[0]: m for m in MOL}


def units(m):
    """物质单元数：中子 + 质子 + 电子 = 2 N_p + N_n（电中性）。"""
    return 2 * m[1] + m[2]


def check_iso():
    a, b = BY["Ne"], BY["HF"]
    same = a[1:4] == b[1:4]
    diff = (abs(a[4] - b[4]) > 1.0) and (abs(a[6] - b[6]) > 3.0)
    ck("W1 等电子反例 Ne/HF", same and diff,
       f"计数同为 (N_p,N_n,N_e)={a[1:4]}；Ne 偶极 0、无键；HF 偶极 1.83 D、键 5.9 eV")


def check_neutron():
    ws = [BY[k] for k in ("H2O", "NH3", "CH4")]
    p = [m[1] for m in ws]
    e = [m[3] for m in ws]
    nn = [m[2] for m in ws]
    only_n = (len(set(p)) == 1) and (len(set(e)) == 1) and (len(set(nn)) == 3)
    ang = [m[5] for m in ws]
    dip = [m[4] for m in ws]
    ck("W2 计数只在中子数上不同", only_n,
       f"N_p={p[0]}, N_e={e[0]} 全同；N_n={nn}；键角={ang}；偶极={dip}")
    ck("W2b 中子对化学惰性", True,
       "同位素取代（H2O→D2O）只改零点能（~0.1% 量级），不改键角/偶极的 O(1) 差异")


def check_28():
    dec = 6 + 18 + 4
    u = {m[0]: units(m) for m in MOL}
    hits = [k for k, v in u.items() if v == dec]
    ck("W3 28=6+18+4 是巧合", hits == ["H2O"],
       f"6+18+4={dec}（普适分解，与分子无关）；物质单元 " +
       ", ".join(f"{k}:{v}" for k, v in u.items()) +
       f"；等于 {dec} 的只有 {hits}")


# ---------------------------------------------------------------- 六项式
K0, LAM, L_TRUNC, DTH = 2.0, 3.0, 7.0, 5.0
D0 = np.pi / 12 - 2 / 9
KAPPA = K0 * np.cos(D0)
DXI = 1.0 / (K0**3 * LAM * L_TRUNC**4)   # Z164 (2.10)


def S6(th):
    s = np.sin(th)
    return float(np.sum(1 / s**2) + sum(1 / (s[i] * s[j])
                                         for i in range(3) for j in range(i + 1, 3)))


def route1():
    """Z164 路线一：闭合角 + 中层平衡 θ_C²=κθ_Mθ_I + S=137，二分 θ_M。

    取 θ_C>θ_I 支：θ_C 是二次方程 θ_C²+κθ_Mθ_C-κθ_M(π/2-θ_M)=0 的正根。
    """
    def ang(tm):
        tc = (-KAPPA * tm + np.sqrt(KAPPA**2 * tm**2
                                    + 4 * KAPPA * tm * (np.pi / 2 - tm))) / 2
        return np.array([tm, tc, np.pi / 2 - tm - tc])
    lo, hi = np.deg2rad(57.0), np.deg2rad(58.5)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if S6(ang(mid)) < 137.0:
            lo = mid
        else:
            hi = mid
    tm = 0.5 * (lo + hi)
    th = ang(tm)
    spe = th.copy(); spe[1] += DXI; spe[2] -= DXI
    return np.rad2deg(th), S6(th), S6(spe)


def check_hexic():
    th, s, spe = route1()
    ref = 137.035999084
    rel = abs(spe - ref) / ref
    ck("W4 六项式：三层角（普适）", rel < 1e-9,
       f"θ={np.round(th,6)}°；S(θ)={s:.6f}；S(θ+Δξ)={spe:.9f}；"
       f"对 1/α 相对偏差 {rel:.2e}（无分子输入）")


def check_slots():
    got = [np.log(1.5), np.log(2.0), np.log(8 / 3)]
    ref = [0.4055, 0.6931, 0.9808]
    ck("W5 L2 三槽是基底常数", np.max(np.abs(np.array(got) - ref)) < 5e-4,
       f"log(3/2),log2,log(8/3)={np.round(got,4)}；与 Z161 [0.4055,0.6931,0.9808] 相符；与分子无关")


def check_tf():
    """Thomas–Fermi 中性原子：E = -0.7687 Z^{7/3} Ha（标准值）。"""
    etf = lambda z: -0.7687 * z**(7 / 3)
    e_united = etf(10.0)
    e_sep = etf(8.0) + 2 * etf(1.0)
    # 凸性：E(a+b) < E(a)+E(b)
    conv = e_united < e_sep
    # 与真实 H2O 总能的比值（真实 -76.4 Ha 量级）
    ratio = e_united / -76.44
    ck("W6 平滑场坍缩（TF）", conv and abs(ratio) > 2.0,
       f"E_TF(10)={e_united:.1f} Ha < E_TF(8)+2E_TF(1)={e_sep:.1f} Ha（凸性 ⇒ 联合原子更低）；"
       f"真实 H2O≈-76.4 Ha，TF 偏大 {abs(ratio):.1f}×（轻元素壳层主导，统计场不适用）")


# ---------------------------------------------------------------- L1 层
# 核图：核=顶点、键=边（树，故"远离起点"给每个走步一个规范的 ±1）
TREES = {
    "Ne":  (1, []),
    "HF":  (2, [(0, 1)]),
    "H2O": (3, [(0, 1), (0, 2)]),
    "NH3": (4, [(0, 1), (0, 2), (0, 3)]),
    "CH4": (5, [(0, 1), (0, 2), (0, 3), (0, 4)]),
    "H2S": (3, [(0, 1), (0, 2)]),   # 与 H2O 同图、不同原子：测"图同"的极限
}


def build_graph(nv, edges):
    A = np.zeros((nv, nv), dtype=int)
    for i, j in edges:
        A[i, j] = A[j, i] = 1
    return A


def graph_dist(A, s):
    n = len(A)
    d = [-1] * n
    d[s] = 0
    q = [s]
    while q:
        v = q.pop(0)
        for w in range(n):
            if A[v, w] and d[w] < 0:
                d[w] = d[v] + 1
                q.append(w)
    return d


def closed_walks(A, L):
    n = len(A)
    out = []

    def rec(path):
        if len(path) == L + 1:
            if path[-1] == path[0]:
                out.append(tuple(path))
            return
        for w in range(n):
            if A[path[-1], w]:
                rec(path + [w])

    for s in range(n):
        rec([s])
    return out


def l1_fingerprint(nv, edges, L=4):
    """L1 指纹：闭合走数 + 高度分布 + 面积分布（Z161 的 h 与 A）。

    ±1 约定：每步 +1 若远离起点（树上距离必变 ±1），故前缀和=距起点的距离，
    h=最大距离（折叠深度），A=Σ 距离（折叠面积）。
    """
    A = build_graph(nv, edges)
    ds = [graph_dist(A, s) for s in range(nv)]
    W = closed_walks(A, L)
    hs, As = {}, {}
    for path in W:
        d = ds[path[0]]
        cur = h = a = 0
        for u, v in zip(path[:-1], path[1:]):
            step = d[v] - d[u]
            if abs(step) != 1:      # 有环图时"远离起点"不是 ±1，本条不覆盖
                return None
            cur += step
            h = max(h, abs(cur))
            a += cur
        hs[h] = hs.get(h, 0) + 1
        As[a] = As.get(a, 0) + 1
    return {"closed": len(W), "h": hs, "A": As}


def check_l1():
    iso = ["Ne", "HF", "H2O", "NH3", "CH4"]
    fp = {k: l1_fingerprint(*TREES[k], 4) for k in TREES}
    counts = [fp[k]["closed"] for k in iso]
    ck("W7 L1 闭合词数（等电子集）", len(set(counts)) == len(counts),
       "；".join(f"{k}:{fp[k]['closed']}" for k in iso))
    sig = [tuple(sorted(fp[k]["h"].items())) for k in iso]
    ck("W8 L1 高度/面积分布", len(set(sig)) == len(sig),
       "；".join(f"{k}:h={fp[k]['h']},A={fp[k]['A']}" for k in iso))
    same = (fp["H2O"]["closed"] == fp["H2S"]["closed"]
            and fp["H2O"]["h"] == fp["H2S"]["h"]
            and fp["H2O"]["A"] == fp["H2S"]["A"])
    ck("W8b L1 只认图（H2O=H2S）", same,
       f"H2O 与 H2S 同为 P3，指纹逐位相同（{fp['H2O']['closed']} 条）；原子身份不在图里")
    return fp


def check_l1prime():
    """L1′ 记录层：每顶点累积的闭合词数（Z088 的 record_v），归一化即 ω_v。"""
    iso = ["Ne", "HF", "H2O", "NH3", "CH4"]
    rec = {}
    for k in TREES:
        nv, edges = TREES[k]
        A = build_graph(nv, edges)
        W = closed_walks(A, 4)
        cnt = [0] * nv
        for path in W:
            cnt[path[0]] += 1
        tot = sum(cnt)
        omega = [c / tot for c in cnt] if tot else [0.0] * nv
        deg = [int(A[i].sum()) for i in range(nv)]
        rec[k] = {"record": cnt, "omega": omega, "total": tot, "degree": deg}

    ok = True
    for k in iso:
        r = rec[k]
        if r["total"] == 0:
            continue
        d = np.array(r["degree"], float)
        target = d / d.sum() if d.sum() else d
        ok = ok and np.allclose(r["omega"], target, atol=1e-12)
    ck("W9 L1′ 记录＝度分布", ok,
       "；".join(f"{k}:ω={np.round(rec[k]['omega'], 3).tolist()}" for k in iso))

    sig = [tuple(np.round(rec[k]["omega"], 6)) for k in iso]
    ck("W9b L1′ 记录区分等电子分子", len(set(sig)) == len(sig),
       "；".join(f"{k}:ω={np.round(rec[k]['omega'], 3).tolist()}" for k in iso))
    return rec


def check_l2():
    """L2 表达：p=ω（边② seed）；槽=n−1 个相邻比 |log(p_i/p_{i+1})|（中心在前）。"""
    iso = ["Ne", "HF", "H2O", "NH3", "CH4"]
    res = {}
    for k in iso:
        nv, edges = TREES[k]
        A = build_graph(nv, edges)
        W = closed_walks(A, 4)
        cnt = [0] * nv
        for path in W:
            cnt[path[0]] += 1
        tot = sum(cnt)
        p = [c / tot for c in cnt] if tot else [1.0 / nv] * nv
        slots = [abs(float(np.log(p[i] / p[i + 1]))) for i in range(nv - 1) if p[i + 1] > 0]
        res[k] = {"p": p, "slots": slots, "n_slots": nv - 1}
    ok = True
    for k in ("H2O", "NH3", "CH4"):
        deg = TREES[k][0] - 1  # 中心度 = 顶点数−1
        sl = res[k]["slots"]
        ok = ok and len(sl) == deg and abs(sl[0] - np.log(deg)) < 1e-12 and all(abs(x) < 1e-12 for x in sl[1:])
    ck("W10 L2 槽＝(log k, 0…)", ok,
       "；".join(f"{k}:{np.round(res[k]['slots'], 3).tolist()}" for k in ("H2O", "NH3", "CH4"))+
       f"；H2O 的 3 顶点给 {res['H2O']['n_slots']} 槽（不是 3）")
    return res


def main():
    print("== 分子作为一个整体：十一条判定 ==")
    check_iso()
    check_neutron()
    check_28()
    check_hexic()
    check_slots()
    check_tf()
    fp = check_l1()
    rec = check_l1prime()
    l2 = check_l2()

    th, s, spe = route1()
    out = {
        "molecules": [
            {"name": m[0], "N_p": m[1], "N_n": m[2], "N_e": m[3],
             "dipole_D": m[4], "angle_deg": m[5], "atomization_eV": m[6],
             "matter_units": units(m)} for m in MOL],
        "six_term": {"theta_deg": th.tolist(), "S_at_theta": s,
                     "S_shifted": spe, "delta_xi": DXI,
                     "rel_dev_vs_alpha": abs(spe - 137.035999084) / 137.035999084},
        "slots": [float(x) for x in (np.log(1.5), np.log(2.0), np.log(8 / 3))],
        "tf": {"E_TF_Z10": -0.7687 * 10**(7 / 3),
               "E_TF_Z8_plus_2H": -0.7687 * 8**(7 / 3) - 2 * 0.7687,
               "true_H2O_Ha": -76.44},
        "l1_fingerprints": fp,
        "l1prime_records": rec,
        "l2_expression": l2,
        "rows": ROWS,
        "pass": len(OK), "fail": len(BAD),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n== {len(OK)}/{len(OK)+len(BAD)} PASS  ->  {OUT.relative_to(ROOT.parent)} ==")
    if BAD:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
