"""
charged_sector.py --- 任务 2：**带电扇区（终端层）** 普查 ＋ 带电局域态（电子候选）

**为什么做这个**
    记录层（L1′）**全是中性的**：闭合词 $\sum_i w_i=0\Rightarrow Q=0$。
    带电的是**非闭合**世界线：计数 $2^n-\binom{n}{n/2}$（$n{=}32$ 时 $3.694\times10^9$，
    **比记录层（$6.01\times10^8$）大 6.1 倍**）。语料把它当死路（$\mathrm{Succ}=\varnothing$），
    **此前一个数都没算过**。

**荷是什么（Z1 §3.1 的散度字典）**
    词 $w$ 的散度 $d_w(v)=\#\{\text{进 }v\}-\#\{\text{出 }v\}$，$\sum_v d_w(v)\equiv0$ 是恒等式。
    闭合 $\iff d_w=0$。从 $u$ 走到 $v$ 的世界线 $\Rightarrow d_w=e_v-e_u$（**偶极**）。
    所以"荷"在 $\Gamma$ 上就是**端点位移**；不带 $\Gamma$ 时就是词的平衡 $b=\sum_i w_i$。

**三部分**
  A **计数普查**（精确整数）：带电计数、与记录层的比、荷分布 $P(b)$ 的结构、带 $\Gamma$ 后带电占比
  B **荷的几何**：荷矢量 $=$ 端点，$d_w=e_v-e_0$ 的逐例核验；带电世界线的空间分布
  C **带电局域态（电子候选）**：把 $+Q,-Q$ 放在相距 $R$ 的两点（**偶极**），
    问最小系统能不能把它**束缚**成局域态 —— 即经典问题「**临界偶极矩**」：
    连续极限下 3D 的临界偶极矩 $p_c\approx0.639$（有限值），2D **任意偶极都束缚**（$p_c=0$）。
    本脚本在格子上测 $p_c(R)=V_c\cdot R$ 是否收敛到常数 —— 这是**电子的经典判据**在 Zero 里的版本。

**等级**：【数值证据】；A 部分是精确整数（【闭式】）。
用法：/usr/bin/python3 charged_sector.py
输出：results/charged_sector.json
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np
from scipy.sparse import csr_matrix, diags, identity, kron
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

RECORDS = []


def rec(sec, obj, readout, x, y, meta=None, note=None, kind="curve", maxpts=64):
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) > maxpts:
        idx = np.unique(np.r_[np.linspace(0, len(x) - 1, maxpts - 1).astype(int), len(x) - 1])
        x, y = x[idx], y[idx]
    RECORDS.append({"id": f"{sec}:{obj}:{readout}", "sector": sec, "object": obj,
                    "readout": readout, "kind": kind,
                    "x": [float(v) for v in x], "y": [float(v) for v in y],
                    "meta": meta or {}, **({"note": note} if note else {})})


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def site_index(coord, L):
    i = 0
    for c in coord:
        i = i * L + int(c)
    return i


# =====================================================================
# A. 计数普查（精确整数）
# =====================================================================
def part_A():
    ns = np.arange(2, 65, 2)
    ch, nu, ratio = [], [], []
    for n in ns:
        n = int(n)
        nu.append(math.comb(n, n // 2))
        ch.append(2 ** n - math.comb(n, n // 2))
        ratio.append((2 ** n - math.comb(n, n // 2)) / math.comb(n, n // 2))
    rec("CHARGED", "全词空间 {±1}^n", "带电计数 2^n - C(n,n/2)", ns, ch,
        {"n32_charged": str(2 ** 32 - math.comb(32, 16)),
         "n32_neutral": str(math.comb(32, 16)),
         "n32_ratio": round((2 ** 32 - math.comb(32, 16)) / math.comb(32, 16), 4)},
        note="n=32：带电 36.9 亿 vs 中性 6.01 亿 = 6.14 倍")
    rec("CHARGED", "全词空间 {±1}^n", "带电/中性 比 vs n", ns, ratio,
        note="闭式 ~ sqrt(pi n/2) - 1（幂律增长）")
    print(f"  A1 n=32: 带电 {2**32 - math.comb(32,16):,} / 中性 {math.comb(32,16):,} "
          f"= {(2**32 - math.comb(32,16))/math.comb(32,16):.4f} 倍", flush=True)

    # 荷分布 P(b)：b = 2k - n
    for n in (16, 32, 64):
        k = np.arange(n + 1)
        b = 2 * k - n
        P = np.array([math.comb(n, int(kk)) for kk in k], float) / 2.0 ** n
        rec("CHARGED", f"荷分布 P(b), n={n}", "P(b) vs b", b, P,
            {"n": n, "variance": n, "note": "b=0 是中性（记录层），b≠0 是带电"})
    # 带电占比随 n
    frac = np.array([1.0 - math.comb(int(n), int(n) // 2) / 2.0 ** int(n) for n in ns])
    rec("CHARGED", "全词空间", "中性占比 C(n,n/2)/2^n vs n", ns, 1 - frac,
        note="~ 1/sqrt(pi n/2) —— 记录层是大 n 的**消失少数**")

    # 带 Γ 后：带电占比 = 1 - n_k(0)/(2D)^k
    for D in (1, 2, 3, 4):
        K = {1: 40, 2: 24, 3: 16, 4: 10}[D]
        ks = np.arange(0, K + 1, 2)          # 只取偶数 k：奇数 k 闭合数恒为 0（二部图），
        closed = np.array([_exact_return(D, int(k)) for k in ks], float)   # 混进去会污染幂律拟合
        tot = np.array([float(2 * D) ** int(k) for k in ks])
        rec("CHARGED", f"Γ=Z^{D} 上的世界线", "中性（闭合）占比 vs k", ks, closed / tot,
            {"D": D}, note="Γ 上的中性占比 ~ k^{-D/2}，比词层掉得更快")
        rec("CHARGED", f"Γ=Z^{D} 上的世界线", "带电占比 vs k", ks, 1 - closed / tot, {"D": D})
    print("  A2 荷分布与带电占比：完成", flush=True)


def _exact_return(D, k):
    """精确整数：n_k(0)（k 奇数为 0）。"""
    if k % 2:
        return 0
    m, tot = k // 2, 0
    if D == 1:
        return math.factorial(2 * m) // (math.factorial(m) ** 2)
    for a1 in range(m + 1):
        f1 = math.factorial(a1) ** 2
        if D == 2:
            tot += math.factorial(2 * m) // (f1 * math.factorial(m - a1) ** 2)
        else:
            for a2 in range(m - a1 + 1):
                f2 = f1 * math.factorial(a2) ** 2
                if D == 3:
                    tot += math.factorial(2 * m) // (f2 * math.factorial(m - a1 - a2) ** 2)
                else:
                    for a3 in range(m - a1 - a2 + 1):
                        a4 = m - a1 - a2 - a3
                        tot += math.factorial(2 * m) // (f2 * math.factorial(a3) ** 2
                                                         * math.factorial(a4) ** 2)
    return tot


# =====================================================================
# B. 荷的几何：d_w = e_v - e_u 的逐例核验
# =====================================================================
def part_B():
    rng = np.random.default_rng(0)
    rows = []
    for D in (1, 2, 3):
        L = 31
        n_ok = 0
        for _ in range(200):
            k = int(rng.integers(4, 40))
            steps = rng.integers(0, 2 * D, size=k)
            pos = np.zeros(D, np.int64)
            div = {}
            start = tuple(pos)
            for s in steps:
                a, sg = int(s) // 2, (1 if s % 2 == 0 else -1)
                nxt = pos.copy()
                nxt[a] = (nxt[a] + sg) % L
                # 有向边计数：出 pos、进 nxt
                div[start] = div.get(start, 0) - 1
                key = tuple(nxt)
                div[key] = div.get(key, 0) + 1
                pos = nxt
                start = key
            # 散度应当 = e_v - e_u
            end = tuple(pos)
            expect = {}
            if end != (0,) * D:
                expect[end] = expect.get(end, 0) + 1
                expect[(0,) * D] = expect.get((0,) * D, 0) - 1
            got = {k2: v for k2, v in div.items() if v != 0}
            if got == expect:
                n_ok += 1
        rows.append({"D": D, "checked": 200, "divergence_equals_endpoint_minus_start": n_ok})
        print(f"  B D={D}: d_w = e_v - e_u 核验 {n_ok}/200", flush=True)
    rec("CHARGED", "世界线散度", "d_w = e_v - e_u 的核验通过率", [0], [0],
        {"runs": rows}, note="Z1 §3.1 的散度字典在随机样本上逐例成立", kind="declared")
    return rows


def _dist_from(i0, D, L):
    """环面上每个格点到 site0 的图距离（BFS）。"""
    from collections import deque
    dist = -np.ones(L ** D, np.int64)
    dist[i0] = 0
    q = deque([i0])
    while q:
        u = q.popleft()
        coord = []
        x = u
        for _ in range(D):
            coord.append(x % L)
            x //= L
        for d in range(D):
            for sg in (1, -1):
                c = list(coord)
                c[d] = (c[d] + sg) % L
                v = site_index(c, L)
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    q.append(int(v))
    return dist


# =====================================================================
# C. 带电局域态：偶极束缚（临界偶极矩）
# =====================================================================
def part_C():
    LD = {1: 401, 2: 81, 3: 31, 4: 11}
    rs = {1: [2, 4, 8, 16, 32], 2: [2, 4, 8, 16], 3: [2, 4, 8], 4: [1, 2, 4]}
    out = []
    for D in (1, 2, 3, 4):
        L = LD[D]
        A = torus_adj(D, L)
        n = L ** D
        deg = np.asarray(A.sum(1)).ravel()
        thr = 20.0 / n
        band = 2.0 * D
        rows = []
        for r in rs[D]:
            if r >= L // 2:
                continue
            coordR = [0] * D
            coordR[0] = r
            i0, iR = site_index([0] * D, L), site_index(coordR, L)
            diag = np.zeros(n)
            diag[i0] = 1.0
            diag[iR] = -1.0
            Dm = diags(diag)
            Vgrid = np.logspace(-2, 2.2, 26)
            ipr_up, ipr_dn, e_up, e_dn = [], [], [], []
            for V in Vgrid:
                H = (A + V * Dm).tocsr()
                try:
                    wu, vu = eigsh(H, k=1, which="LA")
                    wd, vd = eigsh(H, k=1, which="SA")
                except Exception:
                    ipr_up.append(np.nan); ipr_dn.append(np.nan)
                    e_up.append(np.nan); e_dn.append(np.nan)
                    continue
                pu = np.asarray(vu[:, 0]).ravel() ** 2
                pd = np.asarray(vd[:, 0]).ravel() ** 2
                ipr_up.append(float((pu ** 2).sum()))
                ipr_dn.append(float((pd ** 2).sum()))
                e_up.append(float(wu[0]))
                e_dn.append(float(wd[0]))
            ipr_up, ipr_dn = np.array(ipr_up), np.array(ipr_dn)
            e_up, e_dn = np.array(e_up), np.array(e_dn)
            # **首次**越过局域化阈值的 V（不是 IPR 最大处的 V —— 上一版的 bug）
            iu = np.nonzero(np.nan_to_num(ipr_up) > thr)[0]
            idn = np.nonzero(np.nan_to_num(ipr_dn) > thr)[0]
            vu_c = float(Vgrid[iu[0]]) if len(iu) else None
            vd_c = float(Vgrid[idn[0]]) if len(idn) else None
            Vc = min([v for v in (vu_c, vd_c) if v is not None], default=None)
            rec("CHARGE-STATE", f"Z^{D} 偶极 R={r}", "极值态 IPR vs V（上带边）", Vgrid, ipr_up,
                {"D": D, "R": r, "L": L, "IPR_extended=1/N": 1.0 / n, "threshold": thr,
                 "band_edge": band})
            rec("CHARGE-STATE", f"Z^{D} 偶极 R={r}", "极值态 IPR vs V（下带边）", Vgrid, ipr_dn,
                {"D": D, "R": r, "L": L, "IPR_extended=1/N": 1.0 / n, "threshold": thr,
                 "band_edge": band})
            # 阈值处的束缚态：径向轮廓 → 局域长度 ξ（与 z0_bound_states 的中性缺陷对照）
            xi = None; loc_site = None; prof = None
            if len(iu):
                H = (A + float(Vgrid[iu[0]]) * Dm).tocsr()
                wu2, vu2 = eigsh(H, k=1, which="LA")
                psi = np.asarray(vu2[:, 0]).ravel()
                amps = psi ** 2
                loc_site = "site0(+Q)" if amps[i0] > amps[iR] else f"siteR(-Q),R={r}"
                dist = _dist_from(i0, D, L)
                m = (dist >= 1) & (dist <= max(2, L // 2 - 1)) & (amps > 1e-14)
                if m.sum() >= 4:
                    sl = np.polyfit(dist[m], np.log(amps[m]), 1)[0]
                    if sl < 0:
                        xi = -1.0 / sl
                    prof = {"dist": dist[m].astype(int).tolist()[:40],
                            "amp": amps[m].tolist()[:40]}
            rows.append({"R": r, "Vc": None if Vc is None else round(float(Vc), 5),
                         "xi_localization": None if xi is None else round(float(xi), 4),
                         "localized_at": loc_site,
                         "profile": prof,
                         "p_c": None if Vc is None else round(float(Vc) * r, 5),
                         "Vc_top": vu_c, "Vc_bot": vd_c,
                         "E_top_at_Vc": None if not len(iu) else round(float(e_up[iu[0]]), 5),
                         "E_bot_at_Vc": None if not len(idn) else round(float(e_dn[idn[0]]), 5),
                         "E_top_minus_bandedge": None if not len(iu) else round(float(e_up[iu[0]]) - band, 5),
                         "IPR_localized_max": round(float(np.nanmax(ipr_up)), 5),
                         "band_edge": band})
            print(f"  C D={D} R={r:>2} L={L}: Vc={rows[-1]['Vc']}  p_c=Vc*R={rows[-1]['p_c']}  "
                  f"ξ={rows[-1]['xi_localization']}  局域在 {rows[-1]['localized_at']}", flush=True)
        good = [x for x in rows if x["p_c"]]
        if len(good) >= 2:
            RR = np.array([x["R"] for x in good], float)
            PP = np.array([x["p_c"] for x in good], float)
            rec("CHARGE-STATE", f"Z^{D} 偶极", "临界偶极矩 p_c=Vc*R vs R", RR, PP, {"D": D},
                note="收敛到常数 ⟺ 连续极限的经典临界偶极矩（3D ≈0.639）")
        out.append({"D": D, "L": L, "band_edge": band, "rows": rows})
    return out


# =====================================================================
# D. 与记录层的关联：带电尾巴的空间足迹 vs 中性足迹的重叠
# =====================================================================
def part_D(n_walk=4000, seed=3):
    """
    对每条世界线：最后一次回到原点把世界线切成
        中性段（[0,t*]，也就是"已闭合/记录"的部分）与**带电尾巴**（(t*,n]，端点≠原点）。
    测带电尾巴的**新足迹**占比（= 它走出了记录层没覆盖过的地方的比例），以及重叠系数。
    """
    rng = np.random.default_rng(seed)
    out = []
    for D in (1, 2, 3, 4):
        L = {1: 2001, 2: 201, 3: 61, 4: 27}[D]
        n = 400
        rows = []
        for trial in range(n_walk):
            pos = np.zeros(D, np.int64)
            step = rng.integers(0, 2 * D, size=n)
            axis, sign = step // 2, np.where(step % 2 == 0, 1, -1)
            sites = np.zeros(n + 1, np.int64)
            for t in range(n):
                a = int(axis[t])
                pos[a] = (pos[a] + int(sign[t])) % L
                key = 0
                for d in range(D):
                    key = key * L + int(pos[d])
                sites[t + 1] = key
            z = np.nonzero(sites[1:] == 0)[0]
            if len(z) == 0:
                continue                      # 从没闭合过：整条都带电
            tstar = int(z[-1]) + 1
            neu = set(sites[:tstar].tolist())
            chg = set(sites[tstar + 1:].tolist()) - {0}
            if not chg:
                continue
            new = len(chg - neu)
            rows.append((new / len(chg), len(chg & neu) / max(1, len(chg))))
        if not rows:
            continue
        arr = np.array(rows)
        out.append({"D": D, "L": L, "n_steps": n, "walkers_used": len(rows),
                    "charged_tail_new_site_frac": round(float(arr[:, 0].mean()), 4),
                    "charged_tail_overlap_with_record": round(float(arr[:, 1].mean()), 4)})
        print(f"  D D={D}: 带电尾巴的新足迹占比={out[-1]['charged_tail_new_site_frac']}  "
              f"与记录层足迹重叠={out[-1]['charged_tail_overlap_with_record']}  "
              f"({len(rows)} 条)", flush=True)
    rec("CHARGED", "带电尾巴 vs 记录层足迹", "新足迹占比 vs D",
        [r["D"] for r in out], [r["charged_tail_new_site_frac"] for r in out],
        {"detail": out}, note="带电扇区主要在记录层没覆盖过的格点上活动", kind="declared")
    return out


# =====================================================================
def main(argv):
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "definition": "荷 = Z1 散度 d_w；闭合 d_w=0（中性，记录层）；从 u 到 v 的世界线 d_w=e_v-e_u（偶极）",
           "reading": "不吸收（Z0① 字面）"}
    only = argv[0] if argv else None
    if only not in ("C",):
        print("[A] 计数普查")
        part_A()
        print("[B] 荷的几何（散度字典核验）")
        out["divergence_check"] = part_B()
    if only not in ("A", "B"):
        print("[C] 带电局域态：偶极束缚")
        out["dipole"] = part_C()
    if only not in ("A", "B", "C"):
        print("[D] 与记录层的关联")
        out["record_correlation"] = part_D()

    # 结构分类（复用任务 1 的分类器）
    import structures as ST
    for r in RECORDS:
        if r["kind"] == "declared" or len(r["x"]) < 5:
            r["structure"] = {"label": "declared", "params": {}}
        else:
            r["structure"] = ST.classify_curve(r["x"], r["y"])
    out["records"] = RECORDS
    out["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "charged_sector.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n登记 {len(RECORDS)} 条，用时 {out['seconds']}s → results/charged_sector.json")
    print("结构聚类：")
    from collections import Counter
    for k, v in Counter(r["structure"]["label"] for r in RECORDS).most_common():
        print(f"  {k:<26} {v}")


if __name__ == "__main__":
    main(sys.argv[1:])
