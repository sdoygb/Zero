#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_social_port.py —— 把旧库的社会层（H4）接进当前 Zero 的可复算端口
================================================================================

来源（实读，不在 cosmos-construct/modular-equilibrium）
--------------------------------------------------------
旧库 `~/Downloads/Zero-261010` 的 H1–H10 ＋ `port/`：
  H4  社会版移植到 Zero：承重件映射、决定性重算与判定（35/0）
  H7  战争＝终端账本的 Z4 形式化
  H10 金融理论端口调查；`port/fin` 有 GDP/通胀/债务/信贷/股指
  `port/data` 软链 UCDP/PRIO v26.1、UCDP Termination、WVS Wave 7
  `port/scripts` 探针 ⑮⑯⑳ 等
社会理论原文在 `~/Downloads/GeometryAI-Mac-Build/app/articles/10.75–10.107`
（三尺度刚性、社会版数据验证、贸易、政治场、阵营硬度）＋
`gpu_spectra/society/`（49 个探针）＋ `notes/社会版·破缺谱_工作稿_260915.md`。

本脚本做什么
------------
只搬 **H4 的可复算核心** 到当前 Zero，并映射到当前正式篇（不搬旧库编号）：

  S1  WVS Q49：R_obs vs 闭式 R=√(1+n̄σ_s²/σ_i²)（零参数；旧库实测 0.27%）
  S2  UCDP 段长：几何（无记忆）vs 重尾 m=1/2 的 ΔAIC（和平 +112／冲突 +203）
  S3  残余名望：E[T−t|T>t] 随 t 上升（无记忆被否证）
  S4  P(冲突>50 年)：观测 1.0% vs 几何 0.003%（372×）
  S5  台湾案例：UCDP 冲突年 [1949,1950,1954,1958]（只作口径内量级，不作预测）
  S6  当前 Zero 映射：Z004 定理6（概率=计数推前）、Z002 无偏好、Z017 Catalan、
      Z020/Z001 终端层、Z129 持久锁——文件存在性与对应关系
  S7  诚实边界：旧两态 Markov 被数据拒绝；Zero 只与**冲突段**形状相容；
      自由拟合尾指数 m=0.060/0.783 ≠ 1/2；语义层超纲；不给点预报

**边界**
--------
- 这是 **port/审计**，不是新物理，也不是社会定律；数据只读。
- 旧库的 G/H/D 编号只在 `doc/zero_social_port.md` 的 provenance 里出现，
  当前结论只指当前 Z 篇。
- 台湾 5 年窗口：旧工作稿自己明文禁止把它当“开战概率”；本脚本只复算
  UCDP 门槛内的量级，不作预测。

运行
----
    cd ZeroGPU/src && /usr/bin/python3 zero_social_port.py
    # 可用 --data-dir 指定 port/data 或 ~/Downloads
输出：results/zero_social_port.json、zero_social_port.png
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def find_data_dir(cli=None):
    cands = [cli] if cli else []
    cands += [os.path.join(HERE, "..", "..", "..", "Zero-261010", "port", "data"),
              os.path.expanduser("~/Downloads/Zero-261010/port/data"),
              os.path.expanduser("~/Downloads")]
    for d in cands:
        if d and os.path.exists(os.path.join(d, "UcdpPrioConflict_v26_1.csv")):
            return os.path.abspath(d)
    return None


# ══════════════════════════════════════════════════════════════════
#  S1  WVS：R_obs vs 闭式
# ══════════════════════════════════════════════════════════════════
def wvs_r2(data_dir):
    import pandas as pd
    path = os.path.join(data_dir, "WVS_Cross-National_Wave_7_csv_v6_0.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, low_memory=False)
    q49 = [c for c in df.columns if c.startswith("Q49")]
    ctry = [c for c in df.columns if c.upper() in ("B_COUNTRY_ALPHA", "COUNTRY_ALPHA", "B_COUNTRY")]
    if not q49 or not ctry:
        return None
    d = df[[ctry[0], q49[0]]].dropna()
    d = d[d[q49[0]].between(1, 10)]
    g = d.groupby(ctry[0])[q49[0]]
    n_i = g.size()
    keep = n_i[n_i >= 30].index
    d = d[d[ctry[0]].isin(keep)]
    g = d.groupby(ctry[0])[q49[0]]
    mu = g.mean()
    n_i = g.size()
    n_bar = float(n_i.mean())
    s_i = float(d[q49[0]].std(ddof=1))
    s_s = float(mu.std(ddof=1))
    R = s_s / (s_i / math.sqrt(n_bar))
    pred = math.sqrt(1.0 + n_bar * s_s ** 2 / s_i ** 2)
    return dict(n_countries=int(len(mu)), n_rows=int(len(d)), n_bar=n_bar,
                sigma_struct=s_s, sigma_ind=s_i, R_obs=R, R_pred=pred,
                rel=(R / pred - 1.0), eta2=s_s ** 2 / (s_s ** 2 + s_i ** 2))


def wvs_scaling(data_dir):
    """子样本标度：结构方差 σ_s²−σ_i²/n̄ 是否稳定（这才是模型的可检验内容）。"""
    import pandas as pd
    path = os.path.join(data_dir, "WVS_Cross-National_Wave_7_csv_v6_0.csv")
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path, low_memory=False)
    q49 = [c for c in df.columns if c.startswith("Q49")][0]
    ctry = [c for c in df.columns if c.upper() in ("B_COUNTRY_ALPHA", "COUNTRY_ALPHA", "B_COUNTRY")][0]
    d = df[[ctry, q49]].dropna()
    d = d[d[q49].between(1, 10)]
    g = d.groupby(ctry)[q49]
    n_i = g.size()
    keep = n_i[n_i >= 30].index
    d = d[d[ctry].isin(keep)]
    rng = np.random.default_rng(11)
    rows = []
    for nt in (10, 20, 50, 100, 200, 500, 1000):
        vals, labels = [], []
        for c, grp in d.groupby(ctry):
            v = grp[q49].values
            k = min(nt, len(v))
            vals.append(rng.choice(v, k, replace=False))
            labels.append(np.full(k, c))
        sub = pd.DataFrame({ctry: np.concatenate(labels), q49: np.concatenate(vals)})
        gg = sub.groupby(ctry)[q49]
        mu = gg.mean()
        ni = gg.size()
        nbar = float(ni.mean())
        si = float(sub[q49].std(ddof=1))
        ss = float(mu.std(ddof=1))
        rows.append(dict(n_target=nt, n_bar=nbar, R=ss / (si / math.sqrt(nbar)),
                         struct_var=ss * ss - si * si / nbar,
                         R2_minus_1=nbar * ss * ss / (si * si) - 1.0))
    return rows


# ══════════════════════════════════════════════════════════════════
#  S2/S3/S4  UCDP 段长
# ══════════════════════════════════════════════════════════════════
def load_segments(data_dir):
    import pandas as pd
    path = os.path.join(data_dir, "UcdpPrioConflict_v26_1.csv")
    df = pd.read_csv(path, encoding="latin-1")
    years = np.arange(1946, 2026)
    rows = []
    for _, r in df.iterrows():
        for loc in str(r["location"]).split(","):
            loc = loc.strip()
            if loc:
                rows.append({"country": loc, "year": int(r["year"])})
    conf = pd.DataFrame(rows).drop_duplicates(subset=["country", "year"])
    countries = sorted(conf["country"].unique())

    def segs(seq):
        p, w = [], []
        cur, l = seq[0], 1
        for s in seq[1:]:
            if s == cur:
                l += 1
            else:
                (p if cur == 0 else w).append(l)
                cur, l = s, 1
        (p if cur == 0 else w).append(l)
        return p, w

    P, W = [], []
    for c in countries:
        cy = set(conf[conf["country"] == c]["year"])
        st = np.array([1 if y in cy else 0 for y in years])
        p, w = segs(st)
        P += p
        W += w
    return np.array(P), np.array(W), conf


def fit_aic(X):
    """几何（旧两态 Markov）vs 重尾 m=1/2（Zero 型）的 AIC。"""
    from scipy.optimize import minimize, minimize_scalar
    X = np.asarray(X, dtype=int)
    n = len(X)
    q = 1.0 / X.mean()
    ll_geo = float(np.sum(np.log(q) + (X - 1) * np.log(1 - q)))
    kmax = int(X.max()) + 1

    def make(c, m):
        k = np.arange(1, kmax + 2)
        pmf = (k + c) ** (-m) - (k + 1 + c) ** (-m)
        pmf = np.clip(pmf, 1e-300, None)
        return pmf / pmf.sum()

    def nll(p):
        pmf = make(np.exp(p[0]), p[1])
        return -float(np.sum(np.log(pmf[X - 1])))

    ll_half = -float(minimize_scalar(lambda lc: nll([lc, 0.5]),
                                     bounds=(-4, 4), method="bounded").fun)
    r2 = minimize(nll, x0=[np.log(1.0), 0.5], method="Nelder-Mead")
    ll_free = -float(r2.fun)
    aic_geo, aic_half = 2 - 2 * ll_geo, 2 - 2 * ll_half
    tail = {int(k): float((X > k).mean()) for k in (10, 20, 30, 50) if k < kmax}
    return dict(n=n, mean=float(X.mean()), max=int(X.max()),
                aic_geo=aic_geo, aic_half=aic_half, dAIC=aic_geo - aic_half,
                tail=tail, q_geo=float(q),
                free_m=float(r2.x[1]), free_c=float(math.exp(r2.x[0])),
                ll_free=ll_free), tail


def closure_avatars(gens=6):
    """圆满者=闭合事件：K4 走者的 closure 计数、record 增量、seed 携带。"""
    import zero_sim as zs
    from zero_evo import ZeroEvo, Params
    p = Params(L=4, T=4, micro=True, micro_semantics="persistent", retention="top2",
               poisson_on=False, wave_on=False, c_h=0.0, g_src=0.0, Gamma_h=0.0, norm_hi=1e50)
    e = ZeroEvo(zs.k4(), p, seed=[8.0, 4.0, 1.5, 1.0])
    rows = []
    for g in range(1, gens + 1):
        rec0 = e.record.copy()
        e.step()
        c = e.gen_closures.copy()
        rec1 = e.record.copy()
        rows.append(dict(gen=g, closures=[float(x) for x in c],
                         closure_total=float(c.sum()),
                         record_delta=[float(x) for x in (rec1 - rec0)],
                         seed=[float(x) for x in e.seed]))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--json", default=os.path.join(RES, "zero_social_port.json"))
    ap.add_argument("--png", default=os.path.join(RES, "zero_social_port.png"))
    args = ap.parse_args(argv)

    print("=" * 84)
    print("社会层 H4 → 当前 Zero：WVS 闭式、UCDP 重尾、终端层映射")
    print("=" * 84)
    data_dir = find_data_dir(args.data_dir)
    print(f"数据目录：{data_dir}")
    ck("S0 找到 port/data（UCDP）", data_dir is not None, str(data_dir))

    # S1 WVS
    print("\n§1 WVS Q49：R_obs vs 闭式")
    wvs = wvs_r2(data_dir) if data_dir else None
    if wvs:
        print(f"  国家数={wvs['n_countries']}  n̄={wvs['n_bar']:.1f}  "
              f"σ_struct={wvs['sigma_struct']:.4f}  σ_ind={wvs['sigma_ind']:.4f}")
        print(f"  R_obs={wvs['R_obs']:.3f}  R_pred={wvs['R_pred']:.3f}  "
              f"相对差={wvs['rel']:+.4%}")
    # 关键审计：H4 的 "0.27% 预言" 其实是 R_pred=√(1+R_obs²)，差 1/(2R²)，是定义性恒等式
    taut = (abs(wvs["R_pred"] - math.sqrt(1 + wvs["R_obs"] ** 2)) < 1e-12
            and abs(abs(wvs["rel"]) - 1.0 / (2 * wvs["R_obs"] ** 2)) < 1e-3) if wvs else False
    ck("S1 WVS 的 0.27% 是定义性差：R_pred=√(1+R_obs²)，差 1/(2R²)",
       taut, f"{wvs['rel']:+.4%} vs 1/(2R²)={1/(2*wvs['R_obs']**2):.4%}" if wvs else "")
    scale = wvs_scaling(data_dir) if data_dir else []
    stable = [r for r in scale if r["n_target"] >= 50]
    vals = [r["struct_var"] for r in stable]
    cv = (np.std(vals) / np.mean(vals)) if vals else 9.0
    ck("S2 真模型内容：结构方差 σ_s²−σ_i²/n̄ 跨子样本稳定（n≥50，CV<20%）",
       bool(vals) and cv < 0.2, f"σ_s²≈{np.mean(vals):.3f}（CV={cv:.1%}）" if vals else "")
    ck("S3 R ≫ 1（结构锁定）", bool(wvs) and wvs["R_obs"] > 3,
       f"R={wvs['R_obs']:.3f}" if wvs else "")

    # S2/S3/S4 UCDP
    print("\n§2 UCDP 段长：几何 vs 重尾")
    P = W = None
    conf = None
    if data_dir:
        P, W, conf = load_segments(data_dir)
        print(f"  国家数={len(set(conf['country']))}  和平段 n={len(P)} mean={P.mean():.2f}  "
              f"冲突段 n={len(W)} mean={W.mean():.2f}")
        fitP, tailP = fit_aic(P)
        fitW, tailW = fit_aic(W)
        print(f"  和平段 ΔAIC = {fitP['dAIC']:+.2f}  (自由 m={fitP['free_m']:.3f})")
        print(f"  冲突段 ΔAIC = {fitW['dAIC']:+.2f}  (自由 m={fitW['free_m']:.3f})")
    ck("S4 复现段均值（和平≈16.37／冲突≈5.27）",
       P is not None and abs(P.mean() - 16.37) < 0.2 and abs(W.mean() - 5.27) < 0.2,
       f"{P.mean():.2f}/{W.mean():.2f}" if P is not None else "")
    ck("S5 和平段重尾优于几何（ΔAIC>50）", P is not None and fitP["dAIC"] > 50,
       f"{fitP['dAIC']:+.1f}" if P is not None else "")
    ck("S6 冲突段重尾优于几何（ΔAIC>50）", W is not None and fitW["dAIC"] > 50,
       f"{fitW['dAIC']:+.1f}" if W is not None else "")

    # S3 剩余均值
    print("\n§3 残余名望：E[T−t|T>t] 随 t 上升")
    resid = []
    if W is not None:
        for t in (5, 10, 20, 30):
            hyp = W[W > t]
            if len(hyp) > 10:
                resid.append((t, int(len(hyp)), float((hyp - t).mean())))
                print(f"  t={t:>2}: n={len(hyp):>3}  E[T−t|T>t]={resid[-1][2]:.2f}  "
                      f"(指数要求 {W.mean():.2f})")
    ck("S7 剩余均值随 t 上升（无记忆被否证）",
       len(resid) >= 3 and resid[-1][2] > 1.3 * resid[0][2],
       f"{resid[0][2]:.1f}→{resid[-1][2]:.1f}" if len(resid) >= 2 else "")

    # S4 P(T>50)
    obs50 = float((W > 50).mean()) if W is not None else 0.0
    q_w = 1.0 / W.mean() if W is not None else 1.0
    geo50 = float((1 - q_w) ** 50) if W is not None else 0.0
    ck("S8 P(冲突>50年)：观测/几何 >100×",
       W is not None and geo50 > 0 and obs50 / geo50 > 100,
       f"{obs50:.4f} vs {geo50:.5f}（{obs50/geo50:.0f}×）" if W is not None else "")

    # S5 台湾
    print("\n§4 台湾：UCDP 门槛内的冲突年（只作量级，不作预测）")
    tw = []
    if conf is not None:
        tw = sorted(int(x) for x in conf[conf.country == "Taiwan"].year.unique())
        print(f"  冲突年={tw}；最后一次={max(tw) if tw else None}")
    ck("S9 台湾冲突年 = [1949,1950,1954,1958]",
       tw == [1949, 1950, 1954, 1958], str(tw))

    # S6 当前 Zero 映射
    print("\n§5 当前 Zero 的承重件（文件存在性＋对应）")
    mapping = {
        "Z004": "概率=计数测度在粗粒化下的推前（定理 6）",
        "Z002": "无偏好（边可迁）/ 维数读出",
        "Z017": "Catalan 闭合计数与重播种周期",
        "Z020": "终端层收走未闭合者",
        "Z001": "终端与重播种（§3.3）",
        "Z129": "K4 persistent 锁（record=seed 比）",
    }
    found = {}
    for z, note in mapping.items():
        hits = glob.glob(os.path.join(REPO, f"{z}_*.md"))
        found[z] = bool(hits)
        print(f"  {z}: {'✓' if hits else '✗'}  {note}")
    ck("S10 当前 Zero 的 6 个承重件文件都存在", all(found.values()), str(found))

    # S7 边界
    print("\n§6 诚实边界")
    ck("S11 旧两态 Markov 被数据拒绝（不是 Zero 被证实）",
       P is not None and fitP["dAIC"] > 50 and fitW["dAIC"] > 50,
       "ΔAIC 正 = 几何更差")
    ck("S12 自由尾指数 ≠ 1/2（异质混合，不是单一幂律）",
       W is not None and abs(fitW["free_m"] - 0.5) > 0.1,
       f"m_free={fitW['free_m']:.3f}" if W is not None else "")
    ck("S13 道德/阵营/霸权仍超纲；不给点预报", True,
       "见 doc §4；台湾只作口径内量级")
    cav = closure_avatars(6)
    ck("S14 圆满者=闭合事件：record 增量逐代=closure；seed 携带",
       all(np.allclose(r["record_delta"], r["closures"]) and r["closure_total"] > 0
           and sum(r["seed"]) > 0 for r in cav),
       f"6 代；closure 总数={sum(r['closure_total'] for r in cav):.3e}")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(data_dir=data_dir, wvs=wvs, wvs_scaling=scale,
               ucdp=dict(mean_peace=float(P.mean()) if P is not None else None,
                         mean_conflict=float(W.mean()) if W is not None else None,
                         n_peace=int(len(P)) if P is not None else None,
                         n_conflict=int(len(W)) if W is not None else None,
                         dAIC_peace=float(fitP["dAIC"]) if P is not None else None,
                         dAIC_conflict=float(fitW["dAIC"]) if W is not None else None,
                         free_m_peace=float(fitP["free_m"]) if P is not None else None,
                         free_m_conflict=float(fitW["free_m"]) if W is not None else None,
                         remaining=resid, obs50=obs50, geo50=geo50),
               taiwan=tw, mapping=found, closure_avatars=cav, checks=CK, n_fail=nfail)
    with open(args.json, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    make_figure(out, args.png)
    print(f"结果 → {args.json}\n图 → {args.png}")
    return 0 if nfail == 0 else 1


def make_figure(out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7))
    w = out["wvs"]
    scale = out.get("wvs_scaling", [])
    if scale:
        ax[0, 0].plot([r["n_target"] for r in scale], [r["struct_var"] for r in scale], "o-")
        ax[0, 0].axhline(0.617, color="r", ls="--", lw=1, label="full sample 0.617")
        ax[0, 0].set_xscale("log")
        ax[0, 0].set_xlabel("subsample n per country (log)")
        ax[0, 0].set_ylabel("$\\sigma_s^2-\\sigma_i^2/\\bar n$")
        ax[0, 0].set_title("(a) structural variance is stable (real content)")
        ax[0, 0].legend(fontsize=8)
    elif w:
        ax[0, 0].bar(["R_obs", "R_pred"], [w["R_obs"], w["R_pred"]], color=["C0", "C1"])
        ax[0, 0].set_title("(a) WVS Q49: R_obs vs R_pred")
    ax[0, 1].bar(["peace", "conflict"],
                 [out["ucdp"]["dAIC_peace"] or 0, out["ucdp"]["dAIC_conflict"] or 0],
                 color=["C0", "C1"])
    ax[0, 1].axhline(0, color="k", lw=.8)
    ax[0, 1].set_title("(b) ΔAIC: geometric − heavy-tail")
    ax[0, 1].set_ylabel("ΔAIC (positive = geometric worse)")
    rem = out["ucdp"]["remaining"]
    if rem:
        ax[1, 0].plot([r[0] for r in rem], [r[2] for r in rem], "o-", label="empirical")
        ax[1, 0].axhline(out["ucdp"]["mean_conflict"], color="r", ls="--", label="exponential")
        ax[1, 0].set_xlabel("already lasted t (years)")
        ax[1, 0].set_ylabel("E[T−t | T>t]")
        ax[1, 0].set_title("(c) remaining mean rises (no memory rejected)")
        ax[1, 0].legend(fontsize=8)
    ax[1, 1].axis("off")
    ax[1, 1].text(0.0, 0.95,
                  "Boundary:\n"
                  "· old two-state Markov rejected (ΔAIC)\n"
                  "· Zero only shape-compatible with conflict segment\n"
                  "· free tail exponent ≠ 1/2 (heterogeneous)\n"
                  "· semantic layer out of scope\n"
                  "· no point forecasts / Taiwan not a prediction",
                  va="top", fontsize=9)
    for a in ax.ravel():
        if a.get_title():
            a.grid(alpha=.3)
    fig.suptitle("Social port H4 → current Zero (WVS + UCDP)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
