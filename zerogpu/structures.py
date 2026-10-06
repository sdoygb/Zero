"""
structures.py --- 结构分类器（任务 1 的"中性"那一半）

**用途**：给一条数值曲线 (x, y) 只贴**结构标签**，不贴物理名字。
物理名字在文档里另外分开写 —— 让结构自己报名字，而不是我先挑名字再找曲线。

标签集
    constant            常数（无标度）
    logarithmic         y ~ a + b ln x
    power-law           y ~ A x^b
    exponential         y ~ A e^{bx}
    stretched(beta)     y ~ A e^{-c x^beta}（beta=2 即 gaussian-tail）
    saturating          y -> y_inf，差按指数逼近
    delta               分布集中在单个 bin
    parity-alternating  奇偶之一恒为零
    periodic(P)         残差有单一主周期
    gapped              谱密度在一段内部区间恒为零
    identically-zero
    too-short / unclassified

**选择规则（可审计）**：所有候选模型都在**对数空间**算 R²（衰减曲线跨若干数量级，
线性空间的 R² 会被最大值主导）。取 R² 最高者，但若更简单的模型 R² 与之相差
< tol(=2e-3)，则取更简单的。模型复杂度次序：
    constant < logarithmic < power-law < exponential < stretched。
**所有候选连同 R² 一并输出**，便于复核与改判。

等级：工具，无物理断言。
"""
from __future__ import annotations

import numpy as np

TOL = 2e-3
_ORDER = {"constant": 0, "logarithmic": 1, "power-law": 2, "exponential": 3, "stretched": 4}


def _r2(y, yhat):
    y = np.asarray(y, float)
    yhat = np.asarray(yhat, float)
    ss = ((y - y.mean()) ** 2).sum()
    if ss <= 0:
        return 1.0 if np.allclose(y, yhat) else 0.0
    return float(1.0 - ((y - yhat) ** 2).sum() / ss)


def _r2_log(y, yhat):
    """在对数空间比较（要求 y, yhat 同号且非零）。"""
    s = np.sign(y[0]) if y[0] != 0 else 1.0
    ly, lh = np.log(np.abs(y)), np.log(np.abs(yhat) + 1e-300)
    return _r2(ly, lh)


def _lstsq(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(coef[0]), float(coef[1])


# ------------------------------------------------------------------ 曲线分类
def classify_curve(x, y, tol: float = TOL, log_space: bool = True) -> dict:
    """
    对 (x, y) 做结构分类。返回 {'label', 'params', 'candidates', 'n_used', 'dropped_leading_zeros'}。
    """
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 5:
        return {"label": "too-short", "n_used": int(len(x)), "candidates": {}}

    lead = 0
    nz = np.nonzero(y != 0)[0]
    if len(nz) == 0:
        return {"label": "identically-zero", "n_used": int(len(y)), "candidates": {}}
    if nz[0] > 0:
        lead = int(nz[0])
        x, y = x[lead:], y[lead:]

    # 拟合统一在 x>0 的子集上做（幂律/对数需要 x>0；让所有模型用同一子集才可比）
    dropped_nonpos = 0
    if np.any(x <= 0) and int((x > 0).sum()) >= 5:
        keep = x > 0
        dropped_nonpos = int((~keep).sum())
        x, y = x[keep], y[keep]

    # --- 奇偶零结构（只看丢零之后仍存在的零） ---
    zs = np.nonzero(y == 0)[0]
    parity = None
    if len(zs) >= 2:
        if np.all(zs % 2 == 0):
            parity = "even-zeros"
        elif np.all(zs % 2 == 1):
            parity = "odd-zeros"

    out = {"n_used": int(len(x)), "dropped_leading_zeros": lead,
           "dropped_nonpositive_x": dropped_nonpos,
           "x_range": [float(x[0]), float(x[-1])], "candidates": {}}
    if parity:
        out["parity"] = parity

    if len(x) < 5:
        out["label"] = "parity-alternating" if parity else "too-short"
        return out

    cands: dict[str, dict] = {}
    pos = bool(np.all(y > 0))

    # 0) constant
    c0, c1 = _lstsq(np.zeros_like(x), y)
    yhat = np.full_like(y, c0)
    cands["constant"] = {"params": {"a": c0}, "r2": _r2_log(y, yhat) if (pos and log_space) else _r2(y, yhat)}

    if pos:
        ly = np.log(y)
        xp = x > 0
        # 1) logarithmic  y = a + b ln x   （只用 x>0 的点）
        if xp.sum() >= 5:
            xq, yq = x[xp], y[xp]
            a, b = _lstsq(np.log(xq), yq)
            yh = a + b * np.log(xq)
            cands["logarithmic"] = {"params": {"a": a, "b": b, "n_fit": int(xp.sum())},
                                    "r2": _r2_log(yq, yh) if log_space else _r2(yq, yh)}
        # 2) power law  y = A x^b
        if xp.sum() >= 5:
            xq, yq = x[xp], y[xp]
            la, b = _lstsq(np.log(xq), np.log(yq))
            A = float(np.exp(la))
            yh = A * xq ** b
            cands["power-law"] = {"params": {"A": A, "b": b, "n_fit": int(xp.sum())},
                                  "r2": _r2_log(yq, yh) if log_space else _r2(yq, yh)}
        # 3) exponential  y = A e^{bx}
        la, b = _lstsq(x, ly)
        A = float(np.exp(la))
        cands["exponential"] = {"params": {"A": A, "b": b},
                                "r2": _r2_log(y, (A * np.exp(b * x))) if log_space else _r2(y, A * np.exp(b * x))}
        # 4) stretched  y = A e^{-c x^beta}（只在 x>=0 时试）
        if np.all(x >= 0):
            best = None
            for beta in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
                la, bb = _lstsq(x ** beta, ly)
                A2 = float(np.exp(la))
                r2 = _r2_log(y, A2 * np.exp(bb * x ** beta))
                if best is None or r2 > best["r2"]:
                    best = {"params": {"A": A2, "c": -bb, "beta": beta}, "r2": r2}
            cands["stretched"] = best

    # --- 选优 ---
    def score(name, d):
        return d["r2"]

    best_name = max(cands, key=lambda k: score(k, cands[k]))
    best_r2 = cands[best_name]["r2"]
    for name in sorted(cands, key=lambda k: _ORDER[k]):
        if _ORDER[name] < _ORDER[best_name] and cands[name]["r2"] >= best_r2 - tol:
            best_name = name
            break

    out["label"] = best_name
    out["params"] = cands[best_name]["params"]
    out["r2"] = cands[best_name]["r2"]
    out["candidates"] = {k: {"r2": round(v["r2"], 6), **{kk: (round(vv, 8) if isinstance(vv, float) else vv)
                                                         for kk, vv in v["params"].items()}}
                         for k, v in cands.items()}

    # --- 残差周期性 ---
    pm = cands.get("power-law")
    if pm is not None:
        yhat = pm["params"]["A"] * x ** pm["params"]["b"]
        res = y / yhat - 1.0
        if np.all(np.isfinite(res)) and len(res) > 12:
            f = np.abs(np.fft.rfft(res - res.mean())) ** 2
            if f[1:].sum() > 0 and f[1:].max() / f[1:].sum() > 0.30:
                P = len(res) / (np.argmax(f[1:]) + 1)
                out["residual_period"] = round(float(P), 4)
    if parity:
        out["label"] = out["label"] + "+" + parity
    return out


# ---------------------------------------------------------------- 分布分类
def classify_density(centers, dens, tol: float = 0.02, n_eff=None) -> dict:
    """
    对一条**密度曲线** (centers, dens) 做形状分类：
    delta / gapped / continuum×{edge-divergent, edge-vanishing, center-divergent}。
    密度只需相对值（只用比值与指数）。
    """
    c = np.asarray(centers, float).ravel()
    d = np.asarray(dens, float).ravel()
    m = np.isfinite(c) & np.isfinite(d)
    c, d = c[m], d[m]
    if len(c) < 10 or d.sum() <= 0:
        return {"label": "too-short", "n": int(len(c))}
    lo, hi = float(c[0]), float(c[-1])
    width = hi - lo
    out = {"n": n_eff if n_eff is not None else int(len(c)), "band": [lo, hi],
           "bandwidth": float(width)}
    if width <= 0:
        out["label"] = "delta"
        return out
    nb = len(c)
    if d.max() / d.sum() > 0.5:
        out["label"] = "delta"
        out["params"] = {"dominant_fraction": round(float(d.max() / d.sum()), 4)}
        return out
    empty = d <= 0
    cur = best = 0
    for i in range(1, nb - 1):
        cur = cur + 1 if empty[i] else 0
        best = max(best, cur)
    gap_w = best * width / nb
    out["max_interior_gap"] = round(float(gap_w), 6)
    if gap_w > tol * width:
        out["label"] = "gapped"
        out["params"] = {"gap_width": round(float(gap_w), 6),
                         "gap_rel": round(float(gap_w / width), 4)}
        return out

    def edge_alpha(side):
        m = max(4, nb // 10)
        if side == "hi":
            cen, de, delta = c[-m:], d[-m:], hi - c[-m:]
        else:
            cen, de, delta = c[:m], d[:m], c[:m] - lo
        ok = (de > 0) & (delta > 0)
        if ok.sum() < 4:
            return None
        return float(np.polyfit(np.log(delta[ok]), np.log(de[ok]), 1)[0])

    a_hi, a_lo = edge_alpha("hi"), edge_alpha("lo")
    out["edge_exponent"] = {"lo": None if a_lo is None else round(a_lo, 3),
                            "hi": None if a_hi is None else round(a_hi, 3)}
    med = float(np.median(d[d > 0]))
    cen_bin = float(d[nb // 2 - 2: nb // 2 + 3].mean())
    center_div = bool(med > 0 and cen_bin / med > 3.0)
    edge_div = any(a is not None and a < -0.2 for a in (a_lo, a_hi))
    edge_vanish = all(a is not None and a > 0.2 for a in (a_lo, a_hi))
    mods = []
    if edge_div:
        mods.append("edge-divergent")
    if edge_vanish:
        mods.append("edge-vanishing")
    if center_div:
        mods.append("center-divergent(van-Hove)")
    out["label"] = "continuum" + ("/" + "+".join(mods) if mods else "/smooth")
    out["params"] = {"center_over_median": round(cen_bin / med, 3) if med > 0 else None}
    return out


def classify_spectrum(ev, nbins: int = 200, tol: float = 0.02) -> dict:
    """对一组本征值做分布形状分类（先直方图，再交给 classify_density）。"""
    ev = np.asarray(ev, float).ravel()
    ev = ev[np.isfinite(ev)]
    if len(ev) < 8:
        return {"label": "too-short", "n": int(len(ev))}
    lo, hi = float(ev.min()), float(ev.max())
    if hi - lo <= 0:
        return {"label": "delta", "n": int(len(ev)), "band": [lo, hi], "bandwidth": 0.0}
    nbins = int(min(nbins, max(8, len(ev) // 4)))   # 防止 bin 比本征值还多造成的假 gap
    hist, edges = np.histogram(ev, bins=nbins, range=(lo, hi))
    frac_single = float(hist.max()) / len(ev)
    if frac_single > 0.5:
        return {"label": "delta", "n": int(len(ev)), "band": [lo, hi], "bandwidth": hi - lo,
                "params": {"dominant_fraction": round(frac_single, 4)}, "nbins": nbins}
    c = 0.5 * (edges[1:] + edges[:-1])
    out = classify_density(c, hist / max(1, len(ev)), tol=tol, n_eff=int(len(ev)))
    out["nbins"] = nbins
    return out


def summarize(records: list[dict]) -> dict:
    """按结构标签聚类。"""
    from collections import defaultdict
    clusters = defaultdict(list)
    for r in records:
        clusters[r["structure"]["label"]].append(r["id"])
    return {k: sorted(v) for k, v in sorted(clusters.items(), key=lambda kv: -len(kv[1]))}
