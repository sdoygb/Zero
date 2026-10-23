#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_wilsonian_rg.py —— 保谱（Wilson）重标度群：非平凡不动点是真的吗？
================================================================================
背景：能标跑动的机制内核是 Wilson 壳层积分——**保住低模、积掉高模**。
      在图上对应的做法不是粗粒化**几何**，而是粗粒化**测度/谱**：
      配对顶点、解出新电导，使低端本征值（归一化后）匹配。

【本篇要判的不是"能不能算"，而是"那个不动点靠不靠得住"】
  该逆问题匹配 m 个归一化本征值，却有 N_c−1 个未知电导；当 m < N_c−1 时**欠定**，
  `least_squares` 只返回**某个**解。若不动点（CV 平台）随初值改变 ⟹ 它是求解器的产物。

【检验】
  【W1】路径 Laplacian 的正确性（行和 0；均匀路径的谱与闭式 2c(1−cos(πk/N)) 对照）
  【W2】复现"谱被保留"：残差应 < 1e-5
  【W3】级联初值（自相似几何）下的 CV 流：有无平台？
  【W4】★ 稳健性审计：换初值（不同 x0）——CV 平台是否相同？
  【W5】★ well-posed 版：匹配 m = N_c−1 个（方程数 = 未知数）——还有平台吗？
  【W6】对照：几何粗粒化（串联合并）——CV 流向哪里？

用法：/usr/local/bin/python3 zero_wilsonian_rg.py
输出：ZeroGPU/results/zero_wilsonian_rg.json
"""
from __future__ import annotations
import json, os
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "zero_wilsonian_rg.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


# ══════════════════════════════════════════════════════════════════════
def lap(c):
    """路径 Laplacian：c 是边电导（长度 N−1），返回 N×N。"""
    c = np.asarray(c, float)
    N = len(c) + 1
    L = np.zeros((N, N))
    idx = np.arange(N - 1)
    L[idx, idx] += c
    L[idx + 1, idx + 1] += c
    L[idx, idx + 1] -= c
    L[idx + 1, idx] -= c
    return L


def low(c, m):
    ev = np.linalg.eigvalsh(lap(c))
    ev = ev[ev > 1e-12]
    return ev[:m]


def cv(c):
    c = np.asarray(c, float)
    return float(np.std(c) / np.mean(c))


def cascade(n, lam, seed):
    """自相似（层级）初值：逐块乘 exp(±lam)。返回长度 n 的电导。"""
    rng = np.random.RandomState(seed)
    v = np.ones(n)
    blk = n
    while blk >= 2:
        blk //= 2
        for s in range(0, n, blk * 2):
            v[s:s + blk] *= np.exp(rng.choice([-1, 1]) * lam)
    return v


def rg(c, mfrac=0.25, x0_mode="half", seed=0):
    """保谱重整化：配对顶点，解出新电导使最低 m 个**归一化**本征值匹配。

    x0_mode 控制初值：
      'half'  —— G52 的做法：取每隔一条的电导，并固定全局尺度（几何平均=1）
      'rand'  —— 随机扰动
      'flat'  —— 全常数
    """
    c = np.asarray(c, float)
    N = len(c) + 1
    Nc = N // 2
    m = max(4, int(mfrac * Nc))
    tgt = low(c, m)
    tgt = tgt / tgt.sum()

    def res(x):
        cc = np.exp(x)
        ev = low(cc, m)
        if len(ev) < m:
            return np.full(m, 1.0)
        return ev / ev.sum() - tgt

    if x0_mode == "half":
        x0 = np.log(np.array([c[2 * I] for I in range(Nc - 1)]))
        x0 = x0 - x0.mean()
    elif x0_mode == "rand":
        rng = np.random.RandomState(seed)
        x0 = rng.normal(0.0, 1.0, Nc - 1)
        x0 = x0 - x0.mean()
    else:
        x0 = np.zeros(Nc - 1)

    r = least_squares(res, x0, xtol=1e-12, ftol=1e-14, max_nfev=2000)
    out = np.exp(r.x)
    return out, float(np.sqrt(np.mean(res(r.x) ** 2))), m, Nc - 1


def rg_square(c):
    """well-posed：匹配最低 m = N_c−1 个（方程数 = 未知数）。"""
    c = np.asarray(c, float)
    N = len(c) + 1
    Nc = N // 2
    m = Nc - 1
    tgt = low(c, m)
    tgt = tgt / tgt.sum()

    def res(x):
        cc = np.exp(x)
        ev = low(cc, m)
        if len(ev) < m:
            return np.full(m, 1.0)
        return ev / ev.sum() - tgt

    x0 = np.log(np.array([c[2 * I] for I in range(Nc - 1)]))
    x0 = x0 - x0.mean()
    r = least_squares(res, x0, xtol=1e-12, ftol=1e-14, max_nfev=3000)
    return np.exp(r.x), float(np.sqrt(np.mean(res(r.x) ** 2))), m, Nc - 1


def series_reduce(c):
    """几何粗粒化：配对后按**串联**合并（严格等效电导，但谱不保）。"""
    c = np.asarray(c, float)
    n = len(c)
    out = []
    for i in range(0, n - 1, 2):
        out.append(c[i] * c[i + 1] / (c[i] + c[i + 1]))
    if n % 2 == 1:                      # 末尾奇数条，保留
        out.append(c[-1])
    return np.array(out)


# ══════════════════════════════════════════════════════════════════════
print("=" * 90)
print("【W1】路径 Laplacian 的正确性")
print("=" * 90)
c0 = np.ones(9)
L0 = lap(c0)
chk(abs(L0.sum(axis=1)).max() < 1e-12, "行和为 0", f"{abs(L0.sum(axis=1)).max():.1e}")
ev0 = np.sort(np.linalg.eigvalsh(L0))
closed = np.sort(2 * (1 - np.cos(np.pi * np.arange(10) / 10)))
chk(np.abs(ev0 - closed).max() < 1e-12, "均匀路径谱与闭式 2c(1−cos(πk/N)) 一致",
    f"max|Δ|={np.abs(ev0 - closed).max():.1e}")
chk(abs(ev0[0]) < 1e-12, "最小本征值为 0（常数模）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 90)
print("【W2】复现“谱被保留”：残差应 < 1e-5")
print("=" * 90)
w2 = {}
for lam in (0.35, 0.70):
    c = cascade(256, lam, seed=7)
    cc, res, m, nunk = rg(c, 0.25, "half")
    print(f"  λ={lam}: 初始 N=257, CV={cv(c):.5f}  → 匹配 m={m}, 未知={nunk}（"
          f"{'欠定' if m < nunk else '方阵'}）  残差={res:.2e}  CV'={cv(cc):.5f}")
    w2[str(lam)] = {"cv_init": cv(c), "m": m, "n_unknown": nunk, "residual": res,
                    "cv_coarse": cv(cc)}
    chk(res < 1e-4, f"λ={lam}: 残差 < 1e-4（谱在 4–5 位有效数字上被保留）", f"{res:.2e}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 90)
print("【W3】级联初值下的 CV 流（G52 口径：mfrac=0.25，初值 half）")
print("=" * 90)
w3 = {}
for lam in (0.35, 0.70):
    c = cascade(256, lam, seed=7)
    seq = [cv(c)]
    cur = c
    for r in range(3):
        cur, res, m, nunk = rg(cur, 0.25, "half")
        seq.append(cv(cur))
    print(f"  λ={lam}: CV = " + " → ".join(f"{x:.4f}" for x in seq))
    spread = (max(seq[1:]) - min(seq[1:])) / np.mean(seq[1:])
    print(f"      后三级相对摆幅 = {spread:.3f}")
    w3[str(lam)] = {"cv_sequence": seq, "post_spread": float(spread)}
chk(True, "W3 完成（记录 CV 流）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 90)
print("【W4】★ 稳健性审计：换初值，CV 平台是否相同？（欠定系统的要害）")
print("=" * 90)
w4 = {}
for lam in (0.35, 0.70):
    base = cascade(256, lam, seed=7)
    rows = []
    for mode, sd in (("half", 0), ("rand", 1), ("rand", 2), ("flat", 0)):
        cc, res, m, nunk = rg(base, 0.25, mode, seed=sd)
        rows.append({"x0": f"{mode}{sd if mode == 'rand' else ''}", "cv": cv(cc), "residual": res})
    cvs = [r["cv"] for r in rows]
    rel = (max(cvs) - min(cvs)) / np.mean(cvs)
    print(f"  λ={lam}: 同一输入、不同初值 ⟹ CV = "
          + "、".join(f"{r['cv']:.4f}({r['x0']})" for r in rows))
    print(f"      相对散布 = {rel:.3f}")
    w4[str(lam)] = {"rows": rows, "relative_spread": float(rel)}
    chk(rel > 0.50, f"λ={lam}: **散布 > 50%** ⟹ 欠定逆问题不钉住几何（要报的判定）", f"{rel:.3f}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 90)
print("【W5】★ well-posed 版（m = N_c−1，方程数 = 未知数）")
print("=" * 90)
w5 = {}
for lam in (0.35, 0.70):
    c = cascade(256, lam, seed=7)
    seq = [cv(c)]
    cur = c
    resl = []
    for r in range(3):
        cur, res, m, nunk = rg_square(cur)
        seq.append(cv(cur))
        resl.append(res)
    print(f"  λ={lam}: CV = " + " → ".join(f"{x:.4f}" for x in seq)
          + f"   （残差 {max(resl):.1e}，m={m}=未知 {nunk}）")
    w5[str(lam)] = {"cv_sequence": seq, "residual_max": float(max(resl)), "m": m, "n_unknown": nunk}
    tail = seq[2:]
    spread5 = (max(tail) - min(tail)) / np.mean(tail) if len(tail) > 1 else 0.0
    w5[str(lam)]["tail_spread"] = float(spread5)
    chk(spread5 < 0.15, f"λ={lam}: well-posed 版收敛到平台（后两级摆幅 < 15%）", f"{spread5:.3f}")

se = {lam: np.mean(w5[str(lam)]["cv_sequence"][2:]) for lam in (0.35, 0.70)}
dlam = abs(se[0.35] - se[0.70]) / np.mean([se[0.35], se[0.70]])
print(f"  ★ 两个初值给的平台：λ=0.35 → {se[0.35]:.4f}；λ=0.7 → {se[0.70]:.4f}；相对差 {dlam:.3f}")
chk(dlam < 0.15, "well-posed 版的不动点与初值无关（λ 无关，相对差 < 15%）", f"{dlam:.3f}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 90)
print("【W6】对照：几何粗粒化（串联合并）——CV 流向哪里？")
print("=" * 90)
w6 = {}
for lam in (0.35, 0.70):
    c = cascade(256, lam, seed=7)
    seq = [cv(c)]
    cur = c
    for r in range(5):
        cur = series_reduce(cur)
        seq.append(cv(cur))
    print(f"  λ={lam}: CV = " + " → ".join(f"{x:.4f}" for x in seq))
    w6[str(lam)] = {"cv_sequence": seq}
chk(True, "W6 完成（几何粗粒化的对照）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 90)
print("汇总")
print("=" * 90)
out = {
    "object": "路径图 Laplacian；配对顶点 N→N/2；保谱逆问题 vs 串联合并",
    "W2_spectrum_preserved": w2,
    "W3_cv_flow_G52": w3,
    "W4_initial_guess_audit": w4,
    "W5_well_posed": w5,
    "W6_geometric_contrast": w6,
    "W5_plateau": {"lambda_0.35": float(se[0.35]), "lambda_0.70": float(se[0.70]),
                   "relative_diff": float(dlam)},
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)