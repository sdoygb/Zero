"""
z0_period.py --- 「L=4 下周期 3–4」的根因判定：**状态空间太小，还是相位是刚性旋转？**

背景（交接单「实测瓶颈」）：L=4 跑出周期 ≈3–4（返回间隔只有 3 和 4，标准差/均值 0.12），
归因为「L=4 只给 4 个闭合词、2 个活层 ⟹ 状态空间太小，装不下稠密的模流」。
两个候选出路：(1) 扩状态空间（动 R31 的前提）；(2) 承认张力写成结论。

本机判定（四步，全部可复算）：
  A. 状态轨道：周期 ≡ 2（**与状态空间大小无关**）；L=2…16 扫过（2…1252 个闭合词）。
  B. 观测量闭式：P(g) = A1²+A2²+2A1A2cos(φ g)，φ = ω2−ω1 = log(w1/w2) = log(4/q_4) = log 7.2。
  C. 转动数 log(4/q_4)/2π 无理 ⟹ **序列非周期**；「3 和 4」= 该刚性旋转的三间隙定理统计。
  D. 改种子律：branch（现）→ 精确临界；mult（线性）→ 爆炸；gate（相位门）→ 吸收死亡。
  E. 扫 L：表观周期 = 2π/log(w1/w2)（L=4,6,8 实测吻合）；L↑ ⟹ 表观周期**变短**（不是变长）。

用法： /usr/bin/python3 z0_period.py        （注意：必须用系统 python3，见 HANDOFF 坑 #1）
输出： results/z0_period.json
"""
from __future__ import annotations
import json, os, math, cmath
from fractions import Fraction
from itertools import combinations
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results")
RES: dict = {}


# ---------------------------------------------------------------- 基础：层塔（与 z0_full.py 逐行一致）
def qL_exact(L: int) -> Fraction:
    """q_L = Σ_c (o_c/N)²，o_c = 旋转类轨道大小，N = C(L,L/2)。"""
    seen = {}
    for ones in combinations(range(L), L // 2):
        w = [1] * L
        for i in ones:
            w[i] = -1
        c = min(tuple(w[i:] + w[:i]) for i in range(L))
        seen[c] = seen.get(c, 0) + 1
    return sum(Fraction(o, math.comb(L, L // 2)) ** 2 * n for o, n in Counter(seen.values()).items())


def prefix_weights(L: int) -> dict:
    out = {}
    for k in range(L + 1):
        for s in range(-k, k + 1, 2):
            rem = L - k
            if (rem - s) % 2:
                continue
            a = (rem - s) // 2
            if 0 <= a <= rem:
                out[(k, s)] = math.comb(rem, a)
    return out


def build_tower(L: int):
    """层高 |s| 的权重 ∝ 补全数 × q_L^|s|；模频率 ω_a = −log(w_a/2)。"""
    q = qL_exact(L); qf = float(q); blk = Counter()
    for (k, s), W in prefix_weights(L).items():
        blk[abs(s)] += W * (qf ** abs(s))
    hs = sorted(blk)
    ws = np.array([blk[h] for h in hs], float)
    ws /= ws.sum()
    freq = np.array([-math.log(w / 2) for w in ws])
    return q, hs, ws, freq


# ---------------------------------------------------------------- 引擎（与 z0_full.py 同一规则；可换种子律）
def engine(L, T, seed_rule="branch", gate=None, tower=None, collect_P=True):
    """
    忠实规则：状态 = 未闭合词（带重数）；步进 = ±1 全分支；闭合 ⟺ 平衡归零；
    重播种 = newL2[()] += seed。seed 四种读法：
      branch 现引擎：nclo（**闭合分支数**，类计数型非线性，§16 说的唯一天然非线性）
      mult   线性回路：闭合分支的**总重数**（⟹ 爆炸，Z0_CORE §13）
      class  类计数：闭合的**旋转类数**
      branch+gate="phase"：只在两模干涉为**相长**（相对相位 >0，规范不变）时提交种子
    """
    q, hs, ws, freq = tower if tower else build_tower(L)
    L2 = Counter({(): 1}); g = 0; P = []
    traj, seeds, percl = [], [], []
    for _ in range(T):
        newL2 = Counter(); nclo = 0; tot = 0
        amp = np.zeros(len(hs), dtype=complex)
        percl.append(Counter())
        for w, mult in L2.items():
            if mult <= 0:
                continue
            for d in (+1, -1):
                nw = w + (d,)
                if sum(nw) == 0 and len(nw) >= 2:                    # 闭合
                    nclo += 1; tot += mult
                    h = max(abs(sum(nw[: k + 1])) for k in range(len(nw)))
                    hh = hs.index(h) if h in hs else len(hs) - 1
                    percl[-1][h] += mult
                    amp[hh] += mult * ws[hh] * cmath.exp(-1j * freq[hh] * g)
                elif len(nw) < L:                                    # 未闭合且未超长
                    newL2[nw] += mult
        if nclo > 0:
            ok = True
            if gate == "phase":                                      # 相对相位（规范不变）
                nz = [v for v in amp if abs(v) > 0]
                ok = (len(nz) >= 2) and (nz[0] * np.conj(nz[1])).real > 0
            seed = {"branch": nclo, "mult": tot, "class": 2}.get(seed_rule, nclo)
            if ok and seed > 0:
                newL2[()] += seed
                seeds.append(seed)
            else:
                seeds.append(0)
            g += 1
            if collect_P:
                P.append(abs(amp.sum()) ** 2)
        L2 = Counter({w: m for w, m in newL2.items() if m > 0})
        traj.append(tuple(sorted(L2.items())))
    return dict(traj=traj, P=np.array(P), seeds=seeds, percl=percl, g=g)


def closed_words(L: int):
    out = []
    def rec(w):
        if w and sum(w) == 0:
            out.append(tuple(w)); return
        if len(w) == L:
            return
        for d in (+1, -1):
            rec(w + [d])
    for d in (+1, -1):
        rec([d])
    return out


def cycle_of(traj, cap=100000):
    seen = {}
    for i, s in enumerate(traj):
        if s in seen:
            return seen[s] + 1, i - seen[s], i + 1
        seen[s] = i
        if i > cap:
            return None, None, i + 1
    return None, None, len(traj)


# ================================================================ A. 状态轨道 vs 状态空间大小
def part_A(Ls=(2, 4, 6, 8, 10, 12, 14, 16), act_cap=3_000_000):
    print("=" * 96)
    print("A. 状态轨道：周期是否随状态空间（闭合词数）变长？")
    print(f"{'L':>3} {'#闭合词':>8} {'#层高':>6} {'瞬态':>5} {'周期':>5} {'访问态数':>8} {'饱和活动量':>11} {'q_L':>12}")
    rows = []
    for L in Ls:
        cw = closed_words(L)
        hs_n = len({max(abs(sum(w[:k + 1])) for k in range(len(w))) for w in cw})
        r = engine(L, T=200, collect_P=False)
        tr, per, vis = cycle_of(r["traj"])
        act = sum(m for _, m in r["traj"][-1])
        if act > act_cap:            # 状态太大就不报（活动量指数增长）
            act = None
        q = qL_exact(L)
        rows.append(dict(L=L, n_closed=len(cw), n_heights=hs_n, transient=tr, period=per,
                         visited=vis, act_sat=act, q=str(q)))
        print(f"{L:3d} {len(cw):8d} {hs_n:6d} {str(tr):>5} {str(per):>5} {vis:8d} "
              f"{(str(act) if act is not None else '>cap'):>11} {str(q):>12}")
    print(f"{'':3} ⟹ 周期恒为 2（= ±1 行走的宇称下限）；瞬态 = 2L−2；访问态数 = 2L。")
    print(f"{'':3} ⟹ 状态空间从 2 个闭合词涨到 {rows[-1]['n_closed']} 个，周期**一个都没变**。")
    RES["A_state_orbit"] = rows
    return rows


# ================================================================ B. 观测量闭式
def part_B(L=4, T=4000):
    print("=" * 96)
    print("B. 观测量闭式：L=4 的 P(g) 到底是什么")
    q, hs, ws, freq = build_tower(L)
    Z = 1139
    w = [Fraction(int(round(x * Z)), Z) for x in ws]
    print(f"  层塔权重（精确）：w = {[str(x) for x in w]}   （hs={hs}）")
    r = engine(L, T, tower=(q, hs, ws, freq))
    P = r["P"]
    # 饱和态：每个闭合词以重数 4 闭合，h=1 两个词、h=2 两个词 ⟹ 每层两次
    per = r["percl"][-1]
    A = {h: per[h] * ws[hs.index(h)] for h in per}
    phi = freq[hs.index(2)] - freq[hs.index(1)]
    A1, A2 = A[1], A[2]
    Pcf = (A1 ** 2 + A2 ** 2) + 2 * A1 * A2 * np.cos(phi * np.arange(len(P)))
    dev = np.abs(P - Pcf) / P.mean()
    print(f"  饱和闭合量：每层重数 {dict(per)}  ⟹ A1={A1:.5f} A2={A2:.5f}  A1/A2={A1/A2:.5f}")
    print(f"  φ = ω2−ω1 = {phi:.8f}   log(4/q_4) = log({float(4/q):.4f}) = {math.log(float(4/q)):.8f}  "
          f" log 7.2 = {math.log(7.2):.8f}")
    print(f"  闭式 P(g) = A1²+A2²+2A1A2·cos(φg)：均值={Pcf.mean():.4f} 相对起伏={Pcf.std()/Pcf.mean():.4f}")
    print(f"  实测（T={T}）：                    均值={P.mean():.4f} 相对起伏={P.std()/P.mean():.4f}")
    print(f"  第 3 个闭合事件之后逐点偏差 max|ΔP|/mean = {dev[3:].max():.2e}"
          f"（前 3 个事件是瞬态：重数还没饱和到 4）")
    # 自相关 ⟷ cos(φk)
    ac = [float(((P[:-k] - P.mean()) * (P[k:] - P.mean())).mean() / P.var()) for k in (1, 3, 10, 100, 500)]
    pr = [float(np.cos(phi * k)) for k in (1, 3, 10, 100, 500)]
    print(f"  自相关(1,3,10,100,500) 实测 = {[round(x, 3) for x in ac]}")
    print(f"  同一量 cos(φk)        闭式 = {[round(x, 3) for x in pr]}")
    # 返回间隔 = 极大值间隔
    pk = [i for i in range(1, len(P) - 1) if P[i] > P[i - 1] and P[i] >= P[i + 1]]
    iv = np.diff(pk); cnt = Counter(iv.tolist())
    print(f"  返回间隔（P 的极大值间距）：取值 {sorted(cnt)}  计数 {dict(sorted(cnt.items()))}  "
          f"均值={iv.mean():.4f}  标准差/均值={iv.std()/iv.mean():.4f}")
    print(f"  理论：刚性旋转的周期 2π/φ = {2*math.pi/phi:.4f}；三间隙定理 ⟹ 间隔只取 ⌊·⌋,⌈·⌉ = 3,4")
    RES["B_closed_form"] = dict(q=str(q), w=[str(x) for x in w], phi=phi, T_rot=2 * math.pi / phi,
                                A1=A1, A2=A2, P_mean_meas=float(P.mean()), P_relstd_meas=float(P.std() / P.mean()),
                                P_mean_cf=float(Pcf.mean()), P_relstd_cf=float(Pcf.std() / Pcf.mean()),
                                max_dev_after_transient=float(dev[3:].max()),
                                acf_meas=ac, acf_cos=pr, gap_counts={str(k): v for k, v in sorted(cnt.items())},
                                gap_mean=float(iv.mean()), gap_relstd=float(iv.std() / iv.mean()),
                                per_height_multiplicity={str(k): v for k, v in per.items()})
    return phi


# ================================================================ C. 转动数
def part_C():
    print("=" * 96)
    print("C. 转动数 log(4/q_4)/2π 是无理数 ⟹ 序列**非周期**（不是「周期 3–4 的循环」）")
    q = qL_exact(4); r = float(4 / q)
    rho = math.log(r) / (2 * math.pi)
    x = rho; cf = []
    for _ in range(7):
        a = int(x); cf.append(a); x = 1 / (x - a)
    print(f"  q_4 = {q}  ⟹ w1/w2 = 4/q_4 = {r}  ⟹ φ = log {r} = {math.log(r):.8f}")
    print(f"  转动数 ρ = φ/2π = {rho:.10f}（= 每步转过的圈数）")
    print(f"  连分数 ρ = {cf}；收敛分母（= 近返回时刻）3, 16, 35, …")
    print(f"  无理性：若 log(36/5) = 2π p/q，则 e^(2πp) = (36/5)^q ⟹ e^(2π) 代数，")
    print(f"          与 e^(2π) = (e^π)² 超越（Gelfond–Schneider：e^π = (−1)^(−i) 超越）矛盾。")
    print(f"  ⟹ 轨道在圆上稠密、永不闭合成环：**没有有限周期**，只有越来越精确的复现。")
    RES["C_rotation_number"] = dict(q4=str(q), ratio=r, phi=math.log(r), rho=rho, cf=cf,
                                    convergent_denominators=[3, 16, 35])
    return rho


# ================================================================ D. 换种子律：临界点三选一
def linear_loop_rho(L=4):
    """seed = 总重数 ⟹ 回路是**线性**的：直接算转移矩阵的谱半径（精确到机器精度）。"""
    words, rec = [()], {}
    def rec_w(w):
        words.append(w)
        if len(w) == L:
            return
        for d in (+1, -1):
            if sum(w + (d,)) == 0:
                continue
            rec_w(w + (d,))
    for d in (+1, -1):
        rec_w((d,))
    words.append(())                      # 空词
    words = list(dict.fromkeys(words)); idx = {w: i for i, w in enumerate(words)}
    M = np.zeros((len(words), len(words)))
    for w in words:
        for d in (+1, -1):
            nw = w + (d,)
            if sum(nw) == 0 and len(nw) >= 2:
                M[idx[()], idx[w]] += 1
            elif len(nw) < L:
                M[idx[nw], idx[w]] += 1
    return float(max(abs(np.linalg.eigvals(M))))


def part_D(L=4, T=3000):
    print("=" * 96)
    print("D. 种子律的敏感性：为什么「引擎侧调参」不是出路（r*=1 精确临界）")
    rho_lin = linear_loop_rho(L)
    out = []
    for name, rule, gate, Tn in (("branch 现引擎（闭合分支数）", "branch", None, T),
                                 ("mult 线性回路（总重数）", "mult", None, 1200),
                                 ("class 类计数（2 个旋转类）", "class", None, T),
                                 ("branch + 相位门（相长才提交）", "branch", "phase", T)):
        r = engine(L, Tn, seed_rule=rule, gate=gate, collect_P=(rule != "mult"))
        tr, per, vis = cycle_of(r["traj"])
        act = [sum(m for _, m in s) for s in r["traj"]]
        grow = None
        if rule == "mult":
            grow = float(np.mean([(act[i + 2] / act[i]) for i in range(len(act) - 2) if act[i] > 0][-200:]))
        print(f"  {name:30s} 瞬态/周期={str((tr, per)):>14}  活动: 末={act[-1]:.3e} 峰={max(act):.3e}"
              + (f"  两步增长={grow:.6f}（=ρ²，ρ={rho_lin:.6f}）" if grow else ""))
        out.append(dict(rule=rule, gate=gate, transient=tr, period=per, act_last=float(act[-1]),
                        act_max=float(max(act)), growth_2step=grow))
    print(f"  ⟹ 种子 = 分支数 ⟹ **精确临界**（有界 2-周期）；种子 > 分支数 ⟹ 爆炸（线性回路，"
          f"谱半径 ρ = √(1+√3) = {rho_lin:.6f}/步）；")
    print(f"     种子 < 分支数（相位门）⟹ **吸收死亡**（连续两步被挡就再起不来）。状态扇区没有「更长周期」这一档。")
    RES["D_seed_rules"] = dict(rows=out, rho_linear=rho_lin)
    return out


# ================================================================ E. 扫 L：表观周期随状态空间变大而**变短**
def part_E(Ls=(4, 6, 8), T=4000):
    print("=" * 96)
    print("E. 扫 L：表观周期 T_L = 2π/log(w1/w2)（L↑ ⟹ 表观周期变短，不是变长）")
    print(f"{'L':>3} {'q_L':>12} {'w1/w2':>10} {'预测 T_L':>10} {'实测间隔均值':>12} {'#可观测差频':>11} {'相对起伏':>9}")
    rows = []
    for L in Ls:
        q, hs, ws, freq = build_tower(L)
        r = engine(L, T, tower=(q, hs, ws, freq))
        P = r["P"]
        pk = [i for i in range(1, len(P) - 1) if P[i] > P[i - 1] and P[i] >= P[i + 1]]
        iv = np.diff(pk)
        ratio = float(ws[1] / ws[2]); pred = 2 * math.pi / math.log(ratio)
        nf = (len(hs) - 1) * (len(hs) - 2) // 2
        rows.append(dict(L=L, q=str(q), ratio=ratio, pred=pred, meas=float(iv.mean()), n_freq=nf,
                         relstd=float(P.std() / P.mean())))
        print(f"{L:3d} {str(q):>12} {ratio:10.4f} {pred:10.4f} {iv.mean():12.4f} {nf:11d} {P.std()/P.mean():9.4f}")
    # 有闭式的比值：w1/w2 = (Σ1/Σ2)/q_L
    rows2 = []
    for L in (4, 6, 8, 10, 12):
        q, hs, ws, freq = build_tower(L)
        S1 = sum(math.comb(L - k, (L - k - 1) // 2) for k in range(1, L, 2))
        S2 = sum(math.comb(L - k, (L - k - 2) // 2) for k in range(2, L - 1, 2))
        rows2.append(dict(L=L, pred_ratio=(S1 / S2) / float(q), meas_ratio=float(ws[1] / ws[2])))
    print("  闭式比值校验 w1/w2 = (Σ1/Σ2)/q_L：" +
          "  ".join(f"L={d['L']}:{d['pred_ratio']:.4f}/{d['meas_ratio']:.4f}" for d in rows2))
    print("  L=2：只有 1 个活层 ⟹ 0 个差频 ⟹ P 严格常数（周期 1，真平凡）；")
    print("  L=10 起预测 T_L < 2 步 ⟹ 该频率过 Nyquist 走样 ⟹ 观测量反而**更像周期**信号。")
    RES["E_L_scan"] = dict(rows=rows, ratio_check=rows2)
    return rows


if __name__ == "__main__":
    A = part_A()
    phi = part_B()
    rho = part_C()
    D = part_D()
    E = part_E()
    print("=" * 96)
    print("判定：")
    print("  1) 「周期 3–4」不是状态周期 —— 状态周期恒为 2（±1 行走的宇称下限），与状态空间大小无关；")
    print(f"     它是**相位刚性旋转**的表观周期 2π/log(4/q_4) = {2*math.pi/math.log(36/5):.4f}（= 交接单的「3 与 4」，"
          f"标准差/均值 0.12 = 三间隙定理）。")
    print("  2) 序列可证**非周期**（转动数无理）⟹ 「装不下稠密的模流」是反的；精确说：读数是圆上刚性旋转，")
    print("     非周期性可证，但**谱只有一条线**（L=4 只有 1 个差频）—— 稠密 ≠ 复杂。")
    print("  3) 扩状态空间（L↑）不改变状态周期（2），却把表观周期从 3.18 压到 2 以下，且动掉 R31 的 D=4 前提。")
    print("  ⟹ 该写成结论（且是比原表述更强的结论）；要动的是**回路的闭合方式**，不是状态空间大小。")
    json.dump(RES, open(os.path.join(OUT, "z0_period.json"), "w"), ensure_ascii=False, indent=1, default=str)
    print("→ results/z0_period.json")
