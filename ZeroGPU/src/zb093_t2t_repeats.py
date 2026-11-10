#!/usr/bin/env python3
"""命题 4 的实测（第三口径）：T2T 组装里的多拷贝对象 → 拷贝保真度 vs 一重判读的 Boltzmann 上限

前两次口径失败（`_staging/推导/序列可拷贝性_纤维与校验.md` §3.5/§3.6）：
  · 跨物种同源蛋白测的是进化距离；
  · GRCh38 组装把串联重复折叠成单拷贝，没有多拷贝对象。
本轮换 T2T 组装（CHM13v2，含着丝粒/近端着丝粒的完整重复），对象是多拷贝片段本身。

记录层口径（与命题 1、推论 2 一致）：
  · f = "把片段压成 k-mer 指纹"（判读）；
  · g = "由指纹回写片段"（重建，片段即 k-mer 序列本身，可逆）；
  · 冗余 = 位点上存在 ≥2 个可互相校验的拷贝；一重判读的保真度上限 = Boltzmann 因子 e^{-ΔG/kT}。
测：多拷贝片段对的一致度分布——(i) 若出现 ≥99.9% 的拷贝对，说明冗余已被纠错/选择保护；
    (ii) 若一致度停在 Boltzmann 水平，说明一重判读即够。
"""
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_t2t_repeats.json")
SEQ = Path("/tmp/t2t21.txt")


def read_seq():
    return SEQ.read_text().strip().upper()


def periodic_runs(seq, unit=171, tol=0.08, win=2000):
    """找着丝粒型串联重复：滑动窗口内按 unit 周期比对，匹配率 ≥ 1−tol 的窗口。"""
    ok = []
    n = len(seq)
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    for start in range(0, n - win, win):
        w = arr[start:start + win]
        m = np.mean(w[:-unit] == w[unit:])
        if m >= 1 - tol:
            ok.append((start, float(m)))
    return ok


def copy_pairs(seq, unit, tol=0.08, max_pairs=40):
    """在已定位的串联重复区里取相邻拷贝对，逐位统计一致度。"""
    out = []
    n = len(seq)
    i = 0
    while i + 2 * unit <= n and len(out) < max_pairs:
        a = seq[i:i + unit]
        b = seq[i + unit:i + 2 * unit]
        same = sum(1 for x, y in zip(a, b) if x == y)
        ident = same / unit
        if ident >= 1 - tol:              # 只收周期结构内的拷贝
            out.append(round(ident, 4))
            i += unit
        else:
            i += unit
    return out


def scattering(seq, k=60, step=200, tol=0.10, max_hits=60):
    """散在重复：k-mer 抽样，找多处出现（≥2）的片段，取两拷贝做延伸比对。"""
    idx = {}
    n = len(seq)
    for i in range(0, n - k + 1, step):
        km = seq[i:i + k]
        if km in idx:
            idx[km].append(i)
        else:
            idx[km] = [i]
    multi = [(pos[0], pos[1]) for km, pos in idx.items() if len(pos) >= 2][:max_hits]
    idents = []
    for a, b in multi:
        ln = min(500, n - b)
        if ln < 100:
            continue
        same = sum(1 for x, y in zip(seq[a:a + ln], seq[b:b + ln]) if x == y)
        ident = same / ln
        if ident >= 1 - tol:              # 只收"真重复"
            idents.append(round(ident, 4))
    return idents


def main():
    if not SEQ.exists():
        print("[skip] 缺 /tmp/t2t21.txt"); return 1
    seq = read_seq()
    n = len(seq)
    out = {"genome": "T2T-CHM13v2 chr21 (NC_060945.1)", "length": n}

    runs = periodic_runs(seq)
    out["periodic_window_frac"] = len(runs) / max(1, n // 2000)
    print(f"序列 {n} bp；2000 bp 窗口中有周期结构的比例 {out['periodic_window_frac']:.2%}")
    if runs:
        print(f"  例：起点 {runs[0][0]}，匹配率 {runs[0][1]:.3f}；最后一个 {runs[-1][0]}")
        pairs = copy_pairs(seq, unit=171)
        out["copy_pairs"] = pairs
        print(f"  171 bp 单元相邻拷贝一致度（前 {len(pairs)} 对）：{pairs[:8]}…")
        if pairs:
            p = np.asarray(pairs)
            print(f"  一致度 均值 {p.mean():.4f}；≥0.999 的比例 {(p>=0.999).mean():.1%}；<0.99 {(p<0.99).mean():.1%}")

    sc = scattering(seq)
    out["scattered_pairs"] = sc
    print(f"\n散在重复（60-mer 抽样，200 bp 步长）：真重复对 {len(sc)}；一致度样例 {sc[:8]}")
    if sc:
        s = np.asarray(sc)
        print(f"  一致度 均值 {s.mean():.4f}；≥0.999 {(s>=0.999).mean():.1%}")

    print("\n一重判读的 Boltzmann 上限对照：")
    for dG in (1, 3, 5, 10):
        print(f"  ΔG={dG:>2} kcal/mol ⟹ e^(-ΔG/kT) ≈ {np.exp(-dG*4184/(8.314*300)):.1e}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("\n落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())