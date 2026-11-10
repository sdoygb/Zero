#!/usr/bin/env python3
"""序列可拷贝性 · 换口径：同种内拷贝（重复序列）的保真度 vs 一重判读的 Boltzmann 上限

上一轮用跨物种同源蛋白测"位点变异"，口径不合（测的是进化距离，不是复制保真度）。
本轮改用**同一基因组内的重复拷贝**（重复元件/散在重复），这才是"同一基因的一次复制事件"的记录。

测什么：
  · 重复拷贝的**成对一致度**（复制刚发生时的保真度）；
  · 同一族内一致度的**分层**（若有些拷贝完全一致或近完全一致，说明存在纠错/选择；若均匀铺开，说明一重判读即够）。
预测（命题 4）：一重判读下，拷贝间的保守程度应停在 Boltzmann 水平（≈ 1 − e^{-ΔG/kT}，ΔG 取 1–10 kcal/mol 时 10^-1–10^-7），
故**成对一致度不会集中在 99.9% 以上**，除非有第二重判读或选择。

实现：k-mer 哈希找重复片段（k=50），对每个重复族的拷贝做逐位比对，统计一致度分布。
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_repeat_fidelity.json")
GENOMES = {
    "yeast_chrXII": "/tmp/yeast12.fa",        # 1.08 Mb（含 ~100–150 拷贝 rDNA 串联重复）
    "yeast_chrIV": "/tmp/g_NC_001136.fa",     # 1.53 Mb（对照：无大规模串联重复）
    "human_mt": "/tmp/g_NC_012920.fa",        # 16.6 kb（对照：小基因组）
}
K = 50


def read_fasta_text(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def rolling_hashes(seq, k=K, step=50, seed=7):
    """滚动哈希找重复片段的起点：k-mer → 64 位哈希（随机碱基映射 + 多项式滚动），
    与"用字符串键"等价（同一 k-mer 同一哈希），但省内存。"""
    rng = np.random.default_rng(seed)
    code = {c: int(rng.integers(1, 2**62)) for c in "ACGTN"}
    arr = np.fromiter((code.get(c, 0) for c in seq), dtype=np.uint64, count=len(seq))
    # 多项式滚动哈希：H_i = Σ_{t<k} arr[i+t] * B^(k-1-t)，B 取大奇数
    B = np.uint64(1099511628211)
    kpow = np.uint64(1)
    for _ in range(k - 1):
        kpow = np.uint64((int(kpow) * int(B)) & 0xFFFFFFFFFFFFFFFF)
    h = np.zeros(len(seq) - k + 1, dtype=np.uint64)
    acc = np.uint64(0)
    for t in range(k):
        acc = np.uint64((int(acc) * int(B)) & 0xFFFFFFFFFFFFFFFF)
        acc ^= arr[t] if t < len(arr) else np.uint64(0)
    h[0] = acc
    # 递推：去掉首位、加上新位
    for i in range(1, len(h)):
        h[i] = np.uint64(((int(h[i - 1]) - int(arr[i - 1]) * int(kpow)) * int(B) + int(arr[i + k - 1])) & 0xFFFFFFFFFFFFFFFF)
    return h, arr


def find_repeats(seq, k=K, step=50, min_copies=2):
    h, arr = rolling_hashes(seq, k=k, step=step)
    hs = h[::step]
    cnt = Counter(int(x) for x in hs)
    kept = [i for i in range(0, len(h), step) if cnt[int(h[i])] >= min_copies]
    starts = sorted(kept)
    families, cur = [], []
    for a in starts:
        if cur and a - cur[-1] > 4 * step:
            if len(cur) >= 2:
                families.append(cur)
            cur = []
        cur.append(a)
    if len(cur) >= 2:
        families.append(cur)
    return families, seq


def extend_copy(seq, a, b, max_len=5000, allow=0):
    """从两个起点向右延伸，逐位比对，返回 (长度, 一致位数)。允许最多 allow 个错配用于延伸。"""
    i, j = a, b
    n = len(seq)
    same = 0
    ln = 0
    mism = 0
    while i < n and j < n and ln < max_len:
        if seq[i] == seq[j]:
            same += 1
        else:
            mism += 1
            if mism > allow:
                break
        ln += 1
        i += 1
        j += 1
    return ln, same


def analyse(path, name):
    seq = read_fasta_text(path)
    families, sub = find_repeats(seq)
    pairs = []
    for fam in families[:3000]:
        for a, b in [(fam[0], fam[1])] + ([(fam[0], fam[2])] if len(fam) > 2 else []) + \
                    ([(fam[1], fam[2])] if len(fam) > 2 else []):
            ln, same = extend_copy(sub, a, b)
            if ln >= 100:
                pairs.append((a, b, ln, same / ln))
    idents = np.asarray([p[3] for p in pairs])
    res = {
        "genome": name, "length": len(sub), "k": K,
        "n_families": len(families), "n_pairs": len(pairs),
        "identity_mean": float(idents.mean()) if len(idents) else None,
        "identity_frac_ge_999": float((idents >= 0.999).mean()) if len(idents) else None,
        "identity_frac_lt_99": float((idents < 0.99).mean()) if len(idents) else None,
        "identity_quantiles": [float(np.percentile(idents, q)) for q in (5, 25, 50, 75, 95)] if len(idents) else None,
        "examples": [[int(a), int(b), int(ln), round(float(x), 4)] for a, b, ln, x in pairs[:5]],
    }
    print(f"\n{name}（{len(sub)} bp）：重复族 {len(families)}，可比对拷贝对 {len(pairs)}")
    if len(idents):
        print(f"  一致度 均值 {idents.mean():.4f}；5/25/50/75/95 分位 "
              f"{[round(x,4) for x in np.percentile(idents, [5,25,50,75,95])]}")
        print(f"  完全一致(≥99.9%) {res['identity_frac_ge_999']:.1%}；<99% {res['identity_frac_lt_99']:.1%}")
        print(f"  样例（起点A, 起点B, 比对长, 一致度）：{res['examples']}")
        for dG in (1, 3, 5, 10):
            print(f"  一重判读对照：ΔG={dG} kcal/mol ⟹ e^(-ΔG/kT) ≈ {np.exp(-dG*4184/(8.314*300)):.2e}")
    return res


def main():
    rows = {}
    for name, path in GENOMES.items():
        if not Path(path).exists():
            print(f"[skip] {name}: 缺文件 {path}")
            continue
        rows[name] = analyse(path, name)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(rows, f, indent=2, ensure_ascii=False)
    print("\n落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())