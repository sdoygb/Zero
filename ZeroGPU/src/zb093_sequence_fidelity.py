#!/usr/bin/env python3
"""序列可拷贝性 · 命题 4 的实测：一重判读的保真度上限 vs 跨物种位点变异性

命题 4（`_staging/推导/序列可拷贝性_纤维与校验.md`）：
  只有一重判读 f（碱基配对）时每个位的纤维大小 = 1、可纠错扰动集 E = {0}，
  复制误差停在 Boltzmann 因子 e^{-ΔG/kT}；实测错配率 10^-9~10^-10 ⟹ 必有第二重判读。

可测的一侧：同一蛋白在不同物种间的**位点变异性**就是"判读承受过的扰动"实录。
  · 若每个位都能独立地随机出错（无选择/无校验），跨物种的位点变异应"铺平"；
  · 若出现**高度分层的变异性谱**（大量完全保守位 + 少数快变位），则是位点级选择/校验在起作用。

实现：UniProt reviewed 的溶菌酶 C（LYZ）序列，**星形比对**到参照坐标系，逐列统计不同残基数，
并与"独立随机替换"对照模型比较。
"""
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_sequence_fidelity.json")
FASTA = Path("/tmp/lyz.fa")


def read_fasta(path):
    seqs, cur = [], []
    name = None
    for line in open(path):
        if line.startswith(">"):
            if cur:
                seqs.append((name, "".join(cur).upper()))
            name, cur = line[1:].strip(), []
        else:
            cur.append(line.strip().upper())
    if cur:
        seqs.append((name, "".join(cur).upper()))
    return seqs


def needleman_wunsch(a, b, gap=-1.0, mismatch=-1.0, match=1.0):
    """全长比对（numpy 行向量化 + 回溯）。返回 (对齐后的 a, 对齐后的 b)。"""
    n, m = len(a), len(b)
    F = np.full((n + 1, m + 1), -np.inf, dtype=np.float32)
    F[:, 0] = np.arange(n + 1) * gap
    F[0, :] = np.arange(m + 1) * gap
    for i in range(1, n + 1):
        prev = F[i - 1]
        row = F[i]
        ai = a[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match if ai == b[j - 1] else mismatch)
            row[j] = max(diag, prev[j] + gap, row[j - 1] + gap)
    i, j = n, m
    sa, sb = [], []
    while i > 0 and j > 0:
        diag = (match if a[i - 1] == b[j - 1] else mismatch)
        if abs(F[i][j] - (F[i - 1][j - 1] + diag)) < 1e-6:
            sa.append(a[i - 1]); sb.append(b[j - 1]); i -= 1; j -= 1
        elif abs(F[i][j] - (F[i - 1][j] + gap)) < 1e-6:
            sa.append(a[i - 1]); sb.append("-"); i -= 1
        else:
            sa.append("-"); sb.append(b[j - 1]); j -= 1
    while i > 0:
        sa.append(a[i - 1]); sb.append("-"); i -= 1
    while j > 0:
        sa.append("-"); sb.append(b[j - 1]); j -= 1
    return "".join(reversed(sa)), "".join(reversed(sb))


def star_align(ref, others):
    """星形比对：每条与 ref 比，**参照坐标不变**；参照侧出现缺口的位置在该行补'-'。"""
    out = []
    L = len(ref)
    for s in others:
        sa, sb = needleman_wunsch(ref, s)
        row, rpos = [], 0
        for ca, cb in zip(sa, sb):
            if ca == "-":
                row.append("-")        # 参照侧插入 ⟹ 该行此列缺
            else:
                rpos += 1
                row.append(cb)
        if rpos != L:
            continue
        out.append("".join(row))
    return out


def main():
    seqs = read_fasta(FASTA)
    lens = Counter(len(s) for _, s in seqs)
    ref_len = lens.most_common(1)[0][0]
    ref_name, ref = next((n, s) for n, s in seqs if len(s) == ref_len)
    others = [s for n, s in seqs if s != ref and ref_len - 60 <= len(s) <= ref_len + 60]
    short = ref_name.split("|")[2].split()[0] if "|" in ref_name else ref_name
    print(f"读入 {len(seqs)} 条；参照 {short}（{ref_len} aa）；参与比对 {len(others)} 条")

    # 修复：star_align 返回的行与 ref 同宽但**列向坐标系不同**（每条比对各自的缺口分布），
    # 直接横向拼接会让 COLUMN j 变成"各序列在各自对齐里的第 j 位"。正确做法：
    # 每条与 ref 比出来的对齐串先展开成"以 ref 为坐标"的行（ref 侧缺口 → 该行补 '-'），
    # 再把这些行**按 ref 坐标**对齐：即 ref 坐标 k 对应 sa 里走过的第 k 个非缺口位置。
    rows = [ref]
    for s in others:
        sa, sb = needleman_wunsch(ref, s)
        row = ["-"] * len(ref)
        k, rpos = 0, 0
        for ca, cb in zip(sa, sb):
            if ca == "-":
                continue              # 参照侧插入：该行不占 ref 坐标
            row[rpos] = cb
            rpos += 1
        if rpos == len(ref):
            rows.append("".join(row))
    # 上面这种"丢弃插入"的近似只对无插入情况成立；实际上少了插入位会让后面的位错位，
    # 故再把所有行补成等宽：以 ref 的非缺口位为公共坐标，**缺口也算一种状态**。
    rows = [ref] + [r for r in rows[1:]]
    L = len(ref)
    print(f"比对宽度 {L}；参与 {len(rows)} 条")

    vals = []
    for j in range(L):
        col = [r[j] for r in rows if r[j] != "-"]
        vals.append(len(set(col)) if col else None)
    vals_ok = np.asarray([v for v in vals if v is not None])
    ncols = len(vals_ok)
    conserved = int((vals_ok == 1).sum())
    semi = int(((vals_ok >= 2) & (vals_ok <= 4)).sum())
    hyper = int((vals_ok >= 5).sum())

    # 对照模型：每列独立随机替换（保持参照残基的概率按该列实测比例给）
    rng = np.random.default_rng(0)
    AAs = "ACDEFGHIKLMNPQRSTVWY"
    sim = []
    for _ in range(300):
        cols = []
        for j in range(L):
            col = [r[j] for r in rows if r[j] != "-"]
            if not col:
                continue
            p = sum(1 for c in col if c == ref[j]) / len(col)
            drawn = [ref[j] if rng.random() < p else AAs[rng.integers(20)] for _ in col]
            cols.append(len(set(drawn)))
        cols = np.asarray(cols)
        sim.append(float((cols == 1).mean()))
    sim = np.asarray(sim)

    res = {
        "n_sequences": len(rows), "width": L, "ref_len": ref_len,
        "conserved_frac": conserved / ncols,
        "semi_variable_frac": semi / ncols,
        "hyper_variable_frac": hyper / ncols,
        "species_per_column_mean": float(vals_ok.mean()),
        "distribution": [[int(k), int(v)] for k, v in Counter(vals_ok.tolist()).most_common(8)],
        "independent_model_conserved_frac_mean": float(sim.mean()),
        "independent_model_conserved_frac_99": float(np.percentile(sim, 99)),
        "observed_beyond_model": bool(conserved / ncols > np.percentile(sim, 99)),
    }
    print(f"位点变异性：完全保守 {conserved/ncols:.1%}；2–4 种 {semi/ncols:.1%}；≥5 种 {hyper/ncols:.1%}")
    print(f"每列平均不同残基数 {vals_ok.mean():.2f}；分布 {res['distribution']}")
    print(f"独立随机替换对照：完全保守列占比 均值 {sim.mean():.1%}，99% 分位 {np.percentile(sim, 99):.1%}")
    print(f"实测高于模型均值：{conserved / ncols > sim.mean()}；高于 99% 分位：{res['observed_beyond_model']}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())