#!/usr/bin/env python3
"""P2：简并度是否等于写入等级数 N = 1 + floor(λ/κ)？

数据：标准遗传密码表（NCBI 表 1，硬编码，可复核）；
做法：数每个氨基酸的密码子数（= 纤维大小 = 简并度），检验
  (i) 简并度分布是否集中在 {1,2,3,4,6}；
  (ii) 每个简并度 d 是否对应整数级别的 N−1 = d−1（即 d−1 是整数）；
  (iii) 简并度的加权是否与三联体 4³ 相容：Σ_d (#氨基酸 with d) × d = 64。
"""
import json
from collections import Counter
from pathlib import Path

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb095_codon_degeneracy.json")

# NCBI 遗传密码表 1（标准）：64 密码子 → 氨基酸/终止
BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
TABLE = """
TTT F TTC F TTA L TTG L
CTT L CTC L CTA L CTG L
ATT I ATC I ATA I ATG M
GTT V GTC V GTA V GTG V
TCT S TCC S TCA S TCG S
CCT P CCC P CCA P CCG P
ACT T ACC T ACA T ACG T
GCT A GCC A GCA A GCG A
TAT Y TAC Y TAA * TAG *
CAT H CAC H CAA Q CAG Q
AAT N AAC N AAA K AAG K
GAT D GAC D GAA E GAG E
TGT C TGC C TGA * TGG W
CGT R CGC R CGA R CGG R
AGT S AGC S AGA R AGG R
GGT G GGC G GGA G GGG G
"""
pairs = TABLE.split()
code = {pairs[i]: pairs[i + 1] for i in range(0, len(pairs), 2)}

def main():
    assert len(code) == 64, len(code)
    assert set(code) == set(CODONS), "密码子集合不完全"
    degen = Counter(code[c] for c in CODONS)
    aa = {k: v for k, v in degen.items() if k != "*"}
    rows = {
        "n_codons": len(CODONS),
        "n_aa": len(aa),
        "n_stop": degen.get("*", 0),
        "degeneracy_by_aa": dict(sorted(aa.items(), key=lambda kv: (-kv[1], kv[0]))),
        "degeneracy_distribution": dict(sorted(Counter(aa.values()).items())),
        "sum_codons_check": sum(aa.values()) + degen.get("*", 0),
    }
    print(f"密码子 {rows['n_codons']}；氨基酸 {rows['n_aa']}；终止 {rows['n_stop']}")
    print(f"校验：Σ 简并度 = {rows['sum_codons_check']}（应为 64）")
    print("\n简并度分布：")
    for d, n in sorted(rows["degeneracy_distribution"].items()):
        names = [k for k, v in aa.items() if v == d]
        print(f"  d={d}: {n} 个氨基酸  {' '.join(names)}")

    # (i) 是否集中在 {1,2,3,4,6}
    core = {1, 2, 3, 4, 6}
    in_core = sum(n for d, n in rows["degeneracy_distribution"].items() if d in core)
    rows["in_core_set_12346"] = in_core
    rows["all_in_core"] = in_core == rows["n_aa"]
    print(f"\n(i) 简并度全部落在 {{1,2,3,4,6}}：{rows['all_in_core']}")

    # (ii) N = 1 + floor(λ/κ)，故简并度 d 对应 N 级别 d，Δ = d−1 应是非负整数
    rows["d_minus_1_all_integer"] = all(isinstance(d - 1, int) and d - 1 >= 0 for d in aa.values())
    print(f"(ii) 每个 d−1 ≥ 0 且为整数：{rows['d_minus_1_all_integer']}（对数 d 恒为整数，此项必真，仅作记录）")

    # (iii) 4³ 约束
    ok = sum(aa.values()) + degen.get("*", 0) == 64
    rows["triplet_constraint_ok"] = ok
    print(f"(iii) 简并度×氨基酸数 + 终止 = 64：{ok}")

    # 真正的检验：N = 1 + floor(λ/κ) 要求 d ∈ {1+k : k ∈ Z}，任何整数 d 都满足，
    # 故 (ii) 不构成检验；可检验的是"简并度的层级是否只有少数几档"（少量端口数）。
    levels = sorted(rows["degeneracy_distribution"])
    rows["distinct_levels"] = levels
    rows["n_levels"] = len(levels)
    print(f"\n*** 关键读数：简并度只出现 {len(levels)} 档：{levels}")
    print("   N=1+⌊λ/κ⌋ 的模型只有 k 取少数整数时才有意义；若简并度是 1,2,3,4,6 五档，对应 k=0,1,2,3,5 ——")
    print("   即端口带宽只有 5 个可用级别，而 5 是 Zb085 '三笔输入' 之外的一个新数；这一条是本文档需要诚实标注的。")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(rows, f, indent=2, ensure_ascii=False)
    print("\n落盘", OUT)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
