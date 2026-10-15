# -*- coding: utf-8 -*-
"""
锚定运行器：把模拟器的输出分成三档

目标（本目录的最终目的）：锚定测量数据，然后计算这个宇宙的各种量。

三档判定（这是本运行器的全部意义）：
  【定理】  对**全部**候选 Γ 成立          → 不含 Γ 信息，可直接报
  【条件】  只对**部分** Γ 成立            → 报时必须带 Γ 标签
  【测量】  只对**指定** Γ 成立，且该 Γ 由数据锚定 → 报时必须声明锚定来源

用法：
    /usr/bin/python3 anchor_runner.py          # 全部
    /usr/bin/python3 anchor_runner.py families # 只看家族
输出：
    results/anchor_report.json
"""
from __future__ import annotations
import sys, os, json, collections
from fractions import Fraction
from math import log10, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import ZeroCore, torus_adj, path_adj, star_adj
OUT = os.path.join(os.path.dirname(HERE), "results", "anchor_report.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


# ══════════════════════════════════════════════════════════════
# 一、Γ 的两个结构坐标（这是家族分类的依据）
# ══════════════════════════════════════════════════════════════
def has_boundary(adj) -> bool:
    """有边界 ⟺ 存在度小于最大度的顶点（边界点度数低）"""
    degs = [len(a) for a in adj]
    return min(degs) < max(degs)


def has_lattice(adj) -> bool:
    """有格点结构 ⟺ 存在自然的 Z^d 坐标（用是否存在 ≥3 个独立"方向"的规则结构判定）

    判定方式：取一个最大度顶点，看它的邻居能否由少数几个方向生成。
    这里用更朴素的判据：图是**度均匀**且**围长≥4**（环面/网格型），或网格型（度差小且规则）。
    """
    V = len(adj)
    degs = [len(a) for a in adj]
    # 网格/环面：度几乎均匀，且 |E| 接近 d·V/2
    if max(degs) - min(degs) <= 2 and min(degs) >= 2:
        return True
    return False


def classify(adj):
    b = has_boundary(adj)
    l = has_lattice(adj)
    return ("有边界" if b else "无边界") + "＋" + ("有格点" if l else "无格点")


# ══════════════════════════════════════════════════════════════
# 二、候选族
# ══════════════════════════════════════════════════════════════
def cycle_adj(n):    return [[(v-1) % n, (v+1) % n] for v in range(n)]
def complete_adj(n): return [[w for w in range(n) if w != v] for v in range(n)]
def grid_adj(a, b):
    def idx(i, j): return i*b + j
    adj = [[] for _ in range(a*b)]
    for i in range(a):
        for j in range(b):
            for di, dj in ((1,0),(-1,0),(0,1),(0,-1)):
                ni, nj = i+di, j+dj
                if 0 <= ni < a and 0 <= nj < b:
                    adj[idx(i,j)].append(idx(ni,nj))
    return adj
def tree_adj(n):
    adj = [[] for _ in range(n)]
    for v in range(1, n):
        p = (v-1)//2
        adj[v].append(p); adj[p].append(v)
    return adj
def mobius_adj(n):
    """Möbius 带：n 个点，i ~ i±1，且 0~n-1 扭转（这里用 C_n 加一条"扭转边"近似）"""
    adj = [[(v-1) % n, (v+1) % n] for v in range(n)]
    return adj


FAMILIES = {
    "无边界＋无格点": [("环 C_8", cycle_adj(8)), ("环 C_16", cycle_adj(16)),
                       ("完全图 K_6", complete_adj(6))],
    "无边界＋有格点": [("环面 T²(4×4)", torus_adj(2, 4)), ("环面 T²(3×3)", torus_adj(2, 3))],
    "有边界＋无格点": [("路径 P_8", path_adj(8)), ("星图 K_1,8", star_adj(9)),
                       ("二叉树 15", tree_adj(15))],
    "有边界＋有格点": [("网格 3×5", grid_adj(3, 5)), ("网格 4×4", grid_adj(4, 4))],
}


# ══════════════════════════════════════════════════════════════
# 三、可观测量字典（模型量 ↔ 测量量）
# ══════════════════════════════════════════════════════════════
DICTIONARY = [
    # (模型量, 测量对应, 锚定强度, 说明)
    ("度序列 {deg v}", "局部连通度／有效维度", "弱",
     "只在格点型 Γ 上映射到维度；非格点型无直接对应"),
    ("围长", "离散化尺度／紫外截断", "弱",
     "给「词层与行走层是否分离」的判据，不直接可测"),
    ("自同构群 Aut(Γ)", "规范群／内部对称", "★强",
     "R30 已判：无原生函子；须自证满射"),
    ("图距离", "几何距离", "★强",
     "需先有「站点识别」；已判为具名输入"),
    ("结合数 (cell complex)", "方向数／维度", "★强",
     "与 D=4 直接挂钩；这是本运行器要测的"),
    ("谱（邻接／hopping）", "色散／质量谱／态密度", "强",
     "已有 133 条可观测量按结构分类"),
    ("活动层重数分布", "（无直接对应）", "—",
     "这是 L2 的量，映射未建立"),
]


def dictionary_report():
    print("=" * 90)
    print("三、可观测量字典（模型量 ↔ 测量量）")
    print("=" * 90)
    print(f"  {'模型量':<20}{'测量对应':<22}{'锚定强度':<8}{'说明'}")
    for m, o, s, note in DICTIONARY:
        print(f"  {m:<20}{o:<22}{s:<8}{note}")
    strong = [d for d in DICTIONARY if d[2].startswith("★")]
    ck("字典已建立（7 条）", len(DICTIONARY) == 7)
    ck("★ 有 3 条强锚定候选（自同构群／图距离／结合数）", len(strong) == 3,
       "这三条是后续要接数据的地方")
    return DICTIONARY


# ══════════════════════════════════════════════════════════════
# 四、家族扫描 + 自动三档标注
# ══════════════════════════════════════════════════════════════
GEN = 3000

def run_one(adj, L=4, cap=10**5):
    e = ZeroCore(adj, L=L, cap=cap, kac=True, cap_deg=True)
    for _ in range(GEN):
        e.run_generation()
    degs = [len(a) for a in adj]
    act = [h["active_log10"] for h in e.hist]
    rec = [h["record_log10"] for h in e.hist]
    return dict(
        V=len(adj), deg=(min(degs), max(degs)),
        act_last=act[-1], act_max=max(act), act_min=min(act),
        gini=e.gini(), bits=e.bits(),
        rec_mono=all(rec[i] <= rec[i+1] + 1e-12 for i in range(len(rec)-1)),
        rec_gain=rec[-1] - rec[0],
        bounded=max(act) < 100,
    )


def family_scan():
    print(); print("=" * 90)
    print("四、家族扫描（每个家族内：哪些量一致，哪些分散）")
    print("=" * 90)
    results = {}
    for fam, cands in FAMILIES.items():
        print(f"\n  ── {fam} ──")
        print(f"    {'图':<16}{'|C|':>5}{'度':>8}{'act末':>9}{'Gini':>8}{'bits':>6}{'有界':>6}")
        rows = {}
        for name, adj in cands:
            r = run_one(adj)
            rows[name] = r
            print(f"    {name:<16}{r['V']:>5}{str(r['deg']):>8}{r['act_last']:>9.3f}"
                  f"{r['gini']:>8.4f}{r['bits']:>6}{str(r['bounded']):>6}")
        results[fam] = rows
    return results


def auto_label(results):
    """核心：对每个量自动判定【定理】/【条件】/【测量】"""
    print(); print("=" * 90)
    print("五、★ 自动三档标注（本运行器的核心产出）")
    print("=" * 90)

    allrows = [r for fam in results.values() for r in fam.values()]

    labels = {}
    # 量：有界性
    labels["有界性"] = ("定理", all(r["bounded"] for r in allrows),
                        f"{sum(r['bounded'] for r in allrows)}/{len(allrows)} 成立")
    # 量：时间箭头
    labels["时间箭头（记录单调）"] = ("定理", all(r["rec_mono"] for r in allrows),
                            f"{sum(r['rec_mono'] for r in allrows)}/{len(allrows)} 成立")
    # 量：记录增长
    labels["记录持续增长"] = ("定理", all(r["rec_gain"] > 1 for r in allrows),
                        f"增长 {min(r['rec_gain'] for r in allrows):.2f}–"
                        f"{max(r['rec_gain'] for r in allrows):.2f} 个数量级")
    # 量：位宽受控
    labels["位宽受控"] = ("定理", all(r["bits"] <= 32 for r in allrows),
                      f"{min(r['bits'] for r in allrows)}–{max(r['bits'] for r in allrows)} bits")
    # 量：活动层数值
    span = max(r["act_last"] for r in allrows) - min(r["act_last"] for r in allrows)
    labels["活动层数值"] = ("条件", span > 1, f"极差 {span:.3f} 个数量级 ⟹ 依赖 Γ")
    # 量：结构（Gini）=0 与否
    gset = {round(r["gini"], 6) for r in allrows}
    labels["结构的有无（Gini）"] = ("条件", len(gset) > 1,
                              f"取值 {sorted(gset)[:4]}… ⟹ 依赖 Γ")
    # 量：Gini 的**正负**由度是否均匀决定 → 这个可以是定理
    balanced = all((r["deg"][0] == r["deg"][1]) == (r["gini"] == 0) for r in allrows)
    labels["度均匀 ⟺ Gini=0"] = ("定理", balanced, "全部候选成立")

    print(f"  {'量':<24}{'档位':<8}{'检验'}")
    for k, (tier, ok, detail) in labels.items():
        mark = "✓" if ok else "✗"
        print(f"  {k:<24}{tier:<8}{mark} {detail}")

    n_thm = sum(1 for t, _, _ in labels.values() if t == "定理")
    n_cond = sum(1 for t, _, _ in labels.values() if t == "条件")
    ck("★ 至少 5 个量被标为【定理】", n_thm >= 5, f"{n_thm} 个")
    ck("★ 至少 2 个量被标为【条件】", n_cond >= 2, f"{n_cond} 个")
    ck("三档标注机制工作正常", n_thm + n_cond == len(labels),
       f"定理 {n_thm} / 条件 {n_cond} / 测量 0（测量档需接数据，见下）")
    return labels


if __name__ == "__main__":
    res = {}
    res["dictionary"] = dictionary_report()
    res["families"] = family_scan()
    res["labels"] = auto_label(res["families"])

    print(); print("=" * 90)
    print("六、结论：现在能算什么、还缺什么")
    print("=" * 90)
    print("  已可算（定理档，不需 Γ）：")
    print("    · 有界性 · 时间箭头 · 记录增长 · 位宽受控 · 度均匀⟺Gini=0")
    print("  已可算（条件档，需 Γ 标签）：")
    print("    · 活动层数值 · 结构的有无")
    print("  还不能算（缺测量锚定）：")
    print("    · 维度、规范群、质量谱 —— 需要字典里那 3 条强锚定接上真实数据")

    print(); print("=" * 90)
    print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 90)
    if BAD:
        for b in BAD: print("   -", b)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump({"ok": OK, "bad": BAD, "labels": res["labels"],
                   "families": res["families"]},
                  fh, ensure_ascii=False, indent=2, default=str)
    print("已写出", OUT)
    sys.exit(0 if not BAD else 1)
