"""
make_catalogue.py --- 把 results/observable_sweep.json 变成 OBSERVABLES.md

只做**排版**：按结构标签聚类，逐条列出可观测量与拟合参数。
物理命名不在这里生成 —— 见 OBSERVABLES.md §5（后验命名，手写并标注）。
"""
from __future__ import annotations

import json
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "results", "observable_sweep.json")
DST = os.path.join(HERE, "OBSERVABLES.md")


def fmt_params(lab, p):
    if not p:
        return "—"
    if lab.startswith("power-law"):
        return f"b = {p.get('b'):+.4f}"
    if lab.startswith("exponential"):
        return f"b = {p.get('b'):+.4f}"
    if lab.startswith("logarithmic"):
        return f"b = {p.get('b'):+.4f}"
    if lab.startswith("stretched"):
        return f"β = {p.get('beta')}, c = {p.get('c'):+.4f}"
    if lab.startswith("constant"):
        return f"a = {p.get('a'):.4g}"
    return ", ".join(f"{k}={v}" for k, v in list(p.items())[:3])


def main():
    d = json.load(open(SRC))
    recs = d["records"]
    by = defaultdict(list)
    for r in recs:
        by[r["structure"]["label"]].append(r)
    order = sorted(by, key=lambda k: (-len(by[k]), k))

    L = []
    A = L.append
    A("# 按结构分类的可观测量目录（任务 1 产出）\n")
    A(f"> 生成时间 {d['generated']}｜{len(recs)} 条可观测量｜机器 "
      f"[`observable_sweep.py`](observable_sweep.py) ＋ 分类器 [`structures.py`](structures.py)｜"
      f"数据 [`results/observable_sweep.json`](results/observable_sweep.json)\n")
    A("> **等级：全部是【数值证据】，不是【导出】。** 本页只贴**结构标签**；"
      "物理名字在 §5，**后验**给出且与结构分开。\n")

    A("\n---\n\n## §1 生成规则（为什么这是「中性」的）\n")
    A("此前我一直在**挑**可观测量去对 QFT（找极点、找谱密度），那不是中性做法。\n")
    A("这里用机械的叉乘枚举：\n")
    A("```\n可观测量 = (对象 O) × (算子/构造 K) × (读出 R)\n"
      "O ∈ {Γ 上的行走占据 n_k(v)；𝒢_T 闭链图；词层 {±1}^n（无 Γ 对照）；缺陷算符；记录层账本}\n"
      "K ∈ {计数、邻接 A、热核 e^{tA}、谱、商/投影、扰动 A+V|0><0|}\n"
      "R ∈ {点值、谱、分布、熵、参与比、极值、支撑、矩、能隙、关联、标度}\n```\n")
    A("分类器 `structures.classify_curve` 的规则：**所有候选模型都在对数空间算 R²**"
      "（衰减曲线跨若干数量级），取最高者，但若更简单的模型 R² 与之相差 < 2e-3 则取更简单的。"
      "复杂度次序 `constant < logarithmic < power-law < exponential < stretched`。"
      "**所有候选连同 R² 都写进 JSON**，便于复核与改判。\n")

    A("\n## §2 结构聚类\n")
    A("| 结构标签 | 条数 | 涉及扇区 |")
    A("|:--|--:|:--|")
    for k in order:
        secs = sorted({r["sector"] for r in by[k]})
        A(f"| `{k}` | {len(by[k])} | {', '.join(secs)} |")

    A("\n---\n\n## §3 逐结构目录\n")
    for k in order:
        A(f"\n### `{k}`（{len(by[k])} 条）\n")
        A("| 可观测量 | 扇区 | 拟合 | R² | 备注 |")
        A("|:--|:--|:--|--:|:--|")
        for r in sorted(by[k], key=lambda z: z["id"]):
            st = r["structure"]
            r2 = st.get("r2")
            r2s = f"{r2:.5f}" if isinstance(r2, float) else "—"
            note = r.get("note", "")
            m = r.get("meta", {})
            extra = []
            for key in ("T", "D", "L", "n", "N"):
                if key in m:
                    extra.append(f"{key}={m[key]}")
            if extra:
                note = ("；".join(extra) + ("；" + note if note else ""))
            if len(note) > 90:
                note = note[:88] + "…"
            A(f"| `{r['id']}` | {r['sector']} | {fmt_params(k, st.get('params'))} | {r2s} | {note} |")

    # 频谱形状记录单列
    A("\n---\n\n## §4 谱形状（分布类可观测量，逐条）\n")
    A("| 对象 | 标签 | 带 | 带边指数(lo/hi) | 其他 |")
    A("|:--|:--|:--|:--|:--|")
    for r in recs:
        if r["readout"] != "谱密度形状":
            continue
        m = r["meta"]
        sc = m.get("spectrum_class") or {}
        band = sc.get("band")
        band_s = f"[{band[0]:.3f}, {band[1]:.3f}]" if band else "—"
        ee = sc.get("edge_exponent") or {}
        ees = f"{ee.get('lo')} / {ee.get('hi')}" if ee else "—"
        oth = []
        if m.get("bipartite") is not None:
            oth.append(f"二部={m['bipartite']}")
        if m.get("spacing_ratio_r_distinct") is not None:
            oth.append(f"r_distinct={m['spacing_ratio_r_distinct']:.4f}")
        sh = m.get("shape_stats")
        if sh:
            oth.append(f"简并: {sh['distinct']}/{sh['N']} 不同值, 最大重数 {sh['max_multiplicity']}, "
                       f"1%带宽团簇 {sh.get('clusters_at_1pct_bandwidth')}")
        if m.get("lambda_max") is not None:
            oth.append(f"λmax={m['lambda_max']:.3f}")
        A(f"| `{r['id']}` | {sc.get('label')} | {band_s} | {ees} | {'；'.join(oth)} |")

    # ---------------- §5 后验命名 ----------------
    A("\n---\n\n## §5 让结构自己报名字（**后验**，与 §1–§4 分开）\n")
    A("§1–§4 全程只用结构标签，没有任何物理名字参与分类。下表是**事后**把结构对到已知物理概念，"
      "供人读；**它不改变上面的分类，也不构成断言**。\n")
    A("| 结构（§2 的标签） | 出现的对象 | 已知物理里的同名结构 | 等级 |")
    A("|:--|:--|:--|:--|")
    rows5 = [
        ("`power-law` 指数 $-D/2$", "Γ 上行走的回返 $n_k(0)$、参与比 IPR、峰位概率",
         "连续时间随机游走的回返衰减 ⟹ **谱维数** $d_s=D$（扩散）", "【数值证据】＋闭式核验"),
        ("`power-law` 指数 $+1$", "二阶矩 $M_2(k)=\\sum r^2p$",
         "**扩散律** $\\langle x^2\\rangle\\propto t$", "【数值证据】"),
        ("`logarithmic`", "占据熵 $H(k)=-\\sum p\\ln p$",
         "熵的对数增长 $\\tfrac D2\\ln k+\\text{const}$（维数读出）", "【数值证据】"),
        ("`constant`", "超额峰度 $M_4/M_2^2$",
         "**高斯普适性**（中心极限的不动点）", "【数值证据】"),
        ("`exponential`", "闭链图热核 $Z(t)=\\langle e^{t\\lambda}\\rangle$",
         "**Perron/指数增长**，增长率 $=\\lambda_{\\max}$", "【数值证据】"),
        ("`gapped` ＋ $\\langle r\\rangle_{\\rm distinct}\\approx0.39$", "闭链图 $\\mathcal G_T$ 的谱",
         "**Poisson 型能级统计（无能级排斥）** ＋ 团簇/简并谱（非随机矩阵）", "【数值证据】"),
        ("`continuum/edge-divergent`", "$D{=}1$ 邻接谱",
         "**van Hove 奇点**（弧正弦律 $1/\\sqrt{4-\\lambda^2}$）", "【数值证据】"),
        ("`continuum/center-divergent(van-Hove)`", "$D{=}2,3,4$ 与 $\\mathcal G_T$ 的谱",
         "**van Hove 奇点**（态密度对数发散）", "【数值证据】"),
        ("`continuum/edge-vanishing`", "$D{=}3,4$ 邻接谱",
         "带边**幂律消失**（Weyl 律的带边指数）", "【数值证据】"),
        ("`delta` ＋ `parity`", "词层关联 $C(r)=\\langle w_0w_r\\rangle=\\delta_{r0}$、回返奇偶零",
         "**白噪声（无关联）** ／ 二部图奇偶选择定则", "【数值证据】＋精确"),
        ("`constant`（饱和）", "记录层覆盖率 $\\to0.7657$",
         "**饱和/稳态**（非幂律、非指数）", "【数值证据】"),
        ("带边 DOS 指数 $=D/2-1$", "缺陷算符 $A+V|0\\rangle\\langle0|$ 的带边",
         "**晶格 Green 函数在带边的行为**：指数 $<0$ ⟹ 任意弱耦合成键；$>0$ ⟹ **有阈值**",
         "【数值证据】＋理论一致"),
    ]
    for r in rows5:
        A("| " + " | ".join(r) + " |")
    A("\n**注意**：§5 的最后一行的「成键阈值」是[`L2_REGIONS.md`](L2_REGIONS.md) §3.1 的输入；"
      "那里用它把 $D\\le2$ 与 $D\\ge3$ 定性分开。\n")

    with open(DST, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"写出 {DST}（{len(L)} 行，{len(recs)} 条记录）")


if __name__ == "__main__":
    main()
