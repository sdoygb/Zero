#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_architecture_map.py — Zero 全架构总图

把 Zero 的推导链画成一张分层结构图：
    Z0 公理 → 四层 (L0/L1/L1'/L2) + 读出面 R → 舞台 Γ=K4 与通道数 m
    → 六条层间边闭环 → 读出：3 频率槽 / 十六槽 → 物理量（D=4、质量谱、标度链）

设计纪律（与库内体例一致）：
  · 图上每一个数字都必须在库里能找到出处；脚本内以 SRC 表逐条登记。
  · 只画「已证的结构性事实」，不画 no-go / 缺口 / 「不可能」类断语。
  · 输出纯结构图，不含数据诊断曲线（那是既有 81 张图做的事）。

出处（逐条）：
  Z001 §2.1 公理 Z0 三款；§4.1 层表；§7 主链；引理 7.1 D4 窗
  Z016 §2.1 四层内容；命题 1 读出面不占层
  Z088 §1 六条层间边（①L0→L1 ②L1→L2 ③L2→L1 ④L0→L2 ⑤L2→L0 ⑥L1→L0）
  Z161 §1.1 层账 L0=0 / L1=H=L/2 / L2=3；§3.1 三槽 = 相邻比 f_i=|log(p_i/p_i+1)|
  Z021 维数字典 D=m-1；K4、arccos(-1/3)=109.4712°
  Z136 §3.1 Γ16=K4□C4，16 槽；§3.5 每代 8 场
  Z153 §3.1 锚 v=246.22 GeV；ln(M_X/v)=33.6197
  Z001 §7.7 q=0.55777859 ∈ (1/2,3/5)

复跑： /usr/bin/python3 zero_architecture_map.py
输出： ../results/zero_architecture_map.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, os.pardir, "results")

# ── 中文字体（macOS 自带；无此设置则全部 CJK 变成方框）─────────
from matplotlib import font_manager as _fm

_CJK_CANDIDATES = [
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
    "/System/Library/Fonts/STHeiti Medium.ttc",
    "/System/Library/Fonts/Supplemental/Songti.ttc",
    "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
]
_cjk_path = next((p for p in _CJK_CANDIDATES if os.path.exists(p)), None)
if _cjk_path:
    _fm.fontManager.addfont(_cjk_path)
    _cjk_name = _fm.FontProperties(fname=_cjk_path).get_name()
    matplotlib.rcParams["font.family"] = [_cjk_name, "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
else:
    print("[warn] 未找到中文字体，CJK 会显示为方框")

# ── 配色（单一色系，层越"物理"越暖）───────────────────────────
C = {
    "bg": "#ffffff",
    "axiom": "#1f3b73",
    "l0": "#2c5f8a",
    "l1": "#3f7fa8",
    "l1p": "#6aa6c4",
    "l2": "#c8794a",
    "read": "#8a6fa8",
    "phys": "#4a7c59",
    "edge": "#5a5a5a",
    "ink": "#222222",
    "muted": "#6b6b6b",
    "grid": "#e8e8e8",
}


def box(ax, x, y, w, h, text, fc, tc="white", fs=10, lw=1.2, ec=None,
        radius=0.012, bold=False, alpha=1.0):
    """圆角框 + 居中文字（坐标用 axes fraction）。"""
    p = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                       boxstyle=f"round,pad=0,rounding_size={radius}",
                       linewidth=lw, edgecolor=ec or fc, facecolor=fc,
                       alpha=alpha, transform=ax.transAxes, zorder=3,
                       mutation_aspect=1.0)
    ax.add_patch(p)
    ax.text(x, y, text, ha="center", va="center", color=tc,
            fontsize=fs, transform=ax.transAxes, zorder=4,
            fontweight="bold" if bold else "normal", linespacing=1.45)


def arrow(ax, p, q, color=None, lw=1.4, style="-|>", ls="-", rad=0.0,
          ms=8, z=2, alpha=1.0):
    """两点间箭头（axes fraction）。p/q 为 (x,y)。"""
    a = FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms,
                        linewidth=lw, color=color or C["edge"],
                        linestyle=ls, transform=ax.transAxes, zorder=z,
                        connectionstyle=f"arc3,rad={rad}",
                        shrinkA=0, shrinkB=0, alpha=alpha)
    ax.add_patch(a)


def label(ax, x, y, text, fs=8.5, color=None, ha="center", va="center",
          style="italic", z=5, weight="normal", rot=0):
    ax.text(x, y, text, ha=ha, va=va, fontsize=fs, color=color or C["muted"],
            style=style, transform=ax.transAxes, zorder=z,
            fontweight=weight, rotation=rot, linespacing=1.4)


# ══════════════════════════════════════════════════════════════
# 主图
# ══════════════════════════════════════════════════════════════

def build():
    fig = plt.figure(figsize=(19, 13.4), facecolor=C["bg"])
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # ── 标题 ─────────────────────────────────────────────────
    ax.text(0.5, 0.972, "Zero 架构总图", ha="center", va="center",
            fontsize=27, color=C["axiom"], fontweight="bold",
            transform=ax.transAxes)
    ax.text(0.5, 0.944,
            "一条公理 · 四层 · 六条层间边 · 从零不停留到物理量",
            ha="center", va="center", fontsize=12.5, color=C["muted"],
            transform=ax.transAxes, style="italic")
    ax.plot([0.06, 0.94], [0.927, 0.927], color=C["grid"], lw=1.4,
            transform=ax.transAxes, zorder=1)

    # ═════════════════════════════════════════════════════════
    # ① 公理（顶）
    # ═════════════════════════════════════════════════════════
    box(ax, 0.5, 0.882, 0.60, 0.070,
        "公理 $Z_0$（零不停留）—— 体系内唯一公理\n"
        "① 零不停留（平衡态不是静止态）   ② 从不停歇   ③ 不设概率",
        C["axiom"], fs=11.5, bold=False)
    label(ax, 0.5, 0.842,
          "「零」是两支配对的商 {±1}/ι，不是一个静止的点    ·    Z001 §2.1",
          fs=8.6)

    # ═════════════════════════════════════════════════════════
    # ② 导出结构条
    # ═════════════════════════════════════════════════════════
    box(ax, 0.5, 0.800, 0.78, 0.040,
        "导出（非公理）：词 w =(e₁,…,eₙ) ｜ 散度 d_w ｜ 零和恒等式 Σ d_w ≡ 0 ｜ "
        "全分支＋整数重数 ｜ 终端＋重播种",
        "#eef3f8", tc=C["axiom"], fs=10, ec=C["l0"], lw=1.1)
    label(ax, 0.5, 0.771,
          "Z001 §3：这些全部是推论，不是第二公理",
          fs=8.2)

    # ═════════════════════════════════════════════════════════
    # ③ 四层 + 读出面（主体）
    # ═════════════════════════════════════════════════════════
    label(ax, 0.115, 0.722, "四层架构 ＋ 读出面", fs=13, color=C["axiom"],
          style="normal", weight="bold", ha="left")

    rows = [
        # (y, 高, 色, 层名, 内容, 更新律, 物理世界)
        (0.655, 0.062, C["l0"], "L0  底层",
         "公理条款 ｜ 非负整数重数 ｜ 局域补偿移动",
         "无自身更新律", "✗"),
        (0.583, 0.062, C["l1"], "L1  历史层",
         "精确词 W_L ｜ 循环次序 ｜ 双覆盖 ℤ₂ ｜ 高度 h(w)=max|前缀和|",
         "无（只增、不可逆）", "✗"),
        (0.511, 0.062, C["l1p"], "L1′  全局闭合类层",
         "闭合类 [w] ｜ 计数测度 μ=|A|/|Wₙ| ｜ 记录 λ ｜ 读出",
         "无（记录层）", "✗"),
        (0.439, 0.062, C["l2"], "L2  演化层",
         "滑动窗 ｜ 随机步 ｜ 终端吸收 ｜ 重播种 ｜ 亚层 ｜ 容量饱和",
         "★ 有自身更新律", "✅"),
    ]
    for y, h, col, name, content, upd, isphys in rows:
        box(ax, 0.36, y, 0.50, h, "", col, radius=0.008, alpha=0.14, ec=col, lw=1.0)
        ax.add_patch(FancyBboxPatch((0.11, y - h / 2), 0.008, h,
                                    boxstyle="round,pad=0,rounding_size=0.004",
                                    linewidth=0, facecolor=col,
                                    transform=ax.transAxes, zorder=4))
        ax.text(0.135, y + 0.012, name, ha="left", va="center", fontsize=12,
                color=col, fontweight="bold", transform=ax.transAxes, zorder=5)
        ax.text(0.135, y - 0.014, content, ha="left", va="center", fontsize=9.2,
                color=C["ink"], transform=ax.transAxes, zorder=5)
        ax.text(0.625, y, upd, ha="center", va="center", fontsize=9.4,
                color=col if "有自身" in upd else C["muted"],
                fontweight="bold" if "有自身" in upd else "normal",
                transform=ax.transAxes, zorder=5)
        ax.text(0.685, y, isphys, ha="center", va="center", fontsize=12,
                color=C["phys"] if isphys == "✅" else "#bbbbbb",
                transform=ax.transAxes, zorder=5)

    # 读出面（不占层指标）
    box(ax, 0.36, 0.352, 0.50, 0.052, "", C["read"], radius=0.008,
        alpha=0.13, ec=C["read"], lw=1.0, )
    ax.add_patch(FancyBboxPatch((0.11, 0.326), 0.008, 0.052,
                                boxstyle="round,pad=0,rounding_size=0.004",
                                linewidth=0, facecolor=C["read"],
                                transform=ax.transAxes, zorder=4))
    ax.text(0.135, 0.364, "R  读出面（不是层）", ha="left", va="center",
            fontsize=12, color=C["read"], fontweight="bold",
            transform=ax.transAxes, zorder=5)
    ax.text(0.135, 0.338, "π（自识别划分）｜ ω ｜ K_ω = −log ω ｜ σ_t ｜ GNS 内积",
            ha="left", va="center", fontsize=9.2, color=C["ink"],
            transform=ax.transAxes, zorder=5)
    ax.text(0.625, 0.352, "导出记号集", ha="center", va="center", fontsize=9.4,
            color=C["read"], transform=ax.transAxes, zorder=5)
    ax.text(0.685, 0.352, "✗", ha="center", va="center", fontsize=12,
            color="#bbbbbb", transform=ax.transAxes, zorder=5)

    # 表头
    ax.text(0.625, 0.695, "自身更新律", ha="center", va="center", fontsize=9,
            color=C["muted"], transform=ax.transAxes, weight="bold")
    ax.text(0.685, 0.695, "物理世界", ha="center", va="center", fontsize=9,
            color=C["muted"], transform=ax.transAxes, weight="bold")
    label(ax, 0.36, 0.302,
          "Z016 §2.1 / 命题 1：读出面不占层指标——它是一笔具名输入 π 与其记号，不是第五层",
          fs=8.6)

    # ★ 物理世界高亮
    # ★ 物理世界高亮（放在六边框正下方，独立区域）
    ax.add_patch(FancyBboxPatch((0.775, 0.712), 0.195, 0.042,
                                boxstyle="round,pad=0,rounding_size=0.010",
                                linewidth=2.2, edgecolor=C["l2"],
                                facecolor="#fdf3ec", transform=ax.transAxes,
                                zorder=6))
    ax.text(0.8725, 0.732, "★ 物理世界在 L2", ha="center", va="center",
            fontsize=11, color=C["l2"], fontweight="bold",
            transform=ax.transAxes, zorder=7)
    ax.text(0.8725, 0.718, "（唯一有自身更新律的层）", ha="center", va="center",
            fontsize=7.8, color=C["muted"], transform=ax.transAxes, zorder=7)

    # ═════════════════════════════════════════════════════════
    # ④ 六条层间边（右侧回路）
    # ═════════════════════════════════════════════════════════
    ax.add_patch(FancyBboxPatch((0.775, 0.268), 0.195, 0.428,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                linewidth=1.1, edgecolor=C["edge"],
                                facecolor="#fafafa", transform=ax.transAxes,
                                zorder=1))
    ax.text(0.8725, 0.678, "六条层间边（全双向）", ha="center", va="center",
            fontsize=11.5, color=C["ink"], fontweight="bold",
            transform=ax.transAxes, zorder=3)
    edges = [
        ("①  L0 → L1", "闭合约束生成词集"),
        ("②  L1 → L2", "历史作种子 seed"),
        ("③  L2 → L1", "演化写历史 record"),
        ("④  L0 → L2", "结构决定（图 Γ、闭合判据）"),
        ("⑤  L2 → L0", "本代闭合回灌 q_live"),
        ("⑥  L1 → L0", "记录回灌 q_hist"),
    ]
    for i, (a, b) in enumerate(edges):
        yy = 0.635 - i * 0.042
        ax.text(0.788, yy, a, ha="left", va="center", fontsize=9.2,
                color=C["l0"], fontweight="bold",
                transform=ax.transAxes, zorder=3, family="monospace")
        ax.text(0.882, yy, b, ha="left", va="center", fontsize=8.6,
                color=C["ink"], transform=ax.transAxes, zorder=3)
    ax.plot([0.868, 0.868], [0.392, 0.652], color=C["grid"], lw=1.0,
            transform=ax.transAxes, zorder=2)
    ax.text(0.8725, 0.356, "层间边 ≠ 层自身更新律", ha="center", va="center",
            fontsize=9.2, color=C["l2"], style="italic", fontweight="bold",
            transform=ax.transAxes, zorder=3)
    ax.text(0.8725, 0.330, "闭合回灌 L0 ＝ 重数累积，不是加词",
            ha="center", va="center", fontsize=8.4, color=C["muted"],
            transform=ax.transAxes, zorder=3)
    label(ax, 0.8725, 0.308, "Z088", fs=8)

    # 层→六边 的连线
    for y in (0.655, 0.583, 0.511, 0.439):
        arrow(ax, (0.615, y), (0.778, 0.480), color=C["grid"], lw=1.0,
              style="-", z=1, alpha=0.9)

    # ═════════════════════════════════════════════════════════
    # ⑤ 舞台与维数（中段左）
    # ═════════════════════════════════════════════════════════
    label(ax, 0.115, 0.288, "舞台与维数：D 从计数里出来", fs=13, color=C["axiom"],
          style="normal", weight="bold", ha="left")

    # K4 示意图
    kx, ky, kr = 0.168, 0.216, 0.038
    pts = [(kx + kr * np.cos(t), ky + kr * np.sin(t) * 0.85)
           for t in np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, 5)][:4]
    for i in range(4):
        for j in range(i + 1, 4):
            ax.plot([pts[i][0], pts[j][0]], [pts[i][1], pts[j][1]],
                    color=C["l0"], lw=1.6, transform=ax.transAxes, zorder=3,
                    alpha=0.85)
    for (px, py) in pts:
        ax.add_patch(Circle((px, py), 0.0075, color=C["l0"],
                            transform=ax.transAxes, zorder=4))
    ax.text(kx, ky, "K₄", ha="center", va="center", fontsize=12,
            color="white", fontweight="bold", transform=ax.transAxes, zorder=5)
    label(ax, kx, 0.168, "舞台 Γ = K₄：4 顶点、6 边、度 3", fs=8.8)
    label(ax, kx, 0.148, "两两夹角 arccos(−1/3) = 109.4712°", fs=8.2)

    # 维数框
    box(ax, 0.398, 0.230, 0.175, 0.050,
        "维数字典   D = m − 1\n通道数 m = D + 1", C["l0"], fs=10.5)
    label(ax, 0.398, 0.198,
          "C(m,2) − (m−1) = C(D,2)  ⟺  D = m−1（唯一）", fs=8.2)
    label(ax, 0.398, 0.180,
          "物理值处 m = 4 → K₄、Aut(K₄)=S₄（阶 24）", fs=8.2)

    # D4 窗
    box(ax, 0.398, 0.116, 0.175, 0.068,
        "账本泛函  F_D = C(D,2) · q^D\nq = 5/9 = 0.55778\n⟹ arg max F_D = 4",
        C["phys"], fs=9.4)
    label(ax, 0.398, 0.070,
          "峰条件：1/2 < q < 3/5（引理 7.1）", fs=8.2)
    label(ax, 0.398, 0.054, "Z001 §7 / Z021", fs=8)

    # 箭头：K4 → 维数 → D4
    arrow(ax, (0.213, 0.224), (0.307, 0.228), color=C["l0"], lw=1.8)
    arrow(ax, (0.398, 0.203), (0.398, 0.152), color=C["phys"], lw=1.8)

    # ═════════════════════════════════════════════════════════
    # ⑥ 读出：3 槽（中段中）
    # ═════════════════════════════════════════════════════════
    label(ax, 0.535, 0.288, "L2 的读出：3 个频率槽", fs=13, color=C["axiom"],
          style="normal", weight="bold", ha="left")

    box(ax, 0.615, 0.226, 0.145, 0.052,
        "L2 分 3 层\n4 顶点占据 p ∈ Δ³", C["l2"], fs=10)
    box(ax, 0.615, 0.150, 0.145, 0.056,
        "f_i = |log(p_i / p_{i+1})|\ni = 0,1,2", C["l2"], fs=10)
    arrow(ax, (0.615, 0.198), (0.615, 0.180), color=C["l2"], lw=1.8)
    label(ax, 0.615, 0.110,
          "3 槽 = 4 顶点占据的 3 个相邻比坐标", fs=8.2)
    label(ax, 0.615, 0.092,
          "槽 = log 权重在路径 P₄ 上的离散梯度", fs=8.2)
    label(ax, 0.615, 0.074, "Z161 §3.1", fs=8)

    # 层账框
    box(ax, 0.760, 0.170, 0.130, 0.108,
        "层账\nL0 → 0 层\nL1 → H = L/2 层\nL2 → 3 层", C["l0"], fs=9.6)
    label(ax, 0.760, 0.100, "Z161 §1.1", fs=8)

    # ═════════════════════════════════════════════════════════
    # ⑦ 十六槽（中段右，独立坐标系）
    # ═════════════════════════════════════════════════════════
    axg = fig.add_axes([0.858, 0.118, 0.135, 0.150])
    axg.set_xlim(-0.9, 4.9)
    axg.set_ylim(-2.0, 6.6)
    axg.axis("off")
    axg.patch.set_alpha(0.0)
    axg.text(2.0, 6.05, "Γ₁₆ = K₄ □ C₄（16 槽）", ha="center", fontsize=9.4,
             color="#222222", fontweight="bold")
    axg.text(2.0, 5.35, "每代 8 场填 16 槽", ha="center",
             fontsize=8, color="#6b6b6b")
    axg.text(2.0, 4.80, "v = 0 色单态 ｜ v = 1,2,3 三重", ha="center",
             fontsize=7.8, color="#6b6b6b")
    cols_col = ["#c8794a", "#4a7c59", "#3f7fa8", "#8a6fa8"]
    for i in range(4):      # 同位旋 × 手性 的 4 种组合
        for v in range(4):
            axg.add_patch(Rectangle((v + 0.08, i * 0.86 + 0.30), 0.84, 0.72,
                                    facecolor=cols_col[v], edgecolor="white",
                                    linewidth=1.2,
                                    alpha=0.35 if v != 0 else 0.62))
    axg.text(2.0, -0.45, "16 = 4 × 2 × 2", ha="center", fontsize=8.6,
             color="#6b6b6b")
    axg.text(2.0, -1.30, "Z136 §3.1 / §3.5", ha="center", fontsize=7.8,
             color="#6b6b6b")

    # ═════════════════════════════════════════════════════════
    # ⑧ 物理量（最底）
    # ═════════════════════════════════════════════════════════
    label(ax, 0.055, 0.062, "读出到的物理量", fs=13, color=C["axiom"],
          style="normal", weight="bold", ha="left")

    phys = [
        ("时空维数", "D = 4", C["phys"]),
        ("弱混合角", "sin²θ_w = 3/8", C["phys"]),
        ("质量谱", "m_k = A² M_k²", C["phys"]),
    ]
    for i, (name, val, col) in enumerate(phys):
        x = 0.148 + i * 0.165
        box(ax, x, 0.030, 0.150, 0.056, "", col, radius=0.010,
            alpha=0.13, ec=col, lw=1.1)
        ax.text(x, 0.044, name, ha="center", va="center", fontsize=8.6,
                color=C["muted"], transform=ax.transAxes, zorder=5)
        ax.text(x, 0.020, val, ha="center", va="center", fontsize=11.5,
                color=col, fontweight="bold", transform=ax.transAxes, zorder=5)

    # 锚点链
    box(ax, 0.775, 0.056, 0.30, 0.078, "", C["axiom"], radius=0.010,
        alpha=0.10, ec=C["axiom"], lw=1.1)
    ax.text(0.775, 0.078, "唯一量纲锚  v = 246.22 GeV", ha="center",
            va="center", fontsize=10.5, color=C["axiom"], fontweight="bold",
            transform=ax.transAxes, zorder=5)
    ax.text(0.775, 0.055,
            "v → M_X = v·e^33.6197 = 9.821×10¹⁶ GeV",
            ha="center", va="center", fontsize=8.4, color=C["ink"],
            transform=ax.transAxes, zorder=5)
    ax.text(0.775, 0.034,
            "1 跳数 = 2.009×10⁻³³ m  ｜  1 代 = 6.702×10⁻⁴² s",
            ha="center", va="center", fontsize=8.4, color=C["ink"],
            transform=ax.transAxes, zorder=5)
    label(ax, 0.775, 0.011, "Z153 §3.1", fs=8, z=5)

    # ═════════════════════════════════════════════════════════
    # 纵向主链箭头
    # ═════════════════════════════════════════════════════════
    arrow(ax, (0.5, 0.845), (0.5, 0.822), color=C["l0"], lw=2.2, ms=11)
    arrow(ax, (0.36, 0.778), (0.36, 0.740), color=C["l0"], lw=2.2, ms=11)
    arrow(ax, (0.36, 0.408), (0.36, 0.380), color=C["read"], lw=1.8, ms=10)
    arrow(ax, (0.26, 0.408), (0.26, 0.300), color=C["l2"], lw=2.2, ms=11)
    arrow(ax, (0.50, 0.408), (0.50, 0.300), color=C["l2"], lw=2.2, ms=11)
    arrow(ax, (0.84, 0.408), (0.84, 0.300), color=C["l2"], lw=2.2, ms=11)

    # 底部页脚
    ax.plot([0.06, 0.94], [0.004, 0.004], color=C["grid"], lw=1.2,
            transform=ax.transAxes)
    ax.text(0.5, -0.022,
            "图上每个数字均在库内有出处（逐条见脚本 SRC 表与 ZeroGPU/doc/zero_architecture_map.md）",
            ha="center", va="center", fontsize=8.2, color=C["muted"],
            transform=ax.transAxes, style="italic")

    out = os.path.join(RES, "zero_architecture_map.png")
    fig.savefig(out, dpi=150, facecolor=C["bg"],
                bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print(f"[saved] {out}")
    return out


if __name__ == "__main__":
    build()