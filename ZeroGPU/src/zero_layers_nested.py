#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_layers_nested.py — Zero 分层结构图（嵌套版）

用户反馈：第一版画成了"流程图"，不像"一层一层"；L1' 只当了一行表格，
三个亚层、亚亚层都没画出来。

本版改为**嵌套同心结构**：外层包内层，一眼看出"层里还有层"。

结构（逐条有出处，见 SRC 表与 doc）：
  L0   底层        0 层（不加结构）              Z085 / Z129
  L1'  闭合类层    闭合类 [w]、计数测度、记录 λ、读出   Z016 §2.1 / Z001 §4.1
  L1   历史层      H = L/2 层（高度层，静态）     Z129 §2 / Z161 §1.1
                   └ 再分坐标：Z2 对合（Z047）＋ 面积 A_h（Z048）
  L2   演化层      3 个亚层（频率槽 / 代）        Z129 §2 / Z161 §1.1 / §3.1
                   └ 三个亚层各含自己的亚亚层（用户指定归属）：
                       亚层 1 →  4 个亚亚层（K₄ 的 4 个顶点，毁灭时不动）
                       亚层 2 → 18 个亚亚层（3 代 × 6 边，只留最上 2）
                       亚层 3 →  6 个亚亚层（K₄ 的 6 条边，全灭后补回）
                     zero_tau_6_18_4_graph.py：NB=6 边、NG=3 代、NM=18、NT=4
  R    读出面      不占层指标                   Z016 命题 1

层序：按用户指定 L0 → L1′ → L1 → L2。
（注：库内 Z016 §2.1、Z001 §4.1 的层表写作 L0 · L1 · L1′ · L2；本图按用户口径。）

复跑： /usr/bin/python3 zero_layers_nested.py
输出： ../results/zero_layers_nested.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Wedge, Polygon
from matplotlib import font_manager as _fm

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, os.pardir, "results")

# ── 中文字体（macOS 自带）──────────────────────────────────
_CJK = [
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
    "/System/Library/Fonts/STHeiti Medium.ttc",
    "/System/Library/Fonts/Supplemental/Songti.ttc",
    "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
]
_p = next((p for p in _CJK if os.path.exists(p)), None)
if _p:
    _fm.fontManager.addfont(_p)
    matplotlib.rcParams["font.family"] = [
        _fm.FontProperties(fname=_p).get_name(), "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

# ── 配色 ───────────────────────────────────────────────────
C = {
    "l0":   "#1f3b73",
    "l1":   "#2f6fa0",
    "l1p":  "#4f96bd",
    "l2":   "#c8794a",
    "sub":  "#e0a06a",
    "ssub": "#f0c9a8",
    "read": "#8a6fa8",
    "ink":  "#222222",
    "mut":  "#6b6b6b",
    "bg":   "#ffffff",
}


def cjk(fig):
    return fig


# ══════════════════════════════════════════════════════════
def build():
    fig = plt.figure(figsize=(20, 13.5), facecolor=C["bg"])
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.968, "Zero 的分层结构", ha="center", va="center",
            fontsize=28, color=C["l0"], fontweight="bold")
    ax.text(0.5, 0.939,
            "一层包一层：L0 → L1′ → L1（H 层）→ L2（3 亚层：含 4／18／6 个亚亚层）",
            ha="center", va="center", fontsize=13, color=C["mut"], style="italic")
    ax.plot([0.04, 0.96], [0.921, 0.921], color="#e8e8e8", lw=1.4)

    # ══════════════════════════════════════════════════════
    # 左：嵌套同心方框（主体）
    # 层序按用户指定：L0 → L1′ → L1 → L2
    # ══════════════════════════════════════════════════════
    layers = [
        # (x, y, w, h, 边框色, 填充色, 标题, 副题, 层数标注)
        (0.030, 0.040, 0.585, 0.822, C["l0"],  "#eef2f8",
         "L0  底层", "0 层 —— 公理条款、非负整数重数、局域补偿移动（Z085：L0 不加结构）",
         "0 层"),
        (0.072, 0.085, 0.505, 0.737, C["l1p"], "#e3eef6",
         "L1′  全局闭合类层", "闭合类 [w] ｜ 计数测度 μ ｜ 记录 λ ｜ 读出（π, ω, K, σ_t）",
         "记录层"),
        (0.114, 0.135, 0.425, 0.647, C["l1"],  "#e8f0f7",
         "L1  历史层", "H = L/2 层（高度层，静态读出）｜ 再分：ℤ₂ 对合 ＋ 面积 A_h",
         "H = L/2 层"),
    ]
    for (x, y, w, h, ec, fc, title, sub, cnt) in layers:
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0,rounding_size=0.014",
                                    linewidth=2.4, edgecolor=ec, facecolor=fc,
                                    transform=ax.transAxes, zorder=1))
        ax.text(x + 0.022, y + h - 0.026, title, ha="left", va="center",
                fontsize=13.5, color=ec, fontweight="bold",
                transform=ax.transAxes, zorder=4)
        ax.text(x + 0.022, y + h - 0.052, sub, ha="left", va="center",
                fontsize=9, color=C["ink"], transform=ax.transAxes, zorder=4)
        # 层数标注：放在该层框右上角（标题行的最右端）
        ax.text(x + w - 0.022, y + h - 0.026, cnt, ha="right", va="center",
                fontsize=10.5, color=ec, fontweight="bold",
                transform=ax.transAxes, zorder=5)
        ax.text(x + w - 0.018, y + 0.016, "无自身更新律", ha="right",
                va="bottom", fontsize=8.4, color=C["mut"], style="italic",
                transform=ax.transAxes, zorder=4)

    # ── L2（最内层，重点：3 亚层 + 亚亚层）──────────────
    x2, y2, w2, h2 = 0.156, 0.190, 0.345, 0.415
    ax.add_patch(FancyBboxPatch((x2, y2), w2, h2,
                                boxstyle="round,pad=0,rounding_size=0.014",
                                linewidth=3.0, edgecolor=C["l2"],
                                facecolor="#fdf3ec",
                                transform=ax.transAxes, zorder=2))
    ax.text(x2 + 0.020, y2 + h2 - 0.026, "L2  演化层",
            ha="left", va="center", fontsize=13.5, color=C["l2"],
            fontweight="bold", transform=ax.transAxes, zorder=5)
    ax.text(x2 + 0.128, y2 + h2 - 0.026, "★ 物理世界在这里",
            ha="left", va="center", fontsize=9.5, color=C["l2"],
            transform=ax.transAxes, zorder=5)
    ax.text(x2 + 0.020, y2 + h2 - 0.050,
            "唯一有自身更新律的层 ｜ 分 3 个亚层",
            ha="left", va="center", fontsize=9, color=C["ink"],
            transform=ax.transAxes, zorder=5)
    ax.text(x2 + w2 - 0.030, y2 + h2 - 0.026, "3 个亚层",
            ha="right", va="center", fontsize=10.5, color=C["l2"],
            fontweight="bold", transform=ax.transAxes, zorder=5)

    # 3 个亚层（横排三条）：每个亚层各含自己的亚亚层数 4／18／6
    subs = [
        ("亚层 1", "f₀ = |log(p₀/p₁)|", 4,  "顶：K₄ 的 4 个顶点", "#4a7c59"),
        ("亚层 2", "f₁ = |log(p₁/p₂)|", 18, "中：3 代 × 6 边", "#c8794a"),
        ("亚层 3", "f₂ = |log(p₂/p₃)|", 6,  "底：K₄ 的 6 条边", "#a04a4a"),
    ]
    sy = y2 + 0.028
    sw, sh = 0.100, 0.285
    for i, (name, form, nss, sdesc, col) in enumerate(subs):
        sx = x2 + 0.012 + i * 0.110
        ax.add_patch(FancyBboxPatch((sx, sy), sw, sh,
                                    boxstyle="round,pad=0,rounding_size=0.010",
                                    linewidth=1.8, edgecolor=col,
                                    facecolor="#fdf8f2",
                                    transform=ax.transAxes, zorder=4))
        ax.text(sx + sw / 2, sy + sh - 0.020, name, ha="center", va="center",
                fontsize=10, color=C["l2"], fontweight="bold",
                transform=ax.transAxes, zorder=5)
        ax.text(sx + sw / 2, sy + sh - 0.042, form, ha="center", va="center",
                fontsize=7.6, color=C["ink"], transform=ax.transAxes, zorder=5)
        # 亚亚层：该亚层下含 nss 个（4／18／6）
        ax.text(sx + sw / 2, sy + sh - 0.070, f"含 {nss} 个亚亚层",
                ha="center", va="center", fontsize=9.2, color=col,
                fontweight="bold", transform=ax.transAxes, zorder=5)
        ax.text(sx + sw / 2, sy + sh - 0.090, sdesc,
                ha="center", va="center", fontsize=7.2, color=C["mut"],
                transform=ax.transAxes, zorder=5)
        # 按个数画小方块（4／18／6 全部画满，两行布局）
        bx0, bw = sx + 0.012, sw - 0.024
        rows_ = 2
        per = int(np.ceil(nss / rows_))
        for t in range(nss):
            r_, c_ = divmod(t, per)
            yy = sy + 0.148 - r_ * 0.026
            xx = bx0 + bw * c_ / per
            ax.add_patch(Rectangle((xx, yy), bw / per * 0.70, 0.019,
                                   facecolor=C["ssub"], edgecolor=col,
                                   linewidth=0.4, transform=ax.transAxes,
                                   zorder=5))
        ax.text(sx + sw / 2, sy + 0.030,
                "破坏：顶不动／中留 2／底全灭",
                ha="center", va="center", fontsize=6.6, color=C["mut"],
                transform=ax.transAxes, zorder=5)

    ax.text(x2 + w2 / 2, y2 + 0.010,
            "3 个亚层分别含 4 ／ 18 ／ 6 个亚亚层（顶顶点／中3代×6边／底6边）",
            ha="center", va="center", fontsize=8.4, color=C["mut"],
            transform=ax.transAxes, zorder=5)

    # ══════════════════════════════════════════════════════
    # 右上：毁灭-重播种规则（亚亚层的动力学）
    # ══════════════════════════════════════════════════════
    ax.add_patch(FancyBboxPatch((0.640, 0.560), 0.335, 0.335,
                                boxstyle="round,pad=0,rounding_size=0.014",
                                linewidth=1.6, edgecolor=C["l2"],
                                facecolor="#fdf8f2",
                                transform=ax.transAxes, zorder=1))
    ax.text(0.8075, 0.868, "亚亚层上的动力学：毁灭与重播种",
            ha="center", va="center", fontsize=12.5, color=C["l2"],
            fontweight="bold", transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.842, "每 T = 4 代毁灭一次（Z022：周期 T = L = 4）",
            ha="center", va="center", fontsize=9, color=C["mut"],
            transform=ax.transAxes, zorder=4)
    rules = [
        ("亚层 1", "4 个亚亚层", "K₄ 的 4 个顶点 ｜ 毁灭时不动", "#4a7c59"),
        ("亚层 2", "18 个亚亚层", "3 代 × 6 边 ｜ 只留最上 2", "#c8794a"),
        ("亚层 3", "6 个亚亚层", "K₄ 的 6 条边 ｜ 全灭后补回", "#a04a4a"),
    ]
    for i, (nm, act, note, col) in enumerate(rules):
        yy = 0.790 - i * 0.062
        ax.add_patch(FancyBboxPatch((0.662, yy - 0.024), 0.075, 0.048,
                                    boxstyle="round,pad=0,rounding_size=0.008",
                                    linewidth=1.6, edgecolor=col,
                                    facecolor=col, alpha=0.16,
                                    transform=ax.transAxes, zorder=3))
        ax.text(0.6995, yy, nm, ha="center", va="center", fontsize=10,
                color=col, fontweight="bold", transform=ax.transAxes, zorder=4)
        ax.text(0.748, yy + 0.010, act, ha="left", va="center", fontsize=10,
                color=C["ink"], fontweight="bold",
                transform=ax.transAxes, zorder=4)
        ax.text(0.748, yy - 0.012, note, ha="left", va="center", fontsize=8,
                color=C["mut"], transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.604, "份额稳定在 [0.05, 0.198, 0.752]",
            ha="center", va="center", fontsize=9.2, color=C["ink"],
            transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.583, "长跑 10000 代无 NaN ｜ top2 冻结在 [0,15]",
            ha="center", va="center", fontsize=8.6, color=C["mut"],
            transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.567, "6/18/4 是具名输入（不是从底粗粒化导出的）",
            ha="center", va="center", fontsize=8, color=C["mut"],
            style="italic", transform=ax.transAxes, zorder=4)

    # ══════════════════════════════════════════════════════
    # 右中：读出面（不占层）
    # ══════════════════════════════════════════════════════
    ax.add_patch(FancyBboxPatch((0.640, 0.400), 0.335, 0.135,
                                boxstyle="round,pad=0,rounding_size=0.014",
                                linewidth=1.6, edgecolor=C["read"],
                                facecolor="#f4f0f8",
                                transform=ax.transAxes, zorder=1))
    ax.text(0.8075, 0.508, "R  读出面（不是层）",
            ha="center", va="center", fontsize=12.5, color=C["read"],
            fontweight="bold", transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.480,
            "π（自识别划分）｜ ω ｜ K_ω = −log ω ｜ σ_t ｜ GNS 内积",
            ha="center", va="center", fontsize=8.6, color=C["ink"],
            transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.450,
            "读出依赖态（σ_t 显含 ρ），故它是层上的函数，不是第五层",
            ha="center", va="center", fontsize=8.4, color=C["mut"],
            transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.420, "Z016 命题 1 ｜ Z001 推论 4.2",
            ha="center", va="center", fontsize=8, color=C["mut"],
            style="italic", transform=ax.transAxes, zorder=4)

    # ══════════════════════════════════════════════════════
    # 右下：公理（最外层）
    # ══════════════════════════════════════════════════════
    ax.add_patch(FancyBboxPatch((0.640, 0.235), 0.335, 0.135,
                                boxstyle="round,pad=0,rounding_size=0.014",
                                linewidth=2.0, edgecolor=C["l0"],
                                facecolor="#eef2f8",
                                transform=ax.transAxes, zorder=1))
    ax.text(0.8075, 0.344, "公理 Z₀（零不停留）",
            ha="center", va="center", fontsize=12.5, color=C["l0"],
            fontweight="bold", transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.316, "① 零不停留    ② 从不停歇    ③ 不设概率",
            ha="center", va="center", fontsize=9.4, color=C["ink"],
            transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.288,
            "体系内唯一公理；「零」＝商 {±1}/ι",
            ha="center", va="center", fontsize=8.6, color=C["mut"],
            transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.260, "Z001 §2.1",
            ha="center", va="center", fontsize=8, color=C["mut"],
            style="italic", transform=ax.transAxes, zorder=4)

    # ══════════════════════════════════════════════════════
    # 底部：L1 高度层的显式展开（H = L/2 层）
    # ══════════════════════════════════════════════════════
    ax.add_patch(FancyBboxPatch((0.640, 0.055), 0.335, 0.160,
                                boxstyle="round,pad=0,rounding_size=0.014",
                                linewidth=1.6, edgecolor=C["l1"],
                                facecolor="#e8f0f7",
                                transform=ax.transAxes, zorder=1))
    ax.text(0.8075, 0.208, "L1 的 H 层：高度 h = 1 … H（H = L/2）",
            ha="center", va="center", fontsize=11.5, color=C["l1"],
            fontweight="bold", transform=ax.transAxes, zorder=4)
    # 高度层示意：h=1..4，条数 b_h
    bh = [16, 38, 14, 2]
    for i, b in enumerate(bh):
        bx = 0.662 + i * 0.048
        hgt = 0.016 + b / 40 * 0.048
        ax.add_patch(Rectangle((bx, 0.118), 0.032, hgt,
                               facecolor=C["l1"], alpha=0.55,
                               edgecolor=C["l1"], linewidth=1.0,
                               transform=ax.transAxes, zorder=4))
        ax.text(bx + 0.016, 0.114, f"h={i+1}", ha="center", va="top",
                fontsize=7.6, color=C["ink"], transform=ax.transAxes, zorder=5)
        ax.text(bx + 0.016, 0.118 + hgt + 0.006, str(b),
                ha="center", va="bottom", fontsize=7.4, color=C["mut"],
                transform=ax.transAxes, zorder=5)
    ax.text(0.8075, 0.090,
            "L = 8 时 b_h = [16, 38, 14, 2]，总数 C(8,4) = 70",
            ha="center", va="center", fontsize=8.4, color=C["ink"],
            transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.071,
            "每层 b_h 恒偶、顶层 b_H = 2（ℤ₂ 轨道）｜ 面积 A_h 为第二个再分坐标",
            ha="center", va="center", fontsize=8.2, color=C["mut"],
            transform=ax.transAxes, zorder=4)
    ax.text(0.8075, 0.057, "Z161 §1.3 ｜ Z047 ｜ Z048",
            ha="center", va="center", fontsize=7.8, color=C["mut"],
            style="italic", transform=ax.transAxes, zorder=4)

    # 层数总账（放到左下，不压标题）
    ax.add_patch(FancyBboxPatch((0.030, 0.880), 0.585, 0.036,
                                boxstyle="round,pad=0,rounding_size=0.008",
                                linewidth=1.2, edgecolor=C["l0"],
                                facecolor="#eef2f8",
                                transform=ax.transAxes, zorder=3))
    ax.text(0.3225, 0.898,
            "层账：L0 → 0 层　｜　L1′ → 记录层　｜　L1 → H = L/2 层　｜　"
            "L2 → 3 个亚层（分别含 4／18／6 个亚亚层）",
            ha="center", va="center", fontsize=10.5, color=C["l0"],
            fontweight="bold", transform=ax.transAxes, zorder=5)

    out = os.path.join(RES, "zero_layers_nested.png")
    fig.savefig(out, dpi=150, facecolor=C["bg"],
                bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    print(f"[saved] {out}")
    return out


if __name__ == "__main__":
    build()