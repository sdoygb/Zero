#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_layers_spheres.py — Zero 分层结构：同心球（一层一层的球）

用户指定：做成一个一层一层的球。
  外层 → 内层：L0 → L1′ → L1 → L2
  最内的 L2 球一分为三（三个扇区）：
    第 1 份 → 再分 4 份
    第 2 份 → 再分 18 份
    第 3 份 → 再分 6 份

数据出处（逐条）：
  L0=0 层（Z085/Z129）；L1′=闭合类层（Z016 §2.1/Z001 §4.1）；
  L1=H=L/2 层（Z129 §2/Z161 §1.1）；L2=3 亚层（Z161 §3.1 频率槽）；
  亚亚层 4/18/6 = 顶4顶点 / 中3代×6边 / 底6边
    （zero_tau_6_18_4_graph.py: NB,NG,NM,NT = 6,3,18,4）

层序按用户口径 L0 → L1′ → L1 → L2（库内层表 Z016 §2.1 写作 L0·L1·L1′·L2）。

复跑： /usr/bin/python3 zero_layers_spheres.py
输出： ../results/zero_layers_spheres.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, FancyBboxPatch, Polygon
from matplotlib import font_manager as _fm

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, os.pardir, "results")

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

C = {
    "l0":  "#1f3b73",
    "l1p": "#4f96bd",
    "l1":  "#2f6fa0",
    "l2":  "#c8794a",
    "s1":  "#4a7c59",
    "s2":  "#c8794a",
    "s3":  "#a04a4a",
    "ink": "#222222",
    "mut": "#6b6b6b",
    "bg":  "#ffffff",
}


def draw_shell(ax, cx, cy, r, ec, fc, alpha, lw, z=1):
    """画一个球（圆），带一点高光与暗边做出立体感。"""
    ax.add_patch(Circle((cx, cy), r, facecolor=fc, edgecolor=ec,
                        linewidth=lw, alpha=alpha, zorder=z))
    # 高光
    ax.add_patch(Circle((cx - r * 0.30, cy + r * 0.32), r * 0.42,
                        facecolor="white", edgecolor="none",
                        alpha=0.28, zorder=z + 0.1))
    # 暗边（下半弧）
    ax.add_patch(Wedge((cx, cy), r * 0.995, 200, 340, width=r * 0.10,
                       facecolor=ec, edgecolor="none", alpha=0.13,
                       zorder=z + 0.2))


def build():
    fig = plt.figure(figsize=(19, 13), facecolor=C["bg"])
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_aspect("equal")

    ax.text(0.5, 0.968, "Zero 的分层结构：一层一层的球",
            ha="center", va="center", fontsize=27, color=C["l0"],
            fontweight="bold", transform=ax.transAxes)
    ax.text(0.5, 0.938,
            "由外到内  L0 → L1′ → L1 → L2　｜　"
            "最内的 L2 球一分为三，三份各再分 4 ／ 18 ／ 6 份",
            ha="center", va="center", fontsize=12.5, color=C["mut"],
            style="italic", transform=ax.transAxes)
    ax.plot([0.03, 0.97], [0.920, 0.920], color="#e8e8e8", lw=1.4,
            transform=ax.transAxes)

    # ══════════════════════════════════════════════════════
    # 主球：同心壳（画在 axes 坐标，用 aspect equal 保证是圆）
    # ══════════════════════════════════════════════════════
    CX, CY = 0.40, 0.435          # 球心
    shells = [
        # (半径, 边色, 填色, alpha, 层名, 副题, 层数)
        (0.300, C["l0"],  "#dde6f2", 0.55, "L0  底层",
         "公理条款 ｜ 非负整数重数 ｜ 局域补偿移动", "0 层"),
        (0.246, C["l1p"], "#d8e8f4", 0.60, "L1′  全局闭合类层",
         "闭合类 [w] ｜ 计数测度 μ ｜ 记录 λ ｜ 读出", "记录层"),
        (0.192, C["l1"],  "#cfe3f2", 0.62, "L1  历史层",
         "H = L/2 层（高度层）｜ ℤ₂ 对合 ＋ 面积 A_h", "H = L/2 层"),
    ]
    # 三层各自的标注角度（充分分开：L0 左、L1′ 左上、L1 正上）
    shell_angles = [178, 145, 108]     # L0 / L1′ / L1
    for i, (r, ec, fc, al, name, sub, cnt) in enumerate(shells):
        draw_shell(ax, CX, CY, r, ec, fc, al, 2.6, z=1)
        ang = np.deg2rad(shell_angles[i])
        lx = CX + np.cos(ang) * (r + 0.014)
        ly = CY + np.sin(ang) * (r + 0.014)
        # 引线：从壳上拉到标签
        ax.annotate("", xy=(CX + np.cos(ang) * r * 0.99,
                            CY + np.sin(ang) * r * 0.99),
                    xytext=(lx, ly),
                    arrowprops=dict(arrowstyle="-", color=ec, lw=1.3),
                    transform=ax.transAxes, zorder=6)
        ax.text(lx, ly + 0.026, name, ha="center", va="bottom",
                fontsize=11.5, color=ec, fontweight="bold",
                transform=ax.transAxes, zorder=8)
        ax.text(lx, ly + 0.011, sub, ha="center", va="bottom",
                fontsize=7.6, color=C["ink"], transform=ax.transAxes, zorder=8)
        ax.text(lx, ly - 0.003, f"分 {cnt}", ha="center", va="top",
                fontsize=8.8, color=ec, fontweight="bold",
                transform=ax.transAxes, zorder=8)

    # ── 最内层：L2 球，一分为三 ──────────────────────────
    R2 = 0.140
    # 三个扇区：各 120°
    # 每份内再分 4 / 18 / 6 份（用等分扇区表达）
    sectors = [
        (0,   120, C["s1"], 4,  "第 1 份", "4 份（K₄ 顶点）",  62),
        (120, 240, C["s2"], 18, "第 2 份", "18 份（3 代×6 边）", 196),
        (240, 360, C["s3"], 6,  "第 3 份", "6 份（K₄ 的边）",   300),
    ]
    for (a0, a1, col, nsub, sname, sdesc, langle) in sectors:
        # 三份的底色
        ax.add_patch(Wedge((CX, CY), R2, a0, a1, facecolor=col,
                           edgecolor="white", linewidth=2.2, alpha=0.30,
                           transform=ax.transAxes, zorder=3))
        # 该份内的 nsub 个子扇区
        step = (a1 - a0) / nsub
        for t in range(nsub):
            aa = a0 + t * step
            ax.add_patch(Wedge((CX, CY), R2 * 0.94, aa + step * 0.06,
                               aa + step * 0.94, facecolor=col,
                               edgecolor="white", linewidth=0.7,
                               alpha=0.88, transform=ax.transAxes, zorder=4))
        # 三份的标注（引线方向单独指定，避开球壳标签）
        amid = np.deg2rad((a0 + a1) / 2)          # 扇区中角（用于起点）
        lang = np.deg2rad(langle)                 # 标签方向（单独指定）
        sx_ = CX + np.cos(amid) * R2 * 0.99
        sy_ = CY + np.sin(amid) * R2 * 0.99
        lx = CX + np.cos(lang) * (R2 + 0.070)
        ly = CY + np.sin(lang) * (R2 + 0.070)
        ax.annotate("", xy=(sx_, sy_), xytext=(lx, ly),
                    arrowprops=dict(arrowstyle="-", color=col, lw=1.8),
                    transform=ax.transAxes, zorder=6)
        ha = "left" if np.cos(lang) > 0.05 else ("right" if np.cos(lang) < -0.05 else "center")
        dx = 0.014 if ha == "left" else (-0.014 if ha == "right" else 0)
        ax.text(lx + dx, ly + 0.016, f"{sname}：{nsub} 个亚亚层",
                ha=ha, va="center", fontsize=10.5, color=col,
                fontweight="bold", transform=ax.transAxes, zorder=8)
        ax.text(lx + dx, ly - 0.012, sdesc, ha=ha, va="center", fontsize=8.2,
                color=C["mut"], transform=ax.transAxes, zorder=8)

    # L2 球心标注（只放一个标记，其余移到球外，避免重叠）
    ax.text(CX, CY, "L2", ha="center", va="center", fontsize=11,
            color="white", fontweight="bold", transform=ax.transAxes,
            zorder=9,
            bbox=dict(boxstyle="circle,pad=0.34", facecolor=C["l2"],
                      edgecolor="white", linewidth=1.4, alpha=0.95))

    # L2 标题与"物理世界"标注移到球的正下方（球外）
    ax.text(CX, CY - 0.175, "L2  演化层　分 3 个亚层",
            ha="center", va="center", fontsize=12.5, color=C["l2"],
            fontweight="bold", transform=ax.transAxes, zorder=9)
    ax.text(CX, CY - 0.200, "★ 物理世界在这里（唯一有自身更新律的层）",
            ha="center", va="center", fontsize=9.6, color=C["l2"],
            transform=ax.transAxes, zorder=9)

    # ══════════════════════════════════════════════════════
    # 右侧：读出面 + 公理 + 说明
    # ══════════════════════════════════════════════════════
    def panel(x, y, w, h, ec, fc, title, lines, sub=None, fs=12):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0,rounding_size=0.014",
                                    linewidth=1.6, edgecolor=ec, facecolor=fc,
                                    transform=ax.transAxes, zorder=2))
        ax.text(x + w / 2, y + h - 0.026, title, ha="center", va="center",
                fontsize=fs, color=ec, fontweight="bold",
                transform=ax.transAxes, zorder=6)
        for i, ln in enumerate(lines):
            ax.text(x + w / 2, y + h - 0.055 - i * 0.022, ln,
                    ha="center", va="center", fontsize=8.6, color=C["ink"],
                    transform=ax.transAxes, zorder=6)
        if sub:
            ax.text(x + w / 2, y + 0.012, sub, ha="center", va="center",
                    fontsize=7.8, color=C["mut"], style="italic",
                    transform=ax.transAxes, zorder=6)

    panel(0.745, 0.660, 0.240, 0.230, C["l2"], "#fdf8f2",
          "L2 的三份与其亚亚层",
          ["第 1 份 →  4 个亚亚层",
           "　（K₄ 的 4 个顶点；毁灭时不动）",
           "第 2 份 → 18 个亚亚层",
           "　（3 代 × 6 边；只留最上 2）",
           "第 3 份 →  6 个亚亚层",
           "　（K₄ 的 6 条边；全灭后补回）"],
          "每 T = 4 代毁灭一次（Z022：T = L = 4）")

    panel(0.745, 0.430, 0.240, 0.205, C["l0"], "#eef2f8",
          "四层与层数",
          ["L0   底层        0 层",
           "L1′  闭合类层    记录层",
           "L1   历史层      H = L/2 层",
           "L2   演化层      3 个亚层",
           "　　　　　　　　（4 ／ 18 ／ 6 亚亚层）"],
          "层序按用户口径：L0 → L1′ → L1 → L2")

    panel(0.745, 0.235, 0.240, 0.170, "#8a6fa8", "#f4f0f8",
          "R  读出面（不是层）",
          ["π（自识别划分）｜ ω ｜ K_ω = −log ω",
           "σ_t ｜ GNS 内积",
           "读出依赖态，是层上的函数"],
          "Z016 命题 1 ｜ Z001 推论 4.2")

    panel(0.745, 0.055, 0.240, 0.155, C["l0"], "#eef2f8",
          "公理 Z₀（零不停留）",
          ["① 零不停留　② 从不停歇　③ 不设概率",
           "体系内唯一公理",
           "「零」＝商 {±1}/ι"],
          "Z001 §2.1")

    # 球体说明
    ax.text(0.40, 0.055,
            "球面 = 一层；越靠内 = 越靠演化层；层与层之间是包住的关系",
            ha="center", va="center", fontsize=9, color=C["mut"],
            transform=ax.transAxes, zorder=6, style="italic")

    out = os.path.join(RES, "zero_layers_spheres.png")
    fig.savefig(out, dpi=150, facecolor=C["bg"],
                bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    print(f"[saved] {out}")
    return out


if __name__ == "__main__":
    build()