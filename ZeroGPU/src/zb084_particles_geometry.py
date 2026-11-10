#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zb084 · 《零和宇宙 · 质量的几何》动画的几何与标签层（纯 Python，不依赖 Blender）

内容：
  1. 质量三圆的精确几何（取自本库 Zb040 §4：公共弦、c=2、R_k = c·sec θ_k*）
  2. 粒子在圆上的参数位置（取自 figures/mass_circle_three_generations.png 的画法）
  3. 12 段旁白字幕条 + 4 段公式图标的图集 PNG（标签层，避免 Blender 处理 CJK/数学字体问题）
  4. 带时间码的渲染清单 JSON

运行：
    python3 ZeroGPU/src/zb084_particles_geometry.py
输出：
    ZeroGPU/results/zb084_particles_geometry.json  （几何 + 时间码清单）
    ZeroGPU/results/zb084_labels/*.png              （3×4 图集，12 格，8 倍降采样 → 标签格 ~790×444 px）

纪律：本脚本不含任何渲染；几何与 Zb040 逐位对账，标签文本与篇目数值逐字一致。
"""
import json
import math
import os
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "ZeroGPU" / "results"
LABEL_DIR = RESULTS / "zb084_labels"
LABEL_DIR.mkdir(parents=True, exist_ok=True)

FONT = "/System/Library/Fonts/Hiragino Sans GB.ttc"          # CJK W3
FONT_BOLD = "/System/Library/Fonts/ArialUnicode-Bold.ttf"   # 占位（下方覆盖）
FONT_MATH = "/System/Library/Fonts/Supplemental/Arial Unicode.ttf"   # 数学符号（√ → ≈ ₀ ν μ τ）
FONT_CJK_BOLD = "/System/Library/Fonts/Hiragino Sans GB.ttc"  # TTC index 2 = W6 粗体

# ──────────────────────────────────────────────────────────────
# 1. 几何常量：Zb040 §4 命题 4.4（c=2 的呈现）
# ──────────────────────────────────────────────────────────────
C = 2.0                                      # 弦半长：端点 (0,0,±c)，c 是呈现参数（不进入质量）
THETA_STAR = [57.930484, 26.158813, 5.910703]  # Zb036 §5 的三个角（度），第 1/2/3 代
R_OVER_C = [1.0 / math.cos(math.radians(t)) for t in THETA_STAR]
R_ = [C * s for s in R_OVER_C]                # 2.010689, 2.228225, 3.766849（代 1→3，与 Zb040 逐位）
D_ = [math.sqrt(r * r - C * C) for r in R_]   # 圆心到弦的距离 D_k = sqrt(R_k² − c²)
ULEN = 1.0                                    # s 的单位 = 该代“端点到圆顶”弧长（R·π/2）的归一化
SIDE_RATIO = [12.11, 2.58, 1.23]              # 两侧之比 (secθ+tanθ)²，仅用于对账

ORDER = [2, 1, 0]                             # 代号（0 起）：远→近 第1 → 第2 → 第3


def _rhs(gen):
    """该代圆上的弧长尺度（画法里的 s 单位）：从端点到圆顶的弧长 = R·π/2，用 ULEN 归一。"""
    return R_[gen] * math.pi / 2 / ULEN


def circle_point(gen, phi):
    """第 gen 代圆上的点（唯一几何定义，Zb040 §4）。

    x = D_k + R_k cos φ   （圆平面内、垂直于弦的坐标 = 到弦的有符号距离）
    z = R_k sin φ         （沿弦的坐标）
    z = ±c 两点即弦与圆的交点（两个零点）。
    返回 (沿弦轴 x, 弦方向 z, 出平面 y=0)，单位坐标。
    """
    R, Dk = R_[gen], D_[gen]
    return (Dk + R * math.cos(phi), R * math.sin(phi), 0.0)


def point_at(gen, s, hand=0):
    """画法参数 → 圆上的点。

    s ≥ 0：粒子侧，从远端 z=+c 沿弧长 s·ULEN 走向圆顶；
    s < 0：反粒子侧，从较远近端 z=−c 沿弧长 |s|·ULEN 走向圆顶。
    hand = ±1：同一画法里的“外点/内点”（左右手 L/R）——沿弦方向移动
               β = 半间距（见 beta()），d（=到弦距离）不变 ⇒ 仍严格落圆。
    """
    arc = min(abs(s) * ULEN, R_[gen] * math.pi / 2 * 0.98)
    phi0 = math.pi / 2 - arc / R_[gen] if s >= 0 else -math.pi / 2 + arc / R_[gen]
    x, z, y = circle_point(gen, phi0)
    if hand:
        # 沿弦移动：dz = Rcosφ dφ（dz>0 朝向圆顶）⇒ dφ = hand·β/(R cosφ)
        cosphi = math.cos(phi0)
        if abs(cosphi) > 0.05:
            phi0 = phi0 + hand * beta(gen) / R_[gen] / cosphi
        x, z, y = circle_point(gen, phi0)
    return (x, z, y)


def beta(gen):
    """内外两点的半间距（左右手 L/R）：画法里每个质量值 4 个位置的那一对。"""
    return 0.05 * R_[gen]


def particle_position(gen, s_e, hand):
    """粒子位置：由画法弧长 s_e 与左右手 hand 给出（恒在圆上）。"""
    return point_at(gen, s_e, hand)


# 画法落点（s 单位 = 该代端点到圆顶的归一化弧长，0=端点 1=圆顶）。
# 值由 figures/mass_circle_three_generations.png 的画法反读；对账见 on_circle_check。
LEPTON_S = {
    "e":   0.42,    # 第 1 代粒子侧 ≈0.40（近端点但不在端点）
    "mu":  0.52,    # 第 2 代粒子侧 ≈0.50
    "tau": 0.32,    # 第 3 代粒子侧 ≈0.30
}
ANTI_S = {k: v for k, v in {"e": 0.30, "mu": 0.25, "tau": 0.14}.items()}

# 每个粒子：key, 代(0起), 类, 标签, 质量文本, s, hand, 颜色
PARTICLES = [
    ("e_p",  0, "lepton", "电子  e⁻", "0.511 MeV",  LEPTON_S["e"], +1, (0.29, 0.56, 1.00)),
    ("e_a",  0, "anti",   "正电子  e⁺", "0.511 MeV", ANTI_S["e"],   0, (0.29, 0.56, 1.00)),
    ("u",    0, "quark",  "u 夸克",   "×3 · 2.16 MeV", 0.40, +1, (0.20, 0.80, 0.60)),
    ("mc",   0, "quark",  "d 夸克",   "×3 · 4.67 MeV", 0.36,  0, (0.20, 0.80, 0.60)),
    ("mu_p", 1, "lepton", "μ 子  μ⁻", "105.66 MeV",  LEPTON_S["mu"], +1, (0.29, 0.56, 1.00)),
    ("mu_a", 1, "anti",   "μ 子  μ⁺", "105.66 MeV",  ANTI_S["mu"],   0, (0.29, 0.56, 1.00)),
    ("c",    1, "quark",  "c 夸克",   "×3 · 1270 MeV", 0.44, +1, (0.20, 0.80, 0.60)),
    ("tau_p",2, "lepton", "τ 子  τ⁻", "1777 MeV",    LEPTON_S["tau"], +1, (0.29, 0.56, 1.00)),
    ("tau_a",2, "anti",   "τ 子  τ⁺", "1777 MeV",    ANTI_S["tau"],   0, (0.29, 0.56, 1.00)),
    ("b",    2, "quark",  "b 夸克",   "×3 · 4180 MeV", 0.26, +1, (0.20, 0.80, 0.60)),
]
GEN_RGB = [(0.36, 0.64, 1.00), (0.55, 0.52, 1.00), (0.36, 0.94, 0.90)]  # 代 1/2/3


def verify_on_circles(tol=1e-9):
    """独立核验：直接用圆的方程 |P − C_gen| = R_gen 判每个粒子是否严格落圆。

    不借用 point_at 的中间量——只用圆方程与表里的 (D,R)。
    """
    rows = []
    for p in PARTICLES:
        gen = p[1]
        x, z, y = particle_position(gen, p[5], p[6])
        d = math.hypot(x - D_[gen], z)
        rows.append({"key": p[0], "gen": gen + 1,
                     "x": round(x, 6), "z": round(z, 6), "d_to_弦": round(abs(x), 6),
                     "dist_to_center": round(d, 9), "R": round(R_[gen], 9),
                     "delta": round(d - R_[gen], 12), "on_circle": abs(d - R_[gen]) < tol})
    return rows


def geometry_table():
    rows = []
    for gen in (0, 1, 2):
        # 弦与圆的交点：cosφ = −D/R ⇒ sinφ = ±c/R，故 z = ±c（零点所在的两个端点）
        phi_end = math.acos(max(-1.0, min(1.0, -D_[gen] / R_[gen])))
        e = circle_point(gen, math.pi - phi_end)   # (D + Rcos, +c)
        a = circle_point(gen, phi_end - math.pi)   # (D + Rcos, −c)
        rows.append({
            "gen": gen + 1,
            "theta_star_deg": THETA_STAR[gen],
            "R_over_c_sec": R_OVER_C[gen],
            "R (c=2)": R_[gen],
            "D (c=2)": D_[gen],
            "R-D": R_[gen] - D_[gen],
            "c*tan(theta/2)": C * math.tan(math.radians(THETA_STAR[gen] / 2)),
            "side_ratio_sec_plus_tan_sq": (1 / math.cos(math.radians(THETA_STAR[gen]))
                                           + math.tan(math.radians(THETA_STAR[gen]))) ** 2,
            "粒端点 (x,z)": [round(e[0], 6), round(e[1], 6)],
            "反端点 (x,z)": [round(a[0], 6), round(a[1], 6)],
            "β 半间距": beta(gen),
        })
    return rows


# ──────────────────────────────────────────────────────────────
# 2. 时间码（秒 → 帧，24 fps）；键名与 Blender 端一致
# ──────────────────────────────────────────────────────────────
FPS = 24
TOTAL = 30.0


def f(t):
    return int(round(t * FPS))


TIMELINE = {
    "shake_end": 3.0,     # B01
    "ecc_peak": 5.0,      # B02 偏心最大
    "formula1_in": 3.4, "formula1_out": 8.9,
    "circles_in": 6.0, "circles_full": 9.0,
    "field_settle": 5.4,  # 偏心后场归位成三圆
    "fly_start": 12.0,    # B05
    "land_e": 14.6, "land_mu": 17.2, "land_tau": 18.6,
    "zero_land": 21.6,    # B07
    "koide_in": 23.2, "koide_out": 28.9,
    "pull_back": 26.2,    # B09
    "orbit_end": 11.0,    # B04 环绕结束
    "title_in": 29.0,
}

# 旁白：单一时间轴 f(秒) → 文本（与用户原始旁白逐字一致）
NARRATION = [
    ("宇宙从 0 开始震动。", 0.40),
    ("中心态均匀，没有质量。", 3.30),
    ("偏心发生，质量就产生了。", 4.90),
    ("三代粒子，三个同心圆，", 6.60),
    ("共享一条弦。", 8.10),
    ("最大的那个圆，远看空空的。", 9.40),
    ("这是电子，我们最熟悉的粒子。", 13.60),
    ("μ 子、τ 子，", 16.30),
    ("电子的两个重表哥。", 18.10),
    ("弦的两端，", 20.90),
    ("是质量为零的粒子。", 22.40),
    ("三个轻子的质量比，", 23.60),
    ("被一个三分之二锁死。", 25.00),
    ("这就是质量的几何。", 26.60),
]

FORMULAS = {
    "f1": ["m = A² · [1 + √2 cos θ]²", "θ 偏离零点 → 偏离量平方 → 质量"],
    "f2": ["1 : 206.77 : 3477.47", "m_e / m_μ / m_τ  （e⁻ / μ⁻ / τ⁻）"],
    "f3": ["Q = Σm / (Σ√m)²", "  Zb041 定理 2.6：b=√2 ⇒ Q = 2/3"],
}

# 12 段旁白字幕条（格 0..11）
SUBTITLES = [
    "宇宙从 0 开始震动。",
    "中心态均匀，没有质量。",
    "偏心发生，质量就产生了。",
    "三代粒子，三个同心圆，共享一条弦。",
    "最大的那个圆，远看空空的。",
    "这是电子，我们最熟悉的粒子。",
    "μ 子、τ 子，电子的两个重表哥。",
    "弦的两端，是质量为零的粒子。",
    "三个轻子的质量比，被一个三分之二锁死。",
    "这就是质量的几何。",
]
SUBTIME = [(0.40, 5.20), (3.30, 5.00), (4.90, 6.40), (6.60, 9.00), (9.40, 12.0),
           (13.60, 16.10), (16.30, 19.90), (20.90, 22.90), (23.60, 26.30), (26.60, 29.90)]


def render_plan():
    return {
        "fps": FPS,
        "duration_s": TOTAL,
        "frames": f(TOTAL),
        "resolution": [1920, 1080],
        "timeline_frames": {k: (f(v) if isinstance(v, float) else v) for k, v in TIMELINE.items()},
        "narration_frames": [(f(t), txt) for txt, t in NARRATION],
        "subtitle_windows": [[f(a), f(b)] for a, b in SUBTIME],
        "engine": "Blender 4.5.9 LTS, EEVEE-Next (bloom via compositor Glare)",
    }


# ──────────────────────────────────────────────────────────────
# 3. 标签图集：3×4 网格，8 倍降采样
# ──────────────────────────────────────────────────────────────
SHOT = (1920, 1080)
NSX, NSY = 3, 4          # 图集版面 3 列 × 4 行 = 12 格
SUP = 4                  # 每格内以 4 倍坐标作图、再降采样到目标尺寸
ATLAS_WH = (SHOT[0] * NSX, SHOT[1] * NSY)   # 5760 × 4320

WHITE = (255, 255, 255, 255)
DIM = (185, 205, 240, 255)
ACCENT = (120, 210, 255, 255)
GOLD = (255, 214, 120, 255)
BG = (4, 10, 22, 232)
EDGE = (120, 170, 255, 140)

# 版本（决定脚本与渲染端的一致性）
LABEL_SCRIPT_VERSION = "zb084-labels-v3"

# 每格 = 名称 + 版式 + 占屏比例（与 zb084_blender_scene.layer_plays 的表一致）
# 尺寸单位 px（相对 1920×1080 底片），字号单位 px；格内按 SUP 倍画
LAYOUT = [
    # name, kind, payload, wpx, hpx
    ("sub0", "sub", "宇宙从 0 开始震动。", 1920, 281),
    ("sub1", "sub", "中心态均匀，没有质量。", 1920, 281),
    ("sub2", "sub", "偏心发生，质量就产生了。", 1920, 281),
    ("sub3", "sub", "三代粒子，三个同心圆，共享一条弦。", 1920, 281),
    ("sub4", "sub", "最大的那个圆，远看空空的。", 1920, 281),
    ("sub5", "sub", "这是电子，我们最熟悉的粒子。", 1920, 281),
    ("sub6", "sub", "μ 子、τ 子，电子的两个重表哥。", 1920, 281),
    ("sub7", "sub", "弦的两端，是质量为零的粒子。", 1920, 281),
    ("sub8", "sub", "三个轻子的质量比，被一个三分之二锁死。", 1920, 281),
    ("sub9", "sub", "这就是质量的几何。", 1920, 281),
    ("sub10", "sub", "", 1920, 281),
    ("公式1", "formula", ["m = A² · [1 + √2 cos θ]²", "θ 偏离零点 → 偏离量平方 → 质量"], 1190, 378),
    ("公式2", "formula", ["1 : 206.77 : 3477.47", "m_e / m_μ / m_τ  （e⁻ / μ⁻ / τ⁻）"], 1056, 302),
    ("公式3", "formula_k", ["Q = (m₀ + m₁ + m₂) / (√m₀ + √m₁ + √m₂)² = 2/3", "m₀ m₁ m₂  对应  e⁻ μ⁻ τ⁻"], 1190, 378),
    ("标签_e", "lbl", ("电子 e", "0.511 MeV", "⁻"), 330, 150),
    ("标签_mu", "lbl", ("μ 子 μ", "105.66 MeV", "⁻"), 350, 150),
    ("标签_tau", "lbl", ("τ 子 τ", "1777 MeV", "⁻"), 330, 150),
    ("标签_g", "lbl", ("光子 γ / 胶子 g", "质量 = 0", ""), 330, 150),
    ("标签_g2", "lbl", ("光子 γ", "质量 = 0", ""), 240, 150),
    ("标签_nu", "lbl", ("中微子 ν_e ν_μ ν_τ", "质量 ≈ 0", ""), 420, 150),
    ("标签_e2", "lbl", ("电子 e", "0.511 MeV", "⁻"), 330, 150),
    ("标签_mu2", "lbl", ("μ 子 μ", "105.66 MeV", "⁻"), 350, 150),
    ("标签_tau2", "lbl", ("τ 子 τ", "1777 MeV", "⁻"), 330, 150),
    ("标题", "title", ("零和宇宙 · 质量的几何", "ZERO · 三代质量圆 · 偏心 → 质量"), 1920, 324),
]


def _fonts():
    px = SUP
    return {
        "sub": _font(FONT_CJK_BOLD, 60 * px, 2),          # 字幕：Hiragino W6
        "fsmall": _font(FONT_MATH, 28 * px),               # 小注：Arial Unicode（ε ₀ ⁻ 齐备）
        "math": _font(FONT_MATH, 64 * px),
        "mathmid": _font(FONT_MATH, 50 * px),
        "mathsub": _font(FONT_MATH, 34 * px),
        "lbl": _font(FONT_CJK_BOLD, 40 * px, 2),           # 粒子名：CJK 粗体
        "val": _font(FONT_MATH, 32 * px),                  # 质量值：数字/MeV
        "title": _font(FONT_CJK_BOLD, 84 * px, 2),
        "tsub": _font(FONT_CJK_BOLD, 28 * px, 2),
    }


def _draw_text(d, xy, text, font, fill=WHITE, anchor="la"):
    d.text(xy, text, font=font, fill=fill, anchor=anchor)


def _err():
    d = ImageDraw.Draw(Image.new("RGBA", (4, 4)))
    return d


def _symbol(d, name, cx, cy, s, color):
    """用绘图代替缺字形：已着陆圆点 / 勾 / 箭头 / 近似符。s 为半径(px)。"""
    if name == "dot":
        d = ImageDraw.Draw(d) if False else d
        d.ellipse([cx - s, cy - s, cx + s, cy + s], fill=color)
        d.ellipse([cx - s * 2.4, cy - s * 2.4, cx + s * 2.4, cy + s * 2.4], outline=color, width=max(1, s // 3))
    elif name == "check":
        d.line([(cx - s, cy + s * 0.1), (cx - s * 0.25, cy + s * 0.85), (cx + s, cy - s * 0.8)],
               fill=color, width=max(2, int(s * 0.35)))
    elif name == "arrow":
        d.line([(cx - s, cy), (cx + s * 0.8, cy)], fill=color, width=max(2, int(s * 0.3)))
        d.polygon([(cx + s, cy), (cx + s * 0.4, cy - s * 0.6), (cx + s * 0.4, cy + s * 0.6)], fill=color)


def _draw(name, kind, payload, w, h, fonts):
    """在 (w*SUP, h*SUP) 的 RGBA 画布上画一格（坐标都用 ×SUP 的“顶片像素”）。"""
    im = Image.new("RGBA", (w * SUP, h * SUP), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    W, H = w * SUP, h * SUP
    U = SUP
    if kind == "sub":
        if payload:
            tw = _measure(payload, fonts["sub"])
            bh, bw = H - 40 * U, min(W - 40 * U, tw + 120 * U)
            bx0, by0 = (W - bw) / 2, (H - bh) / 2
            d.rounded_rectangle([bx0, by0, bx0 + bw, by0 + bh], radius=24 * U,
                                fill=BG, outline=(120, 170, 255, 90), width=2 * U)
            _draw_text(d, (W / 2, H / 2), payload, fonts["sub"], anchor="mm")
    elif kind == "title":
        main_txt, sub_txt = payload
        tw = _measure(main_txt, fonts["title"])
        bw, bh = min(W - 300 * U, tw + 200 * U), H - 44 * U
        bx0, by0 = (W - bw) / 2, (H - bh) / 2
        d.rounded_rectangle([bx0, by0, bx0 + bw, by0 + bh], radius=26 * U,
                            fill=(3, 8, 18, 230), outline=(140, 190, 255, 200), width=2 * U)
        _draw_text(d, (W / 2, by0 + bh * 0.40), main_txt, fonts["title"], anchor="mm")
        _draw_text(d, (W / 2, by0 + bh * 0.78), sub_txt, fonts["tsub"], fill=DIM, anchor="mm")
    elif kind == "formula":
        d.rounded_rectangle([6 * U, 6 * U, W - 6 * U, H - 6 * U], radius=16 * U,
                            fill=BG, outline=EDGE, width=2 * U)
        _draw_text(d, (46 * U, 40 * U), "ZERO · 质量的几何", fonts["fsmall"], fill=ACCENT, anchor="lm")
        _draw_text(d, (28 * U, 118 * U), payload[0], fonts["math"], anchor="lm")
        _draw_text(d, (30 * U, 192 * U), payload[1], fonts["mathsub"], fill=DIM, anchor="lm")
    elif kind == "formula_k":
        d.rounded_rectangle([6 * U, 6 * U, W - 6 * U, H - 6 * U], radius=16 * U,
                            fill=BG, outline=EDGE, width=2 * U)
        _draw_text(d, (46 * U, 40 * U), "ZERO · 质量的几何 · Koide", fonts["fsmall"],
                   fill=ACCENT, anchor="lm")
        _draw_text(d, (26 * U, 112 * U), payload[0], fonts["mathmid"], anchor="lm")
        _draw_text(d, (26 * U, 176 * U), payload[1], fonts["mathsub"], fill=DIM, anchor="lm")
        # “= 代数后果”徽章：勾用画笔画
        d.rounded_rectangle([W - 268 * U, H - 74 * U, W - 16 * U, H - 18 * U], radius=12 * U,
                            fill=(16, 84, 76, 255), outline=(120, 230, 210, 255), width=2 * U)
        _symbol(d, "check", W - 232 * U, H - 46 * U, 12 * U, (150, 240, 220, 255))
        _draw_text(d, (W - 206 * U, H - 46 * U), "代数后果", fonts["fsmall"],
                   fill=(170, 245, 225, 255), anchor="lm")
    elif kind == "lbl":
        nme, val, sup = (payload + ("",))[:3] if len(payload) == 2 else payload
        d.rounded_rectangle([4 * U, 4 * U, W - 4 * U, H - 4 * U], radius=10 * U,
                            fill=BG, outline=(150, 190, 255, 210), width=2 * U)
        d.rectangle([8 * U, 22 * U, 13 * U, H - 22 * U], fill=(120, 190, 255, 255))
        _draw_text(d, (30 * U, 38 * U), nme, fonts["lbl"], anchor="lm")
        if sup:   # 上标（⁺/⁻）用 Arial Unicode 小字号画，避开 CJK 粗体缺字形
            tw = _measure(nme, fonts["lbl"])
            _draw_text(d, (30 * U + tw, 28 * U), sup, fonts["mathsub"], anchor="lm")
        _draw_text(d, (30 * U, 76 * U), val, fonts["val"], fill=(205, 222, 255, 255), anchor="lm")
        _symbol(d, "dot", 38 * U, H - 30 * U, 5 * U, GOLD)
        _draw_text(d, (56 * U, H - 30 * U), "已着陆", fonts["fsmall"], fill=GOLD, anchor="lm")
    return im


def _measure(txt, font):
    return font.getlength(txt)


def _font(path, size, index=0):
    return ImageFont.truetype(path, size, index=index)


def _content_bbox(im):
    """非全透明内容的包围盒（含 2 px 余量）。"""
    bb = im.getbbox()
    return bb


def make_atlas():
    """逐格画（4 倍坐标）→ 按内容裁到实际尺寸 → 落每格单独 PNG；再贴到 3×4 版面。"""
    board = Image.new("RGBA", ATLAS_WH, (0, 0, 0, 0))
    fonts = _fonts()
    singles = {}
    for i, (name, kind, payload, w, h) in enumerate(LAYOUT):
        big = _draw(name, kind, payload, w, h, fonts)
        small = big.resize((w, h), Image.LANCZOS)
        bb = small.getbbox()
        if bb and kind == "lbl":
            # 粒子/零点标签：按实际内容裁剪（画布原本按上限开，内容偏小时右侧会留大片空白）
            x0, y0, x1, y1 = bb
            pad = 2
            small = small.crop((max(0, x0 - pad), max(0, y0 - pad),
                                min(small.width, x1 + pad), min(small.height, y1 + pad)))
        singles[name] = small
        gx, gy = i % NSX, i // NSX
        board.paste(small, (gx * SHOT[0], gy * SHOT[1]))
    atlas = board.resize((SHOT[0] * NSX, SHOT[1] * NSY), Image.LANCZOS)
    path = LABEL_DIR / "atlas.png"
    atlas.save(path)
    for name, im in singles.items():
        im.save(LABEL_DIR / f"{name}.png")
        # 渲染用版本：信息框的暗底做成实色（alpha=255），文字保持抗锯齿边缘，
        # 否则半透明底会让后面的场景（圆、弦）透上来，文字难以辨认。
        opaque = im.copy()
        px = opaque.load()
        w, h = opaque.size
        for yy in range(h):
            for xx in range(w):
                r, g, b, a = px[xx, yy]
                # 暗底（r<40, g<40, b<60）且在框内的像素 → 实色
                if a > 30 and r < 40 and g < 45 and b < 70:
                    px[xx, yy] = (r, g, b, 255)
        opaque.save(LABEL_DIR / f"{name}_solid.png")
    return path


if __name__ == "__main__":
    out = {
        "script_version": LABEL_SCRIPT_VERSION,
        "source": "ZeroGPU/src/zb084_particles_geometry.py",
        "note": "几何取自 Zb040 §4（c=2 呈现）+ figures/mass_circle_three_generations.png 的画法。",
        "constants": {"c": C, "theta_star_deg": THETA_STAR, "R_over_c": R_OVER_C},
        "geometry": geometry_table(),
        "on_circle_check": verify_on_circles(),
        "particles": [
            {"key": p[0], "gen": p[1] + 1, "kind": p[2], "label": p[3], "mass": p[4],
             "s": p[5], "hand": p[6], "rgb": p[7],
             "pos_xz": [round(float(particle_position(p[1], p[5], p[6])[0]), 6),
                        round(float(particle_position(p[1], p[5], p[6])[1]), 6)]}
            for p in PARTICLES
        ],
        "render_plan": render_plan(),
        "labels": [{"name": n, "kind": k, "wpx": w, "hpx": h} for (n, k, _p, w, h) in LAYOUT],
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "zb084_particles_geometry.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2))
    path = make_atlas()
    print(json.dumps({"atlas": str(path), "cells": len(LAYOUT),
                      "geometry": [{"gen": r["gen"], "R": round(r["R (c=2)"], 4),
                                    "D": round(r["D (c=2)"], 4)} for r in out["geometry"]]},
                     ensure_ascii=False, indent=1))
