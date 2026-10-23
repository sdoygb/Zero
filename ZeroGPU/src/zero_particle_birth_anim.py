#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_particle_birth_anim.py — 「从 0 演化出粒子」动画（无 GUI，纯 Python 逐帧渲染 → ffmpeg）

内容（结构示意，不给数值主张）：空舞台 → seed → K₄ 稳态 → 三个频率槽 → 三代。
口径：三槽 ↔ 三代是结构计数，不是质量数值（Z161 §5）；三槽里只有 log2 属内部 2cos 族，
另两条是输入（Z132 §五）。动画不做数值刻度。

复跑： /usr/bin/python3 zero_particle_birth_anim.py [--frames 240] [--fps 20]
输出： ../results/anim_frames/f_%04d.png
"""
import os
import sys
import math
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, os.pardir, "results")
FRAMES_DIR = os.path.join(RES, "anim_frames")

W, H = 1280, 800
BG = "#070910"
C_GRID = "#2e4a5c"
C_SEED = "#ffd45e"
C_VERT = "#6fb3ea"
C_EDGE = "#e2603c"
C_SLOT = "#5fdcc4"
C_PART = "#eaf2ff"
C_TEXT = "#c8d6e5"

# CJK 字体
_CJK = [
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
    "/System/Library/Fonts/STHeiti Medium.ttc",
    "/System/Library/Fonts/Supplemental/Songti.ttc",
]
_p = next((p for p in _CJK if os.path.exists(p)), None)
if _p:
    fm.fontManager.addfont(_p)
    plt.rcParams["font.family"] = ["Hiragino Sans GB", "DejaVu Sans"]


# ---------------- 数学 ----------------
def norm(v):
    n = np.linalg.norm(v)
    return v / (n if n > 1e-12 else 1.0)


def look_at(eye, target, up=(0.0, 0.0, 1.0)):
    f = norm(np.asarray(target, float) - np.asarray(eye, float))
    r = norm(np.cross(f, np.asarray(up, float)))
    u = np.cross(r, f)
    return np.stack([r, u, -f])          # 世界 → 相机（相机看向 -z）


def project(pts, eye, M, fov_deg=34.0):
    """pts: (N,3) 世界坐标 → (N,2) 像素坐标, (N,) 深度"""
    v = (np.asarray(pts, float) - np.asarray(eye, float)) @ M.T
    z = -v[:, 2]
    fpx = (H / 2.0) / math.tan(math.radians(fov_deg) / 2.0)
    x = W / 2.0 + fpx * v[:, 0] / z
    y = H / 2.0 - fpx * v[:, 1] / z
    return np.stack([x, y], axis=1), z


def ease(x):
    x = min(1.0, max(0.0, x))
    return x * x * (3 - 2 * x)


def clamp01(x):
    return 0.0 if x < 0 else (1.0 if x > 1 else x)


# ---------------- 场景参数 ----------------
HALF = 6.0
NC = 24
TET_C = np.array([0.0, 0.0, 3.4])
TET_S = 2.2
TET_V = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) * (TET_S / math.sqrt(3))
TET_V = TET_V + TET_C
TET_E = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

SLOT_R, SLOT_Z = 2.6, 7.4
SLOT_AZ = [math.pi / 2 + i * 2 * math.pi / 3 for i in range(3)]
SLOT_P = np.array([[SLOT_R * math.cos(a), SLOT_R * math.sin(a), SLOT_Z] for a in SLOT_AZ])
SLOT_TXT = ["0.4055", "0.6931", "0.9808"]

PART_P = np.array([[1.55 * SLOT_R * math.cos(a), 1.55 * SLOT_R * math.sin(a), SLOT_Z + 2.7]
                   for a in SLOT_AZ])
GEN_TXT = ["第 1 代", "第 2 代", "第 3 代"]


# ---------------- 绘制基元 ----------------
def add_ball(items, ax, eye, M, P, r_world, color, alpha, z_extra=0.0):
    """球：按深度排序后画成一串同心圆（自带光晕）"""
    if alpha <= 0.003:
        return
    xy, z = project([P], eye, M)
    xy = xy[0]
    # 世界半径换成像素：用相机距离与焦比
    fpx = (H / 2.0) / math.tan(math.radians(34.0) / 2.0)
    r_px = fpx * r_world / max(z[0], 0.1)
    items.append((float(z[0]) + z_extra, "ball", (xy, r_px, color, alpha)))


def draw_ball(ax, xy, r_px, color, alpha):
    if r_px < 0.4:
        return
    layers = 16
    for i in range(layers):
        f = 1.0 - i / float(layers)
        ax.add_patch(Circle(xy, r_px * f, facecolor=color, edgecolor="none",
                            alpha=alpha * 0.11, zorder=2))
    ax.add_patch(Circle(xy, r_px * 0.42, facecolor="white", edgecolor="none",
                        alpha=alpha * 0.55, zorder=3))
    # 光晕
    for i, halo in enumerate((1.9, 2.6)):
        ax.add_patch(Circle(xy, r_px * halo, facecolor=color, edgecolor="none",
                            alpha=alpha * 0.05 / (i + 1), zorder=1))


def add_line(items, ax, eye, M, A, B, color, alpha, lw=1.0, glow=2):
    if alpha <= 0.003:
        return
    xy, z = project([A, B], eye, M)
    d = float((z[0] + z[1]) / 2.0)
    items.append((d, "line", (xy, color, alpha, lw, glow)))


def draw_line(ax, xy, color, alpha, lw, glow):
    if glow:
        ax.plot(xy[:, 0], xy[:, 1], color=color, lw=lw * 3.2, alpha=alpha * 0.10,
                solid_capstyle="round", zorder=2)
    ax.plot(xy[:, 0], xy[:, 1], color=color, lw=lw, alpha=alpha,
            solid_capstyle="round", zorder=2)


def flush(fig, ax, items):
    items.sort(key=lambda it: it[0], reverse=True)   # 远的先画
    for _, kind, payload in items:
        if kind == "ball":
            draw_ball(ax, *payload)
        else:
            draw_line(ax, *payload)


# ---------------- 单帧 ----------------
def render_frame(f, frames, fps):
    t = f / float(fps)
    dur = frames / float(fps)

    fig = plt.figure(figsize=(W / 100.0, H / 100.0), dpi=100)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(BG)
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")

    az = math.radians(-58.0 + 26.0 * (t / max(dur, 1e-6)))
    el = math.radians(25.0)
    dist = 30.0
    tgt = np.array([0.0, 0.0, 3.6])
    eye = tgt + dist * np.array([math.cos(el) * math.cos(az),
                                 math.cos(el) * math.sin(az),
                                 math.sin(el)])
    M = look_at(eye, tgt)

    items = []

    # ① 舞台：长出来
    gs = ease(clamp01((t - 0.05) / 1.35))
    ga = 0.55 * ease(clamp01((t - 0.05) / 1.1))
    if gs > 0.005:
        half = HALF * (0.02 + 0.98 * gs)
        for k in range(NC + 1):
            v = -half + 2 * half * k / NC
            add_line(items, ax, eye, M, np.array([v, -half, 0]), np.array([v, half, 0]),
                     C_GRID, ga, 0.8, 1)
            add_line(items, ax, eye, M, np.array([-half, v, 0]), np.array([half, v, 0]),
                     C_GRID, ga, 0.8, 1)

    # ② seed
    ss = ease(clamp01((t - 1.15) / 0.95))
    if ss > 0:
        pulse = 1.0 + (0.16 * math.sin(2 * math.pi * (t - 1.2) / 0.9) if ss > 0.99 else 0.0)
        add_ball(items, ax, eye, M, np.array([0, 0, 0.35]), 0.30 * (0.05 + 0.95 * ss) * pulse,
                 C_SEED, ss)

    # ③ K₄：从 seed 位置上长
    ks = ease(clamp01((t - 2.2) / 3.0))
    if ks > 0:
        zlift = 0.35 + (TET_C[2] - 0.35) * ks
        scale = 0.04 + 0.96 * ks
        verts = (TET_V - TET_C) * scale + np.array([0, 0, zlift])
        va = ease(clamp01((t - 2.2) / 2.4))
        ea = ease(clamp01((t - 3.0) / 1.9))
        for i, j in TET_E:
            add_line(items, ax, eye, M, verts[i], verts[j], C_EDGE, ea, 2.6, 2)
        for p in verts:
            add_ball(items, ax, eye, M, p, 0.26 * max(scale, 0.25), C_VERT, max(va, 0.0))

    # ④ 三个频率槽
    sa = ease(clamp01((t - 5.2) / 1.8))
    if sa > 0:
        for i, p in enumerate(SLOT_P):
            k = (0.25 + 0.75 * sa) * (1.0 + 0.10 * math.sin(2 * math.pi * (t - 5.2) / 1.1 + i))
            add_ball(items, ax, eye, M, p, 0.34 * k, C_SLOT, sa)

    # ⑤ 三代
    pa = ease(clamp01((t - 7.0) / 3.5))
    if pa > 0:
        for i in range(3):
            cur = SLOT_P[i] + (PART_P[i] - SLOT_P[i]) * pa
            add_ball(items, ax, eye, M, cur, 0.16 + 0.30 * pa, C_PART,
                     ease(clamp01((t - 7.0) / 1.2)))

    flush(fig, ax, items)

    # ---------------- 文字 ----------------
    ax.text(48, 54, "Zero · 从 0 演化出粒子", color=C_TEXT, fontsize=21, alpha=0.95,
            ha="left", va="center")
    ax.text(48, 84, "结构示意（非数值刻度）", color=C_TEXT, fontsize=12, alpha=0.55,
            ha="left", va="center")

    cap = None
    if t < 1.2:
        cap = "①  空舞台"
    elif t < 2.2:
        cap = "②  seed 点亮：一个零点 / 一条权重"
    elif t < 5.2:
        cap = "③  K₄ 稳态长成：4 顶点 · 6 棱"
    elif t < 7.0:
        cap = "④  三个频率槽亮起（K₄ seed [8, 4, 1.5, 1] 的相邻比）"
    elif t < 10.6:
        cap = "⑤  从三个槽里各自长出一代"
    else:
        cap = "三槽 ↔ 三代：结构计数，不给质量数值"
    if cap:
        ax.text(48, H - 58, cap, color=C_TEXT, fontsize=17, alpha=0.9, ha="left", va="center")

    # 槽读数 / 代标签
    if t >= 5.6:
        a = ease(clamp01((t - 5.6) / 1.2))
        xy, _ = project(SLOT_P, eye, M)
        for i in range(3):
            d = np.array([xy[i, 0] - W / 2.0, xy[i, 1] - H / 2.0])
            d = d / (np.linalg.norm(d) + 1e-9)
            px, py = xy[i] + d * 52.0
            ax.text(px, py, SLOT_TXT[i], color=C_SLOT, fontsize=13,
                    alpha=0.85 * a, ha="center", va="center")
    if t >= 9.4:
        a = ease(clamp01((t - 9.4) / 1.0))
        xy, _ = project(PART_P, eye, M)
        for i in range(3):
            ax.text(xy[i, 0], xy[i, 1] - 52, GEN_TXT[i], color=C_PART, fontsize=14,
                    alpha=0.9 * a, ha="center", va="bottom")

    path = os.path.join(FRAMES_DIR, "f_%04d.png" % f)
    fig.savefig(path, facecolor=BG)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=240)
    ap.add_argument("--fps", type=float, default=20.0)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=-1)
    a = ap.parse_args()
    os.makedirs(FRAMES_DIR, exist_ok=True)
    end = a.frames if a.end < 0 else min(a.end, a.frames)
    import time
    t0 = time.time()
    for f in range(a.start, end):
        p = render_frame(f, a.frames, a.fps)
        if f % 20 == 0:
            print("frame %d/%d  %.1fs" % (f + 1, end, time.time() - t0), flush=True)
    print("done %d frames in %.1fs" % (end - a.start, time.time() - t0))


if __name__ == "__main__":
    main()
