#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_particles_obj.py — 生成「从 0 演化出粒子」场景的 OBJ＋MTL（与 zero_particles_3d.rb 同一套几何）

用途：当 SketchUp 因为授权/登录提示挡住 -RubyStartup、导不出 OBJ 时，用这份等价几何顶上，
让 Blender 侧（zero_particles_anim.py）完全不变。坐标按 SketchUp 的导出习惯乘 25.4（mm）。

分组名与 SketchUp 版一致：stage_bar_all / seed / k4_v0..3 / k4_e0..5 / slot1..3 / particle1..3

复跑：python3 make_particles_obj.py [out.obj]
"""
import sys
import os
import math
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, "results", "zero_particles.obj")
MM = 25.4  # 与 SketchUp 导出（mm）保持一致


def sphere(c, r, nu=16, nv=8):
    verts, faces = [], []
    verts.append((c[0], c[1], c[2] + r))
    for j in range(1, nv):
        th = math.pi * j / nv
        for i in range(nu):
            ph = 2 * math.pi * i / nu
            verts.append((c[0] + r * math.sin(th) * math.cos(ph),
                          c[1] + r * math.sin(th) * math.sin(ph),
                          c[2] + r * math.cos(th)))
    verts.append((c[0], c[1], c[2] - r))
    bottom = len(verts) - 1
    ring = lambda j, i: 1 + (j - 1) * nu + (i % nu)
    for i in range(nu):
        faces.append((0, ring(1, i), ring(1, i + 1)))
    for j in range(1, nv - 1):
        for i in range(nu):
            faces.append((ring(j, i), ring(j, i + 1), ring(j + 1, i + 1), ring(j + 1, i)))
    for i in range(nu):
        faces.append((bottom, ring(nv - 1, i + 1), ring(nv - 1, i)))
    return verts, faces


def tube(a, b, r, seg=10):
    a, b = np.array(a, float), np.array(b, float)
    n = b - a
    L = np.linalg.norm(n)
    n = n / L
    tmp = np.array([0, 0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0, 0])
    u = np.cross(n, tmp)
    u /= np.linalg.norm(u)
    w = np.cross(n, u)
    verts, faces = [], []
    for t in (a, b):
        for i in range(seg):
            ph = 2 * math.pi * i / seg
            p = t + r * (math.cos(ph) * u + math.sin(ph) * w)
            verts.append(tuple(p))
    for i in range(seg):
        i2 = (i + 1) % seg
        faces.append((i, i2, seg + i2, seg + i))
    return verts, faces


def bar(p0, p1, w, h):
    p0, p1 = np.array(p0, float), np.array(p1, float)
    d = p1 - p0
    L = np.linalg.norm(d)
    n = d / L
    side = np.cross(n, np.array([0, 0, 1.0]))
    if np.linalg.norm(side) < 1e-6:
        side = np.array([1.0, 0, 0])
    side /= np.linalg.norm(side)
    up = np.cross(side, n)
    up /= np.linalg.norm(up)
    verts = []
    for t in (p0, p1):
        for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            verts.append(tuple(t + side * (sx * w / 2) + up * (sy * h / 2)))
    faces = [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1), (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)]
    return verts, faces


HALF, NC = 6.0, 24
TET_C = np.array([0.0, 0.0, 3.4])
TET_S = 2.2
s = TET_S / math.sqrt(3.0)
VP = [TET_C + np.array(v, float) * s for v in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
SLOT_R, SLOT_Z = 2.6, 7.4

groups = []   # (name, material, verts, faces)

# 舞台
for k in range(NC + 1):
    v = -HALF + 2 * HALF * k / NC
    groups.append(("stage_bar_x%d" % k, "stage", *bar((v, -HALF, 0), (v, HALF, 0), 0.035, 0.035)))
    groups.append(("stage_bar_y%d" % k, "stage", *bar((-HALF, v, 0), (HALF, v, 0), 0.035, 0.035)))

groups.append(("seed", "seed", *sphere((0, 0, 0.35), 0.30)))
for i, p in enumerate(VP):
    groups.append(("k4_v%d" % i, "k4_vertex", *sphere(tuple(p), 0.26)))
for k, (i, j) in enumerate(EDGES):
    groups.append(("k4_e%d" % k, "k4_edge", *tube(VP[i], VP[j], 0.075)))

for i in range(3):
    a = math.pi / 2 + i * 2 * math.pi / 3
    groups.append(("slot%d" % (i + 1), "slot_%d" % i,
                   *sphere((SLOT_R * math.cos(a), SLOT_R * math.sin(a), SLOT_Z), 0.34)))
    groups.append(("particle%d" % (i + 1), "particle_%d" % i,
                   *sphere((1.55 * SLOT_R * math.cos(a), 1.55 * SLOT_R * math.sin(a), SLOT_Z + 2.7), 0.42)))

COLORS = {
    "stage": (0.27, 0.38, 0.47), "seed": (1.0, 0.84, 0.36),
    "k4_vertex": (0.34, 0.61, 0.84), "k4_edge": (0.84, 0.38, 0.29),
    "slot_0": (0.43, 0.86, 0.78), "slot_1": (0.43, 0.86, 0.78), "slot_2": (0.43, 0.86, 0.78),
    "particle_0": (0.92, 0.95, 1.0), "particle_1": (0.92, 0.95, 1.0), "particle_2": (0.92, 0.95, 1.0),
}

os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
obj = open(OUT, "w")
obj.write("# 从 0 演化出粒子 — 由 make_particles_obj.py 生成（与 zero_particles_3d.rb 同几何）\n")
obj.write("# File units = millimeters\nmtllib %s\n" % os.path.basename(OUT).replace(".obj", ".mtl"))
voff, noff = 1, 1
for name, mname, verts, faces in groups:
    obj.write("g %s\nusemtl %s\n" % (name, mname))
    for v in verts:
        obj.write("v %.5f %.5f %.5f\n" % (v[0] * MM, v[1] * MM, v[2] * MM))
    for f in faces:
        obj.write("f %s\n" % " ".join(str(voff + i) for i in f))
    voff += len(verts)
obj.close()

mtl = open(OUT.replace(".obj", ".mtl"), "w")
mtl.write("# materials for %s\n" % os.path.basename(OUT))
for mname, c in COLORS.items():
    mtl.write("newmtl %s\nKa 0 0 0\nKd %.4f %.4f %.4f\nKs 0.3 0.3 0.3\nd 1.0\n\n" % (mname, *c))
mtl.close()

print("wrote %s  groups=%d  verts=%d" % (OUT, len(groups), voff - 1))
