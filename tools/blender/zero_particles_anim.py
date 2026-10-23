#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_particles_anim.py — 用 Blender 给「从 0 演化出粒子」做动画并渲染（GPU / EEVEE）

几何来源：SketchUp 建好导出的 OBJ（tools/sketchup/zero_particles_3d.rb），
分组名 stage_bar_* / seed / k4_v* / k4_e* / slot* / particle*。
动画在 Blender 里做（SketchUp 只负责建模）：逐帧改缩放/位移，frame_change 处理器算，
所以视口和最终渲染一致。

用法：
  Blender -b --factory-startup -P zero_particles_anim.py -- \
     --obj X.obj --outdir out/ [--frames 240] [--fps 20] [--res 1280x720]
     [--samples 16] [--blend-out X.blend] [--still 120]
"""
import bpy
import sys
import math
import os
import argparse
from mathutils import Vector

argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
ap = argparse.ArgumentParser()
ap.add_argument("--obj", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--frames", type=int, default=240)
ap.add_argument("--fps", type=float, default=20.0)
ap.add_argument("--res", default="1280x720")
ap.add_argument("--samples", type=int, default=16)
ap.add_argument("--blend-out", default="")
ap.add_argument("--still", type=int, default=-1, help="只渲染某一帧（调试用）")
a = ap.parse_args(argv)
W, H = [int(x) for x in a.res.lower().split("x")]
FPS = a.fps
os.makedirs(a.outdir, exist_ok=True)

SC = 0.0254  # OBJ 里是 mm，导入 global_scale=0.001 → m；原建模单位是 inch


def ease(x):
    x = 0.0 if x < 0 else (1.0 if x > 1 else x)
    return x * x * (3 - 2 * x)


# ---------------- 场景 / 导入 ----------------
bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
bpy.ops.wm.obj_import(filepath=a.obj, global_scale=0.001,
                      forward_axis="Y", up_axis="Z",
                      use_split_objects=True, use_split_groups=True)
objs = [o for o in sc.objects if o.type == "MESH"]
print("IMPORT objects=%d" % len(objs), flush=True)


def pick(prefix):
    return [o for o in objs if o.name.split(".")[0].startswith(prefix)]


STAGE = pick("stage_bar")
SEED = pick("seed")
K4V = pick("k4_v")
K4E = pick("k4_e")
SLOTS = sorted(pick("slot"), key=lambda o: o.name)
PARTS = sorted(pick("particle"), key=lambda o: o.name)
print("GROUPS stage=%d seed=%d k4v=%d k4e=%d slots=%d parts=%d"
      % (len(STAGE), len(SEED), len(K4V), len(K4E), len(SLOTS), len(PARTS)), flush=True)

# ---------------- 归一化（缩到半径 1，居中） ----------------
# 先把每个物体的原点挪到它自己的几何中心（OBJ 里顶点是绝对坐标、原点在 0），
# 否则缩放会绕世界原点走，槽/粒子会「飘」。
for o in objs:
    bb = [Vector(c) for c in o.bound_box]
    cen = sum(bb, Vector()) / 8.0
    for v in o.data.vertices:
        v.co -= cen
    o.location = Vector((cen.x * o.scale.x, cen.y * o.scale.y, cen.z * o.scale.z))

mn = Vector((1e9,) * 3)
mx = Vector((-1e9,) * 3)
for o in objs:
    for c in o.bound_box:
        w = o.matrix_world @ Vector(c)
        for i in range(3):
            mn[i] = min(mn[i], w[i])
            mx[i] = max(mx[i], w[i])
center = (mn + mx) / 2.0
radius = max((mx - mn)) / 2.0
S = 1.0 / radius
for o in objs:
    o.location = (o.location - center) * S
    o.scale = tuple(v * S for v in o.scale)
bpy.context.view_layer.update()
print("NORMALIZE radius=%.4f S=%.2f" % (radius, S), flush=True)


def norm_p(x, y, z):
    """建模坐标（inch）→ 归一化世界坐标"""
    p = Vector((x * SC, y * SC, z * SC))
    return (p - center) * S


# 关键位置（归一化后）
TET_C_W = Vector((0, 0, 3.4 * SC))
K4_PIVOT = (TET_C_W - center) * S
SEED_W = Vector((0, 0, 0.35 * SC))
SEED_POS = (SEED_W - center) * S
SLOT_R, SLOT_Z = 2.6, 7.4
SLOT_AZ = [math.pi / 2 + i * 2 * math.pi / 3 for i in range(3)]
SLOT_POS = [norm_p(SLOT_R * math.cos(t), SLOT_R * math.sin(t), SLOT_Z) for t in SLOT_AZ]
PART_POS = [norm_p(1.55 * SLOT_R * math.cos(t), 1.55 * SLOT_R * math.sin(t), SLOT_Z + 2.7)
            for t in SLOT_AZ]
R_VERT = 0.26 * SC * S
R_SLOT = 0.34 * SC * S
R_PART = 0.42 * SC * S

# ---------------- 材质 ----------------
for m in bpy.data.materials:
    if not m.use_nodes:
        continue
    bsdf = next((n for n in m.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
    if bsdf is None:
        continue
    bsdf.inputs["Roughness"].default_value = 0.35
    name = m.name.lower()
    try:
        bsdf.inputs["Emission Strength"].default_value = 0.0
    except Exception:
        pass
    if name.startswith("stage"):
        bsdf.inputs["Base Color"].default_value = (0.10, 0.14, 0.19, 1)
        bsdf.inputs["Roughness"].default_value = 0.6
    if name.startswith("seed"):
        bsdf.inputs["Emission Color"].default_value = (1.0, 0.84, 0.36, 1)
        bsdf.inputs["Emission Strength"].default_value = 1.8
    elif name.startswith("slot"):
        bsdf.inputs["Emission Color"].default_value = (0.20, 0.85, 0.75, 1)
        bsdf.inputs["Emission Strength"].default_value = 0.9
    elif name.startswith("particle"):
        bsdf.inputs["Emission Color"].default_value = (0.85, 0.92, 1.0, 1)
        bsdf.inputs["Emission Strength"].default_value = 1.8

# 球体平滑
for o in SEED + K4V + SLOTS + PARTS:
    for p in o.data.polygons:
        p.use_smooth = True

# ---------------- 父级枢轴（空物体）----------------
def make_empty(name, loc):
    e = bpy.data.objects.new(name, None)
    e.empty_display_size = 0.05
    e.location = loc
    sc.collection.objects.link(e)
    return e


def parent_to(child, empty):
    mw = child.matrix_world.copy()
    child.parent = empty
    child.matrix_parent_inverse = empty.matrix_world.inverted()
    child.matrix_world = mw


e_stage = make_empty("pivot_stage", Vector((0, 0, 0)))
for o in STAGE:
    parent_to(o, e_stage)
e_k4 = make_empty("pivot_k4", K4_PIVOT)
for o in K4V + K4E:
    parent_to(o, e_k4)

# 记录基准
BASE = {}
for o in SEED + SLOTS + PARTS:
    BASE[o.name] = (o.location.copy(), Vector(o.scale))
for o in [e_stage, e_k4]:
    BASE[o.name] = (o.location.copy(), Vector(o.scale))

# ---------------- 灯光 / 世界 / 相机 ----------------
r = 1.0
def area(name, loc, energy, size):
    d = bpy.data.lights.new(name, type="AREA")
    d.energy = energy
    d.size = size
    o = bpy.data.objects.new(name, d)
    o.location = loc
    o.rotation_euler = (Vector((0, 0, 0)) - Vector(loc)).to_track_quat("-Z", "Y").to_euler()
    sc.collection.objects.link(o)
    return o


area("key", Vector((1.2, -1.1, 1.3)), 80.0, 1.6)
area("fill", Vector((-1.4, -0.9, 0.5)), 22.0, 2.2)
area("rim", Vector((-0.4, 1.6, 1.0)), 50.0, 1.2)

world = bpy.data.worlds.new("w")
sc.world = world
world.use_nodes = True
bg = world.node_tree.nodes["Background"]
bg.inputs[0].default_value = (0.02, 0.025, 0.04, 1)
bg.inputs[1].default_value = 0.4

cam_d = bpy.data.cameras.new("cam")
cam_d.lens = 50
cam = bpy.data.objects.new("cam", cam_d)
sc.collection.objects.link(cam)
sc.camera = cam
# 用包围盒自动取景（场景是「竖着长」的，按半径装不下）
bb_mn = Vector((1e9,) * 3)
bb_mx = Vector((-1e9,) * 3)
for o in objs:
    for c in o.bound_box:
        w = (o.matrix_world @ Vector(c))
        for i in range(3):
            bb_mn[i] = min(bb_mn[i], w[i])
            bb_mx[i] = max(bb_mx[i], w[i])
CAM_TGT = (bb_mn + bb_mx) / 2.0 + Vector((0, 0, 0.02))
_size = bb_mx - bb_mn
fov_x = 2 * math.atan(0.5 * cam_d.sensor_width / cam_d.lens)
fov_y = 2 * math.atan(0.5 * cam_d.sensor_width * H / W / cam_d.lens)
DIST = max(_size.x / (2 * math.tan(fov_x / 2)), _size.z / (2 * math.tan(fov_y / 2))) * 1.02
print("FRAME bbox=%s size=%s dist=%.2f" % (tuple(round(v, 2) for v in CAM_TGT),
                                            tuple(round(v, 2) for v in _size), DIST), flush=True)


def set_cam(t, dur):
    az = math.radians(-58 + 26 * (t / max(dur, 1e-6)))
    el = math.radians(20)
    eye = CAM_TGT + DIST * Vector((math.cos(el) * math.cos(az),
                                   math.cos(el) * math.sin(az),
                                   math.sin(el)))
    cam.location = eye
    cam.rotation_euler = (CAM_TGT - eye).to_track_quat("-Z", "Y").to_euler()


# ---------------- 逐帧动画 ----------------
DUR = a.frames / FPS


def animate(frame):
    t = (frame - 1) / FPS

    # ① 舞台
    gs = ease((t - 0.05) / 1.35)
    e_stage.scale = Vector((max(gs, 0.002),) * 3)

    # ② seed
    ss = ease((t - 1.15) / 0.95)
    pulse = 1.0 + (0.16 * math.sin(2 * math.pi * (t - 1.2) / 0.9) if ss > 0.99 else 0.0)
    for o in SEED:
        base_loc, base_scale = BASE[o.name]
        k = (0.05 + 0.95 * ss) * pulse
        o.scale = base_scale * max(k, 0.001)
        o.location = base_loc

    # ③ K₄
    ks = ease((t - 2.2) / 3.0)
    for o in [e_k4]:
        base_loc, base_scale = BASE[o.name]
        o.scale = Vector((max(0.04 + 0.96 * ks, 0.004),) * 3)
        lift = SEED_POS.z + (K4_PIVOT.z - SEED_POS.z) * ks
        o.location = Vector((K4_PIVOT.x, K4_PIVOT.y, lift))

    # ④ 三槽
    sa = ease((t - 5.2) / 1.8)
    for i, o in enumerate(SLOTS):
        base_loc, base_scale = BASE[o.name]
        k = (0.25 + 0.75 * sa) * (1.0 + 0.10 * math.sin(2 * math.pi * (t - 5.2) / 1.1 + i))
        o.scale = base_scale * max(k, 0.001)
        o.location = base_loc

    # ⑤ 三代：从槽飞出去
    pa = ease((t - 7.0) / 3.5)
    for i, o in enumerate(PARTS):
        base_loc, base_scale = BASE[o.name]
        k = 0.12 / 0.42 + (1.0 - 0.12 / 0.42) * pa
        o.scale = base_scale * max(k, 0.001)
        o.location = SLOT_POS[i] + (PART_POS[i] - SLOT_POS[i]) * pa

    set_cam(t, DUR)


bpy.app.handlers.frame_change_pre.append(lambda s, dg: animate(s.frame_current))
sc.frame_start = 1
sc.frame_end = a.frames
sc.render.fps = int(round(FPS))

# ---------------- 渲染设置 ----------------
sc.render.engine = "BLENDER_EEVEE_NEXT"
sc.eevee.taa_render_samples = a.samples
for attr, val in (("use_raytracing", True), ("use_shadows", True)):
    if hasattr(sc.eevee, attr):
        try:
            setattr(sc.eevee, attr, val)
        except Exception:
            pass
try:
    sc.view_settings.view_transform = "Standard"
except Exception:
    pass
sc.render.resolution_x, sc.render.resolution_y = W, H
sc.render.image_settings.file_format = "PNG"
sc.render.filepath = os.path.join(a.outdir, "f_")

if a.blend_out:
    animate(1)
    bpy.ops.wm.save_as_mainfile(filepath=a.blend_out)
    print("SAVED_BLEND %s" % a.blend_out, flush=True)

import time
t0 = time.time()
if a.still >= 0:
    sc.frame_set(a.still)
    animate(a.still)
    sc.render.filepath = os.path.join(a.outdir, "still_%04d" % a.still)
    bpy.ops.render.render(write_still=True)
    print("STILL frame=%d %.2fs -> %s" % (a.still, time.time() - t0, sc.render.filepath), flush=True)
else:
    bpy.ops.render.render(animation=True)
    print("ANIM frames=%d time=%.1fs -> %s" % (a.frames, time.time() - t0, a.outdir), flush=True)
