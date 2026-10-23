#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_render.py — 用 Blender 把 SketchUp 导出的 OBJ 渲染成图（GPU 光栅 / CPU 路径追踪）

用法（Blender 自己的命令行）：
  blender -b -P zero_render.py -- --obj X.obj --out Y.png [--engine EEVEE|CYCLES]
                                 [--samples 64] [--res 1600x1000] [--az 35] [--el 22]

引擎差别（本机实测，RX 570 / Metal 1.2）：
  EEVEE  → Metal 光栅化，用显卡
  CYCLES → 只有 CPU 设备（Blender 的 Cycles Metal 后端不支持这块 AMD 卡），走 CPU 路径追踪
"""
import bpy
import sys
import math
import time
import argparse
from mathutils import Vector

argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
ap = argparse.ArgumentParser()
ap.add_argument("--obj", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--engine", default="EEVEE", choices=["EEVEE", "CYCLES"])
ap.add_argument("--samples", type=int, default=64)
ap.add_argument("--res", default="1600x1000")
ap.add_argument("--az", type=float, default=35.0)
ap.add_argument("--el", type=float, default=22.0)
ap.add_argument("--bg", default="#0a0c12")
ap.add_argument("--glass", action="store_true", help="半透明材质当玻璃渲染")
ap.add_argument("--repeat", type=int, default=1, help="连续渲染多帧测稳态帧时间")
ap.add_argument("--no-rt", action="store_true", help="关掉 EEVEE 光追（更接近实时视口）")
ap.add_argument("--save-blend", default="", help="另存一份 .blend（可直接打开用实时视口看）")
a = ap.parse_args(argv)

W, H = [int(x) for x in a.res.lower().split("x")]


def hex2rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


# ---------- 清场 ----------
bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene

# ---------- 导入 OBJ ----------
t0 = time.time()
bpy.ops.wm.obj_import(filepath=a.obj, global_scale=0.001, forward_axis="Y", up_axis="Z")
meshes = [o for o in sc.objects if o.type == "MESH"]
print("IMPORT meshes=%d %.2fs" % (len(meshes), time.time() - t0), flush=True)

# ---------- 材质：OBJ 的 d（透明度）→ Principled alpha ----------
for m in bpy.data.materials:
    if not m.use_nodes:
        continue
    bsdf = next((n for n in m.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
    if bsdf is None:
        continue
    col = bsdf.inputs["Base Color"].default_value
    alpha = bsdf.inputs["Alpha"].default_value if "Alpha" in bsdf.inputs else 1.0
    # 半透明壳：给一点粗糙 + 透射感
    bsdf.inputs["Roughness"].default_value = 0.35
    if alpha < 0.999 and a.glass:
        # Blender 4.2+ 用 surface_render_method；DITHERED 对多层半透明排序最稳
        if hasattr(m, "surface_render_method"):
            m.surface_render_method = "DITHERED"
        elif hasattr(m, "blend_method"):
            m.blend_method = "BLEND"
        try:
            m.use_transparency_overlap = False
        except Exception:
            pass
        bsdf.inputs["Roughness"].default_value = 0.05
        bsdf.inputs["Alpha"].default_value = 1.0
        for key, val in (("Transmission Weight", 1.0), ("IOR", 1.45)):
            if key in bsdf.inputs:
                bsdf.inputs[key].default_value = val
    # 发光感：让高饱和色的自发光弱一点，暗场里更好看
    try:
        bsdf.inputs["Emission Color"].default_value = (col[0], col[1], col[2], 1.0)
        bsdf.inputs["Emission Strength"].default_value = 0.0
    except Exception:
        pass

for m in list(bpy.data.materials)[:8]:
    bs = next((n for n in m.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None) if m.use_nodes else None
    if bs:
        c = bs.inputs["Base Color"].default_value
        print("MAT %-14s rgba=(%.2f,%.2f,%.2f,%.2f) alpha=%.2f"
              % (m.name, c[0], c[1], c[2], c[3], bs.inputs["Alpha"].default_value), flush=True)
    else:
        print("MAT %-14s (no principled)" % m.name, flush=True)

# ---------- 包围盒 → 相机 ----------
mn = Vector((1e9, 1e9, 1e9))
mx = Vector((-1e9, -1e9, -1e9))
for o in meshes:
    for c in o.bound_box:
        w = o.matrix_world @ Vector(c)
        for i in range(3):
            mn[i] = min(mn[i], w[i])
            mx[i] = max(mx[i], w[i])
center = (mn + mx) / 2.0
radius = ((mx - mn).length) / 2.0
print("BBOX min=%s max=%s radius=%.3f" % (tuple(round(v, 2) for v in mn),
                                          tuple(round(v, 2) for v in mx), radius), flush=True)

# 归一化：把整体缩到半径 1（光照/相机就好算，且与模型原始尺寸无关）
if radius > 1e-9:
    scl = 1.0 / radius
    for o in meshes:
        o.location = (o.location - center) * scl
        o.scale = tuple(v * scl for v in o.scale)
    bpy.context.view_layer.update()
    radius = 1.0
    center = Vector((0, 0, 0))

cam_data = bpy.data.cameras.new("cam")
cam_data.lens = 50
cam = bpy.data.objects.new("cam", cam_data)
sc.collection.objects.link(cam)
az, el = math.radians(a.az), math.radians(a.el)
fov = 2.0 * math.atan(0.5 * cam_data.sensor_width / cam_data.lens)
dist = radius / math.tan(fov / 2.0) * 1.22
cam.location = center + Vector((math.cos(el) * math.cos(az),
                                math.cos(el) * math.sin(az),
                                math.sin(el))) * dist
direction = center - cam.location
cam.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
sc.camera = cam

# ---------- 灯光 ----------
def area(name, loc, energy, size, rot=None):
    d = bpy.data.lights.new(name, type="AREA")
    d.energy = energy
    d.size = size
    o = bpy.data.objects.new(name, d)
    o.location = loc
    if rot:
        o.rotation_euler = rot
    else:
        o.rotation_euler = (center - Vector(loc)).to_track_quat("-Z", "Y").to_euler()
    sc.collection.objects.link(o)
    return o

r = radius
area("key", center + Vector((1.2 * r, -1.1 * r, 1.3 * r)), 150.0, 1.6 * r)
area("fill", center + Vector((-1.4 * r, -0.9 * r, 0.5 * r)), 45.0, 2.2 * r)
area("rim", center + Vector((-0.4 * r, 1.6 * r, 1.0 * r)), 90.0, 1.2 * r)

world = bpy.data.worlds.new("w")
sc.world = world
world.use_nodes = True
bg = world.node_tree.nodes["Background"]
bg.inputs[0].default_value = (*hex2rgb(a.bg), 1.0)
bg.inputs[1].default_value = 0.35

# ---------- 渲染设置 ----------
sc.render.resolution_x, sc.render.resolution_y = W, H
# AgX/Filmic 会把饱和色洗淡，模型示意图用 Standard 更接近 SketchUp 里的颜色
try:
    sc.view_settings.view_transform = "Standard"
    sc.view_settings.look = "None"
except Exception:
    pass
sc.render.resolution_percentage = 100
sc.render.film_transparent = False
sc.render.image_settings.file_format = "PNG"
sc.render.filepath = a.out

if a.engine == "EEVEE":
    sc.render.engine = "BLENDER_EEVEE_NEXT"
    ee = sc.eevee
    ee.taa_render_samples = a.samples
    for attr, val in (("use_raytracing", True), ("use_shadows", True),
                      ("use_volumetric_shadows", True), ("shadow_ray_count", 2),
                      ("use_gtao", True)):
        if hasattr(ee, attr):
            try:
                setattr(ee, attr, val)
            except Exception:
                pass
else:
    sc.render.engine = "CYCLES"
    sc.cycles.device = "CPU"
    sc.cycles.samples = a.samples
    sc.cycles.use_denoising = True
    sc.cycles.max_bounces = 4

if a.save_blend:
    bpy.ops.wm.save_as_mainfile(filepath=a.save_blend)
    print("SAVED_BLEND %s" % a.save_blend, flush=True)

if a.no_rt and a.engine == "EEVEE":
    try:
        sc.eevee.use_raytracing = False
    except Exception:
        pass

for i in range(max(1, a.repeat)):
    t1 = time.time()
    bpy.ops.render.render(write_still=True)
    dt = time.time() - t1
    print("RENDER engine=%s samples=%d res=%dx%d frame=%d time=%.2fs -> %s"
          % (a.engine, a.samples, W, H, i + 1, dt, a.out), flush=True)
