#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zb084 · 《零和宇宙 · 质量的几何》Blender 场景（30 s / 1920×1080 / 24 fps / EEVEE-Next）

前置（纯 Python，无 Blender 依赖）：
    python3 ZeroGPU/src/zb084_particles_geometry.py     几何 JSON + 叠加图集 PNG（见该脚本）

用法：
    blender -b -P ZeroGPU/src/zb084_blender_scene.py -- build    建场景 + 存盘
    blender -b ZeroGPU/results/zb084_mass_geometry.blend \
        -o //zb084_frames/f_#### -F PNG -x 1 -a                   渲帧（全片 720 帧）

设计
  · 单场景：三维几何（AgX）＋ 屏幕层 plane（相机父子、z=-10、COPY 目标=相机）→ 合成器叠加。
  · 屏幕层材质 = 图像 → 加光（叠加混合）；alpha 由图像自身给，淡入淡出由「缩放到 0」实现。
  · 时序唯一来源：zb084_particles_geometry.json（timeline_frames / render_plan）。
  · 排版离线可核：标签图集是纯 Python（Pillow）出的 PNG，Blender 端不碰 CJK 字体。
"""
import json
import math
from pathlib import Path

import bpy
import mathutils

ROOT = Path("/Users/oygb/Downloads/Zero")
RESULTS = ROOT / "ZeroGPU" / "results"
GEOM_JSON = RESULTS / "zb084_particles_geometry.json"
LABELS = RESULTS / "zb084_labels"
BLEND = RESULTS / "zb084_mass_geometry.blend"

FPS = 24
FRAME_END = 720
RES = (1920, 1080)

G = json.loads(GEOM_JSON.read_text())
TL = G["render_plan"]["timeline_frames"]
GEN = {r["gen"]: r for r in G["geometry"]}
PART = {p["key"]: p for p in G["particles"]}
GEN_RGB = {1: (0.36, 0.64, 1.00), 2: (0.62, 0.50, 1.00), 3: (0.36, 0.92, 0.90)}


def F(t):
    return int(round(t * FPS))


def clamp01(x):
    return max(0.0, min(1.0, x))


def ease(t):
    t = clamp01(t)
    return t * t * (3 - 2 * t)


# ──────────────────────────────────────────────────────────────
# 材质工具
# ──────────────────────────────────────────────────────────────
def emit_mat(name, color, strength):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    em = nt.nodes.new("ShaderNodeEmission")
    em.inputs["Color"].default_value = (*color, 1.0)
    em.inputs["Strength"].default_value = strength
    nt.links.new(em.outputs["Emission"], out.inputs["Surface"])
    return m


def link(ob, mat):
    ob.data.materials.append(mat)
    ob.visible_shadow = False
    return ob


# ──────────────────────────────────────────────────────────────
# 清场与渲染设定
# ──────────────────────────────────────────────────────────────
def reset():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    for coll in (bpy.data.images, bpy.data.materials, bpy.data.meshes, bpy.data.fonts):
        for blk in list(coll):
            try:
                coll.remove(blk)
            except Exception:
                pass
    sc = bpy.context.scene
    sc.name = "主"
    sc.render.fps = FPS
    sc.frame_start, sc.frame_end = 1, FRAME_END
    sc.render.resolution_x, sc.render.resolution_y = RES
    sc.render.resolution_percentage = 100
    sc.render.engine = "BLENDER_EEVEE_NEXT"
    sc.eevee.taa_render_samples = 32
    sc.view_settings.view_transform = "AgX"
    sc.view_settings.look = "AgX - Medium High Contrast"
    sc.render.use_compositing = True
    sc.render.use_sequencer = False


def build_world(sc):
    w = bpy.data.worlds.new("Space")
    w.use_nodes = True
    bg = w.node_tree.nodes["Background"]
    bg.inputs["Color"].default_value = (0.003, 0.008, 0.026, 1.0)
    bg.inputs["Strength"].default_value = 1.0
    sc.world = w


# ──────────────────────────────────────────────────────────────
# 几何：弦、三圆、粒子、零点、场、星域
# ──────────────────────────────────────────────────────────────
def curve_obj(name, pts, radius, color, strength):
    cu = bpy.data.curves.new(name, "CURVE")
    cu.dimensions = "3D"
    cu.bevel_depth = radius
    cu.bevel_resolution = 4
    cu.use_fill_caps = True
    cu.bevel_factor_mapping_start = "RESOLUTION"
    cu.bevel_factor_mapping_end = "RESOLUTION"
    sp = cu.splines.new("POLY")
    sp.points.add(len(pts) - 1)
    for i, p in enumerate(pts):
        sp.points[i].co = (p[0], p[1], p[2], 1.0)
    ob = bpy.data.objects.new(name, cu)
    bpy.context.collection.objects.link(ob)
    link(ob, emit_mat(name + "M", color, strength))
    return ob


def circle_pts(gen, n=360):
    R, D = GEN[gen]["R (c=2)"], GEN[gen]["D (c=2)"]
    th = math.radians(GEN[gen]["theta_star_deg"])
    uh = (math.cos(th), -math.sin(th), 0.0)
    pts = []
    for k in range(n + 1):
        phi = 2 * math.pi * k / n
        r = D + R * math.cos(phi)
        s = R * math.sin(phi)
        pts.append((r * uh[0], r * uh[1], s))
    return pts


def P(key):
    """(沿弦轴 x, 出平面 y=0, 弦方向 z)：粒子在弦与圆心所在的平面内。"""
    p = PART[key]
    return (p["pos_xz"][0], 0.0, p["pos_xz"][1])


def build_geometry():
    chord = curve_obj("弦", [(0, 0, -2.0), (0, 0, 2.0)], 0.014, (0.80, 0.88, 1.0), 4.0)
    circles = {g: curve_obj(f"圆{g}", circle_pts(g), 0.012, GEN_RGB[g], 5.0) for g in (1, 2, 3)}
    parts = {}
    for key, p in PART.items():
        rad = 0.038 if p["kind"] == "lepton" else 0.028
        bpy.ops.mesh.primitive_uv_sphere_add(segments=24, ring_count=14, radius=rad,
                                             location=P(key))
        ob = bpy.context.object
        ob.name = "P_" + key
        link(ob, emit_mat("M_" + key, tuple(c * 0.42 for c in p["rgb"]), 1.15 if p["kind"] == "lepton" else 0.8))
        parts[key] = ob
        # 外环：细环，颜色同粒子、亮一些（特写时像行星的光环）
        bpy.ops.mesh.primitive_torus_add(location=P(key), major_radius=rad * 1.75,
                                         minor_radius=rad * 0.035,
                                         major_segments=48, minor_segments=8,
                                         rotation=(math.radians(74), 0, 0))
        ring = bpy.context.object
        ring.name = "R_" + key
        ring.parent = ob
        ring.matrix_parent_inverse = ob.matrix_world.inverted()
        ring.location = ob.location
        link(ring, emit_mat("MR_" + key, p["rgb"], 1.2 if p["kind"] == "lepton" else 0.9))
        parts["ring_" + key] = ring
    return chord, circles, parts


def build_zero():
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=3, radius=0.055, location=(0, 0, 0))
    core = bpy.context.object
    core.name = "零点"
    link(core, emit_mat("M_zero", (1.0, 0.94, 0.80), 9.0))
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2, radius=0.12, location=(0, 0, 0))
    halo = bpy.context.object
    halo.name = "零点光晕"
    hz = emit_mat("M_halo", (1.0, 0.74, 0.42), 0.55)
    hz.blend_method = "BLEND"
    hz.show_transparent_back = True
    hz.use_backface_culling = True
    link(halo, hz)
    return core, halo


def build_field():
    """均匀场：一枚半径为 1 的圆盘。
    B01 它是球状均匀的（竖立 → 透明度低、边缘让 Facing 渐变尽量均匀）；
    B02 偏心后它落下、摊平成三圆所在平面（z 方向压扁），最终成为背景薄雾，不再遮挡几何。
    """
    bpy.ops.mesh.primitive_circle_add(vertices=128, radius=1.0, fill_type="NGON",
                                      location=(0, 0, 0), rotation=(math.pi / 2, 0, 0))
    haze = bpy.context.object
    haze.name = "场"
    m = bpy.data.materials.new("HazeM")
    m.use_nodes = True
    m.blend_method = "BLEND"
    m.show_transparent_back = True
    m.use_backface_culling = False
    nt = m.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    # 径向渐变：中心亮、边缘暗（用 Generated/UV 的半径）
    tc = nt.nodes.new("ShaderNodeTexCoord")
    ln = nt.nodes.new("ShaderNodeVectorMath")
    ln.operation = "LENGTH"
    nt.links.new(tc.outputs["Object"], ln.inputs[0])
    ramp = nt.nodes.new("ShaderNodeValToRGB")
    ramp.color_ramp.elements[0].position = 0.0
    ramp.color_ramp.elements[0].color = (0.16, 0.36, 0.86, 1.0)
    ramp.color_ramp.elements[1].position = 1.0            # 到边缘平滑归零，不出现地平线折痕
    ramp.color_ramp.elements[1].color = (0.02, 0.06, 0.22, 1.0)
    nt.links.new(ln.outputs["Value"], ramp.inputs["Fac"])
    em = nt.nodes.new("ShaderNodeEmission")
    em.inputs["Strength"].default_value = 0.20
    nt.links.new(ramp.outputs["Color"], em.inputs["Color"])
    tr = nt.nodes.new("ShaderNodeBsdfTransparent")
    mix = nt.nodes.new("ShaderNodeMixShader")
    v = nt.nodes.new("ShaderNodeValue")
    v.name = "alpha"
    v.outputs[0].default_value = 0.35
    nt.links.new(v.outputs[0], mix.inputs[0])
    nt.links.new(tr.outputs["BSDF"], mix.inputs[1])
    nt.links.new(em.outputs["Emission"], mix.inputs[2])
    nt.links.new(mix.outputs["Shader"], out.inputs["Surface"])
    haze.data.materials.append(m)
    haze.visible_shadow = False
    return haze


def build_stars():
    import random
    random.seed(84084)
    mat = emit_mat("StarsM", (0.80, 0.86, 1.0), 3.0)
    for i in range(380):
        u, v = random.random(), random.random()
        th, ph = 2 * math.pi * u, math.acos(2 * v - 1)
        r = random.choice([30, 34, 38])
        bpy.ops.mesh.primitive_ico_sphere_add(
            subdivisions=1, radius=random.uniform(0.04, 0.11),
            location=(r * math.sin(ph) * math.cos(th), r * math.sin(ph) * math.sin(th),
                      r * math.cos(ph)))
        ob = bpy.context.object
        ob.name = f"星{i}"
        link(ob, mat)


# ──────────────────────────────────────────────────────────────
# 相机
# ──────────────────────────────────────────────────────────────
def build_camera(sc):
    cd = bpy.data.cameras.new("Cam")
    cd.lens = 40
    cd.clip_start = 0.01
    cd.clip_end = 300
    cam = bpy.data.objects.new("Camera", cd)
    bpy.context.collection.objects.link(cam)
    sc.camera = cam
    tgt = bpy.data.objects.new("CamTarget", None)
    bpy.context.collection.objects.link(tgt)
    con = cam.constraints.new("TRACK_TO")
    con.target = tgt
    con.track_axis = "TRACK_NEGATIVE_Z"
    con.up_axis = "UP_Y"
    return cam, tgt


def camera_poses():
    e, mu, tau = P("e_p"), P("mu_p"), P("tau_p")
    zend = (0.0, 0.0, -2.0)
    N = math.radians

    def near(p, d, az, el):
        """在粒子 p 旁距离 d 取机位：az 绕 z 轴（屏幕的竖轴）、el 抬高。"""
        return (p[0] + d * math.cos(el) * math.cos(az),
                p[1] + d * math.cos(el) * math.sin(az),
                p[2] + d * math.sin(el))

    # 粒子邻域的空镜头侧：朝向一族的 +y（圆的外侧），避免穿进零点/弦
    poses = [
        (F(0.0), (0.0, -3.3, 0.30), (0.0, 0, 0.0), 42),          # B01 光点居中
        (F(2.9), (0.05, -3.1, 0.30), (0.0, 0, 0.0), 42),
        (F(4.4), (0.85, -2.7, 0.40), (0.12, 0, 0.0), 42),        # B02 微推近（偏心）
        (F(6.2), (1.10, -5.0, 0.90), (0.25, 0, 0.0), 38),        # B03 拉远，三圆入画
        (F(8.9), (0.60, -7.4, 1.30), (0.0, 0, 0.0), 34),
        (F(11.0), (-3.0, -6.2, 1.80), (0.0, 0, 0.0), 34),        # B04 环绕
        (F(13.0), (e[0] - 1.6, -4.2, 0.9), (e[0] - 0.2, 0, e[2] - 1.2), 46),   # B05
        (F(14.7), near(e, 1.55, N(-105), N(12)), e, 46),
        (F(16.5), near(e, 1.55, N(-105), N(12)), (e[0], 0, e[2] + 0.4), 46),   # B06
        (F(17.4), near(mu, 1.25, N(-100), N(8)), mu, 46),
        (F(18.7), near(tau, 0.70, N(-70), N(6)), tau, 46),
        (F(19.6), (1.35, -1.85, 0.35), (0.7, 0, -0.6), 48),
        (F(20.3), (0.75, -1.55, -1.45), (0.2, 0, -1.75), 48),               # B07 弦端点（零点）
        (F(21.8), (0.42, -1.20, -1.95), (0.0, 0, -1.98), 48),
        (F(23.6), (2.0, -2.9, 0.20), (0.3, 0, 0.6), 48),         # B08 微拉远，信息栏淡入
        (F(26.3), (2.0, -2.9, 0.20), (0.3, 0, 0.6), 48),         # B09 快拉远
        (F(28.8), (-1.8, -10.4, 2.0), (0.0, 0, 1.1), 31),
        (F(29.9), (-1.7, -10.2, 2.0), (0.0, 0, 1.1), 31),          # B10 定格
    ]
    return poses


def bake_camera(cam, tgt):
    for fr, loc, aim, lens in camera_poses():
        cam.location = loc
        cam.keyframe_insert("location", frame=fr)
        tgt.location = aim
        tgt.keyframe_insert("location", frame=fr)
        cam.data.lens = lens
        cam.data.keyframe_insert("lens", frame=fr)


# ──────────────────────────────────────────────────────────────
# 屏幕层：父=相机、z=-10、只跟随（不旋转），尺寸按当前焦距算屏高
# ──────────────────────────────────────────────────────────────
def screen_plane(name, img_path):
    me = bpy.data.meshes.new(name)
    me.from_pydata([(-.5, -.5, 0), (.5, -.5, 0), (.5, .5, 0), (-.5, .5, 0)], [], [(0, 1, 2, 3)])
    me.update()
    # 顶点顺序 (-,-) (+, -) (+, +) (-, +) ⇒ 与之对应的 UV (0,0)(1,0)(1,1)(0,1)
    uvl = me.uv_layers.new(name="UVMap")
    for j, c in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]):
        uvl.data[j].uv = c
    me.uv_layers.active = uvl
    ob = bpy.data.objects.new(name, me)
    bpy.context.collection.objects.link(ob)
    m = bpy.data.materials.new(name + "M")
    m.use_nodes = True
    m.blend_method = "BLEND"
    m.show_transparent_back = True
    nt = m.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    em = nt.nodes.new("ShaderNodeEmission")
    em.inputs["Strength"].default_value = 1.0
    img = nt.nodes.new("ShaderNodeTexImage")
    img.image = bpy.data.images.load(str(img_path))
    img.image.colorspace_settings.name = "sRGB"
    img.extension = "CLIP"
    img.interpolation = "Linear"
    try:
        img.image.pack()          # 内嵌，避免存盘后相对路径丢失
    except Exception:
        pass
    tc = nt.nodes.new("ShaderNodeTexCoord")
    nt.links.new(tc.outputs["UV"], img.inputs["Vector"])
    nt.links.new(img.outputs["Color"], em.inputs["Color"])
    tr = nt.nodes.new("ShaderNodeBsdfTransparent")
    mix = nt.nodes.new("ShaderNodeMixShader")
    nt.links.new(img.outputs["Alpha"], mix.inputs[0])
    nt.links.new(tr.outputs["BSDF"], mix.inputs[1])
    nt.links.new(em.outputs["Emission"], mix.inputs[2])
    # 混合比例：整体不透明度 × 图像 alpha
    nt.links.new(mix.outputs["Shader"], out.inputs["Surface"])
    ob.data.materials.append(m)
    ob.visible_shadow = False
    return ob, m


# ──────────────────────────────────────────────────────────────
# 帧回调（唯一时序源）
# ──────────────────────────────────────────────────────────────
class Stage:
    rendering = False


ST = Stage()


def elapsed(objects):
    """进度时钟（秒）：由 Camera 的自定义属性 time 驱动，保证 .blend 存盘后仍能求值。"""
    return objects["Camera"]["time"]


def add_clock(cam, sc):
    cam["time"] = 0.0
    cam.id_properties_ui("time").update(min=0.0)
    fc = cam.driver_add("location", 0)          # 占位：真正的时钟驱动见下
    # time 由帧驱动
    cam.driver_remove("location", 0)
    d = cam.driver_add('["time"]')
    drv = d.driver
    drv.type = "SCRIPTED"
    v = drv.variables.new()
    v.name = "f"
    v.type = "SINGLE_PROP"
    v.targets[0].id_type = "SCENE"
    v.targets[0].id = sc
    v.targets[0].data_path = "frame_current"
    drv.expression = "f / 24.0"


def seek(cam, ob, prop, expr, owner=None):
    """给 ob 的自定义属性 prop 加驱动；表达式里可用 clock（= Camera.time）与 owner 的自定义属性。"""
    owner = owner or ob
    if prop not in ob.keys():
        ob[prop] = 0.0
    fc = ob.driver_add('["%s"]' % prop)
    drv = fc.driver
    drv.type = "SCRIPTED"
    for nm in [k for k in owner.keys() if not k.startswith("_")]:
        v = drv.variables.new()
        v.name = nm
        v.type = "SINGLE_PROP"
        v.targets[0].id_type = "OBJECT"
        v.targets[0].id = owner
        v.targets[0].data_path = '["%s"]' % nm
    v = drv.variables.new()
    v.name = "clock"
    v.type = "SINGLE_PROP"
    v.targets[0].id_type = "OBJECT"
    v.targets[0].id = cam
    v.targets[0].data_path = '["time"]'
    drv.expression = expr


# 内置小函数（驱动表达式可调用）
DRIVER_NAMESPACE = {
    "smoothstep": lambda a, b, x: ease((x - a) / (b - a)) if b != a else 1.0,
    "clamp01": clamp01,
    "ease": ease,
}


def init_namespace():
    for f in DRIVER_NAMESPACE:
        if f not in bpy.app.driver_namespace:
            bpy.app.driver_namespace[f] = DRIVER_NAMESPACE[f]


def ss(a, b):
    """帧 → 驱动表达式里的 smoothstep 参数（返回 0..1 的缓动值）。"""
    return "%.6f, %.6f, clock" % (a / FPS, b / FPS)


def clock_at(fr):
    return fr / FPS


def fire(ob, road):
    """road: [(帧, x, y, z, 焦距或None)] —— 直接在存盘前逐帧求值并打关键帧（LINEAR）。"""
    cam = bpy.data.objects["Camera"]
    for (fr, x, y, z, lens) in road:
        ob.location = (x, y, z)
        ob.keyframe_insert("location", frame=fr)
        if lens is not None and ob.type == "CAMERA":
            ob.data.lens = lens
            ob.data.keyframe_insert("lens", frame=fr)


def fire_target(tgt, road):
    for (fr, x, y, z, lens) in road:
        tgt.location = (x, y, z)
        tgt.keyframe_insert("location", frame=fr)


def set_linear(ob):
    ad = ob.animation_data
    if ad and ad.action:
        for fc in ad.action.fcurves:
            for kp in fc.keyframe_points:
                kp.interpolation = "LINEAR"
    return ob


# ──────────────────────────────────────────────────────────────
# 烘焙：每帧用纯 Python 算好状态，逐帧写关键帧（.blend 存盘后无需任何回调/驱动）
# ──────────────────────────────────────────────────────────────
def layer_plays():
    """屏幕层时间表：(入帧, 出帧, 入淡秒, 出淡秒, x, y, 占宽, 占高)；x,y ∈ [-1,1] 相对一整屏。

    排布原则：字幕占底部条带；公式占固定侧栏；粒子/零点标签只在“无人占位”的时段出现，
    同一时刻最多 2 个信息框，避免互相压住。
    """
    F_ = F
    P = {}
    # 字幕（旁白逐句）
    subs = [
        ("sub0", 0.40, 3.30), ("sub1", 3.30, 4.90), ("sub2", 4.90, 6.60),
        ("sub3", 6.60, 8.10), ("sub4", 8.10, 9.40), ("sub5", 9.40, 12.0),
        ("sub6", 13.6, 16.3), ("sub7", 16.3, 20.9), ("sub8", 20.9, 23.6),
        ("sub9", 23.6, 26.6), ("sub10", 26.6, 29.0),
    ]
    for name, a, b in subs:
        P[name] = (F_(a), F_(b), 0.28, 0.24, 0.0, -0.66, 1640.0, 1.0)
    # 统一按 1920×1080 像素坐标排布：x ∈ [120, 1800]，y ∈ [140, 960]（y 向下）
    def at(cx_px, cy_px, wpx):
        return (2.0 * (cx_px / 1920.0) - 1.0, 1.0 - 2.0 * (cy_px / 1080.0), wpx)

    # 公式 1（质量产生）：B02–B03 右上
    x, y, w = at(1400, 250, 960)
    P["公式1"] = (F_(3.5), F_(8.9), 0.6, 0.4, x, y, w, 1.0, w)
    # 公式 2（三代质量比）：μ、τ 两拍之间，左下
    x, y, w = at(520, 880, 720)
    P["公式2"] = (F_(17.0), F_(18.3), 0.4, 0.3, x, y, w, 1.0, w)
    # 公式 3（Koide）：B08 右上
    x, y, w = at(1400, 250, 960)
    P["公式3"] = (F_(23.2), F_(28.6), 0.7, 0.5, x, y, w, 1.0, w)
    # 粒子标签（着陆拍弹出，左下，避开字幕）
    x, y, w = at(300, 560, 470)
    P["标签_e"] = (F_(13.9), F_(16.3), 0.35, 0.28, x, y, w, 1.0, w)
    x, y, w = at(300, 560, 500)
    P["标签_mu"] = (F_(16.9), F_(18.5), 0.30, 0.25, x, y, w, 1.0, w)
    x, y, w = at(1560, 560, 470)
    P["标签_tau"] = (F_(18.3), F_(20.6), 0.30, 0.25, x, y, w, 1.0, w)
    # 零点
    x, y, w = at(280, 760, 350)
    P["标签_g"] = (F_(6.3), F_(9.6), 0.4, 0.35, x, y, w, 1.0, w)
    x, y, w = at(280, 560, 350)
    P["标签_g2"] = (F_(21.0), F_(26.4), 0.3, 0.3, x, y, w, 1.0, w)
    x, y, w = at(1560, 560, 480)
    P["标签_nu"] = (F_(22.2), F_(26.4), 0.3, 0.3, x, y, w, 1.0, w)
    # 全景回收：三代轻子标签（B09 拉远后）
    x, y, w = at(300, 460, 420)
    P["标签_e2"] = (F_(28.7), 9999, 0.4, 0.3, x, y, w, 1.0, w)
    x, y, w = at(1600, 460, 440)
    P["标签_mu2"] = (F_(28.75), 9999, 0.4, 0.3, x, y, w, 1.0, w)
    x, y, w = at(960, 250, 420)
    P["标签_tau2"] = (F_(28.9), 9999, 0.4, 0.3, x, y, w, 1.0, w)
    # 标题：B10 居中
    P["标题"] = (F_(29.0), 9999, 0.6, 0.4, 0.0, -0.80, 1500.0, 1.0)
    return P


def build_screen_layers():
    """建标签 plane（每格一张 PNG），可见性与位置由 bake_screen_layers 逐帧写关键帧。"""
    import os
    subtitles_on = os.environ.get("ZB084_SUBTITLES", "1") != "0"
    plans = layer_plays()
    if not subtitles_on:
        plans = {k: v for k, v in plans.items() if not k.startswith("sub")}
    specs = {}
    for name, t in plans.items():
        solid = LABELS / (name + "_solid.png")
        path = solid if solid.exists() else (LABELS / (name + ".png"))
        if not path.exists():
            raise FileNotFoundError("标签 PNG 缺失: %s（先跑 zb084_particles_geometry.py）" % path)
        ob, _m = screen_plane("L_" + name, path)
        ob.hide_render = True
        img = bpy.data.images[path.name]
        specs[name] = {"obj": ob, "t": t, "px": (img.size[0], img.size[1])}
    return specs


def bake_animation(cam, tgt):
    """一帧一帧求值 → 关键帧（LINEAR）。相机路径另外做平滑（见 bake_camera_smooth）。"""
    sc = bpy.context.scene
    T = TL
    land = {"e_p": T["land_e"], "mu_p": T["land_mu"], "tau_p": T["land_tau"]}
    base_loc = {k: tuple(ob.location) for k, ob in ST.parts.items()}
    for f in range(1, FRAME_END + 1):
        t = f / FPS
        # 零点
        if f < T["ecc_peak"]:
            amp = 0.012 + 0.048 * ease(f / (2.9 * FPS))
        elif f < T["field_settle"]:
            amp = 0.060
        else:
            amp = 0.0
        on = f < T["circles_in"]
        ST.core.location = (amp * math.sin(t * 2 * math.pi * 5.0),
                            amp * 0.3 * math.sin(t * 2 * math.pi * 3.6),
                            amp * math.sin(t * 2 * math.pi * 2.7))
        ST.core.scale = (1, 1, 1)
        ST.core.keyframe_insert("location", frame=f)
        ST.core.keyframe_insert("scale", frame=f)
        ST.halo.scale = (1.0 + 1.5 * ease(f / T["ecc_peak"]),) * 3
        ST.halo.keyframe_insert("scale", frame=f)
        ST.core.hide_render = not on
        ST.halo.hide_render = not on
        ST.core.keyframe_insert("hide_render", frame=f)
        ST.halo.keyframe_insert("hide_render", frame=f)
        ST.haze.keyframe_insert("hide_render", frame=f)

        # 场
        if f < T["ecc_peak"]:
            k = ease(f / T["ecc_peak"])
            ST.haze.scale = (3.0 + 1.6 * k, 3.0 + 1.6 * k, 1.0)
            ST.haze.location = (0.0, 6.0, 0.0)
            ST.haze.rotation_euler = (math.pi / 2, 0, 0)
            alpha = 0.24
        else:
            k = ease((f - T["ecc_peak"]) / (T["field_settle"] - T["ecc_peak"]))
            s_ = 4.4 + 9.0 * k
            ST.haze.scale = (s_, s_, 1.0)
            ST.haze.location = (0.6 * k, 6.0 + 3.0 * k, 0.4 * k)
            ST.haze.rotation_euler = (math.pi / 2, 0, 0)
            alpha = 0.24 - 0.06 * k
        ST.haze.keyframe_insert("scale", frame=f)
        ST.haze.keyframe_insert("location", frame=f)
        ST.haze.keyframe_insert("rotation_euler", frame=f)
        ST.haze_mat.node_tree.nodes["alpha"].outputs[0].default_value = alpha
        ST.haze_mat.node_tree.nodes["alpha"].outputs[0].keyframe_insert("default_value", frame=f)

        # 三圆
        for g, ob in ST.circles.items():
            t0 = T["circles_in"] + (g - 1) * 5
            k = ease((f - t0) / (T["circles_full"] - T["circles_in"] - 5))
            grow = ease((f - T["circles_in"]) / 53)
            ob.data.bevel_factor_end = max(0.001, grow)
            ob.scale = (max(0.02, k), 1, 1)
            ob.hide_render = k <= 0.002
            ob.keyframe_insert("scale", frame=f)
            ob.keyframe_insert("hide_render", frame=f)
            ob.data.keyframe_insert("bevel_factor_end", frame=f)

        # 弦
        ck = ease(f / (2.2 * FPS))
        ST.chord.data.bevel_factor_end = max(0.001, ck)
        ST.chord.data.keyframe_insert("bevel_factor_end", frame=f)
        ST.chord.hide_render = ck <= 0.001
        ST.chord.keyframe_insert("hide_render", frame=f)

        # 粒子
        for key, ob in ST.parts.items():
            if key.startswith("ring_"):
                continue
            p = PART[key]
            base = 1.0 if p["kind"] == "lepton" else 0.82
            if f < T["circles_full"] - 5:
                s_ = 0.001
            else:
                s_ = base
                dt = (f - land.get(key, -9999)) / FPS
                if -0.35 < dt < 0.9:
                    s_ *= 1.0 + 0.5 * math.exp(-((dt - 0.06) ** 2) / (2 * 0.16 ** 2))
            ob.scale = (s_, s_, s_)
            ob.hide_render = s_ < 0.002
            ob.keyframe_insert("scale", frame=f)
            ob.keyframe_insert("hide_render", frame=f)




def bake_labels(cam, f):
    return


def bake_screen_layers(cam):
    """标签：位置/缩放按相机逐帧算好写关键帧。"""
    sc = bpy.context.scene
    for f in range(1, FRAME_END + 1):
        sc.frame_set(f)
        R = cam.matrix_world.to_3x3().normalized()
        right, up, fwd = R.col[0], R.col[1], -R.col[2]
        D = 0.30                              # 标签平面距相机 0.3 m（比所有几何都近 ⇒ 等于屏幕层）
        H = D * 21.6 / cam.data.lens
        W = H * RES[0] / RES[1]
        base = cam.matrix_world.to_translation()
        for name, spec in ST.planes.items():
            ob = spec["obj"]
            a, b, fi, fo, x, y, sx, sy, maxw = (tuple(spec["t"]) + (0,))[:9]
            vis = ease((f - (a - fi)) / fi) if f < a else 1.0
            if f > b - fo:
                vis = min(vis, ease((b - f) / fo))
            if vis <= 0.003:
                ob.hide_render = True
                ob.keyframe_insert("hide_render", frame=f)
                continue
            ob.hide_render = False
            ob.keyframe_insert("hide_render", frame=f)
            ob.location = base + right * (x * W / 2) + up * (y * H / 2) + fwd * D
            ob.rotation_euler = cam.matrix_world.to_euler()
            px_w, px_h = spec["px"]
            if sx >= 100.0:
                # sx 是“目标像素宽”（≥100 视为像素宽）：字幕/公式/标签都用这一支
                pw = sx if not maxw else min(sx, maxw)
                scw = (pw / 1920.0) * W * vis
            else:
                # sx 是相对原图像素的放大倍数（< 100）
                pw = px_w * sx if not maxw else min(px_w * sx, maxw)
                scw = (pw / 1920.0) * W * vis
            sch = scw * (px_h / px_w)
            ob.scale = (scw, sch, 1.0)
            ob.keyframe_insert("location", frame=f)
            ob.keyframe_insert("rotation_euler", frame=f)
            ob.keyframe_insert("scale", frame=f)


def smooth_camera_keys(cam, tgt):
    """相机关键帧：中间镜头用平滑插值，抹掉 Bezier 回冲。"""
    for ob in (cam, tgt):
        ad = ob.animation_data
        if ad and ad.action:
            for fc in ad.action.fcurves:
                for kp in fc.keyframe_points:
                    kp.interpolation = "BEZIER"
                    kp.handle_left_type = "AUTO_CLAMPED"
                    kp.handle_right_type = "AUTO_CLAMPED"
    cd = cam.data
    if cd.animation_data and cd.animation_data.action:
        for fc in cd.animation_data.action.fcurves:
            for kp in fc.keyframe_points:
                kp.interpolation = "BEZIER"
                kp.handle_left_type = "AUTO_CLAMPED"
                kp.handle_right_type = "AUTO_CLAMPED"


def setup_compositor(sc):
    sc.use_nodes = True
    nt = sc.node_tree
    nt.nodes.clear()
    rl = nt.nodes.new("CompositorNodeRLayers")
    gl = nt.nodes.new("CompositorNodeGlare")
    gl.glare_type = "FOG_GLOW"
    gl.quality = "HIGH"
    gl.threshold = 1.2
    gl.size = 7
    gl.mix = -0.35
    comp = nt.nodes.new("CompositorNodeComposite")
    nt.links.new(rl.outputs["Image"], gl.inputs["Image"])
    nt.links.new(gl.outputs["Image"], comp.inputs["Image"])


def main_build():
    reset()
    sc = bpy.context.scene
    build_world(sc)
    build_stars()
    ST.chord, ST.circles, ST.parts = build_geometry()
    ST.core, ST.halo = build_zero()
    ST.haze = build_field()
    ST.haze_mat = ST.haze.data.materials[0]
    cam, tgt = build_camera(sc)
    bake_camera(cam, tgt)
    add_clock(cam, sc)
    setup_compositor(sc)
    ST.planes = build_screen_layers()
    bake_animation(cam, tgt)
    smooth_camera_keys(cam, tgt)
    bake_screen_layers(cam)
    sc.frame_set(1)
    BLEND.parent.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND))
    print("SAVED", BLEND)


if __name__ == "__main__":
    main_build()
