#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_viewer_boot.py — Blender 打开 .blend 后自动切到「实时预览」视口

用法（由 tools/zero_gpu_render.sh --viewer 调用）：
  Blender <file.blend> --python zero_viewer_boot.py

做的事：把 3D 视口切成 Material Preview（EEVEE，走 GPU）、对齐到相机、
关掉大件叠加层、把渲染引擎设成 EEVEE、采样调到视口速度优先。不改几何、不退出。
"""
import bpy

LOG = "/tmp/zero_viewer_boot.log"
out = []


def log(m):
    out.append(str(m))


try:
    sc = bpy.context.scene
    sc.render.engine = "BLENDER_EEVEE_NEXT"
    try:
        sc.eevee.taa_render_samples = 64
    except Exception:
        pass

    w = bpy.context.window_manager.windows[0]
    area = next((a for a in w.screen.areas if a.type == "VIEW_3D"), None)
    if area is None:
        log("no VIEW_3D area")
    else:
        space = area.spaces.active
        space.shading.type = "MATERIAL"          # 实时预览 = EEVEE
        space.shading.use_scene_lights = True
        space.shading.use_scene_world = True
        try:
            space.overlay.show_overlays = False   # 干净点，像渲染
        except Exception:
            pass
        try:
            space.region_3d.view_perspective = "CAMERA"
        except Exception:
            pass
        log("shading=%s engine=%s res=%dx%d" % (space.shading.type, sc.render.engine,
                                                sc.render.resolution_x, sc.render.resolution_y))
        log("hint: 右键拖=转视角 / Z=切换着色 / F12=出终图")

    # 视口采样低一点更跟手（4.2+ 是 taa_render_samples 之外的 viewport 采样）
    try:
        sc.eevee.taa_samples = 8
    except Exception:
        pass
except Exception as e:
    log("ERR %s: %s" % (type(e).__name__, e))

open(LOG, "w").write("\n".join(out))
