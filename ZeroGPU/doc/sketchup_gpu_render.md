# SketchUp 的 GPU 渲染：本机能到什么程度（实测）

## 结论（先给判定）

- **装了**：Blender 4.5.9 LTS 在 `/Applications/Blender.app`，走 **Metal** 后端，用**显卡**做 EEVEE 光栅渲染。
- **能用的那条链**：SketchUp 模型 → 导出 OBJ → Blender 渲染。一条命令：
  `tools/zero_gpu_render.sh <model.skp> [out.png] [EEVEE|CYCLES]`
- **GPU 光栅确实在卡上跑**：1600×1000 / 128 spp，K₄ 模型 **1.6 s**、同心层模型 **3.1 s**。
- **GPU 路径追踪（Cycles）在这台机器上没有**：Blender 的 Cycles 设备列表里**只有 CPU**（实测），
  所以 Cycles 只能 CPU 跑（同机同图 1600×1000 / 48 spp = **34 s**）。
- **SketchUp 里直接调显卡的插件**：厂商口径里能用的（Enscape／D5／Lumion 等）只出 Windows 版，
  V-Ray 的 GPU 模式要 CUDA（Windows／Linux）。**这两条我在这台 Mac 上没法验证**，只能按厂商文档说，
  所以本轮没有装"SketchUp 内的 GPU 插件"，装的是能真正吃到这块卡的 Blender 管线。

## 本机硬件事实（`system_profiler`）

| 项 | 值 |
|:--|:--|
| GPU | AMD Radeon RX 570，8 GB VRAM |
| Metal | **Metal 2**（不是 Metal 3） |
| 系统 | macOS 15.7.3（24G419） |
| CPU／内存 | Intel 2.40 GHz / 16 GB |

Blender 侧实测：`gpu.platform.renderer_get()` = `Metal API`，`vendor` = `AMD Radeon RX 570`；
`cycles.preferences.devices` 在 METAL 后端下只列出 CPU。

## 用法

```bash
# 一条命令（.skp → .obj → png）
tools/zero_gpu_render.sh ZeroGPU/results/zero_stage_3d.skp out.png

# 一条命令连实时预览一起：导 OBJ → 出 png → 存 .blend → 用 Blender 打开 EEVEE 视口
tools/zero_gpu_render.sh --viewer ZeroGPU/results/zero_stage_3d.skp

# 只存交互文件、不渲染
tools/zero_gpu_render.sh --viewer --no-render ZeroGPU/results/zero_layers_3d.skp
```

选项：`--engine EEVEE|CYCLES`、`--samples N`、`--res WxH`、`--az A --el E`、`--viewer`、`--no-render`。

# 只渲染已有的 OBJ，带参数
/Applications/Blender.app/Contents/MacOS/Blender -b --factory-startup \
  -P tools/blender/zero_render.py -- \
  --obj X.obj --out Y.png --engine EEVEE --samples 128 --res 1600x1000 --az 40 --el 18
```

渲染脚本支持的开关：`--engine EEVEE|CYCLES`、`--samples`、`--res WxH`、`--az/--el`（机位）、
`--bg`（背景色）、`--glass`（半透明材质当玻璃渲染）。

## 实时（边转边看）能做到什么

**SketchUp 里没有**：实时渲染插件（Enscape／D5／Lumion／Twinmotion）只出 Windows 版，V-Ray 的 GPU 要 CUDA；
SketchUp 自带的视口是实时的，但那是 OpenGL/Metal 着色显示，不是 PBR 渲染器。

**Blender 的视口是实时渲染器**（EEVEE／Metal，吃显卡）。本机实测：

| 场景 | 设置 | 结果 |
|:--|:--|:--|
| 视口重绘（Material Preview，Material 着色） | 2565×1470 | 首次 **26 fps**，持续 **117 fps** |
| 离线单帧 EEVEE（光追开） | 1280×720 / 8 spp | **0.25 s/帧**（≈4 fps） |
| 离线单帧 EEVEE（光追开） | 1600×1000 / 64 spp | 1.9 s |

`--viewer` 做的事：导 OBJ → 渲染 png → 存 `<模型>_viewer.blend` → 用 Blender 打开，
并跑 `tools/blender/zero_viewer_boot.py` 把视口切成 Material Preview、对齐相机、关叠加层
（日志写 `/tmp/zero_viewer_boot.log`，读到 `shading=MATERIAL` 即已就位）。
实测：`--viewer` 全流程 6 s 内拉起窗口；注意窗口可能开在别的 Space，切过去看。

现成的交互文件：`ZeroGPU/results/zero_layers_3d_viewer.blend`
（双击打开 → 按 `Z` 选 Material Preview 或 `Shift+Z` → 直接转视角／改灯看效果）。

视口截图（8 采样，半透明壳有噪点，属实时预览的正常形态）：
`ZeroGPU/results/render/viewport_realtime.png`

## 实测数字

| 图 | 引擎 | 分辨率／采样 | 时间 |
|:--|:--|:--|:--|
| K₄ ＋棱（`zero_stage_3d`） | EEVEE（GPU） | 1600×1000 / 128 | **1.6 s** |
| 同心层（`zero_layers_3d`） | EEVEE（GPU） | 1600×1000 / 128 | **3.1 s** |
| 同心层 | Cycles（CPU） | 1600×1000 / 48 | 34 s |

## 五个坑（都实际踩到）

1. **Homebrew 装不了 Blender**：`brew install --cask blender` 卡在升级 portable-ruby，
   `ghcr.io:443` SSL 失败。改为直接下 `https://download.blender.org/release/Blender4.5/blender-4.5.9-macos-x64.dmg`，
   `hdiutil attach` → `cp -R` 到 `/Applications`。`www.blender.org` 本机连不上，但 `download.blender.org` 通。
2. **SketchUp 的 `-RubyStartup` 要用 `open -a … --args` 才会执行**，而且反复 `pkill -9` 之后启动会卡死
   （进程在、脚本一行不写）；重启用 SIGTERM ＋ 等 8 s。
3. **SketchUp 导出的 OBJ 单位是毫米**：Blender 导入要 `global_scale=0.001`，否则相机裁切会算错一个数量级。
4. **OBJ 只导面、不导边**：方格舞台的网格线、坐标轴在 Blender 里会消失（要用线就得走 DAE 或在 Blender 里重建）。
5. **颜色对不上／糊成白球**：Blender 默认 AgX 视图变换会把饱和色洗淡 → 用 `Standard`；
   多层半透明壳必须用 `surface_render_method = 'DITHERED'`，用 `BLENDED` 会因为排序糊成一个白球。

## 相关文件

- 导出：`tools/sketchup/zero_export_obj.rb`
- 渲染：`tools/blender/zero_render.py`
- 串联：`tools/zero_gpu_render.sh`
- 产物：`ZeroGPU/results/render/stage_eevee_gpu.png`、`layers_eevee_gpu.png`、`layers_cycles_cpu.png`
