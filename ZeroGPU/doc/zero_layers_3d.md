# 3D 图：用 SketchUp 的 Ruby 接口生成（实测跑通）

## 结论

本机能直接用 SketchUp 2025 Pro（25.0.659）的 **Ruby API** 出 3D 图，不需要手工建模：
脚本建几何 → 存 `.skp` → 同时用 SketchUp 自己的视图把它导出成 png 预览。

两张模型（同一套接口，各一份脚本）：

| 图 | 脚本 | 模型 | 预览 |
|:--|:--|:--|:--|
| 同心层（含 4／18／6 细分） | `tools/sketchup/zero_layers_3d.rb` | `ZeroGPU/results/zero_layers_3d.skp`（2.5 MB） | `zero_layers_3d_view1.png`／`view2.png` |
| 方格舞台 ＋ K₄ | `tools/sketchup/zero_stage_3d.rb` | `ZeroGPU/results/zero_stage_3d.skp`（2.1 MB） | `zero_stage_3d_view1.png`／`view2.png` |

内容为同心层结构 3D 版（外 → 内 `L0 → L1′ → L1 → L2`，最内 L2 一分为三＝3 个 120° 扇面），
层序与配色沿用 2D 版 `ZeroGPU/src/zero_layers_spheres.py`。四个壳各是一个 group（名称 `L0`／`L1'`／`L1`／`L2`），
L2 的三个扇面是独立 group（`L2_sector1..3`），材质带透明度；每个壳带经纬线（1 条子午线系＋5 条纬圈），
另有坐标轴和按各自方位角摆放的引线标签。

## 跑法

```bash
open -a "/Applications/SketchUp 2025/SketchUp.app" \
     --args -RubyStartup "/Users/oygb/Downloads/Zero/tools/sketchup/zero_layers_3d.rb"
```

在已打开的 SketchUp 里改完脚本想立刻重画：Ruby 控制台执行
`load "/Users/oygb/Downloads/Zero/tools/sketchup/zero_layers_3d.rb"`。

整轮耗时 ≈ 40 s（几乎全是 SketchUp 冷启动）；日志写到 `/tmp/zero_layers_3d.log`。

## 五个坑（都实际踩过）

1. **不能直接跑 `SketchUp.app/Contents/MacOS/SketchUp -RubyStartup …`**：进程起来了，脚本一行都不执行
   （探针文件不生成）。必须经 LaunchServices：`open -a "…/SketchUp.app" --args -RubyStartup …`。
   改脚本重跑前要先 `pkill -f "SketchUp 2025/SketchUp.app/Contents/MacOS/SketchUp"`，否则 `open` 只是把旧窗口提到前台，新脚本不跑。
2. **`add_faces_from_mesh(mesh, 4, mat)` 只给正面材质**：球面若法线朝内，看到的是 SketchUp 默认背面灰蓝，
   透明和配色全不显示。要逐个 `face.back_material = mat`。
3. **相机距离按半径手填会爆框**：设完 `camera.set(eye, target, up)` 必须再 `view.zoom_extents` 才把整个模型装进画面。
4. **`Entities` 没有 `add_cylinder`**（2025 版实测抛 `undefined method`）：圆柱/棱要用 `Geom::PolygonMesh` 自己搭侧面，
   或 `add_circle` 出边→`add_face`→`pushpull`。
5. **`ents.grep(Sketchup::Face).each { f.back_material = mat }` 会污染同组里先前建好的面**：K₄ 的棱改成红色后，
   四个顶点球一并变红。只对本次新增的面设背面材质（记下调用前的面数，取 `[before..-1]`）。

## 内容口径

- 同心层：层序 `L0 → L1′ → L1 → L2`（用户口径）；L2 三个扇面各自再分 4／18／6 份
  （4＝K₄ 的 4 个顶点、18＝3 代×6 边、6＝K₄ 的 6 条边，出处 `zero_tau_6_18_4_graph.py:17`）。
  扇面按 120° 等分只是画法，实际三槽不等分（Z161 §3.2 给 $[0.4055,0.6931,0.9808]$）。
- 舞台：方格 24×24 是具名输入（Z021／Z160），K₄ 是局部舞台（Z021 判定）。画法尺度（每格 0.5、四面体半径 2.2）
  不是物理长度。

## 可扩展方向（未做）

- L2 三个扇面各自的 4／18／6 亚亚层细分（2D 版已有口径，3D 里可做成分区鳍片或嵌套薄壳）。
- K₄ 四面体与方格舞台（Z021／Z160）作为独立 .skp。
- 导出 `.dae`/`.glb` 之类给网页或动画用；SketchUp 侧也可 `model.export` 走它自己的导出器。
