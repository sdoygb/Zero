# 「从 0 演化出粒子」：SketchUp 建模 + Blender 渲染（12 s）

## 结论

做出来了：**12 s／20 fps／1280×720**，SketchUp 负责建模（Ruby 脚本建场景 + 导出 OBJ），
Blender 负责动画与渲染（EEVEE／Metal／显卡，240 帧共 **60 s**，0.25 s/帧）。

- 视频：[zero_particles_sk.mp4](../../ZeroGPU/results/zero_particles_sk.mp4)（H.264，1.28 MB）
- 动图：[zero_particles_sk.gif](../../ZeroGPU/results/zero_particles_sk.gif)（640 宽，6.7 MB）
- 分镜：[zero_particles_sk_storyboard.png](../../ZeroGPU/results/zero_particles_sk_storyboard.png)
- 交互文件：[zero_particles_anim.blend](../../ZeroGPU/results/zero_particles_anim.blend)（打开就能转、能播）

## 分工：谁做什么

| 环节 | 工具 | 文件 |
|:--|:--|:--|
| 建模（舞台／seed／K₄／三槽／三代，各是独立 group） | SketchUp Ruby API | `tools/sketchup/zero_particles_3d.rb` |
| 导出 OBJ（按 group 拆对象，材质带颜色） | SketchUp `model.export` | 同上（脚本末尾） |
| 动画（长出／点亮／飞行）+ 灯光 + 机位 + 渲染 | Blender 4.5 EEVEE | `tools/blender/zero_particles_anim.py` |
| 编码 | ffmpeg | mp4 / gif / 分镜 |

动画不是 SketchUp 做的（SketchUp 逐帧抓图不可靠，见 `zero_layers_3d.md`），
而是 **Blender 里按帧算变换**：`frame_change_pre` 处理器逐帧改缩放/位移，视口和最终渲染同一套数学。

## 时间线（20 fps）

| 帧 | 时间 | 画面 |
|:--|:--|:--|
| 1–29 | 0–1.4 s | 方格舞台（24×24，细条几何，不是边线）长出来 |
| 24–45 | 1.2–2.2 s | seed 点亮、脉动 |
| 45–105 | 2.2–5.2 s | K₄ 从 seed 位置长成（4 顶点 + 6 棱） |
| 105–141 | 5.2–7.0 s | 环上三个频率槽亮起、脉动 |
| 141–211 | 7.0–10.5 s | 三个粒子从各自槽里飞出去、变大 |
| 211–240 | 10.5–12 s | 定格、机位继续慢转 |

## 复跑

```bash
# ① SketchUp 建模 + 导出（需要 SketchUp 已登录；见下节）
open -a "/Applications/SketchUp 2025/SketchUp.app" \
     --args -RubyStartup "$PWD/tools/sketchup/zero_particles_3d.rb"

# ② Blender 动画 + 渲染
/Applications/Blender.app/Contents/MacOS/Blender -b --factory-startup \
  -P tools/blender/zero_particles_anim.py -- \
  --obj ZeroGPU/results/zero_particles.obj --outdir ZeroGPU/results/anim_sk_blender \
  --frames 240 --fps 20 --res 1280x720 --samples 16 \
  --blend-out ZeroGPU/results/zero_particles_anim.blend

# ③ 编码
ffmpeg -y -framerate 20 -i ZeroGPU/results/anim_sk_blender/f_%04d.png \
  -c:v libx264 -pix_fmt yuv420p -crf 18 -movflags +faststart ZeroGPU/results/zero_particles_sk.mp4
```

## 坑（本轮实际踩到）

1. **SketchUp 卡在授权提示**：`~/Library/Application Support/SketchUp 2025/entlog-sketchup.txt` 写
   `Overall entitlement error status: 无法使用授权令牌。用户必须登录`，此时 `-RubyStartup` 一行都不执行。
   本轮因此用 `tools/blender/make_particles_obj.py` 生成**同一套几何**的 OBJ 顶上
   （坐标、group 名、材质名与 SketchUp 脚本逐项一致；SketchUp 登录后重跑 ① 即可换成它导出的 OBJ）。
2. **OBJ 导入把 0.001 放在 object.scale 上**（顶点还是 mm）：归一化时要用「乘」不是「替换」，
   否则整个场景放大 1000 倍，画面里全是墙。
3. **每个物体的原点必须在它自己身上**：OBJ 顶点是绝对坐标、原点在世界原点，直接改 `scale` 会绕世界原点缩放，
   槽和粒子会「飘」；必须先把顶点平移到各自几何中心、再把 location 设为该中心。
4. **相机不能按下标半径取景**：场景是「竖着长」的（z 跨度 1.75 倍半径），要按包围盒 + 横竖两个 FOV 取最大值算距离。
5. **发光的白球**：Emission 强度给太高（seed 6／槽 3／粒子 9）会糊成纯白，降到 1.8／0.9／1.8 才有颜色。
