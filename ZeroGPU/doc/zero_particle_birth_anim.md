# 动画：从 0 演化出粒子（12 秒，12 步）

## 结论

做出来了，**12 秒、20 fps、1280×800** 的无界面渲染动画，纯 Python 逐帧出图再编码，
不依赖任何 GUI（SketchUp 那条路逐帧抓图不可靠，见下）。

- 视频：[zero_particle_birth.mp4](../../ZeroGPU/results/zero_particle_birth.mp4)（H.264，0.94 MB）
- 动图：[zero_particle_birth.gif](../../ZeroGPU/results/zero_particle_birth.gif)（720 宽，5.1 MB）
- 分镜：[zero_particle_birth_storyboard.png](../../ZeroGPU/results/zero_particle_birth_storyboard.png)（6 格）
- 脚本：`ZeroGPU/src/zero_particle_birth_anim.py`（帧目录 `ZeroGPU/results/anim_frames/`）

## 内容（结构示意，不做数值刻度）

| 时间 | 画面 | 口径 |
|:--|:--|:--|
| 0–1.4 s | 空舞台长出来（24×24 方格） | 方格舞台是具名输入（Z021／Z160） |
| 1.2–2.2 s | 中心一个 seed 点亮并脉动 | 一个零点／一条权重 |
| 2.2–5.2 s | K₄ 稳态从 seed 位置长成：4 顶点、6 棱 | Z129 稳态＝seed 形状 |
| 5.2–7.0 s | 环上三个频率槽亮起，标 `0.4055 / 0.6931 / 0.9808` | K₄ seed `[8, 4, 1.5, 1]` 的相邻比（Z131）；**三槽里只有 log2 属内部 2cos 族，另两条是输入**（Z132 §五） |
| 7.0–10.6 s | 从每个槽里长出一个球，标「第 1/2/3 代」 | 三槽 ↔ 三代是**结构计数**（Z161 §5） |
| 10.6–12 s | 定格 + 字幕「三槽 ↔ 三代：结构计数，不给质量数值」 | 不给质量数值 |

人物/物体大小、间距、发光都是画法，不是数值；三代球**不按质量缩放**。

## 复跑

```bash
cd ZeroGPU/src && /usr/bin/python3 zero_particle_birth_anim.py --frames 240 --fps 20   # ≈113 s
cd ZeroGPU/results && /usr/local/bin/ffmpeg -y -framerate 20 -i anim_frames/f_%04d.png \
  -c:v libx264 -pix_fmt yuv420p -crf 18 -movflags +faststart zero_particle_birth.mp4
```

渲染 240 帧耗时 113 s（≈0.47 s/帧，单核 matplotlib）。改 `--end`／`--start` 可只渲一段试镜头。

## 坑

1. **SketchUp 逐帧 `view.write_image` 会抓到同一帧**：无 UI 重绘时几何/相机更新不落进输出
   （实测 24 帧同尺寸同内容），所以动画没有走 SketchUp，而是自己把 3D 投到 2D 画。
   `.skp` 的静帧仍然走 SketchUp（见 `zero_layers_3d.md`）；`tools/sketchup/zero_particle_birth.rb` 是那条路的第一版，**未跑通，不作为产物**。
2. **SketchUp 反复 `pkill -9` 之后启动会卡住**（进程在、`-RubyStartup` 脚本不执行，日志一行都不写）。
   要重启就 `pkill -f`（SIGTERM），等 8 秒再 `open -a`。
3. 自己写投影时**球的大小要用 `fpx * r / z` 换算成像素**（不能直接拿世界半径当像素），否则近大远小会错。
4. 中文标签要显式 `addfont` 加载 Hiragino，否则 matplotlib 出方框。
