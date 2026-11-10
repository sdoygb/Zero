# Zb084 · 《零和宇宙 · 质量的几何》动画 · 复现说明

产物：`ZeroGPU/results/Zb084_mass_geometry_30s.mp4`（1920×1080 / 30.0 s / 24 fps / H.264 + AAC）

## 三件套

| 文件 | 作用 | 依赖 |
|---|---|---|
| `ZeroGPU/src/zb084_particles_geometry.py` | 质量三圆几何、粒子落点、标签图（Pillow 出图，24 张 PNG + 3×4 图集） | numpy / Pillow |
| `ZeroGPU/src/zb084_blender_scene.py` | Blender 场景：三圆、公共弦、粒子、场、十镜镜头、标签层、逐帧烘焙 | Blender 4.5.9 |
| `ZeroGPU/src/zb084_audio.py` | BGM / 音效（numpy 合成）+ 旁白（macOS `say` 男声 Eddy）+ ffmpeg 混音 | numpy / say / ffmpeg |

## 重跑

```bash
cd ~/Downloads/Zero

# 1) 几何 + 标签（纯 Python，可离线核验）
python3 ZeroGPU/src/zb084_particles_geometry.py
#    → results/zb084_particles_geometry.json（含 on_circle_check）
#    → results/zb084_labels/*.png

# 2) Blender 场景（逐帧烘焙关键帧后存 .blend）
/Applications/Blender.app/Contents/MacOS/Blender -b -P ZeroGPU/src/zb084_blender_scene.py

# 3) 渲染 720 帧（存盘后无回调依赖，可直接后台批量）
/Applications/Blender.app/Contents/MacOS/Blender -b ZeroGPU/results/zb084_mass_geometry.blend \
  -o //zb084_frames/f_#### -F PNG -a

# 4) 音轨
python3 ZeroGPU/src/zb084_audio.py

# 5) 合成
/usr/local/bin/ffmpeg -y -framerate 24 -i ZeroGPU/results/zb084_frames/f_%04d.png \
  -i ZeroGPU/results/zb084_audio/mix.wav -c:v libx264 -preset slow -crf 17 \
  -pix_fmt yuv420p -movflags +faststart -c:a aac -b:a 192k -shortest \
  ZeroGPU/results/Zb084_mass_geometry_30s.mp4
```

字幕开关：`ZB084_SUBTITLES=0` 时不建字幕层（用户原稿写「默认无字幕」，本片默认带字幕，按需关掉）。

## 画面里的数值取自哪

- 三圆半径 `R_k = c·sec θ_k*`（c=2，θ*=57.930484°/26.158813°/5.910703°）——Zb040 §4 命题 4.4；
  R/c = 1.883424 / 1.114112 / 1.005345，R-D = 0.574804 / 1.245888 / 1.803632（与篇内表逐位一致）。
- 粒子落点：每个粒子逐条核验 `|P − C_gen| = R_gen`（`on_circle_check`，10/10，Δ = 0）。
- `m = A²[1+√2cosθ]²`、`θ0 = 2π/3 + 2/9`、比值 `1 : 206.7703 : 3477.473`——Zb039 命题 5.1（已独立复算）。
- Koide `Q = Σm/(Σ√m)²`、`b=√2 ⇒ Q=2/3`——Zb041 定理 2.4 / 2.6。
  **注意**：用户原稿写的是 `(Σ√m)²/Σm = 2/3`，方向相反；按实测轻子该式给 1.5000，而库内定义给 0.6666605。
  片中已按库内定义上屏，字幕「三分之二」不变（Q=2/3 仍是锁定值）。
- 三代质量 0.511 / 105.66 / 1777 MeV 与 Zb039 比值互洽（偏差 1e-5 / 7e-5）。

## 已知取舍

- 粒子为「示意落点」：画法取自 `figures/mass_circle_three_generations.png`，不是从质量反解的精确位置；
  片中标「已着陆」的是三个轻子，正电子/夸克为背景元素（不标质量）。
- 场（B01 均匀蓝雾）是视觉元素，不是 Zb040 的几何对象；B02 偏心后倒平成背景薄雾。
- 地图集 `atlas.png` 仍在，是遗留产物；渲染实际逐张读 `*_solid.png`（实底版，防止文字后透出线条）。