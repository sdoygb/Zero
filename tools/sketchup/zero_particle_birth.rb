# encoding: UTF-8
# ⚠ 未跑通，不作为产物：SketchUp 逐帧 view.write_image 会抓到同一帧（无 UI 重绘时）。
#   动画的实际产物走 ZeroGPU/src/zero_particle_birth_anim.py（纯 Python 渲染）。
#   本文件保留作参照（几何/时间线的第一版写法）。
# zero_particle_birth.rb — 「从 0 演化出粒子」的 3D 动画（SketchUp Ruby 接口逐帧出图）
#
# 内容（结构示意，不给数值主张）：
#   ① 空舞台（24×24 方格，长出来）→ ② seed 点亮 → ③ K₄ 稳态在 seed 上长成
#   → ④ 三个频率槽亮起（[0.4055, 0.6931, 0.9808]，K₄ seed [8,4,1.5,1] 的相邻比）
#   → ⑤ 从三个槽里各自长出一个「代」（第 1/2/3 代）
# 口径提醒：三槽 ↔ 三代是结构计数，不是质量数值（Z161 §5）；三槽中只有 log2 属内部
# 2cos 族，另两条是输入（Z132 §五）。动画里不做数值刻度，只做结构。
#
# 跑法： open -a "/Applications/SketchUp 2025/SketchUp.app" --args -RubyStartup <本文件>
# 环境变量：ZFRAMES（渲染帧数，默认 240）／ZOUT（输出目录）／ZW、ZH（画面尺寸）

OUT  = ENV['ZOUT'] || "/Users/oygb/Downloads/Zero/ZeroGPU/results/anim_frames"
LOG  = "/tmp/zero_particle_birth.log"
FRAMES = (ENV['ZFRAMES'] || 240).to_i
FPS    = 20.0
DUR    = FRAMES / FPS
WIDTH  = (ENV['ZW'] || 1280).to_i
HEIGHT = (ENV['ZH'] || 800).to_i

def log(m)
  File.open(LOG, "a") { |f| f.puts m }
end
File.delete(LOG) if File.exist?(LOG)
Dir.mkdir(OUT) unless Dir.exist?(OUT)
log "start frames=#{FRAMES} dur=#{DUR} #{WIDTH}x#{HEIGHT}"

model = Sketchup.active_model
model.entities.clear!
model.start_operation("particle birth", true)

def mat(model, name, rgb, alpha)
  m = model.materials[name] || model.materials.add(name)
  m.color = Sketchup::Color.new(*rgb)
  m.alpha = alpha
  m
end

def uv_sphere(ents, c, r, m, nu = 20, nv = 10)
  mesh = Geom::PolygonMesh.new
  top = mesh.add_point(Geom::Point3d.new(c.x, c.y, c.z + r))
  bot = mesh.add_point(Geom::Point3d.new(c.x, c.y, c.z - r))
  ring = []
  (1...nv).each do |j|
    th = Math::PI * j / nv
    row = []
    nu.times do |i|
      ph = 2 * Math::PI * i / nu
      row << mesh.add_point(Geom::Point3d.new(
        c.x + r * Math.sin(th) * Math.cos(ph),
        c.y + r * Math.sin(th) * Math.sin(ph),
        c.z + r * Math.cos(th)))
    end
    ring << row
  end
  nu.times { |i| mesh.add_polygon(top, ring[0][i], ring[0][(i + 1) % nu]) }
  (0...(ring.length - 1)).each do |j|
    nu.times do |i|
      i2 = (i + 1) % nu
      mesh.add_polygon(ring[j][i], ring[j][i2], ring[j + 1][i2], ring[j + 1][i])
    end
  end
  nu.times { |i| mesh.add_polygon(bot, ring[-1][(i + 1) % nu], ring[-1][i]) }
  before = ents.grep(Sketchup::Face).length
  ents.add_faces_from_mesh(mesh, 4, m)
  ents.grep(Sketchup::Face)[before..-1].each { |f| f.back_material = m }
end

def tube(ents, a, b, r, m, seg = 10)
  v = Geom::Vector3d.new(b.x - a.x, b.y - a.y, b.z - a.z)
  return if v.length < 1e-9
  n = Geom::Vector3d.new(v.x, v.y, v.z).normalize
  tmp = n.parallel?(Geom::Vector3d.new(0, 0, 1)) ? Geom::Vector3d.new(1, 0, 0) : Geom::Vector3d.new(0, 0, 1)
  u = n * tmp
  u.normalize!
  mesh = Geom::PolygonMesh.new
  ra = []
  rb = []
  seg.times do |i|
    t = Geom::Transformation.rotation(a, n, 2 * Math::PI * i / seg)
    ra << mesh.add_point(Geom::Point3d.new(a.x + u.x * r, a.y + u.y * r, a.z + u.z * r).transform(t))
    rb << mesh.add_point(Geom::Point3d.new(b.x + u.x * r, b.y + u.y * r, b.z + u.z * r).transform(t))
  end
  seg.times do |i|
    i2 = (i + 1) % seg
    mesh.add_polygon(ra[i], ra[i2], rb[i2], rb[i])
  end
  before = ents.grep(Sketchup::Face).length
  ents.add_faces_from_mesh(mesh, 4, m)
  ents.grep(Sketchup::Face)[before..-1].each { |f| f.back_material = m }
end

# ---------------- 静态几何（先全建好，逐帧只改 alpha / 变换） ----------------
HALF = 6.0
grid = model.entities.add_group
grid.name = "舞台"
NC = 24
(0..NC).each do |k|
  v = -HALF + 2 * HALF * k / NC
  grid.entities.add_line(Geom::Point3d.new(v, -HALF, 0), Geom::Point3d.new(v, HALF, 0))
  grid.entities.add_line(Geom::Point3d.new(-HALF, v, 0), Geom::Point3d.new(HALF, v, 0))
end

seed_m = mat(model, "seed", [255, 214, 92], 1.0)
seed = model.entities.add_group
seed.name = "seed"
uv_sphere(seed.entities, Geom::Point3d.new(0, 0, 0), 0.3, seed_m)
SEED_P = Geom::Point3d.new(0, 0, 0.35)

TET_C = Geom::Point3d.new(0, 0, 3.4)
TET_S = 2.2
s = TET_S / Math.sqrt(3.0)
raw = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
vp = raw.map { |v| Geom::Point3d.new(TET_C.x + v[0] * s, TET_C.y + v[1] * s, TET_C.z + v[2] * s) }
vmat = mat(model, "k4v", [86, 156, 214], 1.0)
emat = mat(model, "k4e", [214, 96, 74], 1.0)
k4v = model.entities.add_group
k4v.name = "K4顶点"
vp.each { |p| uv_sphere(k4v.entities, p, 0.26, vmat) }
k4e = model.entities.add_group
k4e.name = "K4棱"
[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]].each { |i, j| tube(k4e.entities, vp[i], vp[j], 0.07, emat) }

# 三个频率槽：环放在 K₄ 上方
SLOT_R = 2.6
SLOT_Z = 7.4
SLOT_TXT = ["0.4055", "0.6931", "0.9808"]
slot_grps = []
slot_mats = []
3.times do |i|
  a = Math::PI / 2 + i * 2 * Math::PI / 3
  p = Geom::Point3d.new(SLOT_R * Math.cos(a), SLOT_R * Math.sin(a), SLOT_Z)
  m = mat(model, "slot#{i}", [110, 220, 200], 1.0)
  g = model.entities.add_group
  g.name = "槽#{i + 1}"
  uv_sphere(g.entities, Geom::Point3d.new(0, 0, 0), 0.34, m)
  g.transformation = Geom::Transformation.translation(Geom::Vector3d.new(p.x, p.y, p.z))
  slot_grps << g
  slot_mats << m
  model.entities.add_text(SLOT_TXT[i], p,
                          Geom::Vector3d.new(p.x * 0.35, p.y * 0.35, 0.9))
end

# 三个「代」：从各自的槽往外上长
PART_P = []
PART_M = []
part_grps = []
3.times do |i|
  a = Math::PI / 2 + i * 2 * Math::PI / 3
  p = Geom::Point3d.new(1.35 * SLOT_R * Math.cos(a), 1.35 * SLOT_R * Math.sin(a), SLOT_Z + 1.9)
  PART_P << p
  m = mat(model, "part#{i}", [240, 240, 250], 1.0)
  g = model.entities.add_group
  g.name = "代#{i + 1}"
  uv_sphere(g.entities, Geom::Point3d.new(0, 0, 0), 0.42, m)
  part_grps << g
  PART_M << m
end

model.commit_operation
log "geometry built"

# ---------------- 动画辅助 ----------------
def ease(x)
  x = 0.0 if x < 0.0
  x = 1.0 if x > 1.0
  x * x * (3 - 2 * x)
end

def clamp01(x)
  x < 0 ? 0.0 : (x > 1 ? 1.0 : x)
end

view = model.active_view
TGT = Geom::Point3d.new(0, 0, 3.6)
R_CAM = 30.0
EL = 26.0 * Math::PI / 180

caption = nil
cap_text = nil
def set_caption(model, caption, text)
  caption.erase if caption && caption.valid?
  caption = model.entities.add_text(text, Geom::Point3d.new(-6.5, -6.5, -0.2),
                                    Geom::Vector3d.new(0.2, -0.2, -0.3))
  caption
end

t0 = Time.now
(0...FRAMES).each do |f|
  t = f / FPS

  # ① 舞台：0–1.4 s 长出来
  gs = ease(clamp01((t - 0.05) / 1.35))
  grid.transformation = Geom::Transformation.scaling(Geom::Point3d.new(0, 0, 0), 0.02 + 0.98 * gs)

  # ② seed：1.2–2.2 s 出现并脉动
  ss = ease(clamp01((t - 1.15) / 0.95))
  pulse = 1.0 + 0.18 * Math.sin(2 * Math::PI * (t - 1.2) / 0.9) * (ss > 0.99 ? 1 : 0)
  seed.transformation = Geom::Transformation.translation(Geom::Vector3d.new(0, 0, SEED_P.z)) *
                        Geom::Transformation.scaling(Geom::Point3d.new(0, 0, 0), (0.05 + 0.95 * ss) * pulse)
  seed_m.alpha = ss

  # ③ K₄：2.2–5.2 s 从 seed 位置长出
  ks = ease(clamp01((t - 2.2) / 3.0))
  zlift = SEED_P.z + (TET_C.z - SEED_P.z) * ks
  k4v.transformation = Geom::Transformation.translation(Geom::Vector3d.new(0, 0, zlift - TET_C.z)) *
                       Geom::Transformation.scaling(TET_C, 0.04 + 0.96 * ks)
  k4e.transformation = Geom::Transformation.translation(Geom::Vector3d.new(0, 0, zlift - TET_C.z)) *
                       Geom::Transformation.scaling(TET_C, 0.04 + 0.96 * ks)
  k4a = ease(clamp01((t - 2.2) / 2.4))
  vmat.alpha = k4a
  emat.alpha = ease(clamp01((t - 2.9) / 2.0))

  # ④ 三个槽：5.2–7.0 s 亮起
  sa = ease(clamp01((t - 5.2) / 1.8))
  3.times do |i|
    slot_mats[i].alpha = sa
    k = (0.25 + 0.75 * sa) * (1.0 + 0.10 * Math.sin(2 * Math::PI * (t - 5.2) / 1.1 + i))
    slot_grps[i].transformation =
      Geom::Transformation.translation(Geom::Vector3d.new(
        SLOT_R * Math.cos(Math::PI / 2 + i * 2 * Math::PI / 3),
        SLOT_R * Math.sin(Math::PI / 2 + i * 2 * Math::PI / 3), SLOT_Z)) *
      Geom::Transformation.scaling(Geom::Point3d.new(0, 0, 0), k)
  end

  # ⑤ 三个「代」：7.0–10.5 s 从槽里长出来
  pa = ease(clamp01((t - 7.0) / 3.5))
  3.times do |i|
    a = Math::PI / 2 + i * 2 * Math::PI / 3
    from = Geom::Point3d.new(SLOT_R * Math.cos(a), SLOT_R * Math.sin(a), SLOT_Z)
    pp = PART_P[i]
    cur = Geom::Point3d.new(from.x + (pp.x - from.x) * pa,
                            from.y + (pp.y - from.y) * pa,
                            from.z + (pp.z - from.z) * pa)
    rad = 0.12 + 0.30 * pa
    part_grps[i].transformation = Geom::Transformation.translation(Geom::Vector3d.new(cur.x, cur.y, cur.z)) *
                                  Geom::Transformation.scaling(Geom::Point3d.new(0, 0, 0), rad / 0.42)
    PART_M[i].alpha = ease(clamp01((t - 7.0) / 1.2))
  end

  # 相机：缓慢环绕
  az = (-55.0 + 22.0 * (t / DUR)) * Math::PI / 180
  eye = Geom::Point3d.new(TGT.x + R_CAM * Math.cos(EL) * Math.cos(az),
                          TGT.y + R_CAM * Math.cos(EL) * Math.sin(az),
                          TGT.z + R_CAM * Math.sin(EL))
  view.camera.set(eye, TGT, Geom::Vector3d.new(0, 0, 1))

  # 字幕
  newcap = if t < 1.2 then "0 · 空舞台"
           elsif t < 2.2 then "seed 点亮"
           elsif t < 5.2 then "K₄ 稳态长成（4 顶点 / 6 棱）"
           elsif t < 7.0 then "三个频率槽：0.4055 / 0.6931 / 0.9808"
           elsif t < 10.6 then "三槽各自长出一代"
           else "结构示意：口袋里的三代（不给数值刻度）"
           end
  if newcap != cap_text
    caption = set_caption(model, caption, newcap)
    cap_text = newcap
  end

  view.invalidate
  path = File.join(OUT, format("f_%04d.png", f))
  view.write_image(filename: path, width: WIDTH, height: HEIGHT, antialias: true, compression: 0.9)
  log "frame #{f + 1}/#{FRAMES} t=#{format('%.2f', t)}" if (f % 20).zero?
end

log "DONE elapsed=#{(Time.now - t0).round(1)}s"
