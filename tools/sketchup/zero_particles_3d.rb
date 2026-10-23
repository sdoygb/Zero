# encoding: UTF-8
# zero_particles_3d.rb — 「从 0 演化出粒子」的建模底稿（SketchUp 建 → 导出 OBJ → Blender 渲染动画）
#
# 场景分层（每个 group 用名字区分，Blender 侧按前缀取对象做动画）：
#   stage_bar_*   方格舞台（24×24，用细条几何而不是边线，因为 OBJ 不导边）
#   seed          中心零点
#   k4_v0..3      K₄ 四顶点
#   k4_e0..5      K₄ 六棱
#   slot1..3      三个频率槽
#   particle1..3  三代
#
# 跑法：
#   open -a "/Applications/SketchUp 2025/SketchUp.app" \
#        --args -RubyStartup /path/to/zero_particles_3d.rb
# 环境变量：OBJ_OUT（默认 ZeroGPU/results/zero_particles.obj）

OUT_DIR = "/Users/oygb/Downloads/Zero/ZeroGPU/results"
LOG = "/tmp/zero_particles_3d.log"

def log(m)
  File.open(LOG, "a") { |f| f.puts m }
end
File.delete(LOG) if File.exist?(LOG)
Dir.mkdir(OUT_DIR) unless Dir.exist?(OUT_DIR)

model = Sketchup.active_model
model.entities.clear!
model.start_operation("particles 3D", true)

def mat(model, name, rgb, alpha = 1.0)
  m = model.materials[name] || model.materials.add(name)
  m.color = Sketchup::Color.new(*rgb)
  m.alpha = alpha
  m
end

def uv_sphere(ents, c, r, m, nu = 16, nv = 8)
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
  ents
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
  ents
end

# 长方体：给舞台细条用（OBJ 不导边，所以线要用几何做）
def bar(ents, p0, p1, w, h, m)
  v = (p1 - p0)
  n = v.normalize
  up = Geom::Vector3d.new(0, 0, 1)
  side = n * up
  side = Geom::Vector3d.new(1, 0, 0) if side.length < 1e-6
  side.normalize!
  up2 = side * n
  up2.normalize!
  pts = []
  [[-1, -1], [1, -1], [1, 1], [-1, 1]].each do |sx, sy|
    o = side * (sx * w / 2.0) + up2 * (sy * h / 2.0)
    pts << p0 + o
  end
  f = ents.add_face(pts)
  return unless f
  got = f.pushpull(v.length)
  ents.grep(Sketchup::Face).each { |ff| ff.material = m; ff.back_material = m }
  got
end

# ---------- 舞台：24×24，细条 ----------
HALF = 6.0
NC = 24
grid_m = mat(model, "stage", [70, 96, 120])
stage = model.entities.add_group
stage.name = "stage_bar_all"
(0..NC).each do |k|
  v = -HALF + 2 * HALF * k / NC
  bar(stage.entities, Geom::Point3d.new(v, -HALF, 0), Geom::Point3d.new(v, HALF, 0), 0.035, 0.035, grid_m)
  bar(stage.entities, Geom::Point3d.new(-HALF, v, 0), Geom::Point3d.new(HALF, v, 0), 0.035, 0.035, grid_m)
end
log "stage bars built"

# ---------- seed ----------
seed_m = mat(model, "seed", [255, 214, 92])
seed = model.entities.add_group
seed.name = "seed"
uv_sphere(seed.entities, Geom::Point3d.new(0, 0, 0.35), 0.30, seed_m)

# ---------- K₄ ----------
TET_C = Geom::Point3d.new(0, 0, 3.4)
TET_S = 2.2
s = TET_S / Math.sqrt(3.0)
raw = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
vp = raw.map { |v| Geom::Point3d.new(TET_C.x + v[0] * s, TET_C.y + v[1] * s, TET_C.z + v[2] * s) }
vmat = mat(model, "k4_vertex", [86, 156, 214])
emat = mat(model, "k4_edge", [214, 96, 74])
vp.each_with_index do |p, i|
  g = model.entities.add_group
  g.name = "k4_v#{i}"
  uv_sphere(g.entities, p, 0.26, vmat)
end
[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]].each_with_index do |(i, j), k|
  g = model.entities.add_group
  g.name = "k4_e#{k}"
  tube(g.entities, vp[i], vp[j], 0.075, emat)
end
log "k4 built"

# ---------- 三个频率槽 ＋ 三个粒子 ----------
SLOT_R, SLOT_Z = 2.6, 7.4
3.times do |i|
  a = Math::PI / 2 + i * 2 * Math::PI / 3
  p = Geom::Point3d.new(SLOT_R * Math.cos(a), SLOT_R * Math.sin(a), SLOT_Z)
  g = model.entities.add_group
  g.name = "slot#{i + 1}"
  uv_sphere(g.entities, p, 0.34, mat(model, "slot_#{i}", [110, 220, 200]))

  q = Geom::Point3d.new(1.55 * SLOT_R * Math.cos(a), 1.55 * SLOT_R * Math.sin(a), SLOT_Z + 2.7)
  g2 = model.entities.add_group
  g2.name = "particle#{i + 1}"
  uv_sphere(g2.entities, q, 0.42, mat(model, "particle_#{i}", [235, 242, 255]))
end
model.commit_operation
log "scene built: #{model.entities.length} top-level entities"

out = ENV['OBJ_OUT'] || File.join(OUT_DIR, "zero_particles.obj")
begin
  ok = model.export(out, { triangulated_faces: true, doublesided_faces: true,
                           edges: false, texture_maps: true, swap_yz: false })
  log "export #{ok} -> #{out} size=#{(File.size(out) rescue -1)}"
  log "DONE"
rescue => e
  log "EXPORT ERR #{e.class}: #{e.message}"
end
