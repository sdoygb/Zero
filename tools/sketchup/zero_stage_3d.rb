# encoding: UTF-8
# zero_stage_3d.rb — 舞台 3D：方格舞台（24×24）＋ 局部 K₄ 四面体（SketchUp Ruby 接口）
#
# 口径：方格舞台是具名输入（Z160／Z021），K₄ 是局部舞台（Z021 判定）。
# 本图只是结构示意，不含量值主张。
# 跑法见 tools/sketchup/zero_layers_3d.rb 顶部注释（同一个 -RubyStartup 机制）。

OUT_DIR = "/Users/oygb/Downloads/Zero/ZeroGPU/results"
LOG     = "/tmp/zero_stage_3d.log"

def log(msg)
  File.open(LOG, "a") { |f| f.puts msg }
end

File.delete(LOG) if File.exist?(LOG)
Dir.mkdir(OUT_DIR) unless Dir.exist?(OUT_DIR)

model = Sketchup.active_model
model.entities.clear!
model.start_operation("Zero stage 3D", true)

# ---------- 参数 ----------
NC   = 24          # 方格舞台 24×24
CELL = 0.5         # 每格边长（画法尺度，非物理长度）
HALF = NC * CELL / 2.0
TET_C = Geom::Point3d.new(0, 0, 3.2)
TET_S = 2.2        # 正四面体顶点到中心距离
V_R   = 0.22       # 顶点小球半径
E_R   = 0.06       # 棱半径

def material(model, name, rgb, alpha = 1.0)
  m = model.materials[name] || model.materials.add(name)
  m.color = Sketchup::Color.new(*rgb)
  m.alpha = alpha
  m
end

def uv_sphere(ents, c, r, mat, nu = 16, nv = 8)
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
  ents.add_faces_from_mesh(mesh, 4, mat)
  ents.grep(Sketchup::Face).each { |f| f.back_material = mat }
end


# 圆柱管：SketchUp 的 Entities 没有 add_cylinder，用 PolygonMesh 自己搭
def tube(ents, a, b, r, mat, seg = 12)
  v = Geom::Vector3d.new(b.x - a.x, b.y - a.y, b.z - a.z)
  return nil if v.length < 1e-9
  n = Geom::Vector3d.new(v.x, v.y, v.z).normalize
  tmp = n.parallel?(Geom::Vector3d.new(0, 0, 1)) ? Geom::Vector3d.new(1, 0, 0) : Geom::Vector3d.new(0, 0, 1)
  u = n * tmp
  u.normalize!
  mesh = Geom::PolygonMesh.new
  ra = []
  rb = []
  seg.times do |i|
    t = Geom::Transformation.rotation(a, n, 2 * Math::PI * i / seg)
    pa = Geom::Point3d.new(a.x + u.x * r, a.y + u.y * r, a.z + u.z * r).transform(t)
    pb = Geom::Point3d.new(b.x + u.x * r, b.y + u.y * r, b.z + u.z * r).transform(t)
    ra << mesh.add_point(pa)
    rb << mesh.add_point(pb)
  end
  seg.times do |i|
    i2 = (i + 1) % seg
    mesh.add_polygon(ra[i], ra[i2], rb[i2], rb[i])
  end
  before = ents.grep(Sketchup::Face).length
  ents.add_faces_from_mesh(mesh, 4, mat)
  ents.grep(Sketchup::Face)[before..-1].each { |f| f.back_material = mat }
end

# ---------- 方格舞台 ----------
grid = model.entities.add_group
grid.name = "方格舞台 #{NC}x#{NC}"
gmat = material(model, "grid", [110, 120, 130])
(0..NC).each do |k|
  v = -HALF + k * CELL
  grid.entities.add_line(Geom::Point3d.new(v, -HALF, 0), Geom::Point3d.new(v, HALF, 0))
  grid.entities.add_line(Geom::Point3d.new(-HALF, v, 0), Geom::Point3d.new(HALF, v, 0))
end
log "grid #{NC}x#{NC} edges=#{grid.entities.length}"

# ---------- K₄ 四面体 ----------
k4 = model.entities.add_group
k4.name = "K4"
s = TET_S / Math.sqrt(3.0)
raw = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
vpt = raw.map { |v| Geom::Point3d.new(TET_C.x + v[0] * s, TET_C.y + v[1] * s, TET_C.z + v[2] * s) }
vmat = material(model, "k4_vertex", [70, 130, 200])
emat = material(model, "k4_edge", [200, 90, 70])
vpt.each_with_index do |p, i|
  begin
    uv_sphere(k4.entities, p, V_R, vmat)
  rescue => e
    log "vertex #{i} ERR #{e.message}"
  end
end
pairs = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
pairs.each do |i, j|
  a, b = vpt[i], vpt[j]
  begin
    tube(k4.entities, a, b, E_R, emat, 12)
  rescue => e
    log "edge #{i}-#{j} ERR #{e.message}; fallback line"
    k4.entities.add_line(a, b)
  end
end
log "k4 vertices=#{vpt.length} edges=#{pairs.length}"

# ---------- 标签 ----------
model.entities.add_text("方格舞台 24×24（具名输入）",
                        Geom::Point3d.new(HALF, -HALF, 0.1),
                        Geom::Vector3d.new(0.8, -0.8, 0.6))
model.entities.add_text("K₄（局部舞台）",
                        Geom::Point3d.new(0, 0, TET_C.z + TET_S + 0.15),
                        Geom::Vector3d.new(0.3, 0.4, 0.9))
model.commit_operation

# ---------- 视角与出图 ----------
view = model.active_view
cams = [
  [Geom::Point3d.new(16, -16, 12), Geom::Point3d.new(0, 0, 1.4)],
  [Geom::Point3d.new(8, -15, 5), Geom::Point3d.new(0, 0, 2.2)]
]
cams.each_with_index do |(eye, tgt), i|
  view.camera.set(eye, tgt, Geom::Vector3d.new(0, 0, 1))
  view.zoom_extents
  view.invalidate
  path = File.join(OUT_DIR, "zero_stage_3d_view#{i + 1}.png")
  ok = view.write_image(filename: path, width: 1600, height: 1000, antialias: true, compression: 0.9)
  log "image#{i + 1} #{ok} #{path}"
end
view.camera.set(cams[0][0], cams[0][1], Geom::Vector3d.new(0, 0, 1))
view.invalidate

skp = File.join(OUT_DIR, "zero_stage_3d.skp")
begin
  model.save(skp)
  log "saved #{skp} size=#{File.size(skp)}"
rescue => e
  log "SAVE ERR #{e.class}: #{e.message}"
end
log "DONE"
