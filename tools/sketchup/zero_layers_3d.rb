# encoding: UTF-8
begin
  File.open("/tmp/zero_layers_3d.log", "a") { |f| f.puts "SCRIPT STARTED #{Time.now}" }
rescue => e
end
# zero_layers_3d.rb — 用 SketchUp 的 Ruby 接口直接生成 3D 模型（.skp）
#
# 内容：同心层结构（外 → 内）L0 → L1′ → L1 → L2；最内的 L2 球一分为三（3 个 120° 扇面）。
# 与 2D 版 ZeroGPU/src/zero_layers_spheres.py 同一套层序/颜色口径。
#
# 跑法（Mac，走 LaunchServices 才能让 -RubyStartup 生效）：
#   open -a "/Applications/SketchUp 2025/SketchUp.app" --args -RubyStartup tools/sketchup/zero_layers_3d.rb
# 或在 SketchUp 的 Ruby 控制台里： load "/Users/oygb/Downloads/Zero/tools/sketchup/zero_layers_3d.rb"
#
# 输出：ZeroGPU/results/zero_layers_3d.skp ＋ 两张预览 png

OUT_DIR = "/Users/oygb/Downloads/Zero/ZeroGPU/results"
LOG     = "/tmp/zero_layers_3d.log"

def log(msg)
  File.open(LOG, "a") { |f| f.puts msg }
end

File.delete(LOG) if File.exist?(LOG)
Dir.mkdir(OUT_DIR) unless Dir.exist?(OUT_DIR)

model = Sketchup.active_model
model.entities.clear!
model.start_operation("Zero layers 3D", true)

begin
# ---------- 参数 ----------
SHELLS = [
  { key: "L0",  r: 5.0, color: [ 76, 120, 168] },
  { key: "L1p", r: 4.0, color: [ 89, 161,  79] },
  { key: "L1",  r: 3.0, color: [242, 142,  43] },
  { key: "L2",  r: 1.9, color: [225,  87,  89] }
]
SECTOR_COLORS = [[225, 87, 89], [118, 183, 178], [176, 122, 161]]
NU, NV = 48, 24          # 经向/纬向网格
ALPHA_SHELL = 0.22
ALPHA_LUNE  = 0.62

def material(model, name, rgb, alpha)
  m = model.materials[name] || model.materials.add(name)
  m.color = Sketchup::Color.new(*rgb)
  m.alpha = alpha
  m
end

# 球面网格（三角形/四边形），phi 范围可裁剪 → 用来做整球或扇面
def sphere_faces(ents, c, r, mat, phi0 = 0.0, phi1 = 2 * Math::PI, nu = NU, nv = NV)
  mesh = Geom::PolygonMesh.new
  top = mesh.add_point(Geom::Point3d.new(c.x, c.y, c.z + r))
  bot = mesh.add_point(Geom::Point3d.new(c.x, c.y, c.z - r))
  ring = []
  (1...nv).each do |j|
    theta = Math::PI * j / nv
    row = []
    nu.times do |i|
      phi = phi0 + (phi1 - phi0) * i / nu
      row << mesh.add_point(Geom::Point3d.new(
        c.x + r * Math.sin(theta) * Math.cos(phi),
        c.y + r * Math.sin(theta) * Math.sin(phi),
        c.z + r * Math.cos(theta)))
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
  ents.grep(Sketchup::Face).each { |f| f.back_material = mat } if mat
end

# 球面经纬线（线框感）
def graticule(ents, r, segs = 48)
  [0, 30, -30, 60, -60].each do |deg|
    th = deg * Math::PI / 180
    next if th.abs > 1.5
    ents.add_circle(Geom::Point3d.new(0, 0, r * Math.sin(th)),
                    Geom::Vector3d.new(0, 0, 1), r * Math.cos(th).abs, segs)
  end
  6.times do |k|
    a = k * Math::PI / 6
    ents.add_circle(Geom::Point3d.new(0, 0, 0),
                    Geom::Vector3d.new(Math.cos(a), Math.sin(a), 0), r, segs)
  end
end

# ---------- 建模型 ----------
SHELLS.each_with_index do |sh, idx|
  g = model.entities.add_group
  g.name = (sh[:key] == "L1p" ? "L1'" : sh[:key])
  mat = material(model, "shell_#{sh[:key]}", sh[:color], ALPHA_SHELL)
  if sh[:key] == "L2"
    # L2 一分为三：三个 120° 扇面，各给一色
    SUB = [4, 18, 6]   # 三个扇面各自再分的份数：4=K₄ 顶点 / 18=3 代×6 边 / 6=K₄ 边
    3.times do |s|
      sg = g.entities.add_group
      sg.name = "L2_sector#{s + 1}"
      m = material(model, "lune_#{s}", SECTOR_COLORS[s], ALPHA_LUNE)
      a0 = s * 2 * Math::PI / 3
      sphere_faces(sg.entities, Geom::Point3d.new(0, 0, 0), sh[:r], m, a0, a0 + 2 * Math::PI / 3, NU, NV)
      # 亚亚层：把扇面沿经向切成 SUB[s] 条
      n = SUB[s]
      (1...n).each do |k|
        a = a0 + (2 * Math::PI / 3) * k / n
        sg.entities.add_circle(Geom::Point3d.new(0, 0, 0),
                               Geom::Vector3d.new(Math.cos(a), Math.sin(a), 0), sh[:r], 64)
      end
    end
    3.times do |s|
      a = s * 2 * Math::PI / 3
      g.entities.add_circle(Geom::Point3d.new(0, 0, 0),
                            Geom::Vector3d.new(Math.cos(a), Math.sin(a), 0), sh[:r], 64)
    end
    SUB_TXT = ["4 份", "18 份", "6 份"]
    3.times do |s|
      mid = s * 2 * Math::PI / 3 + Math::PI / 3
      d = Geom::Vector3d.new(Math.cos(mid), Math.sin(mid), 0)
      pt = Geom::Point3d.new(d.x * sh[:r] * 1.06, d.y * sh[:r] * 1.06, -0.9)
      model.entities.add_text(SUB_TXT[s], pt, Geom::Vector3d.new(d.x * 0.45, d.y * 0.45, 0))
    end
  else
    sphere_faces(g.entities, Geom::Point3d.new(0, 0, 0), sh[:r], mat)
    graticule(g.entities, sh[:r])
  end
  log "built #{g.name} r=#{sh[:r]}"
end

# 坐标轴
ax = model.entities.add_group
ax.name = "axes"
axmat = material(model, "axis", [90, 90, 90], 1.0)
[[[7, 0, 0], [-7, 0, 0]], [[0, 7, 0], [0, -7, 0]], [[0, 0, 7], [0, 0, -7]]].each do |a, b|
  e = ax.entities.add_line(Geom::Point3d.new(*a), Geom::Point3d.new(*b))
  e.material = axmat if e
end

# 标签（每个壳用自己的方位角，避免叠在一起）
LABEL_AZ = { "L0" => 25, "L1p" => 100, "L1" => 160, "L2" => 250 }
SHELLS.each do |sh|
  az = LABEL_AZ[sh[:key]] * Math::PI / 180
  el = 32 * Math::PI / 180
  d = Geom::Vector3d.new(Math.cos(el) * Math.cos(az), Math.cos(el) * Math.sin(az), Math.sin(el))
  p = Geom::Point3d.new(d.x * sh[:r], d.y * sh[:r], d.z * sh[:r])
  label = (sh[:key] == "L1p" ? "L1′" : sh[:key])
  model.entities.add_text(label, p, Geom::Vector3d.new(d.x * 0.6, d.y * 0.6, d.z * 0.6))
end

model.entities.add_text("L2：3 个 120° 扇面（1→4 份＝K₄ 顶点；2→18 份＝3 代×6 边；3→6 份＝K₄ 边）",
                        Geom::Point3d.new(0, 0, -2.2),
                        Geom::Vector3d.new(0, 0, -0.7))
model.commit_operation

# ---------- 视角与出图 ----------
view = model.active_view
cams = [
  [Geom::Point3d.new(19, -19, 13), Geom::Point3d.new(0, 0, 0)],
  [Geom::Point3d.new(-4, -24, 2), Geom::Point3d.new(0, 0, 0)]
]
cams.each_with_index do |(eye, tgt), i|
  view.camera.set(eye, tgt, Geom::Vector3d.new(0, 0, 1))
  view.zoom_extents
  view.invalidate
  path = File.join(OUT_DIR, "zero_layers_3d_view#{i + 1}.png")
  ok = view.write_image(filename: path, width: 1600, height: 1000, antialias: true, compression: 0.9)
  log "image#{i + 1} #{ok} #{path}"
end
view.camera.set(cams[0][0], cams[0][1], Geom::Vector3d.new(0, 0, 1))
view.invalidate

skp = File.join(OUT_DIR, "zero_layers_3d.skp")
begin
  model.save(skp)
  log "saved #{skp} size=#{File.size(skp)}"
rescue => e
  log "SAVE ERR #{e.class}: #{e.message}"
end
log "DONE"

rescue => e
  File.open(LOG, "a") { |f|
    f.puts "FATAL #{e.class}: #{e.message}"
    f.puts (e.backtrace || []).first(8).join("\n")
  }
end
