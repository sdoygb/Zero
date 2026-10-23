# encoding: UTF-8
# zero_export_obj.rb — 把 SketchUp 模型导出成 OBJ（给 Blender 的 GPU 渲染用）
#
# 用法（SketchUp 已打开目标文件时）：Ruby 控制台 load 本文件，
# 或：
#   open -a "/Applications/SketchUp 2025/SketchUp.app" <file.skp> \
#        --args -RubyStartup /path/to/zero_export_obj.rb
#
# 环境变量：SKP_IN（要打开的 .skp，可空）／OBJ_OUT（输出路径）

LOG = "/tmp/zero_export_obj.log"

def log(m)
  File.open(LOG, "a") { |f| f.puts m }
end

begin
  File.delete(LOG) if File.exist?(LOG)
  in_path = ENV['SKP_IN']
  out = ENV['OBJ_OUT']

  if in_path && !in_path.empty? && File.exist?(in_path)
    Sketchup.open_file(in_path)
    log "opened #{in_path}"
  end

  model = Sketchup.active_model
  if model.nil? || model.entities.length.zero?
    log "EMPTY MODEL — nothing to export"
    raise "empty model"
  end

  p = out || File.join(File.dirname(model.path.to_s), "#{File.basename(model.path.to_s, '.*')}.obj")
  opts = {
    triangulated_faces: true,
    doublesided_faces: true,
    edges: false,
    texture_maps: true,
    swap_yz: false
  }
  ok = model.export(p, opts)
  log "export #{ok} -> #{p} size=#{(File.size(p) rescue -1)}"
  log "DONE"
rescue => e
  File.open(LOG, "a") { |f| f.puts "ERR #{e.class}: #{e.message}" }
end
