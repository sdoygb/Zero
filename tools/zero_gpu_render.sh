#!/usr/bin/env bash
# zero_gpu_render.sh — SketchUp 模型 → OBJ → Blender 渲染 / 实时预览
#
# 用法：
#   tools/zero_gpu_render.sh [选项] <model.skp> [out.png]
#
# 选项：
#   --engine EEVEE|CYCLES   渲染引擎（默认 EEVEE = 走显卡）
#   --samples N             采样（默认 128）
#   --res WxH               分辨率（默认 1600x1000）
#   --az A --el E           机位方位角/仰角
#   --viewer                生成 <model>_viewer.blend 并用 Blender 打开实时预览（GPU 视口）
#   --no-render             只导 OBJ / 存 .blend，不渲染 png
#
# 例：
#   tools/zero_gpu_render.sh ZeroGPU/results/zero_stage_3d.skp
#   tools/zero_gpu_render.sh --viewer ZeroGPU/results/zero_layers_3d.skp
set -euo pipefail

ENGINE="EEVEE"; SAMPLES=128; RES="1600x1000"; AZ="35"; EL="22"; VIEWER=0; DO_RENDER=1
POS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --engine) ENGINE="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --res) RES="$2"; shift 2 ;;
    --az) AZ="$2"; shift 2 ;;
    --el) EL="$2"; shift 2 ;;
    --viewer) VIEWER=1; shift ;;
    --no-render) DO_RENDER=0; shift ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) POS+=("$1"); shift ;;
  esac
done

SKP="${POS[0]:?用法: zero_gpu_render.sh [选项] <model.skp> [out.png]}"
OUT="${POS[1]:-${SKP%.skp}_render.png}"
OBJ="${SKP%.skp}.obj"
BLEND="${SKP%.skp}_viewer.blend"
HERE="$(cd "$(dirname "$0")" && pwd)"
BLENDER="/Applications/Blender.app/Contents/MacOS/Blender"
SKETCHUP="/Applications/SketchUp 2025/SketchUp.app"
[ -x "$BLENDER" ] || { echo "找不到 Blender：$BLENDER"; exit 1; }

# ---------- 1/3 SketchUp → OBJ ----------
if [ ! -f "$OBJ" ]; then
  echo "[1/3] SketchUp 导出 OBJ：$OBJ"
  pkill -f "SketchUp 2025" 2>/dev/null || true
  sleep 8
  rm -f /tmp/zero_export_obj.log
  SKP_IN="$SKP" OBJ_OUT="$OBJ" open -a "$SKETCHUP" "$SKP" \
       --args -RubyStartup "$HERE/sketchup/zero_export_obj.rb"
  for _ in $(seq 1 24); do
    sleep 5
    grep -qE "DONE|ERR|EMPTY" /tmp/zero_export_obj.log 2>/dev/null && break
  done
  grep -q DONE /tmp/zero_export_obj.log 2>/dev/null || {
    echo "导出失败："; cat /tmp/zero_export_obj.log 2>/dev/null; exit 1; }
else
  echo "[1/3] 已有 OBJ，跳过导出：$OBJ"
fi

# ---------- 2/3 Blender 渲染（可选）／存 .blend ----------
ARGS=(--obj "$OBJ" --out "$OUT" --engine "$ENGINE" --samples "$SAMPLES" --res "$RES" --az "$AZ" --el "$EL")
if [ "$VIEWER" = "1" ]; then ARGS+=(--save-blend "$BLEND"); fi
if [ "$DO_RENDER" = "1" ]; then
  echo "[2/3] Blender 渲染（$ENGINE）→ $OUT"
  "$BLENDER" -b --factory-startup -P "$HERE/blender/zero_render.py" -- "${ARGS[@]}" \
    | grep -E "^RENDER|^SAVED_BLEND" || true
else
  echo "[2/3] 只存 .blend（不渲染）"
  "$BLENDER" -b --factory-startup -P "$HERE/blender/zero_render.py" -- "${ARGS[@]}" \
    >/dev/null 2>&1 || true
fi

# ---------- 3/3 打开实时预览 ----------
if [ "$VIEWER" = "1" ]; then
  [ -f "$BLEND" ] || { echo "没生成 $BLEND"; exit 1; }
  echo "[3/3] 打开实时预览（EEVEE 视口，GPU）：$BLEND"
  rm -f /tmp/zero_viewer_boot.log
  "$BLENDER" "$BLEND" --python "$HERE/blender/zero_viewer_boot.py" >/tmp/zero_viewer.log 2>&1 &
  for _ in $(seq 1 20); do sleep 1; [ -f /tmp/zero_viewer_boot.log ] && break; done
  open -a "/Applications/Blender.app" 2>/dev/null || true
  echo "    视口状态：$(tr '\n' ' ' < /tmp/zero_viewer_boot.log 2>/dev/null)"
  echo "    提示：右键拖=转视角，Z=切换着色，F12=出终图"
fi
