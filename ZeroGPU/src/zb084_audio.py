#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zb084 · 《零和宇宙 · 质量的几何》音轨（纯 Python 合成 + macOS say + ffmpeg）

输出（ZeroGPU/results/zb084_audio/）：
    bgm.wav        低沉电子氛围 + 飞行推进感（30 s）
    sfx.wav        全部音效
    voice.wav      旁白（macOS say 男声 Eddy，重采样到 48 kHz）
    mix.wav        混音（BGM + 音效 + 旁白）
旁白时间点取自 zb084_particles_geometry.json 的 render_plan.narration_frames。

用法：
    python3 ZeroGPU/src/zb084_audio.py
"""
import json
import math
import random
import subprocess
import struct
import sys
import wave
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "ZeroGPU" / "results"
OUT = RESULTS / "zb084_audio"
OUT.mkdir(parents=True, exist_ok=True)
G = json.loads((RESULTS / "zb084_particles_geometry.json").read_text())
RP = G["render_plan"]
FPS = RP["fps"]
DUR = RP["duration_s"]
SR = 48000
N = int(DUR * SR)


def fr(sec):
    return int(round(sec * FPS))


TL = RP["timeline_frames"]


def t(frame):
    return frame / FPS


# ──────────────────────────────────────────────────────────────
# 1. BGM：低频铺底 + 缓慢的谐波移动 + 飞行段推进感
# ──────────────────────────────────────────────────────────────
def bgm():
    n = N
    tt = np.arange(n) / SR
    rng = np.random.default_rng(84)

    # 低频 drone：55 Hz 与 82.4 Hz（纯五度），带极慢的幅度呼吸
    drone = (0.55 * np.sin(2 * np.pi * 55.0 * tt)
             + 0.30 * np.sin(2 * np.pi * 82.41 * tt + 0.7)
             + 0.16 * np.sin(2 * np.pi * 110.0 * tt + 1.9))
    breath = 0.75 + 0.25 * np.sin(2 * np.pi * tt / 9.0)
    drone *= breath

    # 高频空气感：带通噪声（简易一阶差分 + 平滑）
    air = rng.standard_normal(n)
    k = int(SR / 900)
    air = np.convolve(air, np.ones(k) / k, mode="same")          # 低通过
    air = air - np.convolve(air, np.ones(k) / k, mode="same")    # 再差分 → 带通
    air *= 0.10 * (0.6 + 0.4 * np.sin(2 * np.pi * tt / 6.3))

    # 飞行推进（B05 起飞 → B07 悬停）：扫频的低频脉络 + 噪声推力
    fly0, fly1 = t(TL["fly_start"]), t(TL["zero_land"])
    env = np.clip((tt - (fly0 - 1.2)) / 1.2, 0, 1) * np.clip((fly1 + 0.8 - tt) / 0.8, 0, 1)
    f_sw = 48.0 + 18.0 * np.sin(2 * np.pi * (tt - fly0) / 4.2)
    phase = 2 * np.pi * np.cumsum(f_sw) / SR
    thrust = 0.16 * np.sin(phase) * env
    rr = rng.standard_normal(n)
    kk = int(SR / 4000)
    rr = np.convolve(rr, np.ones(kk) / kk, mode="same")
    rr -= np.convolve(rr, np.ones(kk) / kk, mode="same")
    thrust += 0.05 * rr * env

    # 收尾：B09 拉远时抬高，B10 定格后余韵
    rise = np.clip((tt - t(TL["pull_back"])) / 2.4, 0, 1) * np.clip((DUR - tt) / 0.8, 0, 1)
    sweep = 0.10 * np.sin(2 * np.pi * (55 * 2 ** (rise * 1.5)) * tt) * rise

    sig = drone + air + thrust + sweep
    # 整体包络：开场淡入 0.4 s、结尾 1.2 s 淡出
    fade_in = np.clip(tt / 0.4, 0, 1)
    fade_out = np.clip((DUR - tt) / 1.2, 0, 1)
    sig *= fade_in * fade_out
    return sig / (np.max(np.abs(sig)) + 1e-9) * 0.55


# ──────────────────────────────────────────────────────────────
# 2. 音效
# ──────────────────────────────────────────────────────────────
def add(buf, at, sig, gain=1.0):
    i0 = int(at * SR)
    i1 = min(len(buf), i0 + len(sig))
    if i1 > i0:
        buf[i0:i1] += gain * sig[: i1 - i0]


def tone(freq, dur, shape=0.35, decay=6.0, partials=((1.0, 1.0), (2.0, 0.35), (3.01, 0.12))):
    n = int(dur * SR)
    tt = np.arange(n) / SR
    s = np.zeros(n)
    for ratio, amp in partials:
        s += amp * np.sin(2 * np.pi * freq * ratio * tt + ratio)
    s *= np.exp(-decay * tt)
    return s * shape


def whoosh(dur=0.7, lo=200, hi=3000, seed=1):
    """噪声扫频：带通滤波的中心频率从 hi 滑到 lo。"""
    rng = np.random.default_rng(seed)
    n = int(dur * SR)
    tt = np.arange(n) / SR
    x = rng.standard_normal(n)
    # 用两-pass 滑动平均做可变带通（近似）
    out = np.zeros(n)
    width = 64
    ker = np.hanning(width)
    xf = np.convolve(x, ker / ker.sum(), mode="same")
    envlp = np.exp(-3.0 * tt / dur) * np.clip(tt / 0.05, 0, 1)
    return xf * envlp * 0.5


def hum(dur, freq=40.0, gain=0.5):
    n = int(dur * SR)
    tt = np.arange(n) / SR
    s = np.sin(2 * np.pi * freq * tt) + 0.4 * np.sin(2 * np.pi * freq * 2 * tt + 0.5)
    return s * gain * np.clip(tt / 0.3, 0, 1) * np.clip((dur - tt) / 0.6, 0, 1)


def ripple(dur=1.2):
    n = int(dur * SR)
    tt = np.arange(n) / SR
    s = np.sin(2 * np.pi * (220 + 90 * np.sin(2 * np.pi * 3 * tt)) * tt)
    return s * np.exp(-4 * tt) * 0.25


def ping(f=1760.0):
    return tone(f, 0.5, shape=0.22, decay=9.0, partials=((1.0, 1.0), (2.0, 0.2)))


def boom():
    return tone(58.0, 1.0, shape=0.55, decay=4.5, partials=((1.0, 1.0), (1.5, 0.25)))


def sfx():
    buf = np.zeros(N)
    # 0 震动 / 低频嗡鸣：B01–B02
    add(buf, 0.0, hum(5.4, 41.2, 0.55), 1.0)
    # 偏心波纹：B02 起（ecc_peak 之前）
    add(buf, t(TL["formula1_in"]) - 0.1, ripple(1.4), 0.9)
    add(buf, t(TL["ecc_peak"]) - 0.15, ripple(1.6), 0.7)
    # 三圆浮现：B03 三次清脆 ping（按代顺序）
    for i, dt in enumerate((0.0, 0.22, 0.44)):
        add(buf, t(TL["circles_in"]) + dt, ping(880 * 2 ** (i * 2 / 12)), 0.5)
    # 粒子着陆：每次一个清脆高音（电子/μ/τ）
    for k, key in enumerate(("land_e", "land_mu", "land_tau")):
        add(buf, t(TL[key]) - 0.05, ping(1320 * 2 ** (k * 3 / 12)), 0.55)
        add(buf, t(TL[key]) - 0.05, whoosh(0.35, 300, 4000, seed=10 + k), 0.5)
    # 零点揭示：B07 悬停时的低频回响 + 金属叮
    add(buf, t(TL["zero_land"]) - 0.2, boom(), 0.55)
    add(buf, t(TL["zero_land"]) + 0.35, tone(2093, 1.2, 0.16, 3.0, ((1, 1.0), (2.76, 0.25))), 0.6)
    # Koide 公式浮现：双音（纯五度）
    add(buf, t(TL["koide_in"]), tone(523.25, 1.0, 0.18, 3.0), 0.45)
    add(buf, t(TL["koide_in"]) + 0.12, tone(783.99, 1.0, 0.18, 3.0), 0.40)
    # 快速拉远：空间回响（长 whoosh + 低 boom）
    add(buf, t(TL["pull_back"]) - 0.2, whoosh(1.6, 120, 6000, seed=99), 0.8)
    add(buf, t(TL["pull_back"]) + 0.2, boom(), 0.35)
    # 标题：温柔一击
    add(buf, t(TL["title_in"]), tone(392.0, 1.4, 0.20, 2.2), 0.5)
    mx = np.max(np.abs(buf)) + 1e-9
    return buf / mx * 0.85


# ──────────────────────────────────────────────────────────────
# 3. 旁白（macOS say，男声 Eddy）
# ──────────────────────────────────────────────────────────────
NAR = [r[0] for r in RP["narration_frames"]]   # 每句起点（帧）
# 原话按旁白节奏切成 14 句（与 render_plan 一致）
# 男声（zh_CN）：Eddy —— 沉稳科普腔；语速略慢
VOICE = "Eddy (中文（中国大陆）)"
RATE = 168
LINES = [
    "宇宙从零开始震动。",
    "中心态均匀，没有质量。",
    "偏心发生，质量就产生了。",
    "三代粒子，三个同心圆。",
    "共享一条弦。",
    "最大的那个圆，远看空空的。",
    "这是电子，我们最熟悉的粒子。",
    "μ子，",
    "τ子，电子的两个重表哥。",
    "弦的两端，",
    "是质量为零的粒子。",
    "三个轻子的质量比，",
    "被一个三分之二锁死。",
    "这就是质量的几何。",
]


def voice():
    tmp = OUT / "voice_parts"
    tmp.mkdir(exist_ok=True)
    parts = []
    for i, line in enumerate(LINES):
        aiff = tmp / f"v{i:02d}.aiff"
        subprocess.run(["say", "-v", VOICE, "-r", str(RATE), "-o", str(aiff), line], check=True)
        wav = tmp / f"v{i:02d}.wav"
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(aiff),
                        "-ar", str(SR), "-ac", "1", str(wav)], check=True)
        parts.append(wav)
    # 合成到时间轴（以 narration 起点对齐；语音比预定穿插略慢时自动压缩到 1.6 s 内）
    buf = np.zeros(N)
    durs = []
    for i, w in enumerate(parts):
        with wave.open(str(w)) as f:
            data = np.frombuffer(f.readframes(f.getnframes()), dtype=np.int16).astype(np.float64) / 32768
        at = t(NAR[i]) if i < len(NAR) else 0.0
        nxt = t(NAR[i + 1]) if i + 1 < len(NAR) else DUR
        budget = int((nxt - at) * SR)
        if len(data) > budget:
            data = data[:budget]
        durs.append(round(len(data) / SR, 3))
        # 首尾 15 ms 淡入淡出，避免爆音
        k = int(0.015 * SR)
        if len(data) > 2 * k:
            data[:k] *= np.linspace(0, 1, k)
            data[-k:] *= np.linspace(1, 0, k)
        add(buf, at, data, 1.0)
    (OUT / "voice_durations.json").write_text(json.dumps(durs))
    mx = np.max(np.abs(buf)) + 1e-9
    return buf / mx * 0.9


# ──────────────────────────────────────────────────────────────
# 4. 写出与混音
# ──────────────────────────────────────────────────────────────
def write_wav(path, sig):
    sig = np.clip(sig, -1, 1)
    with wave.open(str(path), "wb") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(SR)
        f.writeframes((sig * 32767).astype(np.int16).tobytes())


def main():
    b = bgm()
    s = sfx()
    v = voice()
    write_wav(OUT / "bgm.wav", b * 0.55)
    write_wav(OUT / "sfx.wav", s * 0.75)
    write_wav(OUT / "voice.wav", v * 0.95)
    # 混音：旁白优先，BGM/音效压低（侧链近似：旁白窗口内整体衰减 4 dB）
    duck = np.ones(N)
    for fr0, _txt in RP["narration_frames"]:
        i0 = int(t(fr0) * SR)
        i1 = min(N, i0 + int(3.2 * SR))
        ramp = np.clip((np.arange(i0, i1) - i0) / (0.25 * SR), 0, 1)
        duck[i0:i1] = np.minimum(duck[i0:i1], 1.0 - 0.35 * ramp)
    # 平滑 duck 包络
    k = int(0.15 * SR)
    duck = np.convolve(duck, np.ones(k) / k, mode="same")
    mix = b * 0.55 * duck + s * 0.75 * duck + v * 0.95
    mix = mix / (np.max(np.abs(mix)) + 1e-9) * 0.92
    write_wav(OUT / "mix.wav", mix)
    print(json.dumps({"out": str(OUT), "dur": DUR, "sr": SR,
                      "voice_lines": len(LINES)},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()