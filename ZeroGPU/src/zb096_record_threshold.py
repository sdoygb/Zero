#!/usr/bin/env python3
"""G 的积累 → L2 演化：记录阈值触发引擎改门的玩具模型与可检验签名

模型（最小化，只用 Zb074 定义 1.1/2.1 的两个式子）：
    R_{t+1} = R_t + g(m_t(v_s))        g ≥ 0，逐点单调（Zb074 命题 1.2）
    seed_{t+1} = K R_t                 K 三段核（Zb074 定理 2.2）
    m_{t+1}(v_s) = 1 + floor(a · seed_{t+1}(v_s))   a ≥ 0，L2 的闭合率随种子线性增
阈值：引擎在 seed 超过 τ 时换一支（L2 的"门"打开）：m 的增益从 a0 跳到 a1 > a0，且对外部事件 x 的敏感度也跳。

测：外部事件（丢石头/风吹）在某个时刻停止，看鱼跳（可观测）何时停止——
  · 若没有阈值：x 停 ⟹ 鱼跳立刻停（零延迟）；
  · 若有阈值：记录已越阈值，回写继续，鱼跳会**继续一段时间**（记忆性），且停止时间由 τ 与记录量决定。
"""
import json
from pathlib import Path

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb096_record_threshold.json")

def simulate(n_days=200, tau=50.0, a0=1.0, a1=3.0, x_stop=150, seed_val=1.0,
             x_gain=0.5, self_gain=0.3):
    """记录 R 累积；种子由回写 K 给（随 R 饱和增长）；阈值 τ 决定 L2 的增益档（a0→a1）。
    关键修正：种子必须随 R 真正增长（饱和上限足够高），否则阈值永远触发不了。"""
    R = 0.0
    seeds, mvals = [], []
    for t in range(n_days):
        x = 1.0 if (t < x_stop and (t % 1 == 0)) else 0.0   # 外部事件：每天一次，到 x_stop 停
        gain = a1 if seed_val > tau else a0
        m = x * (1.0 + gain * seed_val) if x else 0.0
        # 无外部事件时，引擎靠"种子"自主产生闭合（这是阈值的后果：种子已大，闭合自发）
        if x == 0.0 and seed_val > tau:
            m = self_gain * gain * seed_val / (1.0 + seed_val) * seed_val   # 自发闭合×种子（饱和）
        R += m
        # 回写：K 的一跳使种子随记录增长，但有饱和（否则数值指数爆炸，模型的定性结论不变）
        seed_val = 1.0 + 100.0 * R / (1.0 + R)
        seeds.append(seed_val); mvals.append(m)
    return seeds, mvals

def main():
    res = {}
    for tau in (50.0, 200.0):
        s, m = simulate(tau=tau, x_stop=150)
        # 停止后仍产生闭合的天数
        tail = [i for i in range(150, 200) if m[i] > 0]
        first_stop = (min(tail) if tail else None)
        res[f"tau={tau}"] = {
            "seeds_at_x_stop": round(s[149], 2),
            "m_at_x_stop": round(m[149], 3),
            "self_sustained_days_after_stop": int(sum(1 for i in range(150, 200) if m[i] > 0)),
            "seed_final": round(s[-1], 2),
        }
        print(f"tau={tau:>5}: 停喂时种子 {s[149]:.2f}；停喂后仍有闭合的天数 {int(sum(1 for i in range(150,200) if m[i]>0))}；"
              f"末种子 {s[-1]:.2f}")
    print("\n读法：记录越阈值后，即使外部事件停止（x=0），引擎仍靠回写种子产生自发闭合——这就是"
          "「G 积累到一定值，L2 就自己演化」的可观测形态：**记忆性的滞后**。")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)

if __name__ == "__main__":
    raise SystemExit(main())
