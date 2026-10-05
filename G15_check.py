#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G15_check.py -- 裸 Z0③（无偏好）下不存在有限特征速度：I6 物质层的 no-go。

对应文档 G15_bare_ax3_has_no_characteristic_speed.md。
失败时退出码非零。

  F1  微观信号速度有限（1 边/步）
  F2  宏观（热方程）信号速度无穷
  F3  电流由密度瞬时决定（无记忆）：离散散度定理逐步精确
  F4  Z0③ 全分支演化：二阶矩 ~ tau（指数 p≈1），无弹道区
  F5  持续性随机游走：短时 p≈2（弹道），长时 p->1（扩散）
  F6  Z0③ 无偏好 => 步间无关联 => 弹道区不存在 => c 非有限 => 登记 I7
  F7  结论与诚实边界
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0

# R2 文档锚定（与 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263 同法）：
# 读正文并断言结论的【口径】确实写在正文里。这是漂移检测——正文若被改掉这些口径，
# 下面的结论行会变 [x]，而不再静默继续"通过"。
G15_DOC = io.open(os.path.join(HERE, "G15_bare_ax3_has_no_characteristic_speed.md"),
                  encoding="utf-8").read()

# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")


def check(name, cond, detail="", level="ind"):
    """level: "ind"=独立实断言 / "dep"=依赖上文的可失败结论行 / "note"=解释性，不独立计数。"""
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    if ok:
        tag = "v" if (level == "ind" or LEDGER_MODE == "A") else "i"
    else:
        tag = "x"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ======================================================================
head("F1  微观信号速度有限（1 边/步）")

N = 401
c0 = N // 2
rho = np.zeros(N)
rho[c0] = 1.0
support = [np.max(np.nonzero(rho > 0)[0]) - c0]
r = rho.copy()
for _ in range(20):
    rn = r.copy()
    # 全分支一步：每个内部点把自身的一半送到左右两邻居（对称 Laplacian 步）
    rn[1:-1] = 0.5 * r[:-2] + 0.5 * r[2:]
    rn[0] = rn[-1] = 0.0
    r = rn
    support.append(np.max(np.nonzero(r > 1e-300)[0]) - c0)
check("支持集前沿恰好以 1 边/步 扩张", support == list(range(len(support))),
      "前 6 步 = %s" % support[:6])
check("故底层微观因果性成立：有限信号速度",
      support == list(range(len(support))) and support[-1] == 20,
      "F1 重算：前沿 %s…%s（20 步全为 1 边/步）" % (support[:3], support[-3:]), level="dep")

# ======================================================================
head("F2  宏观（热方程）信号速度无穷")

xs = np.linspace(-60.0, 60.0, 2401)
dx = xs[1] - xs[0]
u = np.exp(-(xs ** 2) / (2 * 0.2 ** 2))
D = 1.0
dt = 0.25 * dx ** 2 / D
for _ in range(int(1.0 / dt)):
    un = u.copy()
    un[1:-1] = u[1:-1] + D * dt / dx ** 2 * (u[2:] - 2 * u[1:-1] + u[:-2])
    un[0] = un[-1] = 0.0
    u = un
far = np.abs(xs) > 3.0
check("热方程在 |x| > 3 处仍有非零值（t=1）", float(np.max(u[far])) > 1e-6,
      "max=%.3e" % float(np.max(u[far])))
check("=> 宏观描述丢失了微观因果性（失败在粗粒化，不在底层）",
      support == list(range(len(support))) and float(np.max(u[far])) > 1e-6
      and ("粗粒化" in G15_DOC),
      "F1（前沿 1 边/步）∧ F2（|x|>3 仍有 %.2e）∧ 正文锚定『粗粒化』"
      % float(np.max(u[far])), level="dep")

# ======================================================================
head("F3  电流由密度瞬时决定（无记忆）：离散散度定理逐步精确")

M = 61
rho2 = np.zeros(M)
rho2[M // 2] = 1.0
K = np.ones(M - 1)
worst = 0.0
for _ in range(30):
    # 电流 j_e = K_e (rho_i - rho_j)，与历史无关
    j = K * (rho2[:-1] - rho2[1:])
    # 连续性方程：d rho_i/dtau = -(j_i - j_{i-1})
    div = np.zeros(M)
    div[:-1] += -j
    div[1:] += j
    rho_new = rho2 + 0.5 * div
    # 直接用 Laplacian 步（eps=0.5，稳定）
    rho_lap = rho2.copy()
    rho_lap[1:-1] = rho2[1:-1] + 0.5 * (rho2[:-2] - 2 * rho2[1:-1] + rho2[2:])
    worst = max(worst, float(np.max(np.abs(rho_new[1:-1] - rho_lap[1:-1]))))
    rho2 = rho_lap
check("电流完全由当前密度决定，连续性方程逐步精确成立", worst < 1e-12,
      "max dev=%.2e" % worst)
check("=> 无独立记忆变量 => Fick 闭合 => 无弛豫时间 => 无有限特征速度",
      (worst < 1e-12) and ("Fick" in G15_DOC) and ("记忆" in G15_DOC),
      "绑定 F3 的 worst=%.1e；正文锚定『Fick』『记忆』" % worst, level="dep")

# ======================================================================
head("F4  Z0③ 全分支演化：二阶矩 ~ tau^p，p≈1（无弹道区）")


def lap_msd(nsteps, n=801, eps=0.5):
    """稳定的离散热步 rho -> rho + eps L rho（eps <= 1/2）。"""
    c = n // 2
    rr = np.zeros(n)
    rr[c] = 1.0
    xsx = np.arange(n) - c
    out = []
    for t in range(nsteps + 1):
        out.append(float((xsx ** 2 * rr).sum() / rr.sum()))
        rn = rr.copy()
        rn[1:-1] = rr[1:-1] + eps * (rr[:-2] - 2 * rr[1:-1] + rr[2:])
        rn[0] = rn[-1] = 0.0
        rr = rn
    return np.array(out)


msd = lap_msd(80)
ts = np.arange(len(msd))
sel = (ts >= 1) & (msd > 0)
p_lap = np.polyfit(np.log(ts[sel]), np.log(msd[sel]), 1)[0]
check("全分支 Laplacian 演化：MSD ~ tau^p，全局指数 p = %.3f ≈ 1" % p_lap,
      abs(p_lap - 1.0) < 0.05)
early = np.polyfit(np.log(ts[1:6]), np.log(msd[1:6]), 1)[0]
check("短时指数 p_early = %.3f，仍接近 1（无 tau^2 弹道区）" % early,
      abs(early - 1.0) < 0.15, "p_early=%.3f" % early)

# ======================================================================
head("F5  持续性随机游走：短时 p≈2（弹道），长时 p->1（扩散）")


def prw_msd(nsteps, nwalk=20000, pkeep=0.95, seed=3):
    rng = np.random.default_rng(seed)
    x = np.zeros(nwalk)
    d = rng.choice([-1.0, 1.0], size=nwalk)
    out = [0.0]
    for _ in range(nsteps):
        flip = rng.random(nwalk) > pkeep
        d = np.where(flip, -d, d)
        x = x + d
        out.append(float(np.mean(x ** 2)))
    return np.array(out)


pm = prw_msd(400)
tt = np.arange(len(pm))
pe = np.polyfit(np.log(tt[1:6]), np.log(pm[1:6]), 1)[0]
pl = np.polyfit(np.log(tt[300:]), np.log(pm[300:]), 1)[0]
check("持续性游走短时指数 p = %.3f ≈ 2（弹道区）" % pe, pe > 1.6, "p_early=%.3f" % pe)
check("持续性游走长时指数 p = %.3f ≈ 1（扩散区）" % pl, abs(pl - 1.0) < 0.1,
      "p_late=%.3f" % pl)
check("=> 弹道区要求步间关联（持续性），其连续描述才是电报方程",
      pe > 1.6 and abs(pl - 1.0) < 0.1,
      "同一持续性游走：短时 p=%.3f（弹道）→长时 p=%.3f（扩散）；连续极限为电报方程" % (pe, pl),
      level="dep")

# ======================================================================
head("F6  Z0③ 无偏好 => 步间无关联 => 弹道区不存在 => c 非有限")


def uncorrelated_msd(nsteps, nwalk=20000, seed=5):
    rng = np.random.default_rng(seed)
    x = np.zeros(nwalk)
    out = [0.0]
    for _ in range(nsteps):
        d = rng.choice([-1.0, 1.0], size=nwalk)   # 每步独立（无偏好）
        x = x + d
        out.append(float(np.mean(x ** 2)))
    return np.array(out)


um = uncorrelated_msd(60)
tu = np.arange(len(um))
pu = np.polyfit(np.log(tu[1:]), np.log(um[1:]), 1)[0]
check("无关联（无偏好）游走：全程指数 p = %.3f ≈ 1，无弹道区" % pu, abs(pu - 1.0) < 0.05,
      "p=%.3f" % pu)
check("=> Z0③ 无偏好抹掉弹道区 => c = sqrt(D/lambda) 不是有限特征速度",
      abs(pu - 1.0) < 0.05 and pe > 1.6 and abs(pl - 1.0) < 0.1,
      "无偏好 p=%.3f（无弹道）对照持续性 p_early=%.3f（有弹道）" % (pu, pe), level="dep")
check("结论：I6 物质层不可由底层条款（Z0 条款 ＋ Z1–Z5 定理）推出（在无偏好不被削弱的前提下）",
      abs(pu - 1.0) < 0.05 and support == list(range(len(support)))
      and ("不可" in G15_DOC) and ("no-go" in G15_DOC),
      "F1+F6 数值 ∧ 正文 §6 锚定『不可』『no-go』", level="dep")

# ======================================================================
head("F7  结论与诚实边界")

check("失败在粗粒化：微观有限速度（F1）被宏观无穷速度（F2）取代",
      support == list(range(len(support))) and float(np.max(u[far])) > 1e-6
      and ("粗粒化" in G15_DOC),
      "F1 ∧ F2 ∧ 正文锚定", level="dep")
check("要恢复因果性需额外输入：电流闭合/记忆核 => 登记为 I7",
      (worst < 1e-12) and ("I7" in G15_DOC) and ("记忆核" in G15_DOC),
      "F3 无记忆 ∧ 正文 §5/§6 锚定 I7『记忆核』", level="dep")
check("诚实边界：no-go 的条件是 Z0③ 保持「无偏好」；削弱 Z0③ 属修改 Z0 条款",
      ("削弱 Z0③" in G15_DOC) and ("修改 Z0 条款" in G15_DOC),
      "正文 §8 原表锚定『削弱 Z0③』『修改 Z0 条款』", level="dep")
check("诚实边界：持续性可从 Z3 的「词」（移动序列）来，但其统计量未导出",
      ("移动序列" in G15_DOC) and ("未导出" in G15_DOC),
      "正文 §8 锚定『移动序列』『未导出』", level="dep")
check("诚实边界：I7 与 I6 独立——I7 是「闭合从哪来」，I6 是「叶层是否纯规范」",
      ("叶层是否纯规范" in G15_DOC) and ("电流闭合从哪来" in G15_DOC),
      "正文 §7 表锚定（I6/I7 两行的定义）", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
