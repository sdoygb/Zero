# Zb073 · 有限能级 BCS 的 d-wave 配对强度：可算估计量与精确对角化参照

> 篇号：Zb073
> 状态：正式篇（推导篇）
> 前置：Zb063（有限能级 BCS 与阈值）、Zb037（能隙方程）、Zb028（张量积）
> 关联脚本：`ZeroGPU/src/zero_dwave_estimate.py`
> 数据：`ZeroGPU/results/osc_4x4_ed.json`
> 登记日期：2026-11-03

## 摘要

把 Zb063 的有限能级 BCS 阈值公式接到方格环面的 d-wave 通道上，得到一个只依赖单粒子谱的可算估计量
$$
\lambda_d(\mu)=\frac1{N_k}\sum_k\frac{g_k^2}{2\,|\varepsilon_k-\mu|},
\qquad g_k=\cos k_x-\cos k_y,
\qquad V_c=\frac1{\lambda_d},
\tag{0.1}
$$
其中 $\lambda_d$ 是 d-wave 通道的配对磁化率，$V_c$ 是该通道的 Thouless 阈值。在 $4\times4$ 与 $6\times6$ 周期环面上按填充 $n$ 求值，并与 OSC 数据集公布的 $4\times4$ 精确对角化 d-wave 关联 $P_d(n)$（$U/t=8$）对照：两种簇尺寸上两者的秩相关均为

$$
\rho_{\rm Sp}(\lambda_d,\,P^{\rm ED}_d)=-0.80 .
\tag{0.2}
$$

即裸带的 BCS 磁化率与相互作用体系的 d-wave 关联在本对照点上方向相反。差别的来源是相互作用自能：$\lambda_d$ 只看见单粒子谱，$P_d$ 还看见 $U/t=8$ 的 Mott 物理。

## §1 问题

OSC 挑战的配对指数取 $A_d(\delta)$ 来自 6×6 Hubbard 模型的变分蒙特卡洛求解。能否用一个**闭式、可秒算**的估计量代替它，是本篇要回答的问题；候选来自库内已有的有限能级 BCS（Zb063）：

$$
1=V\sum_{k=1}^{K}\frac{1}{2E_k},\qquad E_k=\sqrt{\varepsilon_k^2+\Delta^2},
\qquad V_c=\frac{2}{\sum_k 1/|\varepsilon_k|}.
\tag{1.1}
$$

Zb063 把 $V$ 记为**具名条件（C）**、把有限谱上的 BCS 记为有效模型。本篇沿用这两个记账，不改动它们。

## §2 估计量：d-wave 通道

对单带方格紧束缚谱 $\varepsilon_k=-2t(\cos k_x+\cos k_y)$，d-wave 通道的配对核为 $V_{kk'}=V g_k g_{k'}$。线性化能隙方程 $\Delta_k=-\frac1{N_k}\sum_{k'}V_{kk'}\Delta_{k'}/(2E_{k'})$ 给出（在 $T=0$、$\Delta\to0$ 的极限）

$$
1=V\,\lambda_d(\mu),\qquad
\lambda_d(\mu)=\frac1{N_k}\sum_k\frac{g_k^2}{2|\varepsilon_k-\mu|},
\tag{2.1}
$$

即 (0.1)。$\lambda_d$ 的倒数就是该通道的阈值 $V_c$：$\lambda_d$ 越大，越小的配对强度即可打开 d-wave 能隙。

**记账**：$V$ 是具名条件（C）；(2.1) 是线性响应的阈值判据，不是能隙数值；"$\lambda_d$ 大 $\Rightarrow$ $A_d$ 大"是**待检验的假设**，不是本篇的结论——§4 正是检验它。

## §3 方格环面上的算法

1. 取 $N\times N$ 周期环面，$k=(2\pi i/N,\,2\pi j/N)$，$i,j=0,\dots,N-1$；谱 $\varepsilon_k$ 与形状因子 $g_k$ 直接由上式给出，$t=1$。
2. 给定填充 $n\in[0,2]$，用二分法解 $(2/N_k)\sum_k\Theta(\mu-\varepsilon_k)=n$，得 $\mu(n)$。
3. 按 (2.1) 求和；离散谱上 $|\varepsilon_k-\mu|$ 可能为零，取正则化 $|\varepsilon_k-\mu|\to\max(|\varepsilon_k-\mu|,\eta)$，$\eta=10^{-3}t$。
4. 脚本 `ZeroGPU/src/zero_dwave_estimate.py` 输出 $\lambda_d$、$V_c$ 与 ED 的秩相关；结果落在 `ZeroGPU/results/zb073_dwave_estimate_N{N}.json`。

**算法边界**：离散谱上 $\mu$ 只能落在本征值上，$\lambda_d$ 因此对 $n$ 分段常值；$\eta$ 是数值正则化，不是物理展宽。这两点决定了下表的读数只能作**方向性**对照。

## §4 与精确对角化的对照

ED 参照取 OSC 数据集的 `validation/comparison_4x4_ED_vs_VMC.md`（4×4 环面、$U/t=8$、Ku-ku 的精确对角化、主办方归一化），五个填充点的 $P_d$ 收在 `ZeroGPU/results/osc_4x4_ed.json`。

| $n$ | $\mu$（$N=6$） | $\lambda_d$（$N=4$） | $\lambda_d$（$N=6$） | $V_c$（$N=6$） | ED $P_d$ |
|--:|--:|--:|--:|--:|--:|
| 0.500 | $-2.0$ | 125.16 | 0.294 | 3.396 | 0.04808 |
| 0.625 | $-1.0$ | 125.16 | 125.30 | 0.0080 | 0.06012 |
| 0.750 | $0.0$ | 250.13 | 222.48 | 0.0045 | 0.00983 |
| 0.875 | $0.0$ | 250.13 | 222.48 | 0.0045 | 0.01922 |
| 1.000 | $0.0$ | 250.13 | 222.48 | 0.0045 | 0.00620 |

**读数**：
- $\lambda_d$ 随填充单调升向半填充（$n=1$ 处为 van Hove 峰），$N=4$ 与 $N=6$ 方向一致；
- ED 的 $P_d$ 在半填充处是五个点里**最小**的（0.0062），在 $n=0.625$ 处最大（0.0601）；
- 两种簇尺寸上 $\rho_{\rm Sp}=-0.80$：$\lambda_d$ 的排序与 ED 的排序相反。

## §5 记账

| 项 | 类 | 说明 |
|:--|:--|:--|
| 单粒子谱 $\varepsilon_k$、形状因子 $g_k$ | 结构选择（I） | 方格紧束缚、d-wave 通道为选定的舞台 |
| 配对强度 $V$ | **具名条件（C）** | 本篇不取 $V$ 的具体值；只报告阈值 $V_c=1/\lambda_d$ |
| $\mu(n)$ | 由填充定 | 二分法解占据数 |
| $\eta=10^{-3}t$ | 数值正则化 | 非物理展宽 |
| 簇尺寸 $N$ | 结构选择（I） | 4×4 与 6×6 两档 |

## §6 边界

本篇成立的只是 (2.1) 的代数式与 §4 的对照数字；"$\lambda_d$ 可以当作 OSC 的 $A_d$ 估计量"这一假设**在本对照点上不成立**（$\rho_{\rm Sp}=-0.80$）。这不排除换一个单粒子谱（例如已经含自能的谱）后同一公式会给出不同排序——那是另一篇的问题，本篇不作答。

## 自检

- 结论只引脚本与数据文件，不引外部文献；
- 负结果以数字陈述（$\rho_{\rm Sp}=-0.80$、五个对照点），不写方向性断语；
- $V$ 按具名条件（C）记账，未把阈值 $V_c$ 当成能隙；
- 脚本可重跑：`python3 ZeroGPU/src/zero_dwave_estimate.py 4` 与 `... 6`。
