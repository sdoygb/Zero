# zerogpu —— Zero 的 L0–L2 向量化 / GPU 计算层

用 **OpenCL（AMD RX 570）＋ 向量化 numpy** 把 Zero 的 L0–L2 推到文档最远处的数倍距离。
**完整结论见 [`REPORT_L0_L2.md`](REPORT_L0_L2.md)**（含层指标、更正清单与边界）。

---

## 为什么需要这一层：先量交叉点，再分工

本机是 Intel Mac ＋ AMD RX 570：**没有 CUDA，没有 MPS，没有 MLX**，唯一 GPU 通路是 OpenCL。实测三类负载：

| 负载 | GPU | CPU(12 核) | 倍数 |
|:--|--:|--:|--:|
| 朴素 OpenCL matmul 1024³ | 6.5 GFLOP/s | **67.4 GFLOP/s**（NumPy BLAS） | **CPU 快 10×** |
| 位并行 `uint64` 字操作 | **6 604 Mword/s** | 1 108 Mword/s | **6.0×** |
| 大量独立分支整数步进 | **165.3 Gstep/s** | 5.01 Gstep/s | **33×** |

$$
\Longrightarrow\ \textbf{稠密线性代数留给 CPU};\quad \textbf{GPU 只用于位并行枚举与全分支扫描}。
$$

Zero 的 L0（`uint64` 位掩码上的旋转/交换/计数）与 L2（海量互不相干的分支推进）恰好落在 GPU 赢的那一侧。

---

## 文件

| 文件 | 作用 | 产出 |
|:--|:--|:--|
| **`z0_local.py`** | ⭐⭐ **局部异步毁灭 ＋ Kac 寿命** $\tau_v{=}2\lvert E\rvert/\deg v$：正则⟹均匀；非正则⟹活动-度反相关 | `results/z0_local.json` |
| **`z0_core.py`** | ⭐⭐⭐ **Z0 核心（忠实版）**：散度场状态、全分支零权重、终端层、区域间补偿 | `results/z0_core.json` |
| **`l2_selforg.py`** | ⭐⭐ **无参数自反馈**：`CAP2`（D222 保留规则）给出 18.8 个区域 + 生灭循环 + 涌现长度尺度 | `results/l2_selforg.json` |
| **`l2_multi.py`** | ⭐⭐ **多区域**：局部历史反馈下的 L2 相图（1 → 4–7 → 1 区域） | `results/l2_multi.json` |
| **`z0_loop.py`** | ⭐⭐ **零参数自持环**：闭合→记录→毁灭/重播种→回灌；$r^*=1$ 精确临界 | `results/z0_loop.json` |
| **`charged_sector.py`** | ⭐⭐ **任务 2：带电扇区普查 ＋ 带电局域态**（偶极束缚、局域长度） | `results/charged_sector.json` |
| **`r23_check.py`** | ⭐ **任务 3：独立复核 R23 的 $S_D^{\rm evo}=a_L^D$ 与 $F_D$ 的四维窗口** | `results/r23_check.json` |
| **`observable_sweep.py`** | ⭐⭐ **任务 1：中性普查**（对象×算子×读出的机械枚举，只做结构分类） | `results/observable_sweep.json` |
| **`structures.py`** | ⭐ 结构分类器（幂律/指数/对数/delta/带隙/van Hove；**不看物理名字**） | — |
| **`make_catalogue.py`** | 生成 [`OBSERVABLES.md`](OBSERVABLES.md)（按结构聚类的 133 条目录） | `OBSERVABLES.md` |
| **`l2_regions.py`** | ⭐ **L2 分域（亚稳态/慢模）＋ $d_s$ 估计器的环面标定** | `results/l2_regions.json` |
| **`l2_excursions.py`** | ⭐⭐ **远足分解 ＋ 碰撞图**：远足靠"落到同一顶点"相互作用，连通块 = 区域 | `results/l2_excursions.json` |
| **`self_interaction_scaling.py`** | ⭐⭐ **两条世界线是否必然相互作用**（$\iff D\le4$，$D{=}4$ 临界） | `results/self_interaction.json` |
| **`z0_universe.py`** | ⭐⭐ **零参数宇宙生成器**（全 GPU）：从 ∅ 到年龄 30（体积 10.7 亿）**0.46 s** | `results/z0_universe.json` |
| **`z0_mass_spectrum.py`** | ⭐ 区域系综：远足分解精确核验 ＋ **质量谱**（有无极点）＋ 宇宙学标度 | `results/z0_mass_spectrum.json` |
| **`z0_genesis.py`** | ⭐ **零参数最小生成器（宇宙初开）**：从零态出发，全分支，闭合只记录/或吸收。GPU bitset，深度 32 用 0.08 s | `results/z0_genesis.json` |
| **`z0_record.py`** | ⭐ **记录流落到 L0 项链**：覆盖率、单射性、手征结构、沉积时间线（到 $n=32$） | `results/z0_record.json` |
| `z0_mechanism.py` | 验证「被打到 ⟺ 剖面有只访问一次的高度」判据（穷举 $n=8..18$，退出码 0） | — |
| **`z0_record_gpu.py`** | ⭐ 把 expand/canon 搬到 GPU：**12.45 s → 1.08 s（11.5×）**，canon 单项 310× | `results/z0_record_gpu.json` |
| `minimal_audit.py` | 对称群审计（$\lvert\mathrm{Aut}\rvert$ 有限）＋ 独立复核 G89 命题 1 ＋ 重数审计 | `results/minimal_audit.json` |
| `z0_graph.py` | ⭐⭐ **补回 Z1 的图 Γ**：前沿速度/整数重数/色散/谱维数四项核验 | `results/z0_graph.json` |
| `z0_bound_states.py` | ⭐ **粒子机制**：Γ 上缺陷束缚态 ＋ 缺陷间指数衰减相互作用 | `results/z0_bound_states.json` |
| `zcl.py` | OpenCL 引擎：设备选择、内核缓存、缓冲与计时 | — |
| `l0_closure.py` | L0 闭链图：平衡项链枚举（Burnside 核验）＋ 循环对换建图 ＋ 度分布 | — |
| `run_l0_scan.py` | 把闭链图扫到 $T=32$ | `results/l0_scan.json` |
| `dual_graph.py` | **双图对照**：E1 线性词图 vs 闭链循环项链图（谱） | `results/dual_graph.json` |
| `l0_geometry.py` | **真谱维数**（完整谱 / Lanczos 求积）＋ 增长维数到 $T=24$ | `results/l0_geometry.json` |
| `q_L.py` | 类权重 $q_L=\sum_C\omega_C^2$ 推到 $L=32$ | `results/q_L.json` |
| `l2_branch.py` | L2 全分支动力学（首次通过闭合词枚举）到 $n=32$ | `results/l2_branch_n32.json` |
| `l2_capacity.py` | 容量/寿命/记忆 三维扫描（**注意：三维全是塞进去的**） | `results/l2_capacity.json` |
| `analyze_l0.py` | 备用分析（含较慢的完整 §A/§B/§C） | `results/l0_analysis.json` |

> ⚠ **先读 [`UNIVERSE.md`](UNIVERSE.md)**（零参数宇宙生成器跑出来的是什么）、**[`EVAL_four_theories.md`](EVAL_four_theories.md)**（最小系统能否生成四大理论）、**[`EVAL_minimality.md`](EVAL_minimality.md)** 与 **[`RECORD_STREAM.md`](RECORD_STREAM.md)**：它逐项列出前几轮脚本里**被塞进去的 5 样结构**
> （8 个手工种子、寿命、摧毁周期、闭合吸收、容量），其中「闭合即吸收」**与 Z0① 相反**。
> `z0_genesis.py` 才是零参数的最小版本；`z0_record.py` / `z0_mechanism.py` / `RECORD_STREAM.md` 是它的产出分析。

## 交接单

**[`HANDOFF.md`](HANDOFF.md)** —— 本轮做完什么（普查／L2 分域／$D$ 的选择）＋ **还剩三件事**（带电扇区／$\Gamma$ 的来源／洛伦兹锥）＋ 六条坑。**新会话先读这个。**

## 文档（按重要性）

| 文档 | 一句话 |
|:--|:--|
| [`HANDOFF.md`](HANDOFF.md) | **交接单**：本轮做完什么、还剩三件事、六条坑。**新会话先读这个** |
| [`OBSERVABLES.md`](OBSERVABLES.md) | **任务 1 产出**：133 条可观测量按**结构**聚类（§5 才是后验命名） |
| [`L2_REGIONS.md`](L2_REGIONS.md) | **L2 分域 ＋ $D$ 的选择**：远足碰撞图上 $D{=}1$ 裂成两个区域；相互作用必然 $\iff D\le4$ |
| [`Z0_CORE.md`](Z0_CORE.md) | **忠实引擎**：顶点传递⟹不动点；非传递⟹凝聚＋双指数；补偿是全局的 |
| [`MULTI_REGION.md`](MULTI_REGION.md) | **L2 多区域**：局部历史反馈 ⟹ 4–7 个持续区域；$\mathcal G_T$ 上只凝聚不碎裂 |
| [`LOOP.md`](LOOP.md) | **自持环**：$r^*{=}1$ 精确临界、Kac 引理与 1D 首返律双重验证、稳态无极点 |
| [`CHARGED_SECTOR.md`](CHARGED_SECTOR.md) | **带电扇区**：$n{=}32$ 带电 36.94 亿（6.145 倍）；散度字典 600/600；带电局域态与 $\xi$ |
| [`AUDIT_missing_graph.md`](AUDIT_missing_graph.md) | **最重要**：发现最小系统漏了 Z1 的图 Γ，补回后传播/速度/维数/重数全部回来 |
| [`EVAL_focus.md`](EVAL_focus.md) | 我是否在聚焦某些特性；去掉「Γ 均匀」预设后粒子出现 |
| [`UNIVERSE.md`](UNIVERSE.md) | 零参数宇宙生成器跑出来的是什么（物理术语对照） |
| [`EVAL_four_theories.md`](EVAL_four_theories.md) | 最小系统能否生成四大理论；11.5× 加速 |
| [`EVAL_minimality.md`](EVAL_minimality.md) | 前几轮塞进去的 5 样结构 |
| [`RECORD_STREAM.md`](RECORD_STREAM.md) | 记录流落到 L0 项链：覆盖率、判据、手征 |
| [`REPORT_L0_L2.md`](REPORT_L0_L2.md) | 主报告：L0–L2 三层结果 |

## 复现

```bash
cd /Users/oygb/Downloads/lh/zerogpu

python3 l0_closure.py 8 10 12 14 16 18     # 复现文档节点数与平均度（应全部 OK）
python3 run_l0_scan.py 24 26 28 30 32      # 推到 T=32（约 116 s）
python3 dual_graph.py                      # 双图谱对照
python3 q_L.py 32                          # q_L 到 L=32（文档三值应 OK）
python3 l2_branch.py 32                    # 全分支扫描到 n=32
```

**依赖**：**系统 `python3`**（已装 `pyopencl` / `numpy` / `scipy`）。
**不要**用 DSH 自带运行时——它没有 `pyopencl`，会直接失败。

## 纪律

- 本目录全部结论是 **【数值证据】**，不是【导出】。
- 所有断言带**层指标**（L0 / L1 / L1′ / L2）；跨层拼接一处也没有做（`R95` §1 规则 4）。
- 被证伪或被更正的文档断言逐条登记在 `REPORT_L0_L2.md` §6，**不静默修改原文档**。
