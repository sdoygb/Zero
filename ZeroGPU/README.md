# ZeroGPU —— 数值引擎与仿真

**两件事**：① 把零和底层的演化**忠实地**在计算机上跑下去；② 作为**讨论稿工作区**（讨论稿、脚本、结果都在这里，正式稿才进 Zero 根目录）。

---

## 工作流纪律（先读这一节）

| 目录 | 放什么 | 谁决定 |
|:--|:--|:--|
| **`ZeroGPU/`** | 讨论稿、草稿、诊断、脚本、结果文件 | 可自由添加、修改、重写 |
| **`Zero/`（根）** | 只放**确认过的正式稿件** | **必须用户明确下令** |

$$
\boxed{\;\textbf{没用户命令，不在 Zero 根目录新开文章};\ \text{讨论一律在 ZeroGPU 内进行} .\;}
$$

讨论稿**不删**（含每次出错与更正）；正式稿只保留结论与更正条款。

---

## 现状（2026-10-17 清点）

| 项 | 值 |
|:--|:--|
| `src/` 脚本 | **240 个 `.py`** |
| `results/` 结果文件 | **281 个**（`.json`／`.txt`／`.png`） |
| `doc/` 说明 | **39 份**（索引见下） |
| 公共核心自检 | `zero_sim.py`：**96 条判据全过** → `results/zero_sim.json` |
| 动力学引擎自检 | `z0_engine.py`：**24 条断言全过** |
| 统一引擎自检 | `zero_evo.py`：**20 条全过** → `results/zero_evo_selftest.json` |
| 工程侧 | 连续跑数十万代无溢出、有界且有结构（`doc/progress.md`、`doc/longterm_eval.md`） |

**已确认的底座（有文档与脚本支撑）**：

- **图的来源**：$\Gamma=K_4$（4 顶点、度 3、6 边、3 对边对）由前轮链条唯一刻画（`doc/gamma_origin_verdict.md`）；
- **规模＋层对齐**：$\Gamma=$ Heawood 图（$|V|=14$、度 3、围长 6、二部），3-边着色给 3 个轴（`doc/scale_alignment_verdict.md`）；
- **忠实口径**：闭合＝回到起点（`doc/closure_verdict.md`）；
- **容量**：$T=L=4$（用词数 $|W_T|\le2(T+1)$，再由四维窗收到 $\{4\}$；`doc/algebraic_capacity.md`）；
- **饱和律**：$c_T=\binom{T}{\lfloor T/2\rfloor}\le2(T+1)\Rightarrow T\le5$，与分布机的自然容量 4 自动一致（`doc/saturation_law_rediscovered.md`）；
- **维数字典**：$D=m-1$（代数唯一解），格点族无独立字典（`doc/dict_verdict.md`、`doc/lattice_dict_verdict.md`）；
- **耦合**：正则 $D$-单纯形（等角族内唯一），独立身份↔格点轴 $\Rightarrow$ 空间维数 $=3$（`doc/coupling_verdict.md`、`doc/coupling_map.md`、`doc/coupling_recheck.md`）；
- **锚定层**：反问题接口＋三档自动标注（`doc/anchoring_layer.md`）；
- **统一引擎（讨论稿）**：走者／记忆／内容／势／相位／几何／容量在一个状态里；几何回写反转选择律、双稳与手性行波（`doc/zero_evo.md`）；
- **层账（讨论稿）**：L0=0／L1=$H=L/2$ 静态／L2=3 槽；`persistent` 与 `fresh` 闭合口径分叉，L1 高度只在 `fresh` 下逐位（`doc/zero_layers.md`）；
- **第七条边（判定）**：L0 内容→L2 seed 的无量纲耦合 λ；`persistent` 硬 winner-take-all，`fresh` path16 有准定常 8/16 选择、λ≥0.5 振荡；默认 λ=0（`zero_seventh.py`）。

---

## 目录

```
ZeroGPU/
├── README.md              本文件
├── src/                   183 个脚本
│   ├── zero_sim.py        ★ 公共核心：结构（K₄/Heawood/方格/双覆盖）＋图算法
│   │                        （girth/spectrum/aut_order/精确整数幂）＋词层＋耦合几何；
│   │                        96 条判据自检 → results/zero_sim.json。**新脚本从它 import**
│   ├── z0_engine.py       ★ 忠实核心引擎（64 位安全版）＋ 自检（24 条断言）
│   ├── zero_evo.py        ★ 统一多扇区引擎（12 通道）＋ 自检（20 条）
│   ├── zero_evo_phase.py  相图／迟滞／行波／换图／α 扫描
│   ├── zero_evo_speed.py  行波速度 vs 几何响应速度
│   ├── zero_evo_long.py   4 条 10^6 步长跑＋李雅普诺夫＋时空图
│   ├── zero_evo_bistable.py c_h=1.55 的 FID 支复核
│   ├── zero_evo_linear.py 线性化对照（Fiedler 最快）
│   ├── zero_layers.py     L0/L1/L2 层账、口径分叉、L2 槽、脉冲响应
│   ├── zero_seventh.py    第七条边（L0 内容→L2 seed）的 λ 扫描与判据
│   ├── neutrino_lr_oscillation.py  Z019 左右手真质量、相消与三相位震荡
│   ├── closure_reading.py 闭合口径判定（穷举 ＋ 引擎对照）
│   ├── gamma_gap.py       Γ 约束检验 ＋ 差异扫描 ＋ 穷举对照
│   ├── gamma_origin.py    Γ 原生来源核验
│   ├── dict_verdict.py    维数字典裁决
│   ├── lattice_dict.py    格点族字典裁决
│   ├── coupling_verdict.py／coupling_map.py／coupling_recheck.py  耦合方式与映射
│   ├── matching_verdict.py／scale_alignment.py   匹配假设与规模＋层对齐
│   ├── capacity_destruction.py／distribution_engine.py／algebraic_capacity.py  容量线
│   ├── longterm_eval.py   长时间演化评估
│   ├── anchor.py／anchor_runner.py／anchoring_test.py   锚定层
│   ├── field_sim.py／force_from_field.py／field_to_particle.py   场与力
│   ├── curvature_from_defects_check.py／a2_bridge_check.py   曲率与作用量桥核验
│   └── …（其余 150+ 个按篇的判定脚本，命名与篇号对应）
├── doc/                   39 份说明（索引见下）
├── results/               218 个结果文件（脚本输出，JSON/TXT/PNG）
└── logs/                  运行日志
```

---

## 文档索引（`doc/`，39 份）

| 文档 | 内容 |
|:--|:--|
| `zero_sim.md` | ★ 公共核心说明：三层接口 ＋ 判据账 ＋ 两个易错口径 |
| `zero_evo.md` | ★ 统一多扇区引擎：状态／通道／20 条自检／几何反转选择律／相图与迟滞／手性行波 |
| `zero_layers.md` | ★ L0/L1/L2 层账：L1 高度静态读出、L2 三槽、persistent/fresh 口径分叉、脉冲响应 |
| `zero_long_audit.md` | ★ 长跑审计：full 5×10^6 步与 layer 10^7 步；形状冻结、进程不卡 |
| `zero_sublayers.md` | ★ L2 亚亚层 16/8/4 分级毁灭 prototype：底全灭、中留 2、顶不动；周期 4 |
| `zero_faithful_sublayers.md` | ★ 忠实引擎里的 16/8/4：能跑、有周期 4，但 q/n_std 过冲，未复现观测 |
| `zero_reseed_derivation.md` | ★ L1 记录 push-forward 的 16/8/4 reseed：周期 4、未毁层继续被播种；相位平均 0.224/0.360/0.416 |
| `zero_routes_compare.md` | ★ A（L=8 词 top2）vs B（K4 纤维）：A 更接近观测（距离 0.150 vs 0.275），两者都有周期 4 |
| `zero_l4_layers.md` | ★ L=4/D=4 下的 8/4/2：周期 4、未毁顶被 reseed；顶 0.395，低于观测 0.552 |
| `zero_window_readout.md` | ★ L=4 vs L=8：漏掉的是读出窗口≠词长；W8 的 180 个平衡 4-窗口 q_rot=5/9 |
| `zero_window_16.md` | ★ 窗口化 16/8/4：4 site×2×2；周期 4；距观测 0.127（目前最好，保 D=4） |
| `zero_window16_engine.md` | ★ 三拷贝 W + flow/reseed：原始份额 [0.170,0.292,0.537]；q 仍偏高 0.31 |
| `zero_window16_aggregate.md` | ★ 聚合一致 1:2:4 读法：q 0.532、中拷贝独大；暴露“聚合 vs 顶不动”矛盾；原始 0.033 是读法 artifact |
| `zero_window16_residence.md` | ★ B 判定：独立三拷贝 + 驻留权重——周期内 relstd 0.36，周期间 0.138，q 0.456，不成立 |
| `zero_window16_A.md` | ★ A 判定：单 W 读出——份额恒 1/3、q 0.396、无周期 4；A/B 都不成立 |
| `zero_window16_merge.md` | ★ A/B 软约束扫描：无甜点；稳定只在 λ=1（1/3），份额距离仍 0.40 |
| `zero_window16_toppias.md` | ★ C′ 顶偏置：λ=1 份额 [0.235,0.235,0.531]，距离 0.125，周期间稳定；q16=0.091（4 顶点 q 不可比） |
| `progress.md` | ★ 三个实现级根因与「有界＋有结构」的达成 |
| `closure_verdict.md` | ★ 忠实口径判定（闭合＝回到起点） |
| `gamma_origin_verdict.md` | Γ 的原生来源：唯一刻画为 $K_4$ |
| `gamma_gap.md` | Γ 缺口说明书（需要什么／为何给不出／给了能多算什么） |
| `scale_alignment_verdict.md` | 规模＋层对齐：Heawood 图 |
| `matching_verdict.md` | 匹配假设裁决：降为原理 $M$ 并通过对称群检验 |
| `dict_verdict.md` | 维数字典裁决：$D=m-1$（代数唯一解） |
| `lattice_dict_verdict.md` | 格点族字典裁决：格点无独立字典 |
| `coupling_verdict.md` | 耦合方式：正则 $D$-单纯形（等角族内唯一） |
| `coupling_map.md` | 耦合映射：空间维数 $=3$ 被计数强制 |
| `coupling_recheck.md` | 换 Γ 后重检：通过，并得第二项独立证据 |
| `algebraic_capacity.md` | 接通两个 $T$：$T=L=4$ 唯一 |
| `saturation_law_rediscovered.md` | 饱和律早已存在：分布机与之一致 |
| `saturation_mechanism.md` | 饱和的真正机制：容量 1 与被选速度（讨论稿） |
| `capacity_destruction.md` | 触顶即毁灭：机制建成，一致性条件被否证 |
| `distribution_engine.md` | 分布机：$K$ 消除，$T$ 仍缺 |
| `longterm_eval.md` | 长时间演化评估：工程通过；非周期演化被公理排除（定理） |
| `anchoring_layer.md` | 锚定层：反问题、三档标注与两个卡点 |
| `rank_alignment.md` | 把「可区分重数」落实为态的秩（讨论稿） |
| `five_integers.md` | 五个未入账整数（$L$／SPAWN／$N$／$K$／$B$）的账目（讨论稿） |
| `B_determination.md` | $B$ 是否内定：被选速度恰等因果锥速度（讨论稿） |
| `diagnosis.md` | 初版诊断：三重阻塞（含溢出复现） |

---

## 怎么跑

**必须用系统 python**（有 `numpy`）：

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_sim.py            # 公共核心自检：96 条判据 → results/zero_sim.json
/usr/bin/python3 zero_evo.py --selftest # 统一引擎自检：20 条 → results/zero_evo_selftest.json
/usr/bin/python3 zero_evo_phase.py      # 相图/迟滞/行波，约 8 min
/usr/bin/python3 zero_evo_speed.py      # 速度曲线，约 1 min
/usr/bin/python3 zero_evo_long.py       # 10^6 步长跑，约 3–5 min
/usr/bin/python3 zero_evo_bistable.py   # 双稳的 FID 支，约 1 min
/usr/bin/python3 zero_evo_linear.py     # 线性化对照，约 1 s
/usr/bin/python3 zero_seventh.py        # 第七条边 λ 扫描，约 1–2 min
/usr/bin/python3 neutrino_lr_oscillation.py # Z019 左右手相消，约 1 s
/usr/bin/python3 zero_layers.py         # L0/L1/L2 层账与脉冲响应，约 2 min
/usr/bin/python3 zero_long_audit.py --case full --steps 5000000   # 长跑审计，约 15 min
/usr/bin/python3 zero_sublayers.py --steps 200                    # 亚亚层 prototype，约 1 s
/usr/bin/python3 zero_faithful_sublayers.py                        # 忠实引擎 16/8/4 尝试，约 2 s
/usr/bin/python3 zero_reseed_derivation.py                         # L1 push-forward reseed，约 1 s
/usr/bin/python3 zero_routes_compare.py                            # A/B 路线对比，约 1 s
/usr/bin/python3 zero_l4_layers.py                                 # L=4 8/4/2，约 1 s
/usr/bin/python3 zero_window_readout.py                            # W8 的 4-窗口 q_rot，约 1 s
/usr/bin/python3 zero_window_16.py                                 # 窗口化 16/8/4，约 1 s
/usr/bin/python3 zero_window16_engine.py                           # 三拷贝 W micro core，约 1 s
/usr/bin/python3 zero_window16_aggregate.py                        # 聚合一致读法，约 1 s
/usr/bin/python3 zero_window16_residence.py                        # B 驻留权重判定，约 1 s
/usr/bin/python3 zero_window16_A.py                                 # A 单 W 读出判定，约 1 s
/usr/bin/python3 zero_window16_merge.py                             # A/B 软约束扫描，约 1 s
/usr/bin/python3 zero_window16_toppias.py                           # C′ 顶偏置扫描，约 1 s
/usr/bin/python3 z0_engine.py           # 动力学引擎自检：24 条断言
/usr/bin/python3 closure_reading.py     # 闭合口径判定
/usr/bin/python3 anchor.py              # 锚定：反问题 + 三档标注
/usr/bin/python3 gamma_independence.py  # Γ-无关性
/usr/bin/python3 gamma_gap.py           # Γ 缺口
/usr/bin/python3 diagnose.py            # 阻塞诊断（运行模式 + 溢出复现 + 封顶扫描）
```

---

## 引擎的四条纪律

| # | 纪律 |
|--:|:--|
| **1** | 状态是**散度场**：$d=$（终点计数）$-$（起点计数）$+$ 终端层 |
| **2** | 步进是**确定性全分支**（$\mathcal M\leftarrow\mathcal M\cdot A$），整数重数，**不用概率加权** |
| **3** | 闭合 $=$ 起点等于终点；闭合者进**记录层**，未闭合者由**终端层**收走，再从记录层**重播种** |
| **4** | 所有计数用**任意精度整数**（永不溢出）；规模只以**对数与位宽**报告 |

---

## 三个运行模式（已复现）

| 模式 | 触发 | 活动层 | 空间结构 |
|:--|:--|:--|:--|
| **A 不动点** | 传递图 ＋ 均匀种子 | 恒定 | Gini $=0$（结构被抹平） |
| **B 有结构** | **非传递图**（如星图） | 有界 | Gini $=0.89$ |
| **C 自持** | 环上的闭路径系统 | 恒定 | 记录层单调累积（时间箭头） |

**B 模式是唯一「稳定＋有结构」的模式**——它要求度不均匀。

---

## 待做（按优先级，详见对应文档）

| 优先 | 任务 | 类型 |
|--:|:--|:--|
| 1 | 把「中间容量」作为**显式具名输入**引入，测它能否同时给出稳定与结构（`doc/gamma_gap.md`、`doc/saturation_mechanism.md`） | 理论＋实现 |
| 2 | 移植局部异步毁灭（Kac 寿命 $\tau_v=2|E|/\deg v$） | 实现 |
| 3 | 接 GPU：位并行枚举与全分支扫描（实测 GPU 在此类负载上快 $6$–$33$ 倍） | 实现 |

---

*ZeroGPU · 本目录只做一件事：让零和宇宙在计算机上忠实地跑下去。*
