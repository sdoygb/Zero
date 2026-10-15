# ZeroGPU —— 忠实还原零和理论，在计算机上长时间演化

**目标**：把零和底层的演化**忠实地**在计算机上跑下去，跑到"长时间"。

**当前状态**：工程侧已通（可连续跑数万代）；**理论侧有三个阻塞**，详见 [`doc/diagnosis.md`](doc/diagnosis.md)。

---

## ★ 工作流纪律（先读这一节）

```
ZeroGPU/          ← 讨论稿、脚本、结果：**在这里写**
Zero 根目录/       ← 正式稿件：**只有用户明确下令才能新开文章**
```

| 目录 | 放什么 | 谁决定 |
|:--|:--|:--|
| **`ZeroGPU/`** | 讨论稿、草稿、诊断、脚本、结果文件 | 可自由添加、修改、重写 |
| **`Zero/`（根）** | 只放**确认过的正式稿件** | **必须用户明确下令** |

$$
\boxed{\;\textbf{没用户命令，不在 Zero 根目录新开文章};\ \text{讨论一律在 ZeroGPU 内进行} .\;}
$$

**已并入正式稿件的讨论稿**（保留在 `doc/`，作调试轨迹）：

| 讨论稿 | 已并入 |
|:--|:--|
| `diagnosis.md`、`progress.md`、`closure_verdict.md`、`longterm_eval.md` | 根目录 `Z020` |
| `gamma_gap.md`、`gamma_origin_verdict.md`、`dict_verdict.md`、`lattice_dict_verdict.md`、`coupling_verdict.md`、`coupling_map.md`、`scale_alignment_verdict.md`、`coupling_recheck.md`、`matching_verdict.md` | 根目录 `Z021` |
| `capacity_destruction.md`、`saturation_law_rediscovered.md`、`distribution_engine.md` | 根目录 `Z022` |
| `anchoring_layer.md` | 根目录 `Z023` |

$$
\Longrightarrow\quad \text{讨论稿}\,\textbf{不删}\text{（含每次出错与更正）};\ \text{正式稿只保留结论与更正条款} .
$$

---

## 一页结论

$$
\boxed{\;
\underbrace{\text{四个实现级阻塞}}_{\text{已全部修复}};\quad
\underbrace{\text{忠实口径}}_{\text{已钉死}};\quad
\underbrace{\text{维数字典}}_{\text{已裁决}};\quad
\underbrace{\text{格点族字典}}_{\text{已裁决}};\quad
\underbrace{\Gamma\ \text{与耦合方式}}_{\text{未解决·理论}} .
\;}
$$

| 项 | 状态 |
|:--|:--|
| **溢出阻塞** | ✅ **已修复**：任意精度整数；$2\times10^4$ 代后位宽 $\le20$ bits |
| **活动层被丢弃** | ✅ **已修复**：存活走者跨代保留（此前"先收后丢"，活动层恒空） |
| **封顶抹平结构** | ✅ **已修复**：改用**度相关封顶**，Gini 由 $0$ 升到 $0.097$ |
| **吞吐随代数衰减** | ✅ **已修复**：$O(g^2)$ 求和改增量，$100000$ 代仍 $2257$ 代/秒 |
| **长时间演化** | ✅ **已通**：$2000$–$5200$ 代/秒；一天约 $1.7\times10^8$ 代 |
| **有界＋有结构** | ✅ **已达**：$20000$ 代有上界；摆动 $0.38$ 个数量级；Gini $0.097$ |
| **图的来源** | ✅ **已唯一刻画**：$\Gamma=K_4$（4 顶点、度 3、6 边、3 对边对），全部由前轮链条推出。详见 [`doc/gamma_origin_verdict.md`](doc/gamma_origin_verdict.md) |
| **锚定层** | ✅ **已建成**：反问题接口 ＋ 三档自动标注（定理 $5$／条件 $2$） |
| **字典裁决** | ✅ **已裁决**：$D=m-1$ 是代数唯一解；四维窗 $(1/2,3/5)$ 正确。详见 [`doc/dict_verdict.md`](doc/dict_verdict.md) |
| **卡点二** | ✅ **已解除**：正确字典下 $D=4\Rightarrow\deg=4$；已构造 $|C|=26$、围长 $5$ 的例子 |
| **格点族字典** | ✅ **已裁决**：格点无独立字典；字典是通道结构的性质（$D=m-1$）。详见 [`doc/lattice_dict_verdict.md`](doc/lattice_dict_verdict.md) |
| **耦合方式** | ✅ **已给规范构造**：正则 $D$-单纯形，$\rho=-1/D$ 在等角族内唯一。详见 [`doc/coupling_verdict.md`](doc/coupling_verdict.md) |
| **耦合映射** | ✅ **已给出**：独立身份$\leftrightarrow$格点轴 $\Rightarrow$ 空间维数 $=3$（唯一解，两假设同解）；夹角 $=$ 四面体角。详见 [`doc/coupling_map.md`](doc/coupling_map.md) |
| **匹配假设** | ⚠ **已降级并检验**：$\binom D2=D$ 未消掉，降为原理 $M$（两侧同一）；对称群检验一致（$S_4$ 两侧同群）。详见 [`doc/matching_verdict.md`](doc/matching_verdict.md) |
| **规模＋层对齐** | ✅ **同时解决**：$\Gamma=$ **Heawood 图**（$|V|=14$、度 $3$、围长 $6$、二部）；3-边着色给 3 个轴。详见 [`doc/scale_alignment_verdict.md`](doc/scale_alignment_verdict.md) |
| **换 $\Gamma$ 后重检** | ✅ **通过**：身份/关系/维数三项全保留；并得**第二项独立证据**证明"两侧分开"必然。详见 [`doc/coupling_recheck.md`](doc/coupling_recheck.md) |
| **长时间演化评估** | ✅ 一问（工程）通：50 万代、无溢出；⚠ 二问（物理）部分：**非周期演化被公理排除**（定理）。详见 [`doc/longterm_eval.md`](doc/longterm_eval.md) |
| **触顶即毁灭** | ✅ 机制建成、毁灭确实发生；⚠ 有限 $K$ 下周期是 $K$ 的产物；**一致性条件 $|C|=2(T+1)$ 被否证**（命中 $1/5$）⟹ $T$ 目前无确定值。详见 [`doc/capacity_destruction.md`](doc/capacity_destruction.md) |
| **分布机（$K$ 消除）** | ✅ 状态改为"每点分布" ⟹ **自动有界、$K$ 消失**（"$T$ 依赖 $K$"问题消除）；⚠ 但无饱和 ⟹ $T$ 仍无读出量。详见 [`doc/distribution_engine.md`](doc/distribution_engine.md) |
| **饱和律（查旧目录所得）** | ✅ 既有律 $c_T=\binom{T}{\lfloor T/2\rfloor}\le2(T+1)\Rightarrow T\le5$；分布机自然容量 $4$ 与它**自动一致**（$T=1$，落在存活集内）。详见 [`doc/saturation_law_rediscovered.md`](doc/saturation_law_rediscovered.md) |
| **$T$ 的唯一解** | ✅ **$T=L=4$**：饱和判据用**词数** $|W_T|\le2(T+1)\Rightarrow\{2,4\}$；四维窗再收到 $\{4\}$。冲突消解（我上版 $T=5$ 系量纲错配，已作废）。详见 [`doc/algebraic_capacity.md`](doc/algebraic_capacity.md) |
| **遗留** | ⚠ ① $T$ 由**容量＋四维窗**定，**动力学不饱和**（$|W_4|=6\le10$）；② $M$／最小性仍是申报；③ 接口未导出 |
| **忠实口径** | ✅ **已钉死**：闭合 = 回到起点（唯一忠实读法）；详见 [`doc/closure_verdict.md`](doc/closure_verdict.md) |

**三个实现级根因的详情**：[`doc/progress.md`](doc/progress.md)

## 目录

```
ZeroGPU/
├── README.md              本文件
├── src/
│   ├── zero_sim.py        ★ 模拟器公共核心：结构（K₄/Heawood/方格/双覆盖）＋图算法
│   │                        （girth/spectrum/aut_order/精确整数幂）＋词层 ＋耦合几何；
│   │                        92 条判据自检 → results/zero_sim.json。**新脚本从它 import**
│   ├── z0_engine.py       ★ 忠实核心引擎（64 位安全版）＋ 自检
│   ├── closure_reading.py ★ 闭合口径判定（穷举 ＋ 引擎对照）
│   ├── gamma_gap.py       Γ 约束检验 ＋ 差异扫描 ＋ 穷举对照
│   ├── gamma_independence.py  Γ-无关性检验（哪些结论可作定理）
│   ├── dict_verdict.py    维数字典裁决（13 条断言）
│   ├── lattice_dict.py    格点族字典裁决（7 条断言）
│   ├── coupling_verdict.py 耦合方式裁决（16 条断言）
│   ├── coupling_map.py    耦合映射核验（15 条断言）
│   ├── gamma_origin.py    Γ 原生来源核验（16 条断言）
│   ├── matching_verdict.py 匹配假设裁决（10 条断言）
│   ├── scale_alignment.py 规模＋层对齐裁决（15 条断言）
│   ├── coupling_recheck.py 换 Γ 后重检（17 条断言）
│   ├── longterm_eval.py   长时间演化评估（10 条断言）
│   ├── capacity_destruction.py 触顶即毁灭（5 条断言）
│   ├── distribution_engine.py 分布机（7 条断言）
│   ├── algebraic_capacity.py 代数容量机（7 条断言）
│   ├── TL_verdict.py      T 与 L 判定（3 条断言）
│   ├── faithful_reseed.py 归一化重播种检验（4 条断言）
│   ├── anchoring_test.py  锚定可行性（6 条断言）
│   ├── anchor.py          反问题：测量 ⟹ 相容 Γ；正问题：Γ ⟹ 预言/条件
│   ├── anchor_runner.py   家族扫描 ＋ 可观测量字典 ＋ 自动三档标注
│   └── diagnose.py        阻塞诊断（可复现）
├── doc/
│   ├── zero_sim.md        ★ 模拟器公共核心说明（三层接口 ＋ 92 条判据账 ＋ 两个易错口径）
│   ├── algebraic_capacity.md ★ 接通两个 T：T=5 唯一（收窄既有 {5,6}）
│   ├── saturation_law_rediscovered.md 饱和律早已存在：分布机与之一致
│   ├── distribution_engine.md 分布机：K 消除，T 仍缺
│   ├── capacity_destruction.md 触顶即毁灭：机制建成，一致性条件被否证
│   ├── longterm_eval.md   长时间演化评估：工程通过，非周期性被公理排除
│   ├── coupling_recheck.md 换 Γ 后重检：通过，且得第二项独立证据
│   ├── scale_alignment_verdict.md 规模＋层对齐：Heawood 图
│   ├── matching_verdict.md 匹配假设裁决：降为原理 M 并通过对称群检验
│   ├── gamma_origin_verdict.md Γ 的原生来源：唯一刻画为 K_4
│   ├── coupling_map.md    耦合映射：空间维数 = 3 被计数强制
│   ├── coupling_verdict.md 耦合方式：正则 D-单纯形（等角族内唯一）
│   ├── lattice_dict_verdict.md 格点族字典裁决：格点无独立字典
│   ├── dict_verdict.md    维数字典裁决：D=m−1（代数唯一解）
│   ├── anchoring_layer.md 锚定层：反问题、三档标注与两个卡点
│   ├── gamma_gap.md       Γ 缺口说明书（需要什么／为何给不出／给了能多算什么）
│   ├── closure_verdict.md ★ 忠实口径判定（含词/行走层差的新发现）
│   ├── progress.md        ★ 三个实现级根因与"有界＋有结构"的达成
│   └── diagnosis.md       初版诊断：三重阻塞（含溢出复现）
├── results/
│   └── diagnose.json      诊断输出
└── logs/                  运行日志
```

---

## 怎么跑

**必须用系统 python**（有 `numpy`；`pyopencl` 用于后续 GPU 位并行枚举）：

```bash
cd ZeroGPU/src
/usr/bin/python3 z0_engine.py        # 自检：24 条断言（含口径核验）
/usr/bin/python3 closure_reading.py  # 口径判定：8 条断言
/usr/bin/python3 anchor.py            # 锚定：反问题 + 三档标注
/usr/bin/python3 gamma_independence.py  # Γ-无关性：8 条断言
/usr/bin/python3 gamma_gap.py        # Γ 缺口：13 条断言
/usr/bin/python3 diagnose.py         # 诊断：运行模式 + 溢出复现 + 封顶扫描
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

**B 模式是唯一"稳定＋有结构"的模式**——它要求度不均匀。

---

## 下一步（按优先级）

| 优先 | 任务 | 类型 |
|--:|:--|:--|
| ~~1~~ | ~~钉死"忠实"的口径~~ | ✅ **已完成**（闭合 = 回到起点） |
| ~~2~~ | ~~造出 $\Gamma$ 的原生来源~~ | ⚠ **已精确量化，仍需新条款**：三条路线已排除，见 `doc/gamma_gap.md` |
| **3** | 把"中间容量"作为**显式具名输入**引入，测它能否同时给出稳定与结构 | 理论＋实现 |
| **4** | 移植局部异步毁灭（Kac 寿命 $\tau_v=2|E|/\deg v$） | 实现 |
| **5** | 接 GPU：位并行枚举与全分支扫描（实测 GPU 在此类负载上快 $6$–$33$ 倍） | 实现 |

---

*ZeroGPU · 本目录只做一件事：让零和宇宙在计算机上忠实地跑下去。*
