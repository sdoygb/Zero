# ZeroLean — 一个离散物理框架的机器可检查核心

这个仓库把一套离散/格点物理框架（下称 Zero）里**能被 Lean 内核逐行检查**的那部分命题写成 Lean 4 定理。
它不是"把 96 篇都形式化"——那需要 850 条定理头的量级；它给的是**逐篇的机器验证状态**。

## 现状（2026-11-10 机器检查）

| 项 | 值 |
|:--|--:|
| `lake build` | **通过**（9121 jobs） |
| 审计定理（含辅助） | **959 条** |
| 模块 | **75 个** `.lean`，8753 行 |
| 有 kernel 定理的篇 | **42 / 96** |
| `sorryAx` 出现次数 | **0** |
| 公理依赖 | `propext`（136 条只用它）、`Classical.choice`、`Quot.sound`；**60 条零公理** |
| 工具链 | Lean `v4.35.0` + Mathlib（`rev = "v4.35.0"`） |

## 它是什么 / 不是什么

**是**：一批可复算的组合与代数事实，直接从离散出发，不引入连续极限。可复用的部分包括

- **图 Laplacian 的基本事实**：核 = 常向量、半正定、二次型（`Zb012_AuxLaplacian`、`Zb054_Laplacian`）
- **K₄ 与八面体**：线图 = 八面体、$S_4\to S_3$、完美匹配计数（`Zb048_*`、`Zb078_LineGraph`）
- **谱与共扼**：谱事实、共扼对、4d Hodge 星算子（`Zb049_Spectrum`、`Zb076_Conj`、`Zb080_Hodge`）
- **词与计数**：ballot 型双射、runs、composition 计数（`Zb003_*`、`Zb067_*`）
- **有限维复线性代数**：色散恒等式、Hermite 迹、受迫响应矩阵非奇异（`Zb005`、`Zb007`、`Zb014`）
- **手性判据**：自同构群里的镜面元素 $\Rightarrow$ 手征荷为零（`Zb092_Chirality`）

**不是**：

- 不是"这套物理论断已被证明"。Lean 查的是"从给定前提推得对不对"，不是"前提选得对不对"。
- 每篇只形式化了**核心切片**，不是全篇（矩阵见 `_staging/FORMALIZATION_MATRIX_Zb.md`）。
- 物理对象（格点、场、记录的物理解释）都是**建模选择**，只在 Lean 里以公理化形式出现。

## 怎么跑

```bash
cd ZeroLean
lake build                          # 全绿 = 所有定理被内核接受
python3 scripts/check.py            # CI：build + 自动审计 + 三条守则
python3 scripts/check.py --coverage  # 逐篇覆盖
python3 scripts/coverage_matrix.py   # 重生成覆盖矩阵
```

## CI 三条守则

- **R1 禁 `sorryAx`**：证明不许有洞。
- **R2 公理白名单**：只允许 `propext` / `Classical.choice` / `Quot.sound`，以及
  `scripts/axioms-allow.txt` 里登记的具名输入；未登记的公有理 = 未登记的输入。
- **R3 定理名对篇号**：`zbNNN_*` 的 NNN 必须对应真实存在的篇；文件名同理。

审计由 `Zero/Audit.lean` 遍历环境里的 `Zero.*` 定理，逐条打印依赖公理；`check.py` 解析这些行。
**依赖清单 = 机械版的"要素来源表"**：

| 依赖清单 | 判定 |
|:--|:--|
| 只有 `propext` / `Quot.sound` / `Classical.choice` | 导出，不依赖本框架的具名输入 |
| 出现自定义 `axiom` / `opaque` / `constant` | 具名输入，点名列出 |
| 出现 `sorryAx` | 证明有洞（CI 禁止） |

## 文件

| 路径 | 内容 |
|:--|:--|
| `Zero.lean` | 根模块，逐一 import 已接入的 `Zero/Zb*.lean` |
| `Zero/Basic.lean` | K₄ 的核 / 零和超平面 |
| `Zero/Audit.lean` | 遍历 `Zero.*`，打印每条定理的公理依赖 |
| `Zero/AxiomLedger.lean` | 公理账本（具名输入一览） |
| `Zero/ZbNNN_*.lean` | 逐篇定理 |
| `scripts/check.py` | CI：build + 审计 + 守则 R1–R3 |
| `scripts/axioms-allow.txt` | 允许的具名输入清单 |

## 许可

数学内容按 CC BY-NC-SA 4.0，代码按 MIT。

---

问题、纠错、以及"这条其实 Mathlib 已有"的指正都欢迎——尤其是最后一类，能直接减少重复劳动。
