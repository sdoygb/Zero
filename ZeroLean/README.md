# ZeroLean（Zb 新库）— 机器可检查核心

本目录把 `Zb001–Zb065` 中**能被 Lean 内核逐行检查**的那部分命题写成 Lean 4 定理。
它不是「把 65 篇的每一句话都形式化」——那需要 850 条定理头、数周工作量；
它给的是**逐篇的机器验证状态**：哪些篇有 kernel 定理、哪些篇是具名输入、哪些篇依赖数值/拟合。

## 现状（2026-10-31 机器检查）

| 项 | 值 |
|:--|:--|
| `lake build` | **通过**（9075 jobs，约 34 s 增量） |
| 审计定理（含辅助） | **674 条** |
| `zbNNN_*` 具名定理 | **295 条** |
| 有 kernel 定理的篇 | **31 / 68** |
| CI 守则（R1／R2／R3） | **全过** |
| 工具链 | Lean `v4.35.0`（`~/.local/lean4`），Mathlib 已接入 |
| 构建产物 | `ZeroLean/.lake` → `/Volumes/download/lean-build/ZeroLeanZb-lake`（复用旧库的 Mathlib 缓存） |

逐篇状态见 **`_staging/FORMALIZATION_MATRIX_Zb.md`**（机器生成）与
**`_staging/FORMALIZATION_COVERAGE_Zb.md`**（覆盖账）。

## 三分归属（不是「证了/没证」二分）

| 状态 | 含义 | 本轮篇数 |
|:--|:--|--:|
| ① 已证（部分） | 该篇模块里至少有一条 `zbNNN_*` kernel 定理 | 31 |
| ② 具名输入/诠释 | 主结论是建模选择或诠释，已按具名输入记账 | 10 |
| ③ 数值/拟合 | 定量结论依赖实测或拟合参数，不是可脱离数据的定理 | 7 |
| ⏳ 待形式化 | 有可形式化的数学核，本轮未做 | 20 |

**⏳ 不是「已证」**。① 也只覆盖该篇的**部分**命题——每篇的完整定理头是 850 条量级，
本轮做的是核心切片，不是全篇。

## 怎么跑

```bash
export PATH="$HOME/.local/lean4/bin:$PATH"
cd /Users/oygb/Downloads/Zero/ZeroLean

lake build                         # 全绿 = 所有定理被内核接受
python3 scripts/check.py           # CI：build + 自动审计 + 三条守则
python3 scripts/check.py --coverage # 额外打印逐篇覆盖
python3 scripts/coverage_matrix.py  # 重生成 _staging/FORMALIZATION_MATRIX_Zb.md
```

插件侧：`zero_lean` 调用的就是 `scripts/check.py`；改过 `ZeroLean/*.lean` 后重跑它。

## CI 三条守则

- **R1 禁 `sorryAx`**：证明不许有洞。
- **R2 公理白名单**：只允许 `propext` / `Classical.choice` / `Quot.sound`，以及
  `scripts/axioms-allow.txt` 里登记的 Zero 侧具名输入；未登记的公有理 = 未登记的输入。
- **R3 定理名对 Zb 篇号**：`zbNNN_*` 的 NNN 必须对应库里真实存在的 `ZbNNN` 篇；
  文件名 `ZbNNN_*.lean` 同理。

审计由 `Zero/Audit.lean` 自动遍历环境里的 `Zero.*` 定理，逐条打印依赖公理；
`check.py` 解析这些行。**依赖清单 = 机械版的「要素来源表」**：

| 依赖清单 | 判定 |
|:--|:--|
| 只有 `propext` / `Quot.sound` / `Classical.choice` | 导出，不依赖 Zero 侧输入 |
| 出现自定义 `axiom` / `opaque` / `constant` | 具名输入，点名列出 |
| 出现 `sorryAx` | 证明有洞（CI 禁止） |

## Lean 能查什么、不能查什么

**能查**：定理是否真的从前提推出来；有没有跳步；常数与算式的正确性；
反例是否真的是反例；每条结论到底用了哪些前提。

**不能查**：形式化对象与物理意图是否一致；数值仿真的正确性；
「某要素是内生还是外部嫁接」这类建模判断——但它能把输入的清单摆到面前。

一句话：**Lean 管「从前提推得对不对」，不管「前提选得对不对」。**

## 文件

| 路径 | 内容 |
|:--|:--|
| `Zero.lean` | 根模块，逐一 import 已接入的 `Zero/Zb*.lean` |
| `Zero/Basic.lean` | K₄ 的核/零和超平面（Zb001 的底层） |
| `Zero/Audit.lean` | 审计机关：遍历 `Zero.*` 定理并打印依赖公理 |
| `Zero/ZbNNN_*.lean` | 逐篇形式化（31 篇有定理） |
| `scripts/check.py` | CI：build + 审计 + R1/R2/R3 |
| `scripts/coverage_matrix.py` | 生成 `_staging/FORMALIZATION_MATRIX_Zb.md` |
| `scripts/axioms-allow.txt` | 公理白名单 |
| `wip/` | 未接入构建的旧草稿（`Zb024_Boost.lean` 已被 `Zero/Zb024_BoostV2.lean` 取代） |

## 下一批（按可做性排序）

1. **Zb007 的一般算子论版本**：任意 `ρ` 的酉矩阵组装/对角化/`ρ^{it}` 函数演算（2 级块、一般 `n` 对角、2×2 Hermitian 特征多项式/λ±/特征向量/正交/归一化已证）。
2. **Zb060 的矩阵版**：Klein 不等式/最大熵与 pinching 第二定律（对角 KL 非负/等号、2×2 变分已证）。
3. **Zb067 的开放项**：runs 闭式 (4.4) 与扇区一最低层高一般式（`L=4,6`、`L=8` 旋转已证）。

> 诚实边界：本轮**没有**把 65 篇全部证完；`⏳` 的 27 篇是待办，不是已证。
> 若要把「所有文章」都做到 ①，需要按上表逐批推进，且 ③ 类（数值/拟合）篇
> 原则上不应被写成 kernel 定理——它们的正确位置是「具名输入 + 数值核验」。
