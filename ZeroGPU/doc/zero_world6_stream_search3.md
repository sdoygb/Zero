# 试脉冲星/FRB/CMB 极化 + 引擎 record：都不是 L2 沉积流

**性质**：讨论稿/否证，不占正式篇号，不改变根目录 Z 篇的任何判定。
**脚本**：`ZeroGPU/src/zero_world6_stream_search3.py`。
**数据**：`results/zero_world6_stream_search3.json`。
**复跑**：

```bash
cd ZeroGPU/src && /usr/bin/python3 zero_world6_stream_search3.py
```

---

## 摘要

判据同前：间隔两值 $\{s,L=T+2\}$、$s$ 游程 $=T/2-1$。

$$
\boxed{\;\text{脉冲星（0.23）、FRB（0.44）、CMB 极化（游程不对）都不是；}
\text{引擎 record 每代都有闭合（间隔全 1）。}\;}
\tag{0.1}
$$

| 候选 | $n$ | $s$ | $L$ | score | 游程对？ |
|:--|--:|--:|--:|--:|:--|
| 脉冲星 glitch（年） | 14 | 1 | 6 | 0.23 | 偶然对 |
| FRB 爆发（MJD） | 10 | 1 | 74 | 0.44 | ❌ |
| CMB 极化峰 | 9 | 270 | 320 | 1.00 | ❌（游程 ~1） |
| **引擎 gen_closures** | 3000 | 1 | 1 | 0.00 | ❌（每代都有） |

---

## §1 读法

- **脉冲星/FRB/CMB 极化**：间隔不规则或游程不对，都不是 L2 沉积流；
- **引擎 record**：`gen_closures` **每代都非零**（间隔全 1）——不是旧库
  `L2_imprint` 的 $2,2,\dots,T+2$ 模式；
- 这说明：**当前引擎的闭合流与旧库 L2_imprint 的沉积流不是同一个对象**；
  旧库的模式属于旧的闭合模型，当前引擎的 record 是"每代都有闭合"的稠密流。

---

## §2 判定

1. **三条天文候选都失败**：脉冲星/FRB/CMB 极化都不是 L2 沉积流。
2. **引擎 record 也不是**：每代都有闭合，间隔全 1，没有两值游程。
3. **关键发现**：旧库的 $2,2,\dots,T+2$ 模式在当前引擎里**没有出现**；
   要用记录间隔读相位，先得让引擎的闭合流复现旧库的稀疏模式（或找到
   旧库模式对应的引擎变量）。
4. **$k$ 仍未知**。

---

## §3 下一步

- **先修引擎**：让 `gen_closures` 变成稀疏的（旧库模式），或找到引擎里
  对应"沉积事件"的变量（不是每代闭合）；
- 或**承认**：L2 沉积流不是直接天文序列，也不是当前引擎的稠密闭合流；
  需要重建记录层模型；
- 再试候选：Ia 超新星（红移）、CMB 温度谱的峰、脉冲星计时残差。

---

## §4 边界

1. 数据为近似公开值；
2. 脉冲星/FRB 的"事件"选择本身有主观性；
3. 引擎 record 的稠密性是模型实现，不一定是物理；
4. 不是 Z 篇；无新量纲锚。

---

## §5 判据（4/4 PASS）

| # | 判据 | 结果 |
|--:|:--|:--|
| 1 | 脉冲星不是（score 0.23） | PASS |
| 2 | FRB 不是（score 0.44） | PASS |
| 3 | CMB 极化不是（游程不对） | PASS |
| 4 | 引擎 record 也不是（间隔全 1） | PASS |

---

## §6 复跑与文件

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_world6_stream_search3.py     # 4 条判据，约 2 s
```

| 文件 | 内容 |
|:--|:--|
| `src/zero_world6_stream_search3.py` | 脉冲星/FRB/CMB 极化 + 引擎 record |
| `results/zero_world6_stream_search3.json` | 各候选 $s,L,$score、判据 |

---

*ZeroGPU 讨论稿 · zero_world6_stream_search3 · 2026-10-19*
