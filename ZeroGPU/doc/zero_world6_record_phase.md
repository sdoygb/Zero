# 从记录间隔读相位：方法成立，缺的是观测到的 L2 沉积流

**性质**：讨论稿/方法验证，不占正式篇号，不改变根目录 Z 篇的任何判定。
**脚本**：`ZeroGPU/src/zero_world6_record_phase.py`。
**数据**：`results/zero_world6_record_phase.json`。
**复跑**：

```bash
cd ZeroGPU/src && /usr/bin/python3 zero_world6_record_phase.py
```

---

## 摘要

旧库（`L2_imprint_verdict.md`）给出：记录层的沉积流间隔为

$$
\underbrace{2,2,\dots,2}_{T/2-1},\ \mathbf{T+2},\
\underbrace{2,2,\dots,2}_{T/2-1},\ \mathbf{T+2},\dots
\tag{0.1}
$$

由此：

$$
T=\text{long}-2,\qquad
k=\text{现在}-\text{最近一次沉积},\qquad
t_{\rm remain}=\Big(1-\frac{k}{T}\Big)T_{\rm cycle}.
\tag{0.2}
$$

$$
\boxed{\;\text{方法成立（T 可反解、k 可读）；但没有观测到的 L2 沉积流，}k\text{ 仍未知。}\;}
\tag{0.3}
$$

---

## §1 验证

| 检验 | 结果 |
|:--|:--|
| 间隔结构反解 $T=4$ | PASS，$T_{\rm rec}=4$，short $=2.00$ |
| 加噪声（±0.05）仍反解 $T=4$ | PASS |
| 相位 $k=$ 现在 $-$ 最近沉积（示例 $k=1$） | PASS |
| 剩余 $=(1-k/T)\cdot14.51$ Gyr（示例 $k=1$） | PASS，$10.88$ Gyr |
| 没有 L2 沉积流观测 $\Rightarrow k$ 未知 | 方法限制 |

---

## §2 读法

- **$T$ 可读**：长间隔减 2；短间隔为 2；
- **$k$ 可读**：只要有"现在"和"最近一次沉积"两个时间点；
- **剩余**：$(1-k/T)T_{\rm cycle}$；
- **缺的是什么**：**观测到的 L2 沉积流**（记录层的闭合事件序列）。
  我们的引擎有内部 record，但那是模拟的，不是观测。

---

## §3 判定

1. **方法成立**：记录间隔给 $T$ 与 $k$，不需要宇宙学锚、不需要 CP 相位。
2. **$k$ 仍未知**：没有观测到的 L2 沉积流。
3. **这是正确的微观锚方向**（旧库 `L2_imprint` 的路线），比宇宙学/CP 更符合 Zero。
4. **下一步**：找出哪一个观测序列扮演 L2 沉积流（候选：CMB 声学峰序列、
   粒子谱的间隔、或某个记录型观测），用它的间隔结构读 $T$ 与 $k$。

---

## §4 边界

1. 合成流是按旧库规则生成的，不是实测；
2. 加噪声只到 ±0.05；
3. 不是 Z 篇；无新量纲锚。

---

## §5 判据（5/5 PASS，含"方法限制"）

| # | 判据 | 结果 |
|--:|:--|:--|
| 1 | 间隔反解 $T=4$ | PASS |
| 2 | 加噪声仍反解 $T=4$ | PASS |
| 3 | $k=$ 现在 $-$ 最近沉积 | PASS |
| 4 | 剩余公式 | PASS |
| 5 | 无观测流 $\Rightarrow k$ 未知 | 已知限制 |

---

## §6 复跑与文件

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_world6_record_phase.py     # 5 条判据，<1 s
```

| 文件 | 内容 |
|:--|:--|
| `src/zero_world6_record_phase.py` | 沉积流生成、$T/k$ 反解 |
| `results/zero_world6_record_phase.json` | $T$、$k$、剩余、判据 |

---

*ZeroGPU 讨论稿 · zero_world6_record_phase · 2026-10-19*
