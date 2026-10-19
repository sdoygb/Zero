# 给当前相位找锚：φ 依赖“上次毁灭”的标记

**性质**：讨论稿/锚定，不占正式篇号，不改变根目录 Z 篇的任何判定。
**脚本**：`ZeroGPU/src/zero_world6_phase_anchor.py`。
**数据**：`results/zero_world6_phase_anchor.json`。
**复跑**：

```bash
cd ZeroGPU/src && /usr/bin/python3 zero_world6_phase_anchor.py
```

---

## 摘要

周期 $T=14.51$ Gyr（哈勃锚）。当前相位需要一个“周期起点”锚——即“上次毁灭”
发生在什么时候：

$$
\varphi=\frac{t_{\rm age}-t_{\rm event}}{T},\qquad
t_{\rm remain}=T-(t_{\rm age}-t_{\rm event}).
\tag{0.1}
$$

$$
\boxed{\;\text{相位依赖锚：早事件锚 } \varphi\approx0.95\text{、剩余 }0.72\text{ Gyr；}
\text{晚事件锚 } \varphi\approx0.31\text{、剩余 }10\text{ Gyr}.\;}
\tag{0.2}
$$

| 锚（上次毁灭） | $t_0$（Gyr） | $\varphi$ | 剩余 |
|:--|--:|--:|--:|
| 复合 | $0.00038$ | $0.950$ | $0.72$ Gyr |
| 物质-辐射相等 | $0.00005$ | $0.950$ | $0.72$ Gyr |
| **暗能量起始** | $9.3$ | **$0.309$** | **$10.02$ Gyr** |
| 大爆炸（仅对照） | $0$ | $0.950$ | $0.72$ Gyr |

---

## §1 读法

- **早事件锚**（复合/相等/大爆炸）都在 $t\approx0$ 附近，$\varphi\approx0.95$，
  剩余 $\approx0.72$ Gyr（约 7 亿年）——与上一版数值相同，但现在锚是**事件**，
  不是大爆炸；
- **晚事件锚**（暗能量起始 $t\approx9.3$ Gyr）：$\varphi\approx0.31$，
  剩余 $\approx10$ Gyr（约 100 亿年）；
- **相位不是唯一值**：选哪个事件当“上次毁灭”，决定 $\varphi$。

---

## §2 哪个锚更合理？

| 锚 | 作为“毁灭标记”的理由 | 问题 |
|:--|:--|:--|
| 复合 | 早期相变，易测 | 不是“毁灭”，是等离子→原子 |
| 物质-辐射相等 | 早期转折 | 同上 |
| **暗能量起始** | 方程状态的**相变**，像“毁灭/重置” | 仍无 Zero 的直接对应 |
| 大爆炸 | 常规起点 | **Zero 没有大爆炸**，仅对照 |

- 若把“毁灭”理解为**真空/相变**，暗能量起始（$z\sim0.4$）更像一个标记；
  此时 $\varphi\approx0.31$、剩余 $\approx10$ Gyr；
- 若取早期事件，$\varphi\approx0.95$、剩余 $\approx0.72$ Gyr；
- **没有 Zero 内生的锚**：选择是约定，不是推导。

---

## §3 判定

1. **周期已知**：$14.51$ Gyr。
2. **相位可算，但依赖锚**：早事件 $\varphi\approx0.95$，晚事件 $\varphi\approx0.31$。
3. **若取暗能量起始为“上次毁灭”**：当前相位 $\approx0.31$（周期过了约 31%），
   剩余 $\approx10$ Gyr。
4. **若取早期事件**：相位 $\approx0.95$，剩余 $\approx0.72$ Gyr。
5. **没有 Zero 内生的锚**：这是外部约定；要变成预言，需要一个 Zero 能算的标记。

---

## §4 边界

1. 事件时刻取的是 ΛCDM 的近似值（复合 $z\sim1100$、暗能量 $z\sim0.4$）；
2. 周期与锚都是外部输入；
3. 不是 Z 篇；无新量纲锚。

---

## §5 判据（4/4 PASS）

| # | 判据 | 结果 |
|--:|:--|:--|
| 1 | 周期 $14.51$ Gyr | PASS |
| 2 | 复合锚：$\varphi\approx0.95$、剩余 $\approx0.72$ Gyr | PASS |
| 3 | 暗能量锚：$\varphi\approx0.31$、剩余 $\approx10$ Gyr | PASS |
| 4 | 相位依赖锚，无唯一值 | PASS |

---

## §6 复跑与文件

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_world6_phase_anchor.py     # 4 条判据，<1 s
```

| 文件 | 内容 |
|:--|:--|
| `src/zero_world6_phase_anchor.py` | 事件锚与 $\varphi$、剩余 |
| `results/zero_world6_phase_anchor.json` | 全部数值、判据 |

---

*ZeroGPU 讨论稿 · zero_world6_phase_anchor · 2026-10-19*
