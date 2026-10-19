# 用世界模型拟合当前相位：φ≈0.77–0.85，但约束很弱

**性质**：讨论稿/拟合，不占正式篇号，不改变根目录 Z 篇的任何判定。
**脚本**：`ZeroGPU/src/zero_world6_engine_phase.py`。
**数据**：`results/zero_world6_engine_phase.json`、`zero_world6_engine_phase.png`。
**复跑**：

```bash
cd ZeroGPU/src && /usr/bin/python3 zero_world6_engine_phase.py
```

---

## 摘要

用 6 层世界模型（T=4 相对相位调制质量）拟合当前相位：

$$
m_i(\varphi)=m_i^0\exp\!\big(\epsilon\cos(2\pi\varphi+\delta_i)+c\big),
\qquad
\delta_i\in\{0,2\pi/3,4\pi/3\}\ (\text{每代 2 层}).
\tag{0.1}
$$

$m_i^0$ 取 4→6 差映射的预测质量；拟合 $\epsilon,\varphi,c$ 到 6 个观测质量。

$$
\boxed{\;\varphi^*\approx0.77\text{（夸克，}\Delta\chi^2=4.05\text{）},\quad
\varphi^*\approx0.85\text{（轻子，}\Delta\chi^2=2.17\text{）};\quad
\text{剩余}\approx2\text{–}3\ \text{Gyr}.\;}
\tag{0.2}
$$

| 扇区 | $\varphi^*$ | $\chi^2_{\min}$ | $\Delta\chi^2$ | 剩余 |
|:--|--:|--:|--:|--:|
| 6 夸克 | $0.766$ | $0.82$ | $4.05$ | $3.40$ Gyr |
| 6 轻子 | $0.851$ | $1.02$ | $2.17$ | $2.16$ Gyr |

---

## §1 读法

- **两扇区 $\varphi^*$ 接近**（$0.766$ vs $0.851$，差 $0.085$）；
- **但约束很弱**：夸克 $\Delta\chi^2=4.05$（勉强），轻子 $\Delta\chi^2=2.17$（**基本没约束**）；
- 结果依赖 $\delta_i$ 的取法（3 个匹配相位）和调制形式；换 $\delta_i$ 会移 $\varphi^*$；
- 若按面值取：$\varphi\approx0.77$–$0.85$，剩余 $\approx2.2$–$3.4$ Gyr。

---

## §2 与 CP 锚比较

| 锚 | $\varphi$ | 剩余 |
|:--|--:|--:|
| CKM $\gamma$ | $0.182$ | $11.9$ Gyr |
| PMNS $\delta_{CP}$ 0.7π–1.5π | $0.35$–$0.75$ | $3.6$–$9.4$ Gyr |
| **世界模型拟合** | **$0.77$–$0.85$** | **$2.2$–$3.4$ Gyr** |

- 世界模型给的 $\varphi$ 落在 CP 锚的**高端/略超出**；
- 但世界模型的 $\Delta\chi^2$ 只有 2–4，**不比 CP 锚更可靠**；它只是把质量残差拿来定相位，
  而质量残差本身有 ~0.5 dex 的模型误差。

---

## §3 判定

1. **世界模型能拟合出 $\varphi\approx0.8$**，剩余 $\approx2$–$3$ Gyr。
2. **但约束弱**：轻子 $\Delta\chi^2=2.17$ 基本无约束；夸克 $4.05$ 勉强。
3. **不更可靠**：结果依赖 $\delta_i$/调制形式；质量残差 ~0.5 dex。
4. **要更可靠**：需要更多独立观测量（混合角、两扇区联合、轻子与夸克同时），
   或让引擎直接输出相位（不是事后拟合）。

---

## §4 边界

1. $m_i^0$ 来自 4→6 差映射（本身是拟合，~0.9 dex）；
2. $\delta_i$ 的 3 匹配相位是选择；
3. 拟合参数 3 个、数据 6 个，自由度少；
4. 不是 Z 篇；无新量纲锚。

---

## §5 判据（3/4 PASS；FAIL 2 是要报的判定）

| # | 判据 | 结果 |
|--:|:--|:--|
| 1 | 夸克 $\varphi$ 有极小（$\Delta\chi^2>4$） | PASS，$4.05$ |
| 2 | 轻子 $\varphi$ 有极小（$\Delta\chi^2>4$） | **FAIL**，$2.17$ |
| 3 | 两扇区 $\varphi^*$ 一致（$|\Delta\varphi|<0.1$） | PASS，$0.085$ |
| 4 | 拟合可接受（$\chi^2<1.5$） | PASS，$0.82/1.02$ |

---

## §6 复跑与文件

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_world6_engine_phase.py     # 4 条判据，约 10 s
```

| 文件 | 内容 |
|:--|:--|
| `src/zero_world6_engine_phase.py` | 世界模型拟合 $\varphi$ |
| `results/zero_world6_engine_phase.json` | $\varphi^*$、$\chi^2$、判据 |
| `results/zero_world6_engine_phase.png` | $\chi^2(\varphi)$ |

---

*ZeroGPU 讨论稿 · zero_world6_engine_phase · 2026-10-19*
