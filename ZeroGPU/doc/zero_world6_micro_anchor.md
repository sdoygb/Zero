# 从微观世界找相位锚：CP 相位

**性质**：讨论稿/锚定，不占正式篇号，不改变根目录 Z 篇的任何判定。
**脚本**：`ZeroGPU/src/zero_world6_micro_anchor.py`。
**数据**：`results/zero_world6_micro_anchor.json`。
**复跑**：

```bash
cd ZeroGPU/src && /usr/bin/python3 zero_world6_micro_anchor.py
```

---

## 摘要

修正：**不从宇宙学找锚，从微观世界的“相位型”观测量找**。候选是标准模型里的
**CP 相位**（都是 $2\pi$ 周期的相位）：

- CKM 的 $\gamma$（PDG 量级 $65.5^\circ$）；
- PMNS 的 $\delta_{CP}$（全球拟合范围约 $0.7\pi$–$1.5\pi$）。

约定：

$$
\varphi=\frac{\delta}{2\pi}\bmod 1,\qquad
t_{\rm remain}=(1-\varphi)\,T,\quad T=14.51\ {\rm Gyr}.
\tag{0.1}
$$

$$
\boxed{\;\text{微观相位锚给 }\varphi\approx0.18\text{–}0.75,\quad
t_{\rm remain}\approx3.6\text{–}11.9\ {\rm Gyr}.\;}
\tag{0.2}
$$

| 锚 | $\delta$ | $\varphi$ | 剩余 |
|:--|--:|--:|--:|
| CKM $\gamma=65.5^\circ$ | $1.143$ | $0.182$ | $11.87$ Gyr |
| PMNS $\delta_{CP}=0.7\pi$ | $2.199$ | $0.350$ | $9.43$ Gyr |
| PMNS $\delta_{CP}=1.0\pi$ | $3.142$ | $0.500$ | $7.25$ Gyr |
| PMNS $\delta_{CP}=1.2\pi$ | $3.770$ | $0.600$ | $5.80$ Gyr |
| PMNS $\delta_{CP}=1.5\pi$ | $4.712$ | $0.750$ | $3.63$ Gyr |

---

## §1 读法

- **CKM $\gamma$**（测量较准）：$\varphi\approx0.18$，剩余 $\approx11.9$ Gyr；
- **PMNS $\delta_{CP}$**（测量较松，范围 $0.7\pi$–$1.5\pi$）：$\varphi\approx0.35$–$0.75$，
  剩余 $\approx3.6$–$9.4$ Gyr；
- **不同微观锚给不同 $\varphi$**：没有唯一映射。

---

## §2 判定

1. **微观锚可行**：CP 相位是 $2\pi$ 周期量，天然可映射到周期相位。
2. **结果是一段，不是一个数**：$\varphi\approx0.18$–$0.75$，
   剩余 $\approx3.6$–$11.9$ Gyr。
3. **映射是约定**：$\varphi=\delta/(2\pi)$ 不是 Zero 内生的；要变成预言，
   需要 Zero 能算的相位（例如把 $\delta$ 与 6 层的相对相位直接对应）。
4. **比宇宙学锚更合适**：微观相位不依赖大爆炸/宇宙年龄，符合 Zero 没有大爆炸。

---

## §3 边界

1. CP 相位的测量值有不确定度（尤其 PMNS $\delta_{CP}$）；
2. 映射 $\varphi=\delta/(2\pi)$ 是选择；
3. 周期 $T=14.51$ Gyr 仍是哈勃锚；
4. 不是 Z 篇；无新量纲锚。

---

## §4 判据（4/4 PASS）

| # | 判据 | 结果 |
|--:|:--|:--|
| 1 | 微观相位锚给出 $\varphi\in(0,1)$ | PASS |
| 2 | CKM $\gamma$：$\varphi\approx0.18$、剩余 $\approx11.9$ Gyr | PASS |
| 3 | PMNS $\delta_{CP}$：$\varphi$ 0.35–0.75、剩余 3.6–9.4 Gyr | PASS |
| 4 | 不同锚给不同 $\varphi$ | PASS |

---

## §5 复跑与文件

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_world6_micro_anchor.py     # 4 条判据，<1 s
```

| 文件 | 内容 |
|:--|:--|
| `src/zero_world6_micro_anchor.py` | CP 相位锚与 $\varphi$、剩余 |
| `results/zero_world6_micro_anchor.json` | 全部数值、判据 |

---

*ZeroGPU 讨论稿 · zero_world6_micro_anchor · 2026-10-19*
