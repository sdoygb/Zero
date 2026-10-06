# 第一个可测量的数：零参数无量纲预言总账

> **日期**：本轮 · **性质**：**汇总＋独立复算＋更正**（把 `cosmos-construct` 的可用成果收进 Zero 侧）
> **机器**：[`z0_predictions_inventory.py`](z0_predictions_inventory.py)（**12/12**）· [`z0_neutrino_absolute.py`](z0_neutrino_absolute.py)（**8/8**）
> **来源**（`~/Downloads/cosmos-construct/`，3327 文件）：
> `MAPPING_FORMULAS.md`（16 条映射公式总表）· `LEPTON_PHASE_MAPPING.md` · `NEUTRINO_MASS_COMPUTED.md` ·
> `NEUTRINO_ABSOLUTE_MASSES.md` · `Q204_RESULT.md` · `Q244_ALPHA_S.md` · `GUT_PREDICTIONS.md`

---

## §0 一句话

$$
\boxed{\ \textbf{Zero 侧确实有第一个可测量的数 —— 而且是【零参数、无量纲】的：带电轻子质量比。}\ }
$$

$$
\frac{m_\mu}{m_e}=206.770316\ (\text{实测 }206.768283,\ \textbf{偏差 }0.0010\%)\qquad
\frac{m_\tau}{m_e}=3477.4728\ (\text{实测 }3477.2283,\ \textbf{偏差 }0.0070\%)
$$

**输入 0 个自由参数，输出 2 个比值。**

---

## §1 带电轻子（**零参数**）

### §1.1 公式

$$
\delta_0\equiv\frac{\pi}{12}-\frac{k_0}{\Lambda^2}\ \ (k_0=2,\ \Lambda=3)
=0.039577166\ \text{rad}=2.2676046^\circ
$$

$$
m_n\ \propto\ \bigl(1-\cos\delta_n+\sin\delta_n\bigr)^2,\qquad
\delta_n=\delta_0-\frac{2\pi n}{3}\quad(n=0,1,2)
$$

### §1.2 独立复算（本文件）

| 比值 | 框架预言 | 实测 | **偏差** |
|:--|--:|--:|--:|
| $m_\mu/m_e$ | $206.770316$ | $206.768283$ | $\mathbf{0.0010\%}$ |
| $m_\tau/m_e$ | $3477.4728$ | $3477.2283$ | $\mathbf{0.0070\%}$ |
| $m_\tau/m_\mu$ | $16.817662$ | $16.817029$ | $0.0038\%$ |
| **锚 $A^2$** | $313.859230$ MeV | （由 $m_e$ 反解） | 与 `build_166` 拟合差 $0.0019\%$ |

### §1.3 而 Koide 关系**精确**成立

$$
K=\frac{m_e+m_\mu+m_\tau}{(\sqrt{m_e}+\sqrt{m_\mu}+\sqrt{m_\tau})^2}=\frac23
\qquad(\text{框架值偏差 }-5.6\times10^{-14}\%,\ \text{实测 }0.6666605)
$$

**几何来源**：$\sqrt2=2\cos(\theta_{\text{物质}}/2)$，$\theta_{\text{物质}}=90^\circ$。

> **这是第一条能与实验对到 $10^{-5}$ 的零参数预言。**

---

## §2 中微子（**零参数形状 ＋ 一个已测量锚**）

### §2.1 形状（零参数）

$$
m_1=0,\qquad R=\frac{m_2}{m_3}=\frac{73-28\sqrt6}{25}=0.176571488
$$

### §2.2 绝对质量（一个已测量的 $\Delta m^2$ 定标）

$$
\text{[A] 以 }\Delta m^2_{21}\text{ 定标}\ \Longrightarrow\ m_2=\sqrt{\Delta m^2_{21}},\ m_3=m_2/R
\ \Longrightarrow\ \boxed{\ \Delta m^2_{31}\ \textbf{成为预言}\ }
$$

| 定标 | $m_2$ | $m_3$ | $\Sigma m_\nu$ | 预言的另一个 $\Delta m^2$ |
|:--|--:|--:|--:|--:|
| **[A]** $\Delta m^2_{21}$（NuFIT 5.2） | $8.6139$ meV | $48.7844$ meV | $\mathbf{57.3984}$ meV | $\Delta m^2_{31}$ 偏差 $\mathbf{-8.139\%}$ |
| **[B]** $\Delta m^2_{31}$ | $8.8462$ meV | $50.0999$ meV | $58.9461$ meV | $\Delta m^2_{21}$ 偏差 $+5.466\%$ |

### §2.3 本文件补算的两项（原脚本缺）

$$
m_\beta=8.657\ \text{meV}\quad(\text{KATRIN 上限 }450\ \text{meV}\ ✅)\qquad
m_{\beta\beta}\le 3.637\ \text{meV}\ (\text{同相上限})
$$

> ⚠ **注意**：$m_{\beta\beta}\le3.6$ meV **低于**下一代 $0\nu\beta\beta$（nEXO／LEGEND-1000）$\sim10$–$20$ meV 的灵敏度 ⟹ **本预言在可预见的将来不可判**。

### §2.4 可否证性

$$
\Sigma m_\nu=0.0574\ \text{eV}\ <\ \underbrace{0.12\ \text{eV}}_{\text{Planck 2018+BAO}}\ ✅
\qquad <\ \underbrace{0.072\ \text{eV}}_{\text{DESI DR2 2024}}\ ✅
$$

$$
\textbf{正序（NO）＋ 最轻态饱和 }(m_1=0)
$$

---

## §3 角度与群论量（**零参数、精确**）

| # | 量 | 值 | 来源 |
|--:|:--|:--|:--|
| I2 | $(n_a)=(4,3,12)$，$(k_a)=(1,1,5)$，$n_3=n_1n_2=12$ | — | `build_64` |
| I2 | **求和律** $\sum_a k_a/n_a=\frac14+\frac13+\frac5{12}=1$ | **精确** | `Q_LAMBDA` |
| I2 | $\theta_a=(90^\circ,120^\circ,150^\circ)$，$\Sigma\theta=2\pi$ | **精确** | `build_64` |
| I4 | $\lvert X_a\rvert^2=(4,3,1)$ ⟹ $\lvert X_a\rvert=2:\sqrt3:1$ | **精确** | `build_96/97/98` |
| G4 | $\operatorname{tr}(Y^2)/\operatorname{tr}(T_3^2)=5/3$ | **群论必然** | `build_86` |
| G5 | $\sin^2\theta_W(M_X)=1/(1+5/3)=\mathbf{3/8}$ | **精确** | `build_95` |
| I5 | $S=X_1+X_2+X_3=0$（闭合） | **精确** | `build_16` |

---

## §4 ★ 唯一真缺口与已知张力（诚实栏）

### §4.1 `MAPPING_FORMULAS` 的结构性结论

> **所有几何／群论／角度／耦合比／混合角 ⟹ ✅ 全部已覆盖（零参数、无量纲）**
> **所有标度量（$\alpha_X$、$v$、$G_F$）⟹ ❌ 未覆盖**

**但重判之后**（`build_130`／`build_153`／`Q375`）：

| 项 | 判定 |
|:--|:--|
| $\alpha_X,v,G_F$ 的**绝对值** | ✅ **不是缺口** —— 参考能标是**单位**不是参数（`build_130`） |
| $\alpha_s(M_Z)$ | ❌ 几何链**未覆盖**（$M_X$ 由 $\alpha_3=\alpha_2$ 定 ⟹ 循环，`Q244`） |
| **$v/M_X$（一条无量纲比率）** | ❌ ★**唯一真缺口**★；候选机制 ＝ **作用量指数** $v/M_X=e^{-S}$，需 $S=\mathbf{33.6197}$ |

### §4.2 已标明的两处张力

| 张力 | 大小 |
|:--|:--|
| $\sin^2\theta_W$ 的统一链（$M_X$：$1.03\times10^{13}$ GeV 观测反解 vs $9.66\times10^{15}$ GeV 框架） | $\mathbf{-21.9\%}$ |
| 中微子 $\Delta m^2_{31}/\Delta m^2_{21}$ 比值（实测 $33.827$ vs $1/R^2=32.074$） | $\mathbf{+5.466\%}$ |

---

## §5 ★ 本文件更正的两处

### §5.1 `build_171` 的算式错误

```python
p31 = m3A ** 2                 # ✗ 漏了 − m₂²
p31 = m3A ** 2 - m2A ** 2      # ✓ 正确
```

| 项 | 原报告 | **正确值** |
|:--|--:|--:|
| NuFIT 5.2，[A] 预言 $\Delta m^2_{31}$ 的偏差 | $-5.18\%$ | $\mathbf{-8.139\%}$ |
| 结论"两种定标**两边都差 ~5%**" | 对称 | **不对称**：$8.14\%$ vs $5.47\%$ |

（[B] 的 $+5.47\%$ 原本就是对的——因为 $m_1=0$ 时 $\Delta m^2_{21}=m_2^2$。）

### §5.2 本文件自己的一处（记录以免再犯）

我第一次复算时把 $\delta_n$ 写成 $\delta_0+\frac{2\pi n}{3}$，还漏了 $k_3=5$ 与迹比 $5/3$，
得到 $m_\mu/m_e=3477$（错位）与 $\sin^2\theta_W=3/11$。**两处均已在机器里更正并写成断言。**

---

## §6 这对"Zero 是不是物理理论"意味着什么

**我在意见书里说**：
> "给我一个可以放在实验报告里的数，否则你不是物理学。"

**`cosmos-construct` 的回答是**：

| 预言 | 精度 | 零参数？ |
|:--|:--|:--|
| $m_\mu/m_e$ | $\mathbf{0.0010\%}$ | ✅ **零参数** |
| $m_\tau/m_e$ | $\mathbf{0.0070\%}$ | ✅ **零参数** |
| Koide $K=2/3$ | **精确** | ✅ |
| $\sin^2\theta_W(M_X)=3/8$ | **精确** | ✅ |
| $\Sigma m_\nu=57.4$ meV | 通过所有宇宙学上限 | 零参数形状 ＋ 一个锚 |
| $\lvert X_a\rvert=2:\sqrt3:1$ | **精确** | ✅ |

$$
\boxed{\ \textbf{所以意见书的第一击（"没有可测量的数"）【被驳倒】了 —— 但只在这些量上。}\ }
$$

**而仍然成立的**：
- **绝对标度**（$\alpha_X,v,G_F$ 的数值）不可导出 —— 是**单位**问题（`build_130`）；
- **$\alpha_s(M_Z)$** 未覆盖；
- **$v/M_X$** 是唯一真缺口（候选：作用量指数，$S=33.6197$）；
- **两处张力**（$22\%$ 与 $5.5\%$）未被解决。

$$
\boxed{\ \textbf{修正后的判词：这不是"没有数"，而是"有数、但数与数之间有 }5\%\text{–}22\%\text{ 的张力"。}\ }
$$

---

## §7 关联

| 篇／机器 | 内容 |
|:--|:--|
| `z0_predictions_inventory.py`（12/12） | 本清单的机器；带电轻子 ＋ 中微子 ＋ 角度群论 ＋ 缺口张力 |
| `z0_neutrino_absolute.py`（8/8） | 中微子绝对质量的独立复现（含 $m_\beta$、$m_\beta\beta$、更正） |
| `DIMENSION_ORIGIN.md` | 维数是划分的产物（本文与它无关，但同属"可测量性"这一支） |
| `ZERO_LEDGER.md` | 总账 |
| `~/Downloads/cosmos-construct/MAPPING_FORMULAS.md` | 16 条映射公式总表（**外部目录，只读引用**） |
| `~/Downloads/cosmos-construct/LEPTON_PHASE_MAPPING.md` | 带电轻子相角映射（`build_169`，17/17） |
| `~/Downloads/cosmos-construct/NEUTRINO_ABSOLUTE_MASSES.md` | 中微子绝对质量（`build_171`；本文更正其一处算式） |
