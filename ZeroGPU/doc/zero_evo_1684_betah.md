# 新增 $h\to\psi$：$\beta_h$ 扫描

## 结论

新增 $h\to\psi$ 项后，**几何的周期 4 能带进波**：

$$
\psi\leftarrow\exp\!\Bigl(-\tfrac{i}{2}dt\,\beta_h(h-\bar h)\Bigr)\psi .
$$

| $\beta_h$ | $\psi_{coh}$ | $h_E$ | $q_{16}$ | $n_{std}$ | $\psi$ 周期 4 |
|--:|--:|--:|--:|--:|--:|
| $0.1$ | $0.1791$ | $39.6$ | $0.1123$ | $0.8926$ | $0.0501$ |
| $0.3$ | $0.3021$ | $39.6$ | $0.1123$ | $0.8926$ | $0.0105$ |
| $1.0$ | $0.3361$ | $39.6$ | $0.1123$ | $0.8926$ | $0.3562$ |
| $3.0$ | $0.2056$ | $39.6$ | $0.1123$ | $0.8926$ | **$5.939$** |
| $10.0$ | $0.2330$ | $39.6$ | $0.1123$ | $0.8926$ | **$6.877$** |

**长跑（$\beta_h=1$，10000 代）**：$\psi_{coh}=0.2493$、$\psi$ 周期 4 $=0.7411$、$h_E$ 周期 4 $=4.196$、$q_{16}=0.1123$；$\psi_{coh}$ 最大 $0.6515$，无 NaN。

- $\beta_h$ 从 $0.1$ 到 $10$，$\psi$ 周期 4 从 $0.05$ 升到 $6.88$，比 content→phase 的 $\beta$ 强两个量级；
- **长跑 $\beta_h=1$**：10000 代后 $\psi$ 周期 4 仍有 $0.741$，不发散；
- $\psi_{coh}$ 在 $0.18$–$0.34$ 间变，仍有界；
- $h_E$、content $q$、$q_{16}$、$n_{std}$ 不变；
- 无 NaN，4/4 判据。

![β_h 扫描](ZeroGPU/results/zero_evo_1684_betah.png)

## 链条

$$
\underbrace{q_{16}}_{\text{16/8/4 记录}}
\;\to\;
\underbrace{h}_{\text{几何，周期4}}
\;\xrightarrow{\beta_h}\;
\underbrace{\psi}_{\text{波，周期4}}.
$$

- 记录→几何：$\alpha$、$g_{\rm src}$（上一步）；
- 几何→波：**新的 $\beta_h$**（本步）；
- 内容→波：$\beta$（上一步，弱）。

**“记录→几何→波”在探索层面已经接通。**

## 登记：$\beta_h$ 已进入 `Params`

$\beta_h$ 不再是脚本里的子类属性，而是正式参数：

- `zero_evo.Params.beta_h: float = 0.0`；
- `zero_evo.ZeroEvo._wave_step` 里加入
  $$
  \psi\leftarrow\exp\!\Bigl(-\tfrac{i}{2}dt\,\beta_h(h-\bar h)\Bigr)\psi,
  $$
  与已有 $\beta(n-\bar n)$ 同构；
- **默认 $0$，旧行为不变**；`zero_evo.py --selftest` **20/20 PASS**。

这是“登记”，不是“推导”：$\beta_h$ 仍是探索性输入，只是不再藏在脚本子类里。

## 来源：共形红移 $\beta_h=\omega/2$

如果 $h$ 是**共形因子**（$g_{ij}=(1+h)\delta_{ij}$），则固有时

$$
d\tau=(1+h/2)\,dt,
$$

波的相位是 $\omega\tau$，相对相位

$$
\omega\tau-\omega t=\tfrac{\omega}{2}h,
$$

所以

$$
\boxed{\;\beta_h=\frac{\omega}{2}.\;}
$$

如果 $h$ 是 $g_{00}$ 的牛顿势（$d\tau=(1+h)dt$），则 $\beta_h=\omega$；**因子 2 取决于 $h$ 的归一化约定**，需要和 Z093/Z109 的 $h$ 定义对齐。

### 数值检验

从 $\beta_h=0$ 的 2000 代波谱算特征频率：

$$
\omega_{\rm char}=\frac{\sum_k P_k\omega_k}{\sum_k P_k}=1.433,
\qquad
\frac{\omega_{\rm char}}{2}=0.717.
$$

| $\beta_h$ | 来源 | $\psi_{coh}$ | $\psi$ 周期 4 |
|--:|:--|--:|--:|
| $0.717$ | $\omega/2$（共形） | $0.2693$ | $0.1831$ |
| $1.433$ | $\omega$（势） | $0.2471$ | $3.583$ |
| $10$ | 超出 eikonal（$\beta_h\gg\omega$） | $0.2330$ | $6.877$ |

**结论**：

- $\beta_h$ 的候选来源是**共形/引力红移**；
- 在 $\omega/2$ 处，$\psi$ 周期 4 是 $0.183$，不是 $6.88$；
- 扫描里 $\beta_h=3$–$10$ 的强周期 4 **超出 eikonal 范围**（每步相位过大），不能当物理信号；
- 单值 $\beta_h$ 只是**主导模近似**；精确形式是**模依赖**的 $\omega_k$，不是常数。

## 判据（7/7）

1. 所有 $\beta_h$ 无 NaN：PASS；
2. 所有 $\beta_h$ 的 $\psi_{coh}\le1$：PASS；
3. $\beta_h$ 增大时 $\psi$ 周期 4 不消失：PASS（$0.05$–$6.88$）；
4. 所有 $\beta_h$ 的 $h_E$ 仍在：PASS（$39.6$）；
5. 长跑 $\psi$ 有周期 4：PASS（$0.741$）；
6. 长跑 $\psi_{coh}$ 有界：PASS（最大 $0.652$）；
7. 长跑无 NaN：PASS。

## 边界

- **$\beta_h$ 已登记为 `Params.beta_h`（默认 0）**，但它仍是探索性输入，不是从 16/8/4 或 Z 篇推出来的；不能当成“已导出”。
- 形式 $\exp(-\frac{i}{2}dt\beta_h(h-\bar h))$ 与引擎已有的 $\beta(n-\bar n)$ 同构，但作用对象换成了 $h$；这是**类比**，不是定理。
- $\psi_{coh}$ 随 $\beta_h$ 非单调（$0.18\to0.34\to0.21\to0.23$）；相位幅度与相干性不是一回事。
- $\beta_h=10$ 的周期 4 很大（$6.88$），但这是 $\psi_{coh}$ 的 FFT 幅度，不是“物理信号更强”的直接判据；需要看相位/谱形。
- 16/8/4 仍是输入；$q_{16}=0.1123$ 不对目标拟合。

## 下一步

1. **扫更大 $\beta_h$**：$0.1$–$10$，看 $\psi$ 周期 4 是否继续增强、$\psi_{coh}$ 是否失稳。
2. **长跑**：10000 代，看 $h_E$/$psi$ 周期 4 是否持续。
3. **给 $\beta_h$ 找来源**：它现在已登记，但还没有从 16/8/4/Z 篇推出；下一步是解释它的值，而不是再调它。

## 复跑

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_evo_1684_betah.py
```

| 文件 | 内容 |
|:--|:--|
| `ZeroGPU/src/zero_evo_1684_betah.py` | 新增 $h\to\psi$ 的 $\beta_h$ 扫描 |
| `ZeroGPU/results/zero_evo_1684_betah.json` | $\psi_{coh}$、$h_E$、谱 |
| `ZeroGPU/results/zero_evo_1684_betah.png` | $\psi_{coh}$ 与 $h_E$ |

---

*ZeroGPU · 16/8/4 β_h 耦合 · 2026-10-18*
