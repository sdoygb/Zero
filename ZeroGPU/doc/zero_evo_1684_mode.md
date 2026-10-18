# 模依赖共形红移：$\omega_{\rm op}$ 版

## 结论

把 $\beta_h=\omega/2$ 从**常数**升级为**模依赖算子**后，波仍然稳定，周期 4 更强：

| 耦合 | $\psi_{coh}$ | $h_E$ | $q_{16}$ | $\psi$ 周期 4 |
|:--|--:|--:|--:|--:|
| 标量 $\beta_h=\omega_{\rm char}/2=0.717$ | $0.2693$ | $39.625$ | $0.1123$ | $0.1831$ |
| **模依赖 $\omega_{\rm op}$** | $0.1999$ | $39.625$ | $0.1123$ | **$0.2932$** |

- 两种耦合下 content $q$、$q_{16}$、$h_E$、$n_{std}$ 完全一样；
- 模依赖版的 $\psi$ 周期 4 是 $0.293$，比标量版强；
- 无 NaN，$\psi_{coh}$ 有界，4/4 判据。

![模依赖红移](ZeroGPU/results/zero_evo_1684_mode.png)

## 形式

共形红移的相位不是常数 $\omega$，而是算子

$$
\Omega = V\,{\rm diag}(\omega_k)\,V^{\mathsf T},
\qquad
\omega_k=\sqrt{c_w^2\lambda_k+\omega_0^2}.
$$

$h$ 是局域乘法算子 $H={\rm diag}(h-\bar h)$，两者不对易；取对称化

$$
G=\frac12\bigl(\Omega H+H\Omega\bigr),
\qquad
\psi\leftarrow e^{-i\,dt\,G}\psi .
$$

实现上用 $G$ 的本征分解：

$$
G=U_G\,{\rm diag}(g)\,U_G^\dagger,
\qquad
e^{-i\,dt\,G}=U_G\,{\rm diag}(e^{-i\,dt\,g})\,U_G^\dagger .
$$

每步拆成“半相位—自由演化—半相位”：

$$
e^{-i\frac{dt}{2}G}\,e^{-i\,dt\,H_0}\,e^{-i\frac{dt}{2}G}.
$$

## 判据（4/4）

1. 标量 $\beta_h=\omega/2$ 有周期 4：PASS（$0.183$）；
2. 模依赖相位有周期 4：PASS（$0.293$）；
3. 两者 $\psi_{coh}$ 有界：PASS；
4. 两者 $h_E$ 相近：PASS（$39.625$）。

## 意义

- **$\beta_h$ 不再是一个自由常数**：$\Omega$ 由波谱 $\omega_k$ 给；
- 标量 $\omega_{\rm char}/2$ 是它的**主导模近似**；
- 模依赖版保留了各模的不同红移，周期 4 更强（$0.183\to0.293$）；
- $h_E$、$q_{16}$ 不受影响：红移只改相位，不改几何/内容。

## 边界

- $G=\frac12(\Omega H+H\Omega)$ 是**对称化选择**；$\Omega$ 与 $H$ 不对易，顺序（$\Omega H$、$H\Omega$、对称化）会改结果，尚未定。
- 因子 $1/2$ 依赖 $h$ 的归一化（共形因子 vs $g_{00}$ 势）；这里沿用共形约定。
- 这是 **2000 代**短测；长跑/更大 $h$ 未测。
- $\omega_{\rm op}$ 用自由波谱 $\omega_k$，未把 $h$ 的反作用算进去；是**测试场近似**。
- 16/8/4 仍是输入；$q_{16}=0.1123$ 不对目标拟合。

## 下一步

1. **定顺序**：比较 $\Omega H$、$H\Omega$、对称化三种，看 $\psi$ 周期 4/相干性差异；由 $h$ 的物理定义选。
2. **长跑**：10000 代，看模依赖版 $\psi$ 周期 4 $0.293$ 是否持续。
3. **对齐 $h$ 归一化**：查 Z093/Z109，确定因子 $1/2$ 是否要改为 $1$。

## 复跑

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_evo_1684_mode.py
```

| 文件 | 内容 |
|:--|:--|
| `ZeroGPU/src/zero_evo_1684_mode.py` | 模依赖共形红移 |
| `ZeroGPU/results/zero_evo_1684_mode.json` | 标量 vs 模依赖的 $q$/$h_E$/谱 |
| `ZeroGPU/results/zero_evo_1684_mode.png` | 两种耦合的 $\psi_{coh}$ |

---

*ZeroGPU · 模依赖共形红移 · 2026-10-18*
