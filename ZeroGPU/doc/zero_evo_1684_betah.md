# 新增 $h\to\psi$：$\beta_h$ 扫描

## 结论

新增 $h\to\psi$ 项后，**几何的周期 4 能带进波**：

$$
\psi\leftarrow\exp\!\Bigl(-\tfrac{i}{2}dt\,\beta_h(h-\bar h)\Bigr)\psi .
$$

| $\beta_h$ | $\psi_{coh}$ | $h_E$ | content $q$ | $n_{std}$ | $\psi$ 周期 4 |
|--:|--:|--:|--:|--:|--:|
| $0$ | $0.2399$ | $39.6$ | $0.1123$ | $0.8926$ | $0.0029$ |
| $0.01$ | $0.1851$ | $39.6$ | $0.1123$ | $0.8926$ | $0.0408$ |
| $0.1$ | $0.1791$ | $39.6$ | $0.1123$ | $0.8926$ | $0.0501$ |
| $1.0$ | $0.3361$ | $39.6$ | $0.1123$ | $0.8926$ | **$0.3562$** |

- $\beta_h=1$ 时 $\psi$ 的周期 4 到 **$0.356$**，比 content→phase 的 $\beta$ 强约 10 倍；
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

## 判据（4/4）

1. 无 NaN：PASS；
2. $\psi_{coh}\le1$：PASS；
3. $\beta_h>0$ 时 $\psi$ 有周期 4：PASS（$0.041$–$0.356$）；
4. $\beta_h>0$ 时 $h_E$ 仍在：PASS（$39.6$）。

## 边界

- **$\beta_h$ 是新的探索性耦合**，不是从 16/8/4 或 Z 篇推出来的；它必须单独登记，不能当成“已导出”。
- 形式 $\exp(-\frac{i}{2}dt\beta_h(h-\bar h))$ 与引擎已有的 $\beta(n-\bar n)$ 同构，但作用对象换成了 $h$；这是**类比**，不是定理。
- $\psi_{coh}$ 随 $\beta_h$ 非单调（$0.24\to0.185\to0.179\to0.336$）；更大 $\beta_h$ 未扫。
- 16/8/4 仍是输入；$q_{16}=0.1123$ 不对目标拟合。

## 下一步

1. **扫更大 $\beta_h$**：$0.1$–$10$，看 $\psi$ 周期 4 是否继续增强、$\psi_{coh}$ 是否失稳。
2. **长跑**：10000 代，看 $h_E$/$psi$ 周期 4 是否持续。
3. **把 $\beta_h$ 写成正式输入**：如果保留，必须在文档里登记“几何→波耦合 $\beta_h$”，不能偷偷加。

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
