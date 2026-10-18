# 16/8/4 + cycle16：$\beta$ 扫描（content→phase，不是 h→phase）

## 结论

**先纠正一个口径**：引擎的 `beta` 是

$$
\psi \;\leftarrow\; \exp\!\Bigl(-\tfrac{i}{2}dt\,\beta\,(n-\bar n)\Bigr)\psi,
$$

即 **content $n$ → phase $\psi$**，**不是** geometry $h$ → phase $\psi$。所以这轮测的是“内容→波”，不是“几何→波”。

| $\beta$ | content $q$ | $h_E$ | $\psi_{coh}$ | $n_{std}$ | $\psi$ 周期 4 |
|--:|--:|--:|--:|--:|--:|
| $0$ | $0.1123$ | $39.6$ | $0.2399$ | $0.8926$ | $0.0029$ |
| $0.01$ | $0.1123$ | $39.6$ | $0.2424$ | $0.8926$ | $0.0042$ |
| $0.1$ | $0.1123$ | $39.6$ | $0.2540$ | $0.8926$ | $0.0261$ |
| $1.0$ | $0.1123$ | $39.6$ | $0.1656$ | $0.8926$ | $0.0368$ |

- $\beta>0$ 让 $\psi$ 的周期 4 从 $0.0029$ 升到 $0.0368$（仍弱）；
- $\beta=1$ 时 $\psi_{coh}$ 从 $0.24$ 降到 $0.17$（去相干），但有界；
- $h_E$、content $q$、$q_{16}$、$n_{std}$ 不受 $\beta$ 影响；
- 无 NaN，4/4 判据。

![β 扫描](ZeroGPU/results/zero_evo_1684_beta.png)

## 判据（4/4）

1. 无 NaN：PASS；
2. $\psi_{coh}\le1$：PASS；
3. $\beta>0$ 时 $\psi$ 有周期 4：PASS（$0.026$–$0.037$）；
4. $\beta>0$ 时 $h_E$ 仍在：PASS（$39.6$）。

## 意义

- **内容→波**这条边存在，但很弱：$\beta$ 要 $0.1$ 才把 $\psi$ 周期 4 抬到 $2.6\%$；
- **几何→波**这条边**还没接**：引擎没有 $h\to\psi$ 项；所以“记录→几何→波”的链条只完成了“记录→几何”，几何→波仍缺；
- $\beta$ 不是“几何耦合”，不能拿它当 h→psi 的证据。

## 边界

- $\beta$ 是 content→phase 的探索性耦合；不是从 16/8/4 推出来的。
- $\psi$ 周期 4 仍弱（$<0.04$），不能称“波也带周期 4”。
- 没有 $h\to\psi$ 项；若要做几何→波，必须新增一个具名耦合（例如 $\beta_h (h-\bar h)$），并单独登记。
- 16/8/4 仍是输入；$q_{16}=0.1123$ 不对目标拟合。

## 下一步

1. **新增 $h\to\psi$ 项**：在 `_wave_step` 里加
   $$
   \psi\leftarrow\exp\!\Bigl(-\tfrac{i}{2}dt\,\beta_h(h-\bar h)\Bigr)\psi,
   $$
   扫 $\beta_h$，看 $\psi$ 是否被几何带出周期 4；这是新耦合，必须登记。
2. **或接受 content→wave 是唯一已接的波源**：那“几何→波”仍是缺口，不能声称闭环。
3. **不要用 $\beta$ 冒充 $h\to\psi$**：口径已经写清。

## 复跑

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_evo_1684_beta.py
```

| 文件 | 内容 |
|:--|:--|
| `ZeroGPU/src/zero_evo_1684_beta.py` | $\beta$（content→phase）扫描 |
| `ZeroGPU/results/zero_evo_1684_beta.json` | $\psi_{coh}$、$h_E$、谱 |
| `ZeroGPU/results/zero_evo_1684_beta.png` | $\psi_{coh}$ 与 $h_E$ |

---

*ZeroGPU · 16/8/4 β 扫描 · 2026-10-18*
