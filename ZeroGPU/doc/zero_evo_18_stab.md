# 给 γ>0 加份额负反馈：找到稳定窗口 λ≈10

**性质**：讨论稿/原型，不占正式篇号，不改变根目录 Z 篇的任何判定。
**脚本**：`ZeroGPU/src/zero_evo_18_stab.py`。
**数据**：`results/zero_evo_18_stab.json`、`zero_evo_18_stab.png`。
**复跑**：

```bash
cd ZeroGPU/src && /usr/bin/python3 zero_evo_18_stab.py
```

---

## 摘要

上一轮：γ>0 长程漂移。本轮在代权重上加**份额负反馈**：

$$
\text{gen\_eff}_g\propto \text{gen\_w}_g\,
\exp\!\Big(\gamma(\overline{|\Delta n|}_{M_g}-\overline{\cdot})
-\lambda(\text{share}_g-\text{gen\_w}_g)\Big).
\tag{0.1}
$$

扫 $\lambda\in\{0,1,5,10,20\}$，$\gamma=2$，5000 代/点：

$$
\boxed{\;\lambda\approx10\text{ 是稳定窗口：relstd }<0.1\text{、}|{\rm drift}|<10^{-6};\quad
\lambda=0\text{ 不稳，}\lambda=20\text{ 又不稳。}\;}
$$

| $\lambda$ | 份额 | relstd | 漂移 |
|--:|:--|--:|--:|
| 0 | $0.420,0.339,0.241$ | $0.167,0.335,0.187$ | $-2.8\times10^{-6}$ |
| 1 | $0.398,0.366,0.235$ | $0.154,0.283,0.185$ | $-2.1\times10^{-6}$ |
| 5 | $0.365,0.416,0.219$ | $0.098,0.156,0.139$ | $-1.2\times10^{-6}$ |
| **10** | $0.353,0.437,0.210$ | **$0.064,0.096,0.098$** | **$-7.6\times10^{-7}$** |
| 20 | $0.355,0.461,0.184$ | $0.199,0.199,0.130$ | $-5.8\times10^{-7}$ |

---

## §1 读法

- **无反馈不稳**：$\lambda=0$ 的 relstd 到 0.335、漂移 $-2.8\times10^{-6}$；
- **$\lambda=10$ 稳定**：三项 relstd 都 $<0.1$、漂移 $<10^{-6}$；
- **反馈太强也不稳**：$\lambda=20$ 的 relstd 又回到 0.2——不是单调；
- **稳定窗口窄**：$\lambda\approx10$，且稳定后的份额 $[0.353,0.437,0.210]$ 与
  $\gamma=0$ 的 $[0.333,0.472,0.195]$ 不同（反馈把份额拉向 gen_w，但没完全拉回）。

---

## §2 判定

1. **γ>0 可以稳定化**：加份额负反馈 $\lambda\approx10$，relstd $<0.1$、漂移 $<10^{-6}$。
2. **代价是新增一个参数 λ**：它没有来源，是稳定化输入。
3. **稳定是边际的**：relstd 仍在 0.06–0.10，不是严格固定点。
4. **γ=0 仍是最干净的基线**：不漂、不需要 λ，但不含代依赖耦合。

---

## §3 边界

1. 5000 代/点；更长跑要看 $\lambda=10$ 是否真的收敛（relstd 可能继续降，也可能慢漂）。
2. $\lambda$ 未扫更细（8/12/15）；窗口可能更宽或更窄。
3. 反馈形式（线性份额差）是选择；别的形式（饱和、上限）未试。
4. 不是 Z 篇；无新量纲锚。

---

## §4 判据（4/4 PASS）

| # | 判据 | 结果 |
|--:|:--|:--|
| 1 | $\lambda=0$ 不稳 | PASS，relstd 0.335 |
| 2 | 存在稳定窗口 | PASS，$\lambda=10$ |
| 3 | $\lambda=10$ 稳定 | PASS，relstd $<0.1$ |
| 4 | 反馈不破坏有界性 | PASS |

---

## §5 复跑与文件

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_evo_18_stab.py     # 4 条判据，约 5 min
```

| 文件 | 内容 |
|:--|:--|
| `src/zero_evo_18_stab.py` | $\lambda$ 扫描 |
| `results/zero_evo_18_stab.json` | 份额、relstd、漂移、判据 |
| `results/zero_evo_18_stab.png` | gen0 份额轨迹（各 λ） |

---

*ZeroGPU 讨论稿 · zero_evo_18_stab · 2026-10-19*
