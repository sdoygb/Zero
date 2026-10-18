# 16/8/4 + cycle16：site-resolved 映射与 geometry/wave 谱

## 结论

- **聚合量**（$q_{16}$、content $q$、份额）对映射不敏感；
- **site-resolved 量**（$n$ 向量、top1/top2）**随映射改变**；
- **geometry 扇区是零**（$h_E=0$）；**wave 扇区**的周期 4 很弱（$2.9\times10^{-3}$）。

![site-resolved 与谱](ZeroGPU/results/zero_evo_1684_site.png)

## 映射对比（$\alpha=0.1$，2000 代）

| 映射 | $q_{16}$ | content $q$ | $h_E$ | $\psi_{coh}$ | top1/top2 |
|:--|--:|--:|--:|--:|:--|
| `seq` | $0.1123$ | $0.1123$ | $0$ | $0.2409$ | $0/1$ |
| `transpose` | $0.1123$ | $0.1123$ | $0$ | $0.2409$ | $0/4$ |
| `shift` | $0.1123$ | $0.1123$ | $0$ | $0.2409$ | $0/1$ |
| `bitrev` | $0.1123$ | $0.1123$ | $0$ | $0.2409$ | $0/2$ |
| `random` | $0.1123$ | $0.1123$ | $0$ | $0.2409$ | $4/14$ |

- $q_{16}$、content $q$ 完全一样；
- $n$ 向量的最大差 $\lVert n_{\rm map}-n_{\rm seq}\rVert=5.83$——**局域内容不同**；
- top1/top2 随映射变（`seq` 是 $0/1$，`random` 是 $4/14$）。

## 谱

| 量 | 周期 4 bin |
|:--|--:|
| $q_{16}$ | $1.17$（强） |
| $\psi_{coh}$ | $2.9\times10^{-3}$（弱） |
| $h_E$ | $0$（几何扇区无动力学） |

- $h_E=0$ 不是数值问题：默认 $g_{\rm src}=0$、$\beta=0$、$h_0=h_v=0$，几何扇区没有源；
- $\psi$ 有微弱周期 4，但幅度很小；
- **16/8/4 的周期 4 主要在 $q_{16}$ 上，不在 geometry/wave 上**。

## 判据（4/4）

1. $q_{16}$ 不随映射变：PASS（spread $0$）；
2. content $q$ 不随映射变：PASS（spread $0$）；
3. $n$ 向量随映射变（site-resolved 不不变）：PASS（$\Delta=5.83$）；
4. $q_{16}$ 周期 4 仍在：PASS（$1.17$）。

## 边界

- 这是 **2000 代**的短测；长跑 site-resolved 未测。
- `h_E=0` 是默认参数的结果，不是 16/8/4 底层的失败；要测几何谱，必须开 $g_{\rm src}$/$\beta$。
- 映射不唯一是**物理缺口**：class $\to$ 16 顶点的双射没有从 Zero 结构推出来；聚合量不敏感，但局域预测依赖它。
- 16/8/4 仍是输入；$q_{16}=0.1123$ 不对目标拟合。

## 下一步

1. **给映射一个物理来源**：用 K4 邻接/角色或 cycle16 的顶点序，把 class 映射固定下来；否则 site-resolved 预测不可判。
2. **开几何源**：$g_{\rm src}>0$ 或 $\beta>0$，看 $h_E$/$\psi$ 是否出现周期 4；现在 $h_E=0$ 是因为没源。
3. **长跑 site-resolved**：10000 代，看 top1/top2 是否稳定，还是随映射长期漂移。

## 复跑

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_evo_1684_site.py
```

| 文件 | 内容 |
|:--|:--|
| `ZeroGPU/src/zero_evo_1684_site.py` | 五种映射的 site-resolved/谱对比 |
| `ZeroGPU/results/zero_evo_1684_site.json` | $q$、$n$ 差、top 标签、谱 |
| `ZeroGPU/results/zero_evo_1684_site.png` | $q$ 与 geometry/wave |

---

*ZeroGPU · 16/8/4 site-resolved · 2026-10-18*
