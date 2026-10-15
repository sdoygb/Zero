# zero_sim —— 模拟器公共核心（结构 ＋ 算法 ＋ 判据账）

**定位**：`ZeroGPU/src/` 里有 130+ 个判定脚本，但它们各自把公共结构复制了一遍
（`fano_heawood()` 在 `algebraic_capacity.py` / `capacity_destruction.py` /
`saturation_alignment.py` / `distribution_engine.py` 里各写一份，`girth()` 在 `anchor.py`）。
`zero_sim.py` 把这些收拢成一层：**新脚本从它 import，不再复制粘贴**。

**与 `z0_engine.py` 的分工**：`z0_engine.py` 是**动力学**引擎（走者演化、长时间跑）；
`zero_sim.py` 是**结构/组合**核心（图、谱、词层、耦合几何）＋ 判据账。
两者互补：动力学脚本要用图结构时，从 `zero_sim` 取。

---

## 三层接口

| 层 | 函数 | 说明 |
|:--|:--|:--|
| **L0 结构** | `k4()` `heawood()` `cycle(n)` `path(n)` `grid(n)` `double_cover(adj)` `fano_lines()` | `heawood()` = Fano 平面的关联图（点 0–6、线 7–13）；`fano_lines()` 用循环差集 `{j, j+1, j+3}` |
| **L0 算法** | `adj_matrix` `degrees` `is_bipartite` `girth` `spectrum` `aut_order` `stabilizer_edge_action` `adj_power_int` `closed_walks` `closed_walks_at` `return_prob` | `aut_order(adj, fixed={0:0})` 直接给顶点稳定子阶；`adj_power_int` 是**精确整数**矩阵幂（float 在 k≥40 已失真） |
| **L1 词层** | `balanced_words(L)` `layer_of(w)` `layer_distribution(L)` `word_area` `zero_area_counts` `mirror_N_leq` `mirror_shell` `orbit_count` | 层高口径 `h = max(M⁺, M⁻) = max_k |前缀和|`；`layer_distribution` 用 DP（L 任意大），小 L 可与穷举对照 |
| **耦合几何** | `simplex_vectors(D)` `gram(D)` `rho(D)` `peak_argmax(q)` | 正则 D-单纯形规范构造 `u_a = e_a − 1/m·1`（m = D+1）；`peak_argmax` 用精确有理数比较 |

## 判据账（`sim_selftest`，92 条，全过）

```
/usr/bin/python3 zero_sim.py            # 自检 ＋ 写 results/zero_sim.json（退出码 = 失败数>0）
/usr/bin/python3 zero_sim.py --quiet    # 只写文件
```

覆盖（速查）：

| 组 | 判据 |
|:--|:--|
| Heawood 结构 | Fano 7 线/每对点恰一线；\|V\|=14、度 3、\|E\|=21；二部；**围长 6** |
| Heawood 谱 | 谱 = {3, (√2)×6, (−√2)×6, −3}；零化多项式 (A−3)(A+3)(A²−2)=0 |
| Heawood 对称 | \|Aut\| = **336**；顶点稳定子 **24**；稳定子在 3 条入射边上的作用 = **6**（满 S₃） |
| Heawood 输运 | 偶数步精确式 `p(2k) = (9^k + 6·2^k)/(7·9^k)`；偏差率 `(√2/3)^{2k}` ⟹ **质量隙 m = −log(√2/3) = 0.7520386984**、ξ = 1.3297188059；奇数步 = 0（二部图无奇闭道） |
| 词层 | L = 4/8/12/16：`Σ b_h = C(L, L/2)`、`b_h` 恒偶、`b_H = 2`、**Z032(2.1) 镜像闭式与 DP 逐层相符**；L ≤ 12 再与穷举对照；L=4 给 `b = (4, 2)`（层 1 四词、层 2 两词）、层权 2/3 : 1/3 |
| 耦合几何 | D = 1..6：零和、等长（1−1/m）、`ρ = −1/m`（归一化后 −1/D）、**Gram 秩 = D**；D=3 的 ρ = −1/3（四面体角） |
| 峰与窗 | (1/2, 3/5) 含 5/9；`q = 5/9` 时 `F_D = C(D,2)q^D` 唯一峰 D = 4；5/9 不在三维窗 |
| K₄ 与双覆盖 | K₄：度 3、\|E\|=6、围长 3；K₄□K₂：8 顶点 4-正则、谱 {4, 2, 0³, (−2)³} |
| 补偿移动图（Z024/Z072/Z073） | `comp_neighbors` = Johnson 图 J(L,L/2) 的邻接；L=4/6/8/10 **正则、度 = (L/2)²**（4/9/16/25）且**连通**（Z024 不可约）；Z073 结构定理：**上行永远可走**、**下行会被阻**（L=8/10/12 计数 4/32/166，与篇内表逐位相符）；★精化：受阻词并非全在 h=2 —— L=12 起层 3 也出现（L=14：688@h2 + 26@h3） |
| 零面积词（Z048） | `word_area(w) = Σ前缀和` 的**四种写法对平衡词完全等价**（分部求和）；零面积计数 L=4/8/12/16/20 = **2/8/58/526/5448**；L ≡ 2 (mod 4) 恒 0 |
| 3-边着色 | Heawood 的显式 3-边着色：每色 7 边 = 完美匹配、无共端冲突（**这就是 Z021 §8 的三个轴**） |
| K₄ 的 −1 空间 | 维数 3，在 (123) 循环置换下分解为 1 ⊕ ω ⊕ ω²（特征值 {1, ω, ω̄}）—— Z090 |
| Laplacian 谱 | K₄：{0, 4, 4, 4}；Heawood：{0, 3−√2×6, 3+√2×6, 6}（Z002/Z025/Z033/Z102 用的对象） |

## 两个易错口径（脚本会替你把住）

1. **Z032(2.1) 的镜像闭式是"差"**：`b_h = N_≤(L,h) − N_≤(L,h−1)`，且 `N_≤` 的镜像和
   **j 必须取遍 ℤ**（只取 j ≥ 0 会算错——例如 L=4、h=1 需要 j=−1 的 −C(4,0) 项）。
   Z001 的旧式 (7.7) 已判错（步长 h+1 而非 h+2）。
2. **层高口径**：`h(w) = max(M⁺, M⁻) = max_k |前缀和_k|`（不是最大前缀和本身）。
   L=4 的平衡词 6 个：`++−−`、`−−++` 给 h=2，其余 4 个给 h=1。

## 与机器检查（ZeroLean）的关系

同一批结论两侧都有：Lean 侧是 `Zero/Z065_Heawood.lean`（显式 Fano 关联图、度/二部/无 C₄/
6-环/边数/圈秩）、`Zero/Z047_Involution.lean` ＋ `Zero/Z049_Orbits.lean`（层内对合的
二元素轨道）、`Zero/Z021_Dict.lean` ＋ `Z021_Simplex.lean` ＋ `Z021_Coupling.lean`（字典、单纯形、计数）。
**判据互补**：Lean 给"对任意参数成立"的证明；`zero_sim.py` 给"具体数值对得上 + 谱/自同构/闭式"
的实测账。两边不一致就是有 bug——先查 `zero_sim` 的口径（数值侧改起来便宜）。