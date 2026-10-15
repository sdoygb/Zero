# ZeroLean — Zero 推导的机器可检查核心（Lean 4）

这里放**能被内核逐行检查**的那部分 Zero 推导。Python 侧的核验是"数值验算"（跑一遍看数字对不对）；
Lean 侧给的是**证明**：结论由前提经内核检查推出，中间没有"显然"、没有跳步、没有隐藏假设。

## 怎么跑

```bash
export PATH="$HOME/.local/lean4/bin:$PATH"     # 工具链装在 ~/.local/lean4（手动装的官方 release，不需要 elan）
cd /Users/oygb/Downloads/Zero/ZeroLean
lake build                                      # 全绿 = 所有定理被内核接受
```

`#print axioms` 的输出随构建一起打印（那是本项目最有用的东西，见下）。

## 目录

| 文件 | 内容 |
| --- | --- |
| `Zero/Basic.lean` | 第一片：**K₄ 谱** —— 核空间 = 常向量；λ = 4 的特征空间 = 零和超平面 |
| `Zero/AxiomLedger.lean` | 要素来源表的**机械版**演示：自定义 `axiom` 如何出现在依赖清单里 |
| `Zero.lean` | 根模块（import 上面两个） |
| `lakefile.toml` / `lean-toolchain` | Lake 工程（Lean v4.35.0） |

## 为什么从 K₄ 谱开始

Z115 / Z116 的波数 `k = √λ_max(K₄) = 2` **全部**建立在这条谱事实之上；Python 侧只算了特征值。
把它变成定理之后，形式上的边界也就清楚了：

- 定理只说了 **K₄**（4 顶点、度 3）的谱；
- 它**没有**说任何别的图 —— 所以"拿 K₄ 的谱给方格图上的波定波数"在形式上就缺一步
  （这正是 Z117 §一E5 那条判定）。

## `#print axioms` = 机械版「要素来源表」

`lake build` 打印：

```
'Zero.K4L_eq_zero_iff'      depends on axioms: [propext, Quot.sound]
'Zero.K4L_eigen_four_iff'   depends on axioms: [propext, Quot.sound]
'Zero.eigen_four_zero_sum'  depends on axioms: [propext, Quot.sound]
'Zero.AuditDemo.needs_input' depends on axioms: [k_is_two]
'Zero.AuditDemo.no_input'    depends on axioms: [propext, Quot.sound]
```

读法（**这就是 Zero 审计要的那张表**）：

| 依赖清单 | 判定 |
| --- | --- |
| 只有 `propext` / `Quot.sound` / `Classical.choice`（Lean 基石） | 【导出】：不依赖 Zero 的任何输入 |
| 出现自定义 `axiom` / `opaque` / `constant` 的名字 | 【具名输入】：**点名列出**，一个都跑不掉 |
| 出现 `sorryAx` | 证明有洞（CI 应当禁止） |

所以「这条结论是不是全 Zero 内生」在 Lean 里不再靠人读散文判断：**看依赖清单里有没有 Zero 侧的公理**。
工程做法：把 Zero 的每条输入写成一个 `axiom`（或放进 `structure` 的字段/类型类假设），
然后 CI 规定"除基石公理外，只允许出现白名单里的 Zero 公理"。

## Lean 能查什么、不能查什么（对 Zero 的判断）

**能查**：定理是否真的从前提推出来；有没有跳步/循环论证；"显然"是否成立；
常数与算式的正确性（`decide` / `omega` / `norm_num`）；**每条结论到底用了哪些前提**（上面的依赖清单）；
反例是否真的是反例（有限反例最合适）。

**不能查**：形式化对象与物理意图是否一致（把"图上波的波数"编码成什么，是人决定的）；
数值仿真的正确性（那是 Python 的活，Lean 只能查它的数学结论）；
"某要素是内生的还是外部嫁接"这类建模判断 —— 但它能把**输入的清单**摆到你面前，让判断有据可依。

一句话：**Lean 管"从前提推得对不对"，不管"前提选得对不对"。** Zero 现在的缺口（能量标度从哪来、舞台图是什么）属于后者，Lean 不会替你补；但它能把这两类问题彻底分开。

## 网络与镜像（本机实测，2026-10）

| 源 | 结果 |
| --- | --- |
| GitHub release 资产（直连 `github.com/.../releases/download/...`） | ✗ 超时（API 与网页能通，资产不通）|
| `https://ghproxy.net/<原始 github 链接>` | ✓ **最快 ~18 MB/s**（576 MB 用 60 s）|
| `https://gh-proxy.com/...` | ✓ 可用但慢 ~110 KB/s |
| `hub.gitmirror.com` / `gh.llkk.cc` / `github.moeyy.xyz` | ✗ 不可用 |
| `api.github.com`（查版本）| ✓ |
| `leanprover-community.github.io`（Mathlib olean 缓存源）| ✓ 可达 |
| `reservoir.lean-lang.org` | ✓ 可达 |

手动装工具链（不用 elan，因为 elan 拉工具链同样走被卡的 release 通道）：

```bash
curl -L -o lean.tar.zst "https://ghproxy.net/https://github.com/leanprover/lean4/releases/download/v4.35.0/lean-4.35.0-darwin.tar.zst"
mkdir -p ~/.local && tar --zstd -xf lean.tar.zst -C /tmp/leantmp && mv /tmp/leantmp/* ~/.local/lean4
```

## 下一步（按价值排序）

1. **接 Mathlib**（`require mathlib` + `lake exe cache get`，缓存源可达）：有了它才能写 `Matrix (Fin n) (Fin n) ℤ`、
   一般图 Laplacian、特征值与重数、`Real.cos`/`Real.sqrt`/`π` —— 也就是把"λ(q) = 4 − 2cos qx − 2cos qy"、
   "L1 ≠ 欧氏程函"、"亚 Nyquist 混叠"这些**承重**命题形式化。
2. **排形式化队列**（从 Zero 现行结论里挑，先小后大）：
   - K₄ 谱的重数 3（维数/线性代数，Mathlib）；
   - 路径图/方格图的 Laplacian 谱与色散（Mathlib）；
   - `C(D,2) = D ⟺ D ∈ {0,3}`（Z021 的维数字典；Mathlib 一行 `nlinarith` 级别）；
   - 亚 Nyquist：采样后 4 与 4−2π 不可区分（复指数 + 整数采样，Mathlib）；
   - Z004 里那条"Gleason 假设在有限格上不成立"的**有限反例**（有限维线性代数，Mathlib）。
3. **CI 守则**：禁止 `sorryAx`；`#print axioms` 白名单（只允许基石公理 + 明确登记的 Zero 公理）；
   每条定理与它对应的 Z 篇号建立映射（文件名/定理名带篇号）。
4. **接进 dsh-zerovault 插件**：加一个 `zero_lean` 工具（跑 `lake build`，把错误与依赖清单回给模型），
   并在 post-execute 里对 `ZeroLean/**/*.lean` 的改动自动构建 —— 这样"写一篇 → 机器检查一遍"就成了默认动作。

---

## 2026-10-14 更新：Mathlib 已接入 + 形式化队列 + CI

### 现状

- **Mathlib 已接入**（`require mathlib`），源码与 oleans 共 **7.2 GB**，放在**外接卷**
  `/Volumes/download/lean-build/ZeroLean-lake`；库内 `.lake` 是指向它的符号链接（主盘只剩 23 GB，所以构建产物走外接盘）。
  Mathlib 的 olean 缓存已自动下完（9034 个文件），工具链 `~/.local/lean4`（Lean 4.35.0，与 Mathlib master 要求的版本一致）。
- **形式化队列五项全部完成**（48 条定理，全部 `lake build` 通过）：

| 文件 | 定理 | 内容 |
| --- | --- | --- |
| `Zero/Basic.lean` | `K4L_eq_zero_iff` / `K4L_eigen_four_iff` | K₄ 核 = 常向量；λ=4 特征空间 = 零和超平面 |
| `Zero/Z104_Dispersion.lean` | `z104_grid_dispersion` | 方格图色散 `λ(q) = 4 − 2cos qx − 2cos qy` |
| `Zero/Z021_Dimension.lean` | `z021_choose_two_eq_self` / `z021_dim_forced_three` | `C(D,2) = D ⟺ D ∈ {0,3}` |
| `Zero/Z117_Aliasing.lean` | `z117_aliasing` / `z117_beyond_nyquist` | 亚 Nyquist：采样后 4 与 4−2π 不可区分；`π < 4` |
| `Zero/Z117_Metric.lean` | `z117_euclid_ne_l1` | L1 ≠ 欧氏程函（离轴即分道，共线行才重合）|
| `Zero/Z117_CrossGraph.lean` | `z117_grid_lambda_at_pi` 等 | 方格图色散上界 8 ≠ K₄ 的 λ_max 4（跨图取值不可行）|
| `Zero/Z004_Gleason.lean` | `z004_born_one_forces_parallel` | Born {0,1}-赋值必集中于一条射线（Gleason 前提在二维失效的最小片段）|
| `Zero/AxiomLedger.lean` | `needs_input` / `no_input` | 要素来源表的机械版演示 |

### CI：`scripts/check.py`

```bash
python3 scripts/check.py            # lake build + 自动审计 + 三条守则
python3 scripts/check.py --json     # 机器可读（插件 zero_lean 用的就是这个）
```

三条守则（任一违反 → 非零退出）：

- **R1 禁 `sorryAx`** —— 证明不许有洞；
- **R2 公理白名单** —— `scripts/axioms-allow.txt`（**写全名**，含命名空间）；未登记的公有理 = 未登记的具名输入，会被点名列出；
- **R3 定理名对 Z 篇号** —— `zNNN_*` 的 NNN 必须对应库里真实存在的 Z 篇。

审计**自动遍历环境里所有 `Zero.*` 定理**（`Zero/Audit.lean` 用 `collectAxioms`），
不依赖作者手写 `#print axioms` —— 早期版本靠手写，结果漏掉了故意植入的违规样例，
改成自动后立刻多抓到一条先前漏审的定理。三类违规都做了负测试（植入 `sorry`、未登记公理、假篇号 → 全部拦下，exit 1）。

### 磁盘布局（主盘紧张）

| 内容 | 位置 | 占用 |
| --- | --- | --- |
| Lean 工具链 | `~/.local/lean4` | 2.8 GB（主盘）|
| Mathlib + 依赖 + oleans | `/Volumes/download/lean-build/ZeroLean-lake`（外接 APFS）| 7.2 GB |
| 源码 | `<vault>/ZeroLean/` | 几百 KB |
