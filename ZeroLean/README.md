# ZeroLean — Zero 推导的机器可检查核心（Lean 4）

这里放**能被内核逐行检查**的那部分 Zero 推导。Python 侧的核验是「数值验算」（跑一遍看数字对不对）；
Lean 侧给的是**证明**：结论由前提经内核检查推出，中间没有「显然」、没有跳步、没有隐藏假设。

## 现状（2026-10-17 机器检查）

| 项 | 值 |
|:--|:--|
| `lake build` | **通过** |
| 审计 `Zero.*` 定理数 | **807 条** |
| CI 守则（R1／R2／R3） | **全过** |
| 源码 | `Zero/` 下 **88 个 `.lean`**，`Zero.lean` 逐一 import |
| 公理白名单 | 4 个名字：`propext`／`Classical.choice`／`Quot.sound`／`Zero.AuditDemo.k_is_two` |
| 工具链 | Lean v4.35.0（`lean-toolchain`），装在 `~/.local/lean4` |
| Mathlib | 已接入（`lakefile.toml` 的 `require mathlib`）；构建产物在外部卷，`ZeroLean/.lake` 是指向它的符号链接 |

机器可读记录：`.lean-report.json`（`{"ok": true, "theorems": 807, "problems": []}`）。

## 怎么跑

```bash
export PATH="$HOME/.local/lean4/bin:$PATH"     # 工具链在 ~/.local/lean4（手动装的官方 release，不需要 elan）
cd /Users/oygb/Downloads/Zero/ZeroLean
lake build                                      # 全绿 = 所有定理被内核接受
python3 scripts/check.py                        # CI：build ＋ 自动审计 ＋ 三条守则
python3 scripts/check.py --json                 # 机器可读（插件用的就是这个）
```

`#print axioms` 的输出随审计一起给出（那是本项目最有用的东西，见下）。

## 目录

| 路径 | 内容 |
|:--|:--|
| `Zero.lean` | 根模块，逐一 import `Zero/` 下的全部文件 |
| `Zero/Basic.lean` | K₄ 谱：核空间＝常向量；$\lambda=4$ 的特征空间＝零和超平面 |
| `Zero/AxiomLedger.lean` | 要素来源表的**机械版**演示：自定义 `axiom` 如何出现在依赖清单里 |
| `Zero/Z<篇号>_*.lean` | 逐篇形式化（88 个文件，807 条定理；文件名带篇号，CI 据此对篇号） |
| `Zero/Audit.lean` | 审计机关：`collectAxioms` 遍历环境里全部 `Zero.*` 定理 |
| `scripts/check.py` | CI：`lake build` ＋ 自动审计 ＋ 三条守则 |
| `scripts/axioms-allow.txt` | 公理白名单（写全名，含命名空间） |
| `scripts/coverage_matrix.py` | 生成形式化覆盖矩阵（读 `Zero/*.lean` ＋ 本地 `_staging` 覆盖账） |
| `wip/` | **未接入构建**的草稿与探针（见 [`wip/README.md`](wip/README.md)） |
| `lakefile.toml` / `lake-manifest.json` / `lean-toolchain` | Lake 工程与版本锁定 |

## CI：`scripts/check.py` 的三条守则

任一违反 → 非零退出：

- **R1 禁 `sorryAx`** —— 证明不许有洞；
- **R2 公理白名单** —— `scripts/axioms-allow.txt`（**写全名**，含命名空间）；未登记的公有理＝未登记的具名输入，会被点名列出；
- **R3 定理名对 Z 篇号** —— `zNNN_*` 的 NNN 必须对应库里真实存在的 Z 篇。

审计**自动遍历环境里所有 `Zero.*` 定理**（`Zero/Audit.lean` 用 `collectAxioms`），不依赖作者手写 `#print axioms`——早期版本靠手写，漏掉了故意植入的违规样例，改成自动后立刻多抓到一条先前漏审的定理。三类违规都做过负测试（植入 `sorry`、未登记公理、假篇号 → 全部拦下，exit 1）。

## `#print axioms` = 机械版「要素来源表」

审计给出的依赖清单读法（**这就是 Zero 审计要的那张表**）：

| 依赖清单 | 判定 |
|:--|:--|
| 只有 `propext`／`Quot.sound`／`Classical.choice`（Lean 基石） |【导出】：不依赖 Zero 的任何输入 |
| 出现自定义 `axiom`／`opaque`／`constant` 的名字 |【具名输入】：**点名列出**，一个都跑不掉 |
| 出现 `sorryAx` | 证明有洞（CI 禁止） |

所以「这条结论是不是全 Zero 内生」在 Lean 里不靠人读散文判断：**看依赖清单里有没有 Zero 侧的公理**。工程做法：把 Zero 的每条输入写成一个 `axiom`（或放进 `structure` 的字段／类型类假设），CI 规定「除基石公理外，只允许出现白名单里的 Zero 公理」。

## Lean 能查什么、不能查什么

**能查**：定理是否真的从前提推出来；有没有跳步／循环论证；「显然」是否成立；常数与算式的正确性（`decide`／`omega`／`norm_num`）；**每条结论到底用了哪些前提**（依赖清单）；反例是否真的是反例（有限反例最合适）。

**不能查**：形式化对象与物理意图是否一致（把「图上波的波数」编码成什么是人决定的）；数值仿真的正确性（那是 Python 的活）；「某要素是内生的还是外部嫁接」这类建模判断——但它能把**输入的清单**摆到面前，让判断有据可依。

一句话：**Lean 管「从前提推得对不对」，不管「前提选得对不对」。**

## 安装与网络（本机实测，2026-10）

| 源 | 结果 |
|:--|:--|
| GitHub release 资产（直连 `github.com/.../releases/download/...`） | ✗ 超时（API 与网页能通，资产不通） |
| `https://ghproxy.net/<原始 github 链接>` | ✓ 最快 ~18 MB/s |
| `https://gh-proxy.com/...` | ✓ 可用但慢 ~110 KB/s |
| `hub.gitmirror.com`／`gh.llkk.cc`／`github.moeyy.xyz` | ✗ 不可用 |
| `api.github.com`（查版本） | ✓ |
| `leanprover-community.github.io`（Mathlib olean 缓存源） | ✓ 可达 |
| `reservoir.lean-lang.org` | ✓ 可达 |

工具链手动装（不用 elan，因为 elan 拉工具链同样走被卡的 release 通道）：

```bash
curl -L -o lean.tar.zst "https://ghproxy.net/https://github.com/leanprover/lean4/releases/download/v4.35.0/lean-4.35.0-darwin.tar.zst"
mkdir -p ~/.local && tar --zstd -xf lean.tar.zst -C /tmp/leantmp && mv /tmp/leantmp/* ~/.local/lean4
```

## 磁盘布局（主盘紧张）

| 内容 | 位置 |
|:--|:--|
| Lean 工具链 | `~/.local/lean4`（主盘） |
| Mathlib＋依赖＋oleans（7.2 GB） | 外部卷的 `ZeroLean-lake`；`ZeroLean/.lake` 是指向它的符号链接 |
| 源码 | `ZeroLean/`（几百 KB） |

## 下一步

- **按覆盖矩阵推进缺口**：`python3 scripts/coverage_matrix.py` 读出逐篇覆盖（已证／具名输入／不可形式化）与缺口清单，按价值排序补证明。
- **新增篇目随篇形式化**：新 `Zero/Z<篇号>_*.lean` 直接进 `Zero.lean`，CI 自动覆盖。
