# 旧库盘点：cosmos-construct × modular-equilibrium（ZeroGPU 讨论稿）

**性质**：讨论稿／盘点稿，不占正式篇号；只记录**实读**到的内容，不替旧库发言。
**对象**：
- `~/Downloads/cosmos-construct`（3327 文件，119 MB；git 103 commits，仅 `main`）
- `~/Downloads/modular-equilibrium`（847 文件，16 MB；无 git）

**结论先说**：

$$
\boxed{\;\text{这两个目录是\textbf{物理／数学}旧库，不是社会／政治／军事研究};\quad
\text{唯一成体系的\textbf{人文}内容是佛教／华严设计注记}.\;}
\tag{0.1}
$$

---

## §1 实际内容（按文件与主题）

### 1.1 `cosmos-construct`

| 项 | 数 |
|:--|--:|
| 顶层条目 | 468 |
| `*.md` | 461 |
| `build/*.py` | 429 |
| `design2/*.py` | 76 |
| `pdf_export/` | 391 |
| `Q*.md`（Q1–Q860） | 数百 |

**主题**：$\alpha$ 公式与六项式、$E_6$／$SU(3)^3$／GUT、质量与 Koide、中微子、GR 六环节、宇宙学（KK 塔、$\Omega_k$、热历史）、谱作用量、量子信息玩具（K2/K4/K5）。典型标题：`ALPHA_FORMULA.md`、`GR_FROM_GEOMETRY.md`、`Q620_STRUCTURAL_UNIFICATION.md`、`Q860_THERMO_AUDIT.md`。

### 1.2 `modular-equilibrium`

| 项 | 数 |
|:--|--:|
| `AXIOMS.md` | 252 KB |
| `derivations/D*.md` | 259（D100–D259 等） |
| `simulations/` | 48 项（含 `zero_sum_*`） |
| `verify/*.py` | 263 |

**主题**：模流／Tomita–Takesaki、模 Hamilton、支持与面积律、UV 条件、谱作用量、GR 恢复门、层读出。典型标题：`D216_history_phase_state_for_u3.md`、`D248_interlayer_readout_dynamics.md`、`D129_area_law_support_no_go.md`。

---

## §2 社会／政治／军事：没有

**三条独立检查**（都不依赖"我记得"）：

1. **文件名**：两个目录的文件名（含 461+259 个 md、692 个 py）**没有一个**含「社会／政治／军事／战争／文明／制度／国家／权力／宗教／伦理／经济／法律／战略／治理／城邦／帝国」。
2. **标题行**：扫 `*.md` 的 `#` 标题，只有两处命中，且都是物理用法：
   - `Q496` 的「**战略回撤**」（方法论上的撤退）；
   - `Q349` 的「膨胀**历史**奇异」。
3. **正文词频**：社会词命中最高的文件是 `Q_REGISTRY.md`（83 次）与 `AXIOMS.md`（84 次），逐行看**全是物理用法**：
   - 「**民主**组合」= 表示论中的界对称组合（不是政治民主）；
   - 「**历史**」= 历史相位／记录历史；
   - 「**冲突**」= 两条结论不一致（不是军事冲突）；
   - 「**制度**」= 方法论上的常规动作（如"把退化检验制度化"）。

> 结论：**"cosmos-construct / modular-equilibrium 有大量社会、人文、政治、军事研究"这一前提，对这两个目录不成立。** 后续全目录检索发现：社会研究实际在 **`GeometryAI-Mac-Build/app/articles/10.75–10.107`**（30+ 篇）、`gpu_spectra/society/`（49+ 探针）、`notes/社会版·破缺谱_工作稿_260915.md`（2137 行），以及旧库 **`Zero-261010/H1–H10`＋`port/`**；已接 [`zero_social_port.md`](zero_social_port.md)（13/13）。

---

## §3 实际存在的「人文」子集：佛教／华严

| 文件 | 内容 | 机器核验 |
|:--|:--|:--|
| `NEW_UNIVERSE_DESIGN_BUDDHIST.md` | 刹那生灭／缘起／无我／二谛／无始无作 → 五条设计约束 | `design2/k2_buddhist_abhidharma.py`（11/0） |
| `NEW_UNIVERSE_DESIGN_HUAYAN.md` | 华严卷7/9/45 的嵌套、因陀罗网、十方入一、数量级递推 | `design2/k5_huayan_nesting.py`（全过） |
| `NEW_UNIVERSE_DESIGN_HUAYAN_CRITICAL.md` | 三译本层判定（早层／译名层／后世诠释层） | 文本底账 |
| `design2/k4_indras_net.py` | 因陀罗网的精确形式：GHZ 一与整体、非一与单珠 | 11/0 |

**三条硬结果**（旧库自己算的，本库已复算前两条）：

- **K5**：GHZ 任意非空真子集约化态秩恒 2、$S(A)=1$ bit，与 $m,N$ 无关；W 态 $S(A)=H_2(1/N)<1$，随 $N$ 变。
- **K4**：GHZ 的 $I(A{:}B)=I(A{:}C)=1$ bit，但两体 concurrence $=0$；「一即一切」的精确版是**一与整体**，不是一与任一单珠。
- **K2**：刹那 ⇒ 离散时间到连续时间的收敛阶 $1.995$（理论 2）。

---

## §4 能接什么、不能接什么

| 旧库内容 | 与 Zero 的关系 | 处置 |
|:--|:--|:--|
| $\alpha$／$E_6$／GUT／质量／GR／宇宙学 | 与 Zero 的正式篇大量重叠 | **不重复导入**；只在需要时对照 |
| 模流／面积律／支持 | 与 Z005/Z044/Z143 同族 | **不重复导入** |
| 佛教／华严（人文子集） | 与 Z008/Z036（代数关联、无信号、无双分划）、Z161/Z164（层与亚层）有结构对应 | **已接**：`zero_huayan_net.py`＋`zero_huayan_net.md`（10/10） |
| 社会／政治／军事 | **不存在** | 不伪造；见 `zero_social_bridge_proposal.md`（新线提案） |

---

## §5 复跑与文件

```bash
cd /Users/oygb/Downloads/cosmos-construct && /usr/bin/python3 design2/k2_buddhist_abhidharma.py
cd /Users/oygb/Downloads/cosmos-construct && /usr/bin/python3 design2/k4_indras_net.py
cd /Users/oygb/Downloads/cosmos-construct && /usr/bin/python3 design2/k5_huayan_nesting.py
cd /Users/oygb/Downloads/Zero/ZeroGPU/src && /usr/bin/python3 zero_huayan_net.py
```

| 文件 | 内容 |
|:--|:--|
| `ZeroGPU/src/zero_huayan_net.py` | 华严／因陀罗网 → Zero 读出：GHZ 局部身份、一与整体、无信号、嵌套计数（10/10） |
| `ZeroGPU/doc/zero_huayan_net.md` | 映射与边界 |
| `ZeroGPU/doc/zero_social_bridge_proposal.md` | 社会／政治／军事新线提案（未导入，只给映射与判据） |

*ZeroGPU · 讨论稿 · 2026-10-22*
