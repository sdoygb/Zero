# `verify/` —— **借入参照层**（不是本理论的核验层）

> **【状态归属｜[`STATUS.md`](../STATUS.md)】** 本目录只核验借入的旧理论 `D210–D259`，其“未解输入”按 **D259 路线本地状态**读取。它们既不推翻主 `Z/G` 路线的条件恢复，也不能被主路线的 E1–E4 自动结清；两条路线尚未证明等价。全项目当前状态见 [`STATUS.md`](../STATUS.md)。

**来源**：本目录 51 个文件（50 个 `d2xx_*.py` ＋ `latex_utils.py`）逐字节复制自
`../modular-equilibrium/verify/`。它们核验的是**旧理论**的 `D210`–`D259`，不是本理论的 G 系列。

**为什么在 `lh/` 里也有一份**：本语料的 `D2xx_*.md` 本身也是从那边的 `derivations/` 复制来的，
但复制时**只搬了正文、没搬核验脚本**，导致 50 篇 D 文档的
核验行 `**核验**：verify/d2xx_*.py —— **N 通过 / 0 不符**` 在本语料里**无从核验**。
本目录把那一层补齐，使这些声明**在语料内可复跑**。

**一处必要的改动**：这 50 个脚本原本用

```python
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
```

读取**旧理论的** `AXIOMS.md` 与 `derivations/`。搬到 `lh/verify/` 后该路径会解析到 `lh/`，
那些文件在这里并不存在（**不应**把旧理论的公理与推导搬进本语料）。故每个脚本的 ROOT 改为：

```python
# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
```

即：**脚本在 `lh/`，被核验的语料留在原仓库**。这是"引用"而非"吞并"。

**实跑**（50 个全部通过，累计 **1340 项「通过」**，退出码全 0）：

```bash
cd /Users/oygb/Downloads/lh
for f in verify/d2*.py; do python3 "$f" | tail -1; done
```

**它们不进入 `G10` 账本的合计**：账本只统计 `G*_check.py`（`ledger_sync.py` 的 glob），
且本层的结论等级是**旧理论的**，与本理论的主张必须分开记账。
本层在 `G10` §6 末尾单列说明。
