#!/usr/bin/env python3
"""Zero 记忆写入判官（dsh-memory 外部钩子）。

为什么存在
----------
记忆条目是**散文**，插件只能校验格式与凭据，无法判断内容真假。
Zero 的推导有两个反复出现的污染源，都是"写的时候没问题、读的时候已经错了"：

  ① **参数数值**被抄进记忆。数值会随修订漂移，且同一个数在不同语境里可以
     合法地不同（`1/c_ℓ` 按代不同：784.493/784.501/784.549；`137.905` 是
     Ba-138 质量不是 α⁻¹）。记忆里留一个旧值 ⟹ 下游把它当现值。
  ② **档位断言**（某量属 D/C/F/I、或"已导出/已证/开放"）被写成记忆事实。
     档位是**写作时点**的记录；后来的篇把它改判后，记忆里的旧档位仍然在。
     实例：Z167 曾把 `sin²θ_W(M_U)=3/8` 记作 D，而 Z151 §7 自陈唯一非 Zero
     输入是"迹范数相等"，按定义应记 C。

所以本判官只做两件事：**拦下裸数值、要求档位断言带指针**。
其余一律放行——它必须**失败即放行**，不能因为自己出错就堵死记忆写入。

协议（与插件 README 一致）
--------------------------
stdin  : payload JSON（tool / path / target / layer / added[] / removed[] / …）
stdout : {"decision": "allow"|"deny"|"ask", "reason": "..."}
exit 0 = 放行或按 stdout 决策；exit 2 = 硬拒绝（stderr 为原因）。
任何异常 ⟹ exit 0 + allow（放行），并在 note 里说明，绝不阻塞。

用法
----
settings.json（~/.dsh/memory/settings.json）：绝对路径，库改名／搬库后必须同步改
（用 `tools/zero_vault_migrate.py --apply` 一步改完并备份）：

    "writeGuard": "full",
    "writeHookCommand": "python3 <库根>/tools/zero_mem_judge.py"

自测：

    python3 tools/zero_mem_judge.py --selftest
"""

from __future__ import annotations

import json
import re
import sys

# ── 结构常数名（出现这些名 + 数值 ⟹ 疑似把参数值抄进记忆）──────────────────
CONST_NAMES = [
    r"S_e", r"1/c_\s*\u2113", r"1/c_\\ell", r"c_\s*\u2113", r"c_\\ell",
    r"\\delta_0", r"\bdelta0\b", r"\\Delta\\Theta", r"\\Lambda_H",
    r"\\alpha\^\{-1\}", r"\balpha_inv\b", r"\bm_p\b", r"\bm_e\b",
    r"\bM_e\b", r"\\tau_\{\\rm step\}", r"\\chi_T", r"\\chi_L",
]

# 结构常数的**已知数值**：即使没写名字，出现这些数也视为"在抄参数值"。
# 与 CONST_NAMES 分开是刻意的——名字是"在说什么"，值是"抄了多少"。
KNOWN_VALUES = [
    r"\b137\.03\d*", r"\b784\.5\d*", r"\b434\.7\d*", r"\b313\.8\d*",
    r"\b1836\.\d+", r"\b0\.0395\d*", r"\b0\.7078\d*", r"\b0\.0403\d*",
    r"\b9\.11\d*e-22", r"\b2\.73\d*e-13",
]

# 裸数值：≥6 位有效数字的浮点，或 4 位以上小数的数
NUM_RE = re.compile(r"(?<![\w.])\d+\.\d{4,}(?![\w])")

# 档位断言词
TIER_RE = re.compile(
    r"【已证】|【导出】|【条件】|【输入】|"
    r"\bD\s*档|\bC\s*档|\bF\s*档|\bI\s*档|"
    r"判定为导出|判定为输入|属\s*[\u4e00-\u9fa5]{0,3}档"
)

# 作用域：只在 Zero 项目层的记忆文件上生效
ZERO_SCOPE = "--Users-oygb-Downloads-Zero--"

# 指针：篇号、脚本名、章节
POINTER_RE = re.compile(r"Z\d{3}|§\s*\d|[\w\-]+\.py|[\w\-]+\.json|STATUS")


def _norm(s: str) -> str:
    return s.replace(" ", "")


def judge(payload: dict) -> tuple[str, str]:
    """返回 (decision, reason)。decision ∈ allow/deny/ask。"""
    # ── 作用域门：只对 Zero 项目层生效 ─────────────────────────────────
    # 本判官的两条规则都是 Zero 特有的（结构常数名、D/C/F/I 四档）；
    # 全局 topics 与别的项目可能合法地写数值/别的档位，不能一刀切。
    path = str(payload.get("path") or "")
    if ZERO_SCOPE not in path:
        return "allow", ""

    added = payload.get("added") or []
    if not added:
        return "allow", ""

    text = "\n".join(added)

    # ── 规则 1：结构常数名与数值同段 ⟹ ask（指向权威位置，不硬拒）──────────
    # 用 ask 而不是 deny：确有必须记录数字的场合（如"当年记为 807 条"），
    # 这类应当由人/模型自己确认，而不是被一刀切挡掉。
    for line in added:
        has_name = any(re.search(pat, line) for pat in CONST_NAMES)
        has_num = bool(NUM_RE.search(line)) or bool(
            re.search(r"(?<![\w.])\d+\.\d{2,}(?![\w])", line)
        )
        has_known_value = any(re.search(v, line) for v in KNOWN_VALUES)
        if (has_name and has_num) or has_known_value:
            return "ask", (
                "dsh-memory/Zero 判官：这一条把**结构常数的数值**写进记忆了。"
                "数值只应住在 Lean／ZeroGPU/src／ZeroGPU/results／STATUS，"
                "记忆里请只记「值住在哪（Z篇号§节 或 脚本名）＋ 为什么是它」。"
                "同一个数在不同语境里可以合法地不同（1/c_ℓ 按代不同；137.905 是 Ba-138 质量），"
                "所以记忆里留旧值会被下游当现值读。"
                f" → 命中行：{line.strip()[:90]}"
            )

    # ── 规则 2：档位断言没有指针 ⟹ ask ─────────────────────────────────
    # 档位是写作时点的记录（引用纪律），脱离了篇号指针就会陈旧。
    for line in added:
        if TIER_RE.search(line) and not POINTER_RE.search(line):
            return "ask", (
                "dsh-memory/Zero 判官：这一条断言了**档位**（已证／导出／条件／输入／D-C-F-I），"
                "但没给指针。档位是写作时点的记录，篇被改判后记忆里的旧档位仍在"
                "（实例：3/8 曾记作 D，实为 C）。请补 Z篇号§节 指针，或改写成「判定见 Zxxx」。"
                f" → 命中行：{line.strip()[:90]}"
            )

    return "allow", ""


def _selftest() -> int:
    cases = [
        # (说明, added, 期望决策)
        ("裸常数名+数值", ["- [2026-10-28] S_e = 137.036，与实测一致"], "ask"),
        ("精确小数无常数名", ["- [2026-10-28] 脚本给出 0.03957717 的偏差"], "ask"),
        ("任意长小数（非已知常数）", ["- [2026-10-28] 拟合优度 0.98712345"], "allow"),
        ("常数名+指针无数值", ["- [2026-10-28] S_e 的取值见 Z164 §2"], "allow"),
        ("档位断言无指针", ["- [2026-10-28] 3/8 属 D 档"], "ask"),
        ("档位断言带指针", ["- [2026-10-28] 3/8 属 C 档（Z167 §2.3）"], "allow"),
        ("普通条目", ["- [2026-10-28] 拆分片后检索成本下降一个量级"], "allow"),
        ("常数名+数值+指针仍拦", ["- [2026-10-28] S_e = 137.036（Z164）"], "ask"),
        ("空 added", [], "allow"),
    ]
    bad = 0
    for desc, added, want in cases:
        got, reason = judge({"added": added, "path": "/x/" + ZERO_SCOPE + "/t.md"})
        ok = got == want
        if not ok:
            bad += 1
        print(f"  {'ok  ' if ok else 'FAIL'}  {desc}: got={got} want={want}")
        if not ok and reason:
            print(f"         reason={reason[:100]}")
    print(f"\n{'全部通过' if bad == 0 else f'{bad} 项不符'}")
    return 0 if bad == 0 else 1


def main() -> int:
    if "--selftest" in sys.argv:
        return _selftest()

    try:
        raw = sys.stdin.read()
        payload = json.loads(raw) if raw.strip() else {}
        decision, reason = judge(payload)
    except Exception as exc:  # 失败必须放行，不能堵死记忆写入
        print(json.dumps({"decision": "allow",
                          "note": f"zero_mem_judge 异常，已放行：{type(exc).__name__}"},
                         ensure_ascii=False))
        return 0

    out = {"decision": decision}
    if reason:
        out["reason"] = reason
    print(json.dumps(out, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
