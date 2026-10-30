#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_vault_migrate.py —— 库目录改名／换库的外部配置迁移（默认 dry-run）

为什么需要它
------------
Zero 库的“根路径”散落在**三处外部配置**与**一批脚本字面量**里。只做 `mv`，
下列东西会静默失联（不报错，只是搜不到／守卫失效／结果写错库）：

  ① `~/.dsh/profiles/desktop/cordis.patch.yml` —— dsh-zerovault-v3／v4 的
     `vault:` 字段（各一行）。失联现象：插件的检索门与库状态指向空目录。
  ② `~/.dsh/memory/settings.json` —— `writeHookCommand`（记忆写入判官）。
     失联现象：记忆写入的数值／档位守卫**不再执行**，且不会有人告诉你。
  ③ 库内脚本的硬编码绝对前缀（如 `OUT = "/Users/oygb/Downloads/Zero1030/ZeroGPU/results/x.json"`）。
     失联现象：脚本在新库里跑，结果写回**另一个**库，两库互相污染。

本工具只做“改配置”，**不搬目录、不动正文**。搬目录由人执行（见
`_staging/重写/改名预案.md`）。

用法
----
    python3 tools/zero_vault_migrate.py --plan                 # 只报告现状（默认）
    python3 tools/zero_vault_migrate.py --apply --new-root /Users/oygb/Downloads/Zero
    python3 tools/zero_vault_migrate.py --verify --new-root …  # 只读核验

纪律
----
• 默认 dry-run：不加 `--apply` 绝不写盘。
• 每次写盘前把目标文件复制成 `<name>.bak-<UTC时间戳>`，并打印备份路径。
• 幂等：已经是目标值就不再写。
• 任何一步失败都**不继续**后续步骤（先改坏一处再改另一处会更难查）。
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import time
from pathlib import Path

HOME = Path.home()
CORDIS = HOME / ".dsh" / "profiles" / "desktop" / "cordis.patch.yml"
MEMSET = HOME / ".dsh" / "memory" / "settings.json"

VAULT_LINE = re.compile(r"^(\s*vault:\s*)'([^']*)'\s*$")


# ─────────────────────────── 探测 ───────────────────────────

def probe_cordis(path: Path):
    """返回 [(行号, 旧值)]，只认 `vault: '…'` 这种形态（v3/v4 各一行）。"""
    if not path.exists():
        return None, []
    hits = []
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        m = VAULT_LINE.match(line)
        if m:
            hits.append((i, m.group(2)))
    return path, hits


def probe_memset(path: Path):
    if not path.exists():
        return None, None
    try:
        cfg = json.loads(path.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return path, "<无法解析 JSON>"
    return path, cfg.get("writeHookCommand")


HARDCODE_RE = re.compile(r"/Users/[^\"'\s]*/Zero(?=[/\"'])")


def probe_hardcodes(tree: Path, limit: int = 0):
    """扫 tree 下 .py 里的硬编码 `/…/Zero` 前缀（只看代码，不碰 .md 正文）。"""
    rows = []
    if not tree.exists():
        return rows
    for p in sorted(tree.rglob("*.py")):
        if any(part in {".git", "__pycache__", "node_modules"} for part in p.parts):
            continue
        try:
            txt = p.read_text(encoding="utf-8")
        except Exception:  # noqa: BLE001
            continue
        n = len(HARDCODE_RE.findall(txt))
        if n:
            rows.append((p, n))
    rows.sort(key=lambda r: -r[1])
    return rows[:limit] if limit else rows


# ─────────────────────────── 写盘 ───────────────────────────

def _backup(path: Path) -> Path:
    stamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime())
    bak = path.with_name(path.name + f".bak-{stamp}")
    shutil.copy2(path, bak)
    return bak


def apply_cordis(path: Path, new_root: str, dry: bool) -> int:
    if not path.exists():
        print(f"  ✗ 缺 {path}（跳过）")
        return 0
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    changed = []
    for i, line in enumerate(lines):
        m = VAULT_LINE.match(line.rstrip("\n"))
        if not m:
            continue
        old = m.group(2)
        if old == new_root:
            continue
        indent, key = m.group(1), m.group(1)[: len(m.group(1)) - len("vault: ")]
        lines[i] = f"{indent}'{new_root}'\n"
        changed.append((i + 1, old))
    if not changed:
        print("  = cordis.patch.yml 已是目标值（无改动）")
        return 0
    for ln, old in changed:
        print(f"  · 行 {ln}: {old}  →  {new_root}")
    if dry:
        print("  （dry-run：未写）")
        return len(changed)
    bak = _backup(path)
    print(f"  备份 → {bak}")
    path.write_text("".join(lines), encoding="utf-8")
    return len(changed)


def apply_memset(path: Path, judge: Path, dry: bool) -> int:
    if not path.exists():
        print(f"  ✗ 缺 {path}（跳过）")
        return 0
    cfg = json.loads(path.read_text(encoding="utf-8"))
    old = cfg.get("writeHookCommand")
    new = f"python3 {judge}"
    if old == new:
        print("  = settings.json 的 writeHookCommand 已是目标值")
        return 0
    print(f"  · writeHookCommand: {old}\n      → {new}")
    if dry:
        print("  （dry-run：未写）")
        return 1
    bak = _backup(path)
    print(f"  备份 → {bak}")
    cfg["writeHookCommand"] = new
    path.write_text(json.dumps(cfg, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 1


def apply_hardcodes(tree: Path, new_prefix: str, dry: bool) -> int:
    """把 `/…/Zero` 前缀改成 `new_prefix`。用于搬库后统一脚本里的绝对路径。"""
    rows = probe_hardcodes(tree)
    if not rows:
        print("  = 无硬编码前缀")
        return 0
    total = sum(n for _, n in rows)
    print(f"  命中 {len(rows)} 个 .py、共 {total} 处")
    for p, n in rows[:10]:
        print(f"    · {p.relative_to(tree)} ×{n}")
    if len(rows) > 10:
        print(f"    …其余 {len(rows) - 10} 个")
    if dry:
        print("  （dry-run：未写）")
        return total
    written = 0
    for p, _n in rows:
        txt = p.read_text(encoding="utf-8")
        new_txt = HARDCODE_RE.sub(new_prefix, txt)
        if new_txt != txt:
            p.write_text(new_txt, encoding="utf-8")
            written += 1
    print(f"  已改写 {written} 个文件")
    return written


# ─────────────────────────── 核验 ───────────────────────────

def verify(new_root: Path, judge: Path, cordis: Path = None, memset: Path = None) -> int:
    cordis = Path(cordis) if cordis else CORDIS
    memset = Path(memset) if memset else MEMSET
    fails = []
    root_s = str(new_root)

    p, hits = probe_cordis(cordis)
    if p is None:
        fails.append(f"缺 {CORDIS}")
    else:
        bad = [(ln, v) for ln, v in hits if v != root_s]
        if bad:
            fails.append(f"cordis.patch.yml 仍有 {len(bad)} 行指向别处：{bad}")
        elif not hits:
            fails.append("cordis.patch.yml 里找不到任何 vault: 行（结构可能变了）")
        else:
            print(f"  ok    cordis.patch.yml 的 {len(hits)} 行都指向 {root_s}")

    p, cmd = probe_memset(memset)
    if p is None:
        fails.append(f"缺 {MEMSET}")
    elif cmd is None:
        fails.append("settings.json 无 writeHookCommand")
    elif str(judge) not in str(cmd):
        fails.append(f"writeHookCommand 未指向新判官：{cmd}")
    else:
        print(f"  ok    writeHookCommand → {cmd}")

    if new_root.is_dir():
        if not (new_root / "tools" / "zfind.py").exists():
            fails.append(f"{new_root} 下没有 tools/zfind.py")
        else:
            print(f"  ok    {new_root}/tools/zfind.py 在")
    else:
        fails.append(f"{new_root} 不是目录")

    if fails:
        print("\n核验 FAIL：")
        for f in fails:
            print("  ✗ " + f)
        return 1
    print("\n核验 PASS")
    return 0


# ─────────────────────────── 入口 ───────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description="Zero 库根路径的外部配置迁移（默认 dry-run）")
    ap.add_argument("--new-root", required=True,
                    help="改名／搬库后的新库绝对路径（如 /Users/oygb/Downloads/Zero）")
    ap.add_argument("--old-root",
                    help="旧库路径（用于 --fix-tree 的定位；省略则从 cordis 现状读）")
    ap.add_argument("--apply", action="store_true", help="真的写盘（默认只报告）")
    ap.add_argument("--verify", action="store_true", help="只读核验，不改任何东西")
    ap.add_argument("--only-fix-tree", action="store_true",
                    help="只执行 --fix-tree（跳过 ①② 外部配置）；搬库后修旧树时用")
    ap.add_argument("--fix-tree", metavar="DIR", default=None,
                    help="改写 DIR 下 .py 的硬编码 /…/Zero 前缀（默认改成 --new-root）")
    ap.add_argument("--fix-tree-to", default=None,
                    help="--fix-tree 的目标前缀；省略则等于 --new-root。"
                         "注意：搬库后旧树里的前缀应改成**旧树的新路径**")
    ap.add_argument("--cordis", default=str(CORDIS), help="（测试用）指定 cordis.patch.yml")
    ap.add_argument("--memset", default=str(MEMSET), help="（测试用）指定 settings.json")
    a = ap.parse_args()

    new_root = Path(a.new_root).expanduser().resolve()
    cordis = Path(a.cordis).expanduser()
    memset = Path(a.memset).expanduser()
    dry = not a.apply

    if a.verify:
        return verify(new_root, new_root / "tools" / "zero_mem_judge.py",
                      cordis=cordis, memset=memset)

    print(f"目标库根：{new_root}（{'dry-run' if dry else 'APPLY'}）\n")

    # 硬前置：写盘前必须确认目标库真的存在且是 Zero 库。
    # 没有这一条，一次误用 --apply 就会把插件指向不存在的目录（实测踩过）。
    if not dry and not a.only_fix_tree:
        if not new_root.is_dir():
            print(f"✗ 目标库不存在：{new_root}\n  拒绝写盘。先建好库（tools/zfind.py 要在），再回来。")
            return 2
        if not (new_root / "tools" / "zfind.py").exists():
            print(f"✗ {new_root} 下没有 tools/zfind.py，不像 Zero 库。拒绝写盘。")
            return 2
        if not (new_root / "tools" / "zero_mem_judge.py").exists():
            print(f"✗ {new_root} 下没有 tools/zero_mem_judge.py。拒绝写盘。")
            return 2

    _, hits = probe_cordis(cordis)
    old_root = a.old_root or (hits[0][1] if hits else "")
    print(f"旧库根：{old_root or '（未探到）'}\n")

    n1 = n2 = 0
    if a.only_fix_tree:
        print("① 插件配置：跳过（--only-fix-tree）")
        print("② 记忆判官：跳过（--only-fix-tree）")
    else:
        print(f"① 插件配置 {cordis.name}")
        n1 = apply_cordis(cordis, str(new_root), dry)
        print(f"\n② 记忆判官 {memset.name}")
        n2 = apply_memset(memset, new_root / "tools" / "zero_mem_judge.py", dry)

    n3 = 0
    if a.fix_tree:
        prefix = a.fix_tree_to or str(new_root)
        print(f"\n③ 库内 .py 的硬编码前缀（扫 {a.fix_tree} → 改为 {prefix}）")
        n3 = apply_hardcodes(Path(a.fix_tree).expanduser().resolve(), prefix, dry)

    print("\n汇总：" + ("无需改动" if n1 + n2 + n3 == 0 else f"待改 {n1 + n2 + n3} 处"))
    if dry:
        print("加 --apply 才会写盘；写盘前每个文件会备份成 <name>.bak-<时间戳>。")
    else:
        print("\n下一步：重启会话（插件在启动时读路径），再跑 "
              "`python3 tools/zero_vault_migrate.py --verify --new-root …`")
    return 0


if __name__ == "__main__":
    sys.exit(main())
