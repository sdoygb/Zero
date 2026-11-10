#!/usr/bin/env python3
"""
ZeroLean CI 守则（新库 Zb 命名口径）。

三条规则（任一违反即非零退出）：
  R1 禁 sorryAx —— 证明里不许有洞。
  R2 axiom 白名单 —— 除 Lean 基石公理与 scripts/axioms-allow.txt 里登记的 Zero 输入外，
     不允许出现别的公有理；出现即列出（那就是"具名输入"，必须显式登记）。
  R3 定理名对得上 Zb 篇号 —— 形如 `zbNNN_*` 的定理，NNN 必须对应库里真实存在的 Zb 篇；
     基础文件（Basic.lean / AxiomLedger.lean）里的定理允许不带篇号。

用法：
  python3 scripts/check.py                 # 跑 lake build 并审计
  python3 scripts/check.py --no-build      # 只审计上次构建产物
  python3 scripts/check.py --json          # 机器可读输出
  python3 scripts/check.py --coverage      # 额外打印 65 篇的模块覆盖情况（不改变退出码）
"""
import argparse, json, os, re, subprocess, sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
VAULT = ROOT.parent
BASE_MODULES = {"Basic", "AxiomLedger"}
LEAN_BIN = pathlib.Path.home() / ".local" / "lean4" / "bin"

def run_build():
    env = dict(os.environ)
    env["PATH"] = f"{LEAN_BIN}:{env.get('PATH','')}"
    p = subprocess.run(["lake", "build"], cwd=ROOT, env=env, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr

LEDGER = re.compile(r"LEDGER\|([^|]+)\|(.*)")

def run_audit():
    env = dict(os.environ)
    env["PATH"] = f"{LEAN_BIN}:{env.get('PATH','')}"
    p = subprocess.run(["lake", "env", "lean", "Zero/Audit.lean"], cwd=ROOT, env=env,
                       capture_output=True, text=True)
    return p.stdout + p.stderr

def parse_ledger(out):
    res = {}
    for m in LEDGER.finditer(out):
        name, ax = m.group(1).strip(), m.group(2).strip()
        res[name] = [a.strip() for a in ax.split(",") if a.strip()]
    return res

def load_allow():
    f = ROOT / "scripts" / "axioms-allow.txt"
    allow = set()
    if f.exists():
        for line in f.read_text(encoding="utf-8").splitlines():
            line = line.split("#")[0].strip()
            if line:
                allow.add(line)
    return allow

def vault_notes():
    notes = {}
    for p in VAULT.glob("Zb[0-9][0-9][0-9]_*.md"):
        m = re.match(r"Zb(\d{3})_", p.name)
        if m:
            notes[int(m.group(1))] = p.name
    return notes

def module_coverage(notes):
    """篇号 -> 该篇模块里的 theorem/lemma 名（按文件名 ZbNNN_*.lean）。"""
    cov = {n: [] for n in notes}
    for p in sorted((ROOT / "Zero").glob("Zb[0-9][0-9][0-9]_*.lean")):
        m = re.match(r"Zb(\d{3})_", p.name)
        if not m:
            continue
        n = int(m.group(1))
        txt = p.read_text(encoding="utf-8", errors="ignore")
        names = re.findall(r"^\s*(?:theorem|lemma)\s+(zb\d{3}_[A-Za-z0-9_']*)", txt, re.M)
        cov.setdefault(n, []).extend(names)
    return cov

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--coverage", action="store_true")
    a = ap.parse_args()

    out, rc = "", 0
    if not a.no_build:
        rc, out = run_build()
        if rc != 0:
            print("lake build 失败 —— 审计不可信（oleans 可能是旧的）：")
            print("\n".join(l for l in out.splitlines() if "error" in l)[:3000])
            return 1
    out = run_audit()
    allow, notes = load_allow(), vault_notes()
    ledger = parse_ledger(out)

    problems, rows = [], []
    for name, ax in sorted(ledger.items()):
        short = name.split(".")[-1]
        bad_sorry = [x for x in ax if x == "sorryAx"]
        bad_ax = [x for x in ax if x not in allow and x != "sorryAx"]
        m = re.match(r"zb(\d{3})_", short)
        r3 = None
        if m and int(m.group(1)) not in notes:
            r3 = f"定理 {short} 声称属于 Zb{m.group(1)}，但库里没有这一篇"
        if bad_sorry:
            problems.append(f"R1 {short}: 依赖 sorryAx（证明有洞）")
        if bad_ax:
            problems.append(f"R2 {short}: 未登记的公有理 {bad_ax}")
        if r3:
            problems.append(f"R3 {r3}")
        rows.append({"theorem": name, "axioms": ax,
                     "note": f"Zb{m.group(1)}" if m else "-",
                     "ok": not (bad_sorry or bad_ax or r3)})

    for p in sorted((ROOT / "Zero").glob("Zb[0-9][0-9][0-9]_*.lean")):
        mm = re.match(r"Zb(\d{3})_", p.name)
        if mm and int(mm.group(1)) not in notes:
            problems.append(f"R3 文件 {p.name} 声称属于 Zb{mm.group(1)}，但库里没有这一篇")

    cov = module_coverage(notes)
    covered = sum(1 for n in notes if cov.get(n))
    if a.coverage:
        print(f"覆盖：{covered}/{len(notes)} 篇有 Lean 模块；"
              f"{sum(len(v) for v in cov.values())} 条定理/引理")
        for n in sorted(notes):
            names = cov.get(n, [])
            print(f"  Zb{n:03d}  {'✓ ' + str(len(names)) if names else '—'}  {notes[n]}")

    if a.json:
        print(json.dumps({"build_ok": rc == 0, "theorems": rows, "problems": problems,
                          "allowed": sorted(allow), "notes": len(notes),
                          "covered_articles": covered,
                          "module_theorems": sum(len(v) for v in cov.values())},
                         ensure_ascii=False, indent=2))
    else:
        print(f"lake build: {'OK' if rc == 0 else f'FAILED (rc={rc})'}")
        print(f"审计定理 {len(rows)} 条 ｜ 白名单 {len(allow)} 个名字 ｜ 库内 Zb 篇 {len(notes)} 篇 ｜ 有模块 {covered} 篇")
        pad = max((len(r['theorem'].split('.')[-1]) for r in rows), default=10)
        for r in rows:
            print(f"  {'✓' if r['ok'] else '✗'} {r['theorem'].split('.')[-1]:<{pad}} [{r['note']}] ← {', '.join(r['axioms']) or '（无）'}")
        if problems:
            print("\n守则违反：")
            for p in problems:
                print("  ✗", p)
        else:
            print("\n守则全过：无 sorryAx、无未登记公理、定理名与 Zb 篇号一致。")
    return 1 if (problems or rc != 0) else 0

if __name__ == "__main__":
    sys.exit(main())
