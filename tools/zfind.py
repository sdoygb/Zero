#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zfind.py — Zero 库的检索器 ＋ 否定结论门 ＋ 写作机检。

为什么存在：库有 143 篇正式稿、500+ 篇文档/暂存。靠 grep 扫一两篇就下
"库里没有／不能导出"的结论，已经出过错（Z102 层间耦合、D87 连续极限）。
本工具把"查得够不够"变成可检查的条件，而不是靠记性。

子命令
  preflight                 开工卡 ＋ 库况 ＋ 台账（每会话第一条命令）
  search "词1 词2" --family obj|syn|pos|hist
                            查库；打印「必读（标题命中）」「正文命中」「片段＋行号」
  prior "词"                库里关于该词的**既有判定**行（定理／判定／开放／具名输入）
  read Z117 Z120            登记已读（结论门据此判定）
  conclude --claim "…"      结论门：族不全／必读未读／已读<3 篇 ⟹ 拒绝（退出码 1）
  lint <篇>                 写作机检：断语／工作痕迹／悬空篇号／开放无背靠
  groups [--audit]          同义词族；--audit 列出库内 0 命中的死词
  selftest                  工具自检（改过本文件必跑）

查询族（family）的含义
  obj  = 对象名（要判定的东西）        syn = 同义／别名（换一套词再查）
  pos  = 正面形式（导出／定理／构造／已证）  hist = 既有判定（开放／判定／具名输入）
  否定性结论要求四族齐全；正面结论至少三族且必须含 obj＋syn。

退出码：0 = 通过；1 = 拒绝／有命中；2 = 用法错误。
"""
import argparse
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LEDGER = ROOT / "_staging" / "检索台账"
GROUPS_PATH = ROOT / "tools" / "zfind_groups.json"

# 不参与检索的路径片段（快照会把同一篇复制多份，污染召回）
EXCLUDE = (
    ".git", "tools", "__pycache__", ".obsidian", ".DS_Store",
    "zenodo_tools", "zhihu_publish", "github_tools",
    "_snapshot", "移出_", "归档_",
)
TIER_WEIGHT = {"篇": 1.25, "文档": 1.0, "暂存": 0.85, "旧稿": 0.6}

# 判定型标记词：只用来找"既有判定"，不算对象名，也不触发必读（否则「导出」命中 50 篇标题）
MARKER = {"导出", "定理", "构造", "已证", "判定", "开放", "具名输入", "存在",
          "不能", "不可能", "排除", "缺", "没有", "证明", "定义"}

NEG_RE = re.compile(r"(不能|无法|不可|不可能|不存在|没有|排除|排除掉|不行|未导出|未证|缺|做不到|不成立)")
JUDGE_RE = re.compile(r"(定理|判定|开放|具名输入|已证|导出|构造|不能|不可能|排除|不存在|未)")
ASCII_TERM = re.compile(r"^[A-Za-z][A-Za-z0-9_.^+\-]*$")


# ─────────────────────────── 语料 ───────────────────────────

def _excluded(part: str) -> bool:
    return part in EXCLUDE or part.startswith(("移出_", "归档_", "_snapshot", "_backup"))


def tier_of(rel: Path):
    parts = rel.parts
    if len(parts) == 1:
        return "篇" if re.match(r"Z\d{3}_", parts[0]) else None
    top = parts[0]
    if top == "_staging":
        return "暂存"
    if top == "paper":
        return "旧稿"
    if top == "QM":
        return "文档"
    if top == "ZeroGPU" and ((len(parts) > 1 and parts[1] == "doc") or parts[-1] == "README.md"):
        return "文档"
    if top == "ZeroLean" and parts[-1] == "README.md":
        return "文档"
    return None


def load_docs():
    docs = []
    for p in sorted(ROOT.rglob("*.md")):
        try:
            rel = p.relative_to(ROOT)
        except ValueError:
            continue
        if p.name == "README.md":
            continue
        if any(_excluded(part) for part in rel.parts):
            continue
        tier = tier_of(rel)
        if not tier:
            continue
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        low = text.lower()
        heads = [ln.strip() for ln in text.split("\n") if ln.lstrip().startswith("#")]
        fname = p.name[:-3]
        docs.append({
            "name": p.name,
            "key": fname,
            "rel": str(rel),
            "tier": tier,
            "text": text,
            "low": low,
            "heads": "\n".join(heads),
            "heads_low": "\n".join(heads).lower(),
            "fname_low": fname.lower(),
        })
    return docs

_LAZY = {"docs": None}


def docs():
    """惰性加载全部可检索文档；audit/lint 不需要它（省 ~1.5s）。"""
    if _LAZY["docs"] is None:
        _LAZY["docs"] = load_docs()
    return _LAZY["docs"]



def count(hay: str, term: str) -> int:
    """出现次数。纯 ASCII 词加字母边界（避免 GR 命中 graph）；CJK 直接数。"""
    t = term.lower()
    if ASCII_TERM.match(term):
        pat = r"(?<![A-Za-z])" + re.escape(t) + r"(?![A-Za-z])"
        return len(re.findall(pat, hay, re.IGNORECASE))
    return hay.count(t)


def split_terms(q: str):
    parts = re.split(r"[\s,，、;；|]+", q.strip())
    out = []
    for p in parts:
        p = p.strip(" \t\"'“”‘’()（）[]【】")
        if p:
            out.append(p)
    return out


def load_groups():
    try:
        data = json.loads(GROUPS_PATH.read_text(encoding="utf-8"))
        return [list(g) for g in data.get("groups", [])]
    except (OSError, json.JSONDecodeError):
        return []


def expand(terms, groups=None):
    groups = groups if groups is not None else load_groups()
    out = []
    for t in terms:
        for g in groups:
            if t in g or any(x.lower() == t.lower() for x in g):
                for x in g:
                    if x not in out:
                        out.append(x)
                break
        if t not in out:
            out.append(t)
    return out


# ─────────────────────────── 排名 ───────────────────────────

def rank(terms, expanded):
    n = len(docs())
    # 判定型标记词（导出/定理/判定/开放…）只用来找既有判定，不参与打分：
    # 否则 pos/hist 查询会被「导出」这种高频词带成"命中 143 篇"。
    scored = [t for t in expanded if t not in MARKER] or expanded
    idf, df = {}, {}
    for t in scored:
        df[t] = sum(1 for d in docs()
                    if count(d["fname_low"] + "\n" + d["heads_low"], t) or count(d["low"], t))
        idf[t] = math.log((n + 1) / (df[t] + 1)) + 1.0
    rows = []
    for d in docs():
        tb = sum(count(d["fname_low"] + "\n" + d["heads_low"], t) for t in scored)
        bb = {t: count(d["low"], t) for t in scored}
        total = sum(bb.values())
        if not total:
            continue
        score = 0.0
        for t in scored:
            c = bb[t]
            if c:
                score += idf[t] * (1.0 + math.log1p(c))
            ttl = count(d["fname_low"] + "\n" + d["heads_low"], t)
            if ttl:
                score += idf[t] * 5.0 * min(ttl, 3)
        score *= TIER_WEIGHT.get(d["tier"], 1.0)
        obj_terms = [t for t in terms if t not in MARKER] or terms
        primary_title = any(count(d["fname_low"] + "\n" + d["heads_low"], p) > 0 for p in obj_terms)
        rows.append({
            "doc": d["name"], "tier": d["tier"], "score": score,
            "primary_title": primary_title,
            "per": {t: bb[t] for t in scored if bb[t]},
            "tb": tb,
        })
    rows.sort(key=lambda r: -r["score"])
    return rows, df, idf


def doc_by_name(name):
    for d in docs():
        if d["name"] == name:
            return d
    return None


def snippet(d, terms, width=70):
    for t in terms:
        i = d["low"].find(t.lower())
        if i >= 0:
            line = d["text"].count("\n", 0, i) + 1
            s, e = max(0, i - width), min(len(d["text"]), i + width)
            frag = d["text"][s:e].replace("\n", " ")
            return f"L{line}: …{frag}…"
    return ""


# ─────────────────────────── 台账 ───────────────────────────

def L() -> Path:
    return Path(os.environ.get("ZFIND_LEDGER", str(DEFAULT_LEDGER)))


def window_path() -> Path:
    return L() / "current.json"


def load_window():
    p = window_path()
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def save_window(w):
    L().mkdir(parents=True, exist_ok=True)
    window_path().write_text(json.dumps(w, ensure_ascii=False, indent=1), encoding="utf-8")


def new_window():
    return {"id": time.strftime("%Y-%m-%dT%H:%M:%S"), "queries": [], "reads": [], "concluded": None}


def append_jsonl(name, obj):
    L().mkdir(parents=True, exist_ok=True)
    with open(L() / name, "a", encoding="utf-8") as f:
        f.write(json.dumps(obj, ensure_ascii=False) + "\n")


def resolve_doc(s):
    s = s.strip()
    if s.endswith(".md"):
        s = s[:-3]
    hits = [d for d in docs() if d["key"] == s]
    if hits:
        return hits[0]
    if re.fullmatch(r"Z\d{3}", s):          # 篇号：优先正式篇
        form = [d for d in docs() if d["tier"] == "篇" and d["key"].startswith(s + "_")]
        if len(form) == 1:
            return form[0]
        if len(form) > 1:
            raise SystemExit(f"「{s}」在正式篇里有 {len(form)} 篇同名号："
                             + ", ".join(d["key"][:30] for d in form))
    hits = [d for d in docs() if d["key"].startswith(s + "_")]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        cand = ", ".join(d["key"][:30] for d in hits[:6])
        raise SystemExit(f"「{s}」有歧义：{cand}")
    return None


# ─────────────────────────── 命令 ───────────────────────────

def cmd_search(a):
    terms = split_terms(a.query)
    if not terms:
        print("用法：zfind.py search \"词1 词2\" --family obj")
        return 2
    expanded = expand(terms)
    rows, df, _ = rank(terms, expanded)
    dead = [t for t in expanded if df.get(t, 0) == 0]
    must = [r for r in rows if r["primary_title"] and r["tier"] == "篇"]
    refs = [r for r in rows if r["primary_title"] and r["tier"] != "篇"]
    formal = [r for r in rows if r["tier"] == "篇"]

    w = load_window()
    if not w or w.get("concluded"):
        w = new_window()
    w["queries"].append({
        "family": a.family,
        "terms": terms,
        "expanded": expanded,
        "hits": [{"doc": r["doc"], "score": round(r["score"], 3),
                  "primary_title": r["primary_title"]} for r in rows[:400]],
    })
    save_window(w)

    hits = []
    for r in rows[:a.top]:
        d = doc_by_name(r["doc"])
        hits.append({
            "doc": r["doc"],
            "key": d["key"],
            "tier": r["tier"],
            "score": round(r["score"], 1),
            "per": " ".join(f"{t}×{c}" for t, c in sorted(r["per"].items(), key=lambda kv: -kv[1])[:4]),
            "snippet": snippet(d, terms),
        })
    data = {
        "ok": True,
        "family": a.family,
        "terms": terms,
        "expanded": expanded,
        "dead": dead,
        "total": len(rows),
        "formal": len(formal),
        "searchable": len(docs()),
        "must_read": [r["doc"] for r in must],
        "refs": [r["doc"] for r in refs],
        "hits": hits,
    }
    if getattr(a, "json", False):
        print(json.dumps(data, ensure_ascii=False))
        return 0

    print(f"━━ 查询族 {a.family}｜词：{' '.join(terms)}")
    print(f"   同义扩展：{'、'.join(x for x in expanded if x not in terms) or '（无）'}")
    print(f"   可搜 {len(docs())} 篇（正式 {sum(1 for d in docs() if d['tier'] == '篇')}）；本题命中 {len(rows)} 篇")
    if dead:
        print(f"   ⚠ 库内 0 命中：{'、'.join(dead)}（拼写？该换说法？）")
    if must:
        print(f"   ★ 必读（正式篇标题命中，{len(must)} 篇）：{'  '.join(r['doc'].split('_')[0] for r in must)}")
    else:
        print("   ★ 必读（正式篇标题命中）：无 —— 没有任何正式篇的标题含这个词，别据此说库里没有")
    if refs:
        print(f"   参考（暂存/文档标题命中，{len(refs)} 篇）："
              + "；".join(r["doc"][:36] for r in refs[:6]))
    for r in rows[:a.top]:
        d = doc_by_name(r["doc"])
        per = " ".join(f"{t}×{c}" for t, c in sorted(r["per"].items(), key=lambda kv: -kv[1])[:4])
        print(f"   [{r['tier']}] {d['key'][:44]:<44} score={r['score']:6.1f}  {per}")
        sn = snippet(d, terms)
        if sn:
            print(f"        {sn}")
    print(f"   正式篇命中 {len(formal)} 篇；正文命中明细见上（只列前 {a.top}，--top 可调）")
    return 0


def cmd_prior(a):
    terms = split_terms(a.query)
    expanded = expand(terms)
    out = []
    for d in docs():
        found = []
        heading = ""
        lines = d["text"].split("\n")
        for i, ln in enumerate(lines, 1):
            if ln.lstrip().startswith("#"):
                heading = ln.strip("# ").strip()
            low = ln.lower()
            if any(t.lower() in low for t in expanded) and JUDGE_RE.search(ln):
                found.append((i, heading, ln.strip()))
        if found:
            out.append((d, found))
    out.sort(key=lambda x: (TIER_WEIGHT.get(x[0]["tier"], 0), len(x[1])), reverse=True)

    items = []
    for d, found in out[:a.top]:
        for i, heading, ln in found[:3]:
            items.append({"doc": d["name"], "key": d["key"], "tier": d["tier"],
                          "heading": heading[:60], "line": i, "text": ln[:150]})
    data = {"ok": True, "terms": terms, "expanded": expanded,
            "docs_hit": len(out), "items": items}
    if getattr(a, "json", False):
        print(json.dumps(data, ensure_ascii=False))
        return 0

    print(f"━━ 既有判定｜词：{' '.join(terms)}（扩展 {len(expanded)} 词）｜命中 {len(out)} 篇")
    for d, found in out[:a.top]:
        print(f"\n   [{d['tier']}] {d['key']}")
        for i, heading, ln in found[:3]:
            frag = ln[:150]
            print(f"      L{i:<5d} {('§' + heading[:28] + ' ') if heading else ''}{frag}")
    if not out:
        print("   （无）—— 没有篇把该词与『定理/判定/开放/具名输入』写在同一行；这不等于库里没有，")
        print("   只说明该词的判定散落在正文里，需换同义词再查或直接读标题命中篇。")
    return 0


def cmd_read(a):
    if not a.docs:
        print("用法：zfind.py read Z117 Z120 …")
        return 2
    w = load_window() or new_window()
    added = []
    for s in a.docs:
        d = resolve_doc(s)
        if not d:
            print(f"   ✗ 找不到篇：{s}")
            continue
        if d["name"] not in w["reads"]:
            w["reads"].append(d["name"])
            added.append(d["key"])
    save_window(w)
    data = {"ok": True, "added": added, "reads": w["reads"], "count": len(w["reads"])}
    if getattr(a, "json", False):
        print(json.dumps(data, ensure_ascii=False))
        return 0
    print(f"   登记已读 {len(added)} 篇：{'、'.join(added) or '（无新增）'}；本窗口共 {len(w['reads'])} 篇")
    return 0


def gate_check(w, claim, reads=None):
    """返回 [(ok, 说明)]。规则见 AGENTS.md §1。"""
    reads = set(w.get("reads") or []) if reads is None else set(reads)
    qs = w.get("queries") or []
    fams = {q["family"] for q in qs}
    neg = bool(NEG_RE.search(claim))
    checks = []

    missing = {"obj", "syn"} - fams
    if neg:
        missing |= {"pos", "hist"} - fams
    ok = not missing and len(fams) >= 3
    checks.append((ok, f"查询族：已查 {sorted(fams) or '无'}；"
                       f"{'缺 ' + '、'.join(sorted(missing)) if missing else '齐全'}"
                       + (f"（否定性结论要 obj/syn/pos/hist 四族）" if neg else "")))

    # 只对 obj／syn 族要求「标题命中全读」——pos/hist 是找既有判定的，不强制读
    unread = []
    for q in qs:
        if q.get("family") not in ("obj", "syn"):
            continue
        for h in q.get("hits", []):
            if not h.get("primary_title") or h["doc"] in reads:
                continue
            d = doc_by_name(h["doc"])
            if d and d["tier"] == "篇":
                unread.append(h["doc"].split("_")[0])
    unread = sorted(set(unread))
    checks.append((not unread, f"标题命中必读：{'全部已读' if not unread else '未读 ' + '、'.join(unread)}"))

    # 全部命中里得分前 3 的篇必须读过
    allhits = {}
    for q in qs:
        if q.get("family") not in ("obj", "syn"):
            continue
        for h in q.get("hits", []):
            allhits[h["doc"]] = max(allhits.get(h["doc"], 0), h.get("score", 0))
    top3 = [k for k, _ in sorted(allhits.items(), key=lambda kv: -kv[1])[:3]]
    miss3 = [k.split("_")[0] for k in top3 if k not in reads]
    checks.append((not miss3, f"最高分 3 篇：{'已读' if not miss3 else '未读 ' + '、'.join(miss3)}"))

    formal = [r for r in reads if (doc_by_name(r) or {}).get("tier") == "篇"]
    checks.append((len(formal) >= 3, f"已读正式篇：{len(formal)} 篇（要求 ≥3）"))
    return checks


def short_name(name):
    d = doc_by_name(name)
    if not d:
        return name
    m = re.match(r"(Z\d{3})_", d["key"])
    return m.group(1) if m else d["key"][:24]


def render_card(w, claim):
    qs = w.get("queries") or []
    fams = {}
    for q in qs:
        fams.setdefault(q["family"], []).append(q["terms"])
    lines = [
        f"## 结论卡 · {claim}",
        f"- 判定性质：{'否定性' if NEG_RE.search(claim) else '正面'}",
        "- 查询族：" + "；".join(f"{k}「{' '.join(v[0])}」" for k, v in fams.items()),
        f"- 已读篇（{len(w.get('reads') or [])}）：" + "、".join(sorted(short_name(x) for x in w.get("reads") or [])),
        "- 必读标题命中：已全部读过（见 search 输出）",
        "- 正文允许写法：正面部分照写；未证的连线写「开放（背靠 Zxxx 定理 n）」，点出背靠哪条定理",
        f"- 台账编号：{w.get('id')}",
    ]
    return "\n".join(lines)


def cmd_conclude(a):
    w = load_window()
    if not w or not w.get("queries"):
        if getattr(a, "json", False):
            print(json.dumps({"ok": False, "claim": a.claim, "checks": [],
                              "card": "", "error": "没有检索记录；先 search 四族再 conclude"}, ensure_ascii=False))
            return 2
        print("✗ 没有检索记录。先 search（obj/syn/pos/hist 四族），再 conclude。")
        return 2
    checks = gate_check(w, a.claim)
    ok = all(c[0] for c in checks)
    card = render_card(w, a.claim) if ok else ""
    if getattr(a, "json", False):
        print(json.dumps({"ok": ok, "claim": a.claim,
                          "checks": [{"ok": c[0], "msg": c[1]} for c in checks],
                          "card": card, "error": "" if ok else "有 FAIL 项"}, ensure_ascii=False))
        if ok:
            append_jsonl("conclusions.jsonl", {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
                                               "claim": a.claim, "reads": w.get("reads"),
                                               "families": sorted({q["family"] for q in w["queries"]})})
            w["concluded"] = {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"), "claim": a.claim}
            save_window(w)
        return 0 if ok else 1

    print(f"━━ 结论门｜{a.claim}")
    for c_ok, msg in checks:
        print(f"   {'PASS' if c_ok else 'FAIL'}  {msg}")
    if not ok:
        print("\n✗ 拒绝出结论。把上面的 FAIL 逐条补上（换同义词查／读标题命中篇），再跑一次。")
        return 1
    append_jsonl("conclusions.jsonl", {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
                                       "claim": a.claim, "reads": w.get("reads"),
                                       "families": sorted({q["family"] for q in w["queries"]})})
    w["concluded"] = {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"), "claim": a.claim}
    save_window(w)
    print("\nPASS —— 结论卡（贴进工作笔记/正文依据）：\n")
    print(card)
    print("\n提醒：否定性结论还需 ≥2 个独立角度复算（代数/实跑/统计/Jacobian）。")
    return 0


def _style_check():
    p = ROOT / "_staging" / "style_check.py"
    spec = importlib.util.spec_from_file_location("style_check", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cmd_lint(a):
    p = Path(a.file)
    if not p.exists():
        if getattr(a, "json", False):
            print(json.dumps({"ok": False, "file": a.file, "issues": [],
                              "error": "文件不存在"}, ensure_ascii=False))
        else:
            print(f"✗ 文件不存在：{p}")
        return 2
    text = p.read_text(encoding="utf-8")
    issues = []

    try:
        sc = _style_check()
        for i, w, frag in sc.scan(p):
            issues.append((i, f"写作纪律[{w}]", frag))
    except Exception as e:  # noqa: BLE001
        issues.append((0, "style_check 无法运行", str(e)))

    existing = {int(q.name[1:4]) for q in ROOT.glob("Z[0-9][0-9][0-9]_*.md")}
    nmax = max(existing) if existing else 0
    archived = set(range(1, nmax + 1)) - existing
    for m in re.finditer(r"Z(\d{3})", text):
        n = int(m.group(1))
        if n in existing:
            continue
        line = text.count("\n", 0, m.start()) + 1
        tag = "缺号（已归档/撤号）" if n in archived else "不存在的篇号"
        issues.append((line, f"悬空引用 {m.group(0)}：{tag}", ""))

    lines = text.split("\n")
    for i, ln in enumerate(lines):
        if "开放" in ln and not re.search(r"Z\d{3}", ln):
            win = "\n".join(lines[max(0, i - 2):i + 3])
            if not re.search(r"Z\d{3}", win):
                issues.append((i + 1, "「开放」未点背靠定理", ln.strip()[:110]))

    issues.sort()
    data = {"ok": not issues, "file": p.name,
            "issues": [{"line": i, "tag": tag, "frag": frag} for i, tag, frag in issues]}
    if getattr(a, "json", False):
        print(json.dumps(data, ensure_ascii=False))
        return 1 if issues else 0
    print(f"━━ lint｜{p.name}｜{len(issues)} 处")
    for i, tag, frag in issues:
        print(f"   L{i:<5d} [{tag}]  {frag}")
    if not issues:
        print("   干净。")
    return 1 if issues else 0


def cmd_groups(a):
    groups = load_groups()
    print(f"同义词族 {len(groups)} 组：")
    for g in groups:
        print("   " + "、".join(g))
    if a.audit:
        dead = []
        for g in groups:
            for t in g:
                df = sum(1 for d in docs() if count(d["fname_low"] + "\n" + d["heads_low"], t) or count(d["low"], t))
                if df == 0:
                    dead.append((g[0], t))
        print(f"\n库内 0 命中的词：{len(dead)} 个")
        for head, t in dead:
            print(f"   [{head}组] {t}")
    return 0


def audit_data():
    """快照：篇数/缺号/可搜篇数/README 漂移/体例基线/检索窗口（不加载正文，~0.3s）。"""
    notes = sorted(ROOT.glob("Z[0-9][0-9][0-9]_*.md"))
    nums = sorted(int(q.name[1:4]) for q in notes)
    nmax = nums[-1] if nums else 0
    have = set(nums)
    missing = [n for n in range(1, nmax + 1) if n not in have]
    searchable = 0
    for q in ROOT.rglob("*.md"):
        if q.name == "README.md":
            continue
        rel = q.relative_to(ROOT)
        if any(_excluded(part) for part in rel.parts):
            continue
        if tier_of(rel):
            searchable += 1
    dirty, hits = [], 0
    try:
        sc = _style_check()
        for q in notes:
            hh = sc.scan(q)
            if hh:
                dirty.append(q.name.split("_")[0])
                hits += len(hh)
    except Exception:  # noqa: BLE001
        pass
    drift = []
    rd = ROOT / "README.md"
    if rd.exists():
        rt = rd.read_text(encoding="utf-8")
        m1 = re.search(r"正式稿 \*\*(\d+) 篇\*\*（Z(\d{3})–Z(\d{3})）", rt)
        m2 = re.search(r"下一篇编号 Z(\d{3})", rt)
        if m1 and (int(m1.group(1)) != len(notes) or int(m1.group(3)) != nmax):
            drift.append(f"README 写 {m1.group(1)} 篇 Z{m1.group(2)}–Z{m1.group(3)}，实测 {len(notes)} 篇 Z001–Z{nmax:03d}")
        if m2 and int(m2.group(1)) != nmax + 1:
            drift.append(f"README 下一篇写 Z{m2.group(1)}，应为 Z{nmax + 1:03d}")
    w = load_window()
    window = None
    if w and w.get("queries"):
        window = {"families": sorted({q["family"] for q in w["queries"]}),
                  "reads": len(w.get("reads") or []),
                  "concluded": bool(w.get("concluded"))}
    return {"ok": True, "formal": len(notes), "max": nmax, "next": nmax + 1,
            "missing": missing, "searchable": searchable,
            "drift": drift, "style": {"dirty": dirty, "hits": hits},
            "window": window}


def cmd_audit(a):
    data = audit_data()
    if getattr(a, "json", False):
        print(json.dumps(data, ensure_ascii=False))
        return 0
    print(f"正式 {data['formal']} 篇（Z001–Z{data['max']:03d}，缺 {len(data['missing'])} 个）；可搜 {data['searchable']} 篇")
    print("README 漂移：" + ("；".join(data["drift"]) if data["drift"] else "无"))
    print(f"体例机检：{len(data['style']['dirty'])} 篇 / {data['style']['hits']} 处")
    return 0


def cmd_preflight(a):
    print("═" * 72)
    print(" Zero 开工卡（每会话第一条命令；正文口径以 AGENTS.md 为准）")
    print("═" * 72)
    agents = ROOT / "AGENTS.md"
    if agents.exists():
        t = agents.read_text(encoding="utf-8")
        m = re.search(r"<!-- CARD:BEGIN -->(.*?)<!-- CARD:END -->", t, re.S)
        print((m.group(1).strip() if m else t)[:4000])
    else:
        print("⚠ 缺 AGENTS.md（开工卡不在库里）")

    print("\n" + "─" * 72)
    formal = [d for d in docs() if d["tier"] == "篇"]
    nums = sorted(int(d["key"][1:4]) for d in formal)
    nmax = nums[-1] if nums else 0
    missing = [n for n in range(1, nmax + 1) if n not in set(nums)]
    print(f"库况：正式 {len(formal)} 篇（Z{nums[0]:03d}–Z{nmax:03d}，缺 {len(missing)} 个）；"
          f"可检索 {len(docs())} 篇")

    readme = ROOT / "README.md"
    if readme.exists():
        rt = readme.read_text(encoding="utf-8")
        m1 = re.search(r"正式稿 \*\*(\d+) 篇\*\*（Z(\d{3})–Z(\d{3})）", rt)
        m2 = re.search(r"下一篇编号 Z(\d{3})", rt)
        drift = []
        if m1 and (int(m1.group(1)) != len(formal) or int(m1.group(3)) != nmax):
            drift.append(f"篇数/区间写 {m1.group(1)} 篇 Z{m1.group(2)}–Z{m1.group(3)}，实测 {len(formal)} 篇 Z001–Z{nmax:03d}")
        if m2 and int(m2.group(1)) != nmax + 1:
            drift.append(f"下一篇写 Z{m2.group(1)}，实测应为 Z{nmax + 1:03d}")
        print("README 漂移：" + ("；".join(drift) if drift else "无"))

    try:
        sc = _style_check()
        dirty, hits = [], 0
        for q in sorted(ROOT.glob("Z*.md")):
            hh = sc.scan(q)
            if hh:
                dirty.append(q.name.split("_")[0])
                hits += len(hh)
        print(f"体例机检：{len(dirty)} 篇 / {hits} 处命中（要求 0）"
              + (f"；{'、'.join(dirty[:6])}" if dirty else ""))
    except Exception as e:  # noqa: BLE001
        print(f"体例机检：无法运行（{e}）")

    w = load_window()
    print("\n" + "─" * 72)
    if w and w.get("queries"):
        fams = sorted({q["family"] for q in w["queries"]})
        print(f"当前检索窗口：{w['id']}｜族 {fams}｜已读 {len(w.get('reads') or [])} 篇"
              f"｜{'已出结论' if w.get('concluded') else '未出结论'}")
    else:
        print("当前检索窗口：无")
    concl = L() / "conclusions.jsonl"
    if concl.exists():
        n = sum(1 for _ in concl.open(encoding="utf-8"))
        print(f"累计过门结论：{n} 条（{concl}）")
    print("\n开工三步：① 本命令 ② 判定类任务先进检索门（search×4 族）③ 落盘前 lint")
    print("═" * 72)
    return 0


# ─────────────────────────── 自检 ───────────────────────────

def cmd_selftest(a):
    fails = []

    def chk(cond, msg):
        print(f"   {'ok  ' if cond else 'FAIL'}  {msg}")
        if not cond:
            fails.append(msg)

    formal = [d for d in docs() if d["tier"] == "篇"]
    print("1) 语料")
    chk(len(formal) >= 140, f"正式篇 {len(formal)}（≥140）")
    chk(len(docs()) >= 200, f"可检索 {len(docs())}（≥200）")

    def title_hits(q):
        terms = split_terms(q)
        rows, _, _ = rank(terms, expand(terms))
        return [r["doc"].split("_")[0] for r in rows if r["primary_title"]]

    print("2) 已知标题命中")
    h = title_hits("双缝")
    chk("Z117" in h, f"双缝 → 标题命中含 Z117（实际 {h[:8]}）")
    h = title_hits("精细结构")
    chk("Z164" in h, f"精细结构 → 含 Z164（实际 {h[:8]}）")
    h = title_hits("Koide")
    chk("Z019" in h or "Z082" in h, f"Koide → 含 Z019/Z082（实际 {h[:8]}）")
    h = title_hits("量纲锚")
    chk("Z153" in h, f"量纲锚 → 含 Z153（实际 {h[:8]}）")

    print("3) 同义扩展")
    e = expand(["Weinberg"])
    chk("弱混合角" in e, f"Weinberg → 扩展出弱混合角（{e}）")
    rows, _, _ = rank(["Weinberg"], e)
    chk(any("Z151" in r["doc"] for r in rows[:5]), "Weinberg → 前 5 命中含 Z151")

    print("4) prior 视图")
    out = []
    for d in docs():
        for i, ln in enumerate(d["text"].split("\n"), 1):
            if "禁闭" in ln and JUDGE_RE.search(ln):
                out.append(d["key"])
                break
    chk(len(out) >= 3, f"禁闭 的既有判定行覆盖 {len(out)} 篇（≥3）")

    print("5) 结论门逻辑")
    d117 = resolve_doc("Z117")
    d120 = resolve_doc("Z120")
    d014 = resolve_doc("Z014")
    chk(all([d117, d120, d014]), "测试篇存在")

    def mk(fams, reads, docs):
        return {"id": "test", "queries": [
            {"family": f, "terms": ["双缝"], "expanded": ["双缝"],
             "hits": [{"doc": d["name"], "score": 10.0, "primary_title": True} for d in docs]}
            for f in fams], "reads": [d["name"] for d in reads]}

    claim = "双缝条纹无法在 L2 上产生"
    w1 = mk(["obj"], [d117], [d117, d120, d014])
    c1 = gate_check(w1, claim)
    chk(not all(ok for ok, _ in c1), "只查 1 族 ⟹ 拒绝")
    w2 = mk(["obj", "syn", "pos", "hist"], [d117, d120, d014], [d117, d120, d014])
    c2 = gate_check(w2, claim)
    chk(all(ok for ok, _ in c2), "四族全 ＋ 必读全读 ⟹ 放行")
    w3 = mk(["obj", "syn", "pos", "hist"], [d117], [d117, d120, d014])
    c3 = gate_check(w3, claim)
    chk(not all(ok for ok, _ in c3), "必读未读 ⟹ 拒绝")
    card = render_card(w2, claim)
    chk("结论卡" in card and "开放" in card, "结论卡含标题与写法提醒")

    print("6) lint")
    with tempfile.TemporaryDirectory() as td:
        clean = Path(td) / "Z999_干净.md"
        clean.write_text("# 标题\n\n已证：$a=b$（定理 1）。开放（背靠 Z014 定理 22）。\n", encoding="utf-8")
        dirty = Path(td) / "Z998_脏.md"
        dirty.write_text("# 标题\n\n本轮把结论收窄了。这不可能导出。开放。\n", encoding="utf-8")

        def run_lint(p):
            r = subprocess.run([sys.executable, str(Path(__file__).resolve()), "lint", str(p)],
                               capture_output=True, text=True, cwd=str(ROOT))
            return r.returncode

        chk(run_lint(clean) == 0, "干净稿 → 0 命中")
        chk(run_lint(dirty) == 1, "脏稿 → 有命中")

    print("7) 死词审计")
    groups = load_groups()
    dead = []
    for g in groups:
        for t in g:
            df = sum(1 for d in docs() if count(d["low"], t))
            if df == 0:
                dead.append(t)
    chk(not dead, f"同义词族无死词（{len(dead)}）" + (f"：{dead}" if dead else ""))

    print(f"\n{'全部通过' if not fails else '失败 %d 项' % len(fails)}")
    return 0 if not fails else 1


# ─────────────────────────── 入口 ───────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Zero 库检索/结论门/写作机检")
    sub = ap.add_subparsers(dest="cmd")

    s = sub.add_parser("search")
    s.add_argument("query")
    s.add_argument("--family", default="obj", choices=["obj", "syn", "pos", "hist"])
    s.add_argument("--top", type=int, default=10)
    s.add_argument("--json", action="store_true")
    s.set_defaults(fn=cmd_search)

    s = sub.add_parser("prior")
    s.add_argument("query")
    s.add_argument("--top", type=int, default=10)
    s.add_argument("--json", action="store_true")
    s.set_defaults(fn=cmd_prior)

    s = sub.add_parser("read")
    s.add_argument("docs", nargs="*")
    s.add_argument("--json", action="store_true")
    s.set_defaults(fn=cmd_read)

    s = sub.add_parser("conclude")
    s.add_argument("--claim", required=True)
    s.add_argument("--json", action="store_true")
    s.set_defaults(fn=cmd_conclude)

    s = sub.add_parser("lint")
    s.add_argument("file")
    s.add_argument("--json", action="store_true")
    s.set_defaults(fn=cmd_lint)

    s = sub.add_parser("audit")
    s.add_argument("--json", action="store_true")
    s.set_defaults(fn=cmd_audit)

    s = sub.add_parser("groups")
    s.add_argument("--audit", action="store_true")
    s.set_defaults(fn=cmd_groups)

    s = sub.add_parser("preflight")
    s.set_defaults(fn=cmd_preflight)

    s = sub.add_parser("selftest")
    s.set_defaults(fn=cmd_selftest)

    a = ap.parse_args()
    if not getattr(a, "fn", None):
        ap.print_help()
        return 2
    return a.fn(a)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        traceback.print_exc()
        sys.exit(2)
