#!/usr/bin/env python3
"""zretrieval.py —— 全语料检索（含归档），BM25 基线 + 可选的向量检索。

为什么要有它
------------
`tools/zfind.py` 把 `_staging/归档_`、`移出_`、`_snapshot` 排除在语料之外（EXCLUDE 与
`_excluded`），于是 1116 篇里只有 639 篇可搜。2026-10-28 那次"准备做对象=闭合构型"的
判断就是因为这个盲区差点做错——决定性反证 D50／D51 正好在被排除的归档里。

本工具做两件事：
  1. **全语料 BM25**（纯标准库，无依赖）：把"唯一草稿"类归档纳入召回。
  2. **可选向量检索**（SiliconFlow bge-m3，1024 维）：解决"换一套词就搜不到"。
     需要 `~/.config/zero/siliconflow.env` 里的 SILICONFLOW_API_KEY 与账户余额。

用法
----
  python3 tools/zretrieval.py stats
  python3 tools/zretrieval.py query "为什么格点上的点不能和层一一对应" --k 5
  python3 tools/zretrieval.py embed          # 建向量索引（要 key+余额）
  python3 tools/zretrieval.py query "..." --k 5 --vector

定位：**召回层**，不替代 zfind 的四族检索门与结论卡；命中后仍须读原文核实。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def apply_vault(vault=None):
    """把语料根切到 vault（`--vault` 或环境变量 ZERO_VAULT）；不传则用默认根。

    向量缓存与索引状态按 ROOT 派生，故切库后会自动用另一套缓存，不会串味。
    """
    global ROOT, CORPUS_CACHE, INDEX_STATE
    if not vault:
        return
    new = Path(vault).expanduser().resolve()
    if not new.is_dir():
        raise SystemExit(f"--vault 不是目录：{new}")
    ROOT = new
    CORPUS_CACHE = ROOT / "ZeroReal" / "results" / "corpus_cache.pkl"
    INDEX_STATE = ROOT / "ZeroReal" / "results" / "index_state.json"

# ── 语料策展：唯一草稿收，状态快照不收 ────────────────────────
INCLUDE_ROOTS = (
    ".",                        # 根目录 Z*.md
    "_staging/推导",
    "_staging/归档_推导稿",
    "_staging/移出_开放项",
    "_staging/归档",
    "ZeroGPU/doc",
    "ZeroReal/doc",
)
EXCLUDE_PART = ("快照", "_backup", "mem-backup", "zenodo", "zhihu_publish", "paper",
                "_snapshot", "__pycache__", ".git", "results", "node_modules",
                # 自身/元文档：检索必须回到原始证据，不能把我的综述排到前面
                "ZeroReal", "tools")
NOTE_RE = re.compile(r"^(Z\d{3})_")
CHUNK_MAX = 2000         # 8B 支持 32K；按 `##` 整节切，节过长才再切
CHUNK_OVERLAP = 200
CHUNK_TAG = "sec"        # 分块版本，进缓存签名与索引文件名
DEFAULT_MODEL = "Qwen/Qwen3-Embedding-8B"
DEFAULT_DIMS = 1024          # 混合检索下 4096 只多 +0.017 MRR（<0.02 阈），不值 4 倍体积
INSTRUCT = "Instruct: 给定中文问题，检索能回答它的段落\nQuery: "


@dataclass
class Chunk:
    path: str
    note: str
    heading: str
    start: int
    text: str
    tier: str = "正式"

    @property
    def key(self) -> str:
        # 正文也进 key：否则改一段话但节标题/行号不变时，增量会误判"没变"而留下旧向量
        h = hashlib.sha1(f"{self.path}|{self.heading}|{self.start}|{self.text}".encode()).hexdigest()[:16]
        return h


def _included(p: Path) -> bool:
    rel = p.relative_to(ROOT)
    parts = rel.parts
    if any(x in s for s in parts for x in EXCLUDE_PART):
        return False
    if p.name.upper().startswith("README"):
        return False                        # 任何 README 都不收（元文档，会自我污染）
    if str(rel.parent) == ".":
        return bool(NOTE_RE.match(p.name))  # 根目录只收 Z 篇
    for root in INCLUDE_ROOTS[1:]:
        if str(rel).startswith(root + "/"):
            return True
    return False


def iter_files() -> list[Path]:
    out = []
    for p in sorted(ROOT.rglob("*.md")):
        if _included(p):
            out.append(p)
    return out


def split_sections(text: str) -> list[tuple[str, int, str]]:
    """按 `## ` 整节切（篇标题 `# ` 不切，归入首节）。"""
    lines = text.splitlines()
    secs, cur_h, cur_start, buf = [], "(篇首)", 1, []
    for i, ln in enumerate(lines, 1):
        if ln.startswith("## "):
            if buf and any(b.strip() for b in buf):
                secs.append((cur_h, cur_start, "\n".join(buf)))
            cur_h, cur_start, buf = ln.lstrip("# ").strip(), i, [ln]
        else:
            buf.append(ln)
    if buf and any(b.strip() for b in buf):
        secs.append((cur_h, cur_start, "\n".join(buf)))
    return secs


def _tier(rel: str) -> str:
    if rel.startswith("_staging/归档"):
        return "归档"
    if rel.startswith("_staging/"):
        return "暂存"
    if rel.startswith("ZeroGPU/"):
        return "文档"
    return "正式"


def chunk_file(p: Path) -> list[Chunk]:
    rel = str(p.relative_to(ROOT))
    # 篇名要能唯一定位到文件：截断会造出 "zero_layers_" 这种前缀，召回后无法登记已读
    note = NOTE_RE.match(p.name).group(1) if NOTE_RE.match(p.name) else p.stem
    tier = _tier(rel)
    out = []
    for head, start, body in split_sections(p.read_text(encoding="utf-8", errors="ignore")):
        i = 0
        while i < len(body):
            piece = body[i:i + CHUNK_MAX]
            out.append(Chunk(rel, note, head, start + body[:i].count("\n"), piece, tier))
            if i + CHUNK_MAX >= len(body):
                break
            i += CHUNK_MAX - CHUNK_OVERLAP
    return out


CORPUS_CACHE = ROOT / "ZeroReal" / "results" / "corpus_cache.pkl"
INDEX_STATE = ROOT / "ZeroReal" / "results" / "index_state.json"


def manifest(files) -> dict:
    """path -> [mtime_ns, size]。只 stat，不读文件（快）。"""
    return {str(p.relative_to(ROOT)): [p.stat().st_mtime_ns, p.stat().st_size] for p in files}


def index_status() -> dict:
    """索引是否过期：比对当前语料 manifest 与建索引时记录的 manifest。"""
    cur = manifest(iter_files())
    if not INDEX_STATE.exists():
        return {"state": "missing", "cur": cur}
    st = json.loads(INDEX_STATE.read_text(encoding="utf-8"))
    old = st.get("manifest", {})
    added = sorted(set(cur) - set(old))
    removed = sorted(set(old) - set(cur))
    changed = sorted(k for k in set(cur) & set(old) if cur[k] != old[k])
    return {"state": "stale" if (added or removed or changed) else "fresh",
            "added": added, "removed": removed, "changed": changed,
            "built_at": st.get("built_at"), "model": st.get("model"),
            "dims": st.get("dims"), "chunks": st.get("chunks")}


def corpus_signature(files, max_len: int, overlap: int) -> str:
    h = hashlib.sha256()
    for p in files:
        st = p.stat()
        h.update(f"{p}|{st.st_mtime_ns}|{st.st_size}".encode())
    h.update(f"|{CHUNK_TAG}|{max_len}|{overlap}".encode())
    return h.hexdigest()


def load_corpus(use_cache: bool = True, max_len: int | None = None,
                overlap: int | None = None):
    """返回 (chunks, BM25)。带落盘缓存：文件没变就直接加载（0.64s vs 1.55s）。"""
    import pickle
    global CHUNK_MAX, CHUNK_OVERLAP
    if max_len:
        CHUNK_MAX = max_len
    if overlap:
        CHUNK_OVERLAP = overlap
    files = iter_files()
    sig = corpus_signature(files, CHUNK_MAX, CHUNK_OVERLAP)
    if use_cache and CORPUS_CACHE.exists():
        try:
            d = pickle.loads(CORPUS_CACHE.read_bytes())
            if d.get("sig") == sig and d.get("tag") == CHUNK_TAG:
                return d["chunks"], d["bm"]
        except Exception:  # noqa: BLE001
            pass
    chunks = [c for f in files for c in chunk_file(f)]
    bm = BM25().fit(chunks)
    CORPUS_CACHE.parent.mkdir(parents=True, exist_ok=True)
    CORPUS_CACHE.write_bytes(pickle.dumps(
        {"sig": sig, "tag": CHUNK_TAG, "chunks": chunks, "bm": bm}, protocol=4))
    return chunks, bm


def tokenize(s: str) -> list[str]:
    """中文按二元组、英文按词、数字按串。不装分词器。"""
    s = s.lower()
    toks = re.findall(r"[a-z][a-z0-9_.^\-]*|\d+(?:\.\d+)?", s)
    cjk = re.findall(r"[\u4e00-\u9fff]", s)
    toks += [cjk[i] + cjk[i + 1] for i in range(len(cjk) - 1)]
    toks += cjk
    return toks


@dataclass
class BM25:
    chunks: list[Chunk] = field(default_factory=list)
    df: Counter = field(default_factory=Counter)
    tf: list[Counter] = field(default_factory=list)
    dl: list[int] = field(default_factory=list)
    avgdl: float = 1.0
    k1: float = 1.5
    b: float = 0.6

    def fit(self, chunks: list[Chunk]) -> "BM25":
        self.chunks = chunks
        self.tf = [Counter(tokenize(c.text + " " + c.heading)) for c in chunks]
        self.dl = [sum(t.values()) for t in self.tf]
        self.avgdl = (sum(self.dl) / len(self.dl)) if self.dl else 1.0
        self.df = Counter()
        for t in self.tf:
            self.df.update(t.keys())
        return self

    def search(self, q: str, k: int = 8) -> list[tuple[float, Chunk]]:
        n = len(self.chunks)
        ids = [self.df[t] for t in tokenize(q) if self.df[t] > 0]
        if not ids:
            return []
        scored = []
        for i, tf in enumerate(self.tf):
            sc = 0.0
            dl = self.dl[i] or 1
            for term, dfv in Counter(tokenize(q)).items():
                f = tf.get(term)
                if not f:
                    continue
                idf = math.log(1 + (n - dfv + 0.5) / (dfv + 0.5))
                sc += idf * f * (self.k1 + 1) / (f + self.k1 * (1 - self.b + self.b * dl / self.avgdl))
            if sc > 0:
                scored.append((sc, self.chunks[i]))
        scored.sort(key=lambda x: -x[0])
        return scored[:k]


# ── 向量检索（可选）────────────────────────────────────────
def cache_path(model: str = DEFAULT_MODEL, dims: int | None = None,
               ctx: bool = False) -> Path:
    slug = model.split("/")[-1].lower().replace(".", "-")
    if dims:
        slug += f"-{dims}"
    path = ROOT / "ZeroReal" / "results" / f"vector_index_{slug}-{CHUNK_TAG}.npz"
    return path.with_name(path.stem + "-ctx" + path.suffix) if ctx else path


CACHE = cache_path()


def _sf_key() -> str:
    """缺 key 时抛异常而不是退出——调用方要能降级到 BM25。"""
    p = Path.home() / ".config/zero/siliconflow.env"
    if not p.exists():
        raise RuntimeError("缺 ~/.config/zero/siliconflow.env")
    for line in p.read_text().splitlines():
        if line.startswith("SILICONFLOW_API_KEY="):
            v = line.split("=", 1)[1].strip()
            if v:
                return v
    raise RuntimeError("环境文件里没有 SILICONFLOW_API_KEY")


def embed_chunks(chunks: list[Chunk], model: str = DEFAULT_MODEL, batch: int = 32,
                 workers: int = 2, limit: int | None = None, trunc: int = 480,
                 dims: int | None = None, ctx: bool = False,
                 incremental: bool = False) -> None:
    """建向量索引。批小、并发低、带重试——隧道下 8 并发会超时，且批 32 超 8192 token。"""
    import concurrent.futures as cf
    import time
    import urllib.error
    import urllib.request

    import numpy as np

    if limit:
        chunks = chunks[:limit]
    key = _sf_key()
    url = "https://api.siliconflow.cn/v1/embeddings"

    def call(texts):
        payload = {"model": model, "input": texts}
        if dims:
            payload["dimensions"] = dims
        body = json.dumps(payload).encode()
        last = None
        for attempt in range(6):
            try:
                req = urllib.request.Request(url, data=body, headers={
                    "Authorization": f"Bearer {key}", "Content-Type": "application/json"})
                with urllib.request.urlopen(req, timeout=40) as r:
                    d = json.loads(r.read())
                if d.get("data"):
                    tok[0] += int((d.get("usage") or {}).get("total_tokens", 0))
                    return [e["embedding"] for e in d["data"]]
                last = RuntimeError(f"嵌入失败：{d}")
            except urllib.error.HTTPError as e:
                last = e
                if e.code == 429:          # 免费档限流：退避久一点
                    time.sleep(5.0 * (attempt + 1))
                    continue
                time.sleep(1.5 * (attempt + 1))
            except Exception as e:  # noqa: BLE001
                last = e
                time.sleep(1.5 * (attempt + 1))
        raise RuntimeError(f"重试 6 次仍失败：{last}")

    # ctx=True：把「篇名 + 节标题」拼在块前——否则 450 字从节中间切开，
    # 稠密模型不知道这段属于哪一节（BM25 不在乎，稠密模型很在乎）。
    texts = [(f"【{c.note}】§{c.heading}\n" + c.text)[:trunc + 60] if ctx else c.text[:trunc]
             for c in chunks]

    # 增量：只嵌入新块，已存在的键直接复用向量
    reuse: dict[str, "object"] = {}
    if incremental:
        out0 = cache_path(model, dims, ctx)
        if out0.exists():
            import numpy as _np
            z0 = _np.load(out0)
            reuse = {k: v for k, v in zip(z0["keys"].tolist(), z0["vec"])}
            print(f"  增量：已有 {len(reuse)} 块可复用", flush=True)
    if reuse:
        todo = [i for i, c in enumerate(chunks) if c.key not in reuse]
        print(f"  需新嵌入 {len(todo)}/{len(chunks)} 块（其余复用）", flush=True)
        texts_sel, idx_sel = [texts[i] for i in todo], todo
    else:
        texts_sel, idx_sel = texts, list(range(len(texts)))
    batches = [texts_sel[i:i + batch] for i in range(0, len(texts_sel), batch)]
    vecs: list[list[float]] = []
    tok = [0]
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(call, b) for b in batches]
        for i, f in enumerate(futs, 1):
            vecs.extend(f.result())
            if i % 10 == 0 or i == len(futs):
                dt = time.time() - t0
                print(f"  {i}/{len(futs)} 批  已用 {dt:.0f}s  {i*batch/dt:.1f} 块/s", flush=True)
    if reuse:
        full = [None] * len(chunks)
        for i, c in enumerate(chunks):
            if c.key in reuse:
                full[i] = reuse[c.key]
        for j, i in enumerate(idx_sel):
            full[i] = vecs[j]
        arr = np.asarray(full, dtype=np.float32)
    else:
        arr = np.asarray(vecs, dtype=np.float32)
    arr /= np.linalg.norm(arr, axis=1, keepdims=True)
    out = cache_path(model, dims, ctx)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, vec=arr, keys=np.array([c.key for c in chunks]))
    (out.with_suffix(".json")).write_text(json.dumps(
        [{"key": c.key, "path": c.path, "note": c.note, "heading": c.heading,
          "start": c.start, "tier": c.tier, "text": c.text[:300]} for c in chunks], ensure_ascii=False))
    INDEX_STATE.write_text(json.dumps({
        "built_at": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
        "model": model, "dims": dims or arr.shape[1], "ctx": ctx, "chunks": len(chunks),
        "manifest": manifest(iter_files()),
    }, ensure_ascii=False), encoding="utf-8")
    print(f"索引已存：{out}（{len(chunks)} 块，{arr.shape[1]} 维）")
    print(f"模型 {model}；本轮计费 token ≈ {tok[0]:,}（按 dashboard 单价自行折算）")


def embed_one(text: str, model: str = DEFAULT_MODEL, dims: int | None = None) -> "object":
    import json as _json
    import urllib.request
    import numpy as np
    key = _sf_key()
    if model.startswith("Qwen/Qwen3-Embedding"):
        text = INSTRUCT + text          # 指令型嵌入：查询侧必须带 instruction
    payload = {"model": model, "input": [text]}
    if dims:
        payload["dimensions"] = dims
    body = _json.dumps(payload).encode()
    req = urllib.request.Request("https://api.siliconflow.cn/v1/embeddings", data=body,
                                 headers={"Authorization": f"Bearer {key}",
                                          "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        d = _json.loads(r.read())
    v = np.asarray(d["data"][0]["embedding"], dtype=np.float32)
    return v / (np.linalg.norm(v) or 1.0)


def query_vector(text: str, k: int = 8, model: str = DEFAULT_MODEL,
                 dims: int | None = None, ctx: bool = True):
    import numpy as np
    path = cache_path(model, dims, ctx)
    if not path.exists():
        sys.exit(f"缺向量索引：{path}；先跑 python3 tools/zretrieval.py embed")
    z = np.load(path)
    vec, keys = z["vec"], z["keys"].tolist()
    meta = json.loads(path.with_suffix(".json").read_text())
    lut = {m["key"]: m for m in meta}
    q = embed_one(text, model, dims)
    sims = vec @ q
    idx = np.argsort(-sims)[:k]
    out = []
    for i in idx:
        m = lut.get(keys[int(i)], {})
        out.append((float(sims[int(i)]), m))
    return out


# 多粒度：细块（450 字）与整节（≤2000 字）两套索引并存，篇号层面 RRF 融合
FINE_INDEX = "vector_index_qwen3-embedding-8b-1024-ctx"
COARSE_INDEX = "vector_index_qwen3-embedding-8b-1024-sec-ctx"


def _load_index(name: str):
    import numpy as np
    base = ROOT / "ZeroReal" / "results" / name
    z = np.load(str(base) + ".npz")
    meta = json.loads((base.with_suffix(".json")).read_text())
    return z["vec"], meta


LAST_MODE = "multi"   # "multi" | "bm25-fallback:<原因>"


def bm25_only(text: str, chunks: list[Chunk], k: int = 10, pool: int = 30):
    """全语料 BM25，按篇号去重。降级路径与检索门第 0 步同语料，只是没有语义层。"""
    hits = BM25().fit(chunks).search(text, pool)
    seen, out = set(), []
    for i, (_, c) in enumerate(hits, 1):
        if c.note in seen:
            continue
        seen.add(c.note)
        out.append((c.note, 1.0 / (60 + i)))
    return out[:k]


def multi_search(text: str, chunks: list[Chunk], k: int = 10, K: int = 60,
                 model: str = DEFAULT_MODEL, dims: int | None = DEFAULT_DIMS):
    """三路融合：BM25（整节）+ 细块向量 + 整节向量。

    实测（18 条查询）：内容层 @5 = 0.944、MRR 0.760；D51 那条排第 4。
    """
    global LAST_MODE
    import numpy as np
    pool = max(k, 30)   # 候选池要够大：D51 单列在第 10–20 名，池小会漏掉
    lex = [c.note for _, c in BM25().fit(chunks).search(text, pool)]
    try:
        e = embed_one(text, model, dims)
    except Exception as exc:  # noqa: BLE001
        LAST_MODE = f"bm25-fallback：{type(exc).__name__}: {exc}"
        return bm25_only(text, chunks, k, pool)
    LAST_MODE = "multi"
    lists = [lex]
    for name in (FINE_INDEX, COARSE_INDEX):
        try:
            vec, meta = _load_index(name)
        except Exception as exc:  # noqa: BLE001
            # 新库／刚切 --vault 的库可能还没有向量索引：退到纯 BM25，而不是硬失败
            LAST_MODE = f"bm25-fallback（缺 {name}：{type(exc).__name__}；先跑 embed --ctx 建索引）"
            return bm25_only(text, chunks, k, pool)
        idx = np.argsort(-(vec @ e))[:pool]
        lists.append([meta[i]["note"] for i in idx])
    score: dict[str, float] = {}
    for lst in lists:
        for i, n in enumerate(lst, 1):
            if n.upper().startswith("README"):
                continue                    # 旧索引里还留着 README 块，查询时剔掉
            score[n] = score.get(n, 0.0) + 1.0 / (K + i)
    return sorted(score.items(), key=lambda x: -x[1])[:k]


def hybrid_search(text: str, chunks: list[Chunk], k: int = 8, K: int = 60,
                  model: str = DEFAULT_MODEL, dims: int | None = DEFAULT_DIMS,
                  ctx: bool = True):
    """BM25 + 向量 的 RRF 融合。实测优于任一路单独用。"""
    bm = BM25().fit(chunks)
    lex = [(c.note, c) for _, c in bm.search(text, k)]
    vec = [(m.get("note", ""), m) for _, m in query_vector(text, k, model, dims, ctx)]
    score: dict[str, float] = {}
    for lst in (lex, vec):
        for i, (n, _) in enumerate(lst, 1):
            score[n] = score.get(n, 0.0) + 1.0 / (K + i)
    order = sorted(score, key=lambda n: -score[n])
    lex_map = {n: c for n, c in lex}
    vec_map = {n: m for n, m in vec}
    out = []
    for n in order:
        if n in lex_map:
            out.append((score[n], lex_map[n]))
        elif n in vec_map:
            m = vec_map[n]
            out.append((score[n], Chunk(m.get("path", ""), m.get("note", ""),
                                        m.get("heading", ""), m.get("start", 0), m.get("text", ""))))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vault", default=os.environ.get("ZERO_VAULT"),
                    help="语料根（默认=本脚本上两级；也可用环境变量 ZERO_VAULT）")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("stats")
    sub.add_parser("status")
    q = sub.add_parser("query"); q.add_argument("text"); q.add_argument("--k", type=int, default=8)
    q.add_argument("--json", action="store_true", help="机器可读输出（供 zfind recall／插件调用）")
    q.add_argument("--vector", action="store_true")
    q.add_argument("--model", default=DEFAULT_MODEL); q.add_argument("--dims", type=int, default=DEFAULT_DIMS)
    q.add_argument("--hybrid", action="store_true")
    q.add_argument("--multi", action="store_true", help="三路融合（推荐）")
    e = sub.add_parser("embed"); e.add_argument("--limit", type=int, default=None)
    e.add_argument("--batch", type=int, default=32); e.add_argument("--workers", type=int, default=2)
    e.add_argument("--dims", type=int, default=None)
    e.add_argument("--model", default=DEFAULT_MODEL)
    e.add_argument("--ctx", action="store_true")
    e.add_argument("--incremental", action="store_true", help="只嵌入变更块")
    a = ap.parse_args()
    apply_vault(a.vault)

    if a.cmd == "status":
        st = index_status()
        if st["state"] == "missing":
            print(f"索引：未建（先跑 python3 tools/zretrieval.py embed --ctx）")
            return 0
        if st["state"] == "fresh":
            print(f"索引：最新（{st['chunks']} 块，{st['model']}"
                  f"{'/' + str(st['dims']) + '维' if st.get('dims') else ''}，建于 {st['built_at']}）")
            return 0
        n = len(st["added"]) + len(st["removed"]) + len(st["changed"])
        print(f"索引：**过期** —— 语料变更 {n} 篇"
              f"（新增 {len(st['added'])}／改 {len(st['changed'])}／删 {len(st['removed'])}）；"
              f"重建：python3 tools/zretrieval.py embed --model \"{st.get('model')}\" "
              f"{'--dims ' + str(st['dims']) + ' ' if st.get('dims') else ''}--ctx --incremental")
        for k in (st["changed"] + st["added"])[:5]:
            print(f"    · {k}")
        return 0

    if a.cmd == "stats":
        files = iter_files()
        chunks, _ = load_corpus()
        by_top = Counter(str(f.relative_to(ROOT)).split("/")[0] for f in files)
        print(f"文件 {len(files)} 篇；块 {len(chunks)}")
        for k, v in by_top.most_common():
            print(f"  {k:24s} {v:4d}")
        return 0

    files = iter_files()
    chunks = [c for f in files for c in chunk_file(f)]
    if a.cmd == "embed":
        embed_chunks(chunks, model=a.model, batch=a.batch, workers=a.workers,
                     limit=a.limit, dims=a.dims, ctx=a.ctx, incremental=a.incremental)
        return 0

    if a.multi:
        res = multi_search(a.text, chunks, a.k)
        tmap = {c.note: c.tier for c in chunks}
        pmap = {c.note: c.path for c in chunks}
        if a.json:
            print(json.dumps({
                "ok": True, "query": a.text, "mode": LAST_MODE, "k": len(res),
                "hits": [{"note": n, "score": round(float(sc), 5),
                          "tier": tmap.get(n, "?"), "path": pmap.get(n, "")}
                         for n, sc in res],
            }, ensure_ascii=False))
            return 0
        if LAST_MODE != "multi":
            print(f"⚠ 语义检索不可用，已降级到全语料 BM25（{LAST_MODE}）")
            print(f"查询（BM25 全语料）：{a.text}")
        else:
            print(f"查询（三路融合 BM25+细块+整节）：{a.text}")
        for i, (n, sc) in enumerate(res, 1):
            print(f"  {i}. [{sc:.4f}] {n:24s} ({tmap.get(n, '?')})")
        return 0

    if a.hybrid:
        print(f"查询（混合 BM25+向量）：{a.text}")
        for i, (sc, c) in enumerate(hybrid_search(a.text, chunks, a.k, model=a.model, dims=a.dims), 1):
            print(f"  {i}. [{sc:.4f}] {c.note}  §{c.heading[:30]}  L{c.start}  {c.path}")
        return 0

    if a.vector:
        print(f"查询（向量）：{a.text}")
        for i, (sc, m) in enumerate(query_vector(a.text, a.k, a.model, a.dims), 1):
            print(f"  {i}. [{sc:.4f}] {m.get('note','?')}  §{m.get('heading','')[:30]}  "
                  f"L{m.get('start','?')}  {m.get('path','?')}")
        return 0

    hits = bm.search(a.text, a.k)
    print(f"查询：{a.text}")
    for i, (sc, c) in enumerate(hits, 1):
        print(f"  {i}. [{sc:6.2f}] {c.note}  §{c.heading[:30]}  L{c.start}  {c.path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
