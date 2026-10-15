#!/usr/bin/env python3
"""Batch convert Markdown articles (Chinese, LaTeX math) to PDF.

md -> (preprocess math shorthands) -> pandoc -> xelatex -> PDF.
Output: <src_dir>/pdf_export/<article>.pdf ; logs in pdf_export/_logs/
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

TEXBIN = "/Library/TeX/texbin"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LUA_FILTER = os.path.join(SCRIPT_DIR, "strip_heading_math.lua")
PREPROC = os.path.join(SCRIPT_DIR, "preprocess_math.py")
CJK_FONT = os.environ.get("GEO_CJK_FONT", "Heiti SC")


EXTRA_CODES = [0x02B0, 0x02B3, 0x02E1, 0x0302, 0x1E21, 0x2032, 0x2033, 0x207A, 0x207D, 0x207E, 0x207F, 0x208A, 0x208B, 0x2102, 0x210B, 0x210D, 0x210F, 0x2111, 0x2112, 0x2113, 0x2115, 0x2118, 0x211D, 0x2124, 0x212C, 0x212D, 0x2130, 0x2153, 0x2154, 0x2197, 0x21A6, 0x21D4, 0x2200, 0x2202, 0x2203, 0x2205, 0x2207, 0x2208, 0x2209, 0x220F, 0x2211, 0x2213, 0x2216, 0x2218, 0x221D, 0x2223, 0x2225, 0x2227, 0x2228, 0x2229, 0x222A, 0x222B, 0x222E, 0x223C, 0x2245, 0x2254, 0x2260, 0x2261, 0x2264, 0x2265, 0x226A, 0x226B, 0x2272, 0x2273, 0x2282, 0x2286, 0x228A, 0x228B, 0x2295, 0x2297, 0x22A5, 0x22C9, 0x22EF, 0x2308, 0x2309, 0x230A, 0x230B, 0x2460, 0x2461, 0x2462, 0x2463, 0x2464, 0x2465, 0x2466, 0x2467, 0x2468, 0x2469, 0x246A, 0x246B, 0x246C, 0x246D, 0x246E, 0x246F, 0x2470, 0x2471, 0x2472, 0x2473, 0x250C, 0x2510, 0x2514, 0x2518, 0x251C, 0x2524, 0x252C, 0x2534, 0x2550, 0x2551, 0x255A, 0x255D, 0x25A1, 0x25BA, 0x25BC, 0x25C4, 0x25CF, 0x2605, 0x26A0, 0x27C2, 0x27F7, 0x27F8, 0x27F9, 0x27FA, 0x5185, 0x5B9A, 0x5B9E, 0x5E38, 0x6570, 0x68D2, 0x6D4B, 0x73AF, 0x771F, 0x89C2, 0x951A, 0xFF0C, 0xFF1A, 0x1D42E, 0x1D49C, 0x1D4A2, 0x1D4A6, 0x1D4A9, 0x1D4AC, 0x1D4AE, 0x1D4B0, 0x1D4B1, 0x1D4B2, 0x1D4B3, 0x1D4B5, 0x1D530, 0x1D532, 0x1D6FD, 0x1D705, 0x2500, 0x2501, 0x2502, 0x2503, 0x2504, 0x2505, 0x2506, 0x2507, 0x2508, 0x2509, 0x250A, 0x250B, 0x250C, 0x250D, 0x250E, 0x250F, 0x2510, 0x2511, 0x2512, 0x2513, 0x2514, 0x2515, 0x2516, 0x2517, 0x2518, 0x2519, 0x251A, 0x251B, 0x251C, 0x251D, 0x251E, 0x251F, 0x2520, 0x2521, 0x2522, 0x2523, 0x2524, 0x2525, 0x2526, 0x2527, 0x2528, 0x2529, 0x252A, 0x252B, 0x252C, 0x252D, 0x252E, 0x252F, 0x2530, 0x2531, 0x2532, 0x2533, 0x2534, 0x2535, 0x2536, 0x2537, 0x2538, 0x2539, 0x253A, 0x253B, 0x253C, 0x253D, 0x253E, 0x253F, 0x2540, 0x2541, 0x2542, 0x2543, 0x2544, 0x2545, 0x2546, 0x2547, 0x2548, 0x2549, 0x254A, 0x254B, 0x254C, 0x254D, 0x254E, 0x254F, 0x2550, 0x2551, 0x2552, 0x2553, 0x2554, 0x2555, 0x2556, 0x2557, 0x2558, 0x2559, 0x255A, 0x255B, 0x255C, 0x255D, 0x255E, 0x255F, 0x2560, 0x2561, 0x2562, 0x2563, 0x2564, 0x2565, 0x2566, 0x2567, 0x2568, 0x2569, 0x256A, 0x256B, 0x256C, 0x256D, 0x256E, 0x256F, 0x2570, 0x2571, 0x2572, 0x2573, 0x2574, 0x2575, 0x2576, 0x2577, 0x2578, 0x2579, 0x257A, 0x257B, 0x257C, 0x257D, 0x257E, 0x257F]
DECLS = ", ".join(
    f'"{c:04X}'
    for c in sorted(set(
        list(range(0x0391, 0x03A9 + 1))
        + list(range(0x03B1, 0x03C9 + 1))
        + [0x00A7, 0x00B0, 0x00B2, 0x00B3, 0x00B9, 0x00D7, 0x2013, 0x207B,
           0x2080, 0x2081, 0x2082, 0x2083, 0x2084, 0x2085, 0x2086, 0x2087,
           0x2070, 0x2074, 0x2075, 0x2076, 0x2077, 0x2078, 0x2079, 0x2088, 0x2089, 0x2110, 0x2133, 0x2190, 0x2191, 0x2192, 0x2193,
           0x21D2, 0x2194, 0x220E, 0x2212, 0x221A, 0x2248, 0x2500, 0x2501, 0x2502,
           0x253C, 0x254B, 0x3001, 0x3002, 0xFF0C, 0xFF0E, 0x2705, 0x2713, 0x27E8, 0x27E9, 0x1D49E] + EXTRA_CODES
    ))
)

HEADER_TEX = r"""
% --- extra packages ---
\usepackage{amsmath, amssymb, amsthm}
\usepackage{mathtools}
\usepackage[all]{xy}
\usepackage{fontspec}
\newfontfamily\symfallback{STIX Two Math}
\setCJKfallbackfamilyfont{\CJKrmdefault}{STIX Two Math}
\AtBeginDocument{%
\makeatletter
\xeCJKsetup{fallback=true}
\xeCJKDeclareCharClass{CJK}{__DECLS__}
\makeatother
}
""".replace("__DECLS__", DECLS)


def find_pandoc():
    for cand in ("/usr/local/bin/pandoc", "pandoc"):
        if os.path.exists(cand):
            return cand
    return "pandoc"


def convert_one(md_path, out_pdf, log_path, header, work_dir):
    env = dict(os.environ)
    env["PATH"] = f"{TEXBIN}:" + env.get("PATH", "")
    base = os.path.basename(md_path)
    tmp_md = os.path.join(work_dir, ".pp_" + base)
    subprocess.run([sys.executable, PREPROC, md_path, tmp_md],
                   capture_output=True, check=False)
    cmd = [
        find_pandoc(), tmp_md, "-o", out_pdf,
        "--pdf-engine=xelatex",
        "--from=markdown-yaml_metadata_block+tex_math_single_backslash+tex_math_double_backslash",
        "--lua-filter", LUA_FILTER,
        "-H", header,
        "-V", f"CJKmainfont={CJK_FONT}",
        "-V", f"CJKsansfont={CJK_FONT}",
        "-V", "CJKmonofont=Menlo",
        "-V", "geometry:margin=2cm",
        "-V", "colorlinks=true",
        "-V", "urlcolor=blue",
        "-V", "linkcolor=blue",
        "-V", "toc-depth=2",
        "--standalone",
    ]
    start = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900, env=env)
        ok = proc.returncode == 0 and os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 1000
        with open(log_path, "w", encoding="utf-8", errors="replace") as f:
            f.write(f"cmd: pandoc {os.path.basename(tmp_md)} -> {os.path.basename(out_pdf)}\n")
            f.write(f"returncode: {proc.returncode}\n")
            f.write("--- stdout tail ---\n")
            f.write(proc.stdout[-4000:])
            f.write("\n--- stderr tail ---\n")
            f.write(proc.stderr[-8000:])
        return ok, time.time() - start
    except subprocess.TimeoutExpired:
        with open(log_path, "w", encoding="utf-8") as f:
            f.write("TIMEOUT after 900s\n")
        return False, 900.0
    except Exception as e:  # noqa: BLE001
        with open(log_path, "w", encoding="utf-8") as f:
            f.write(f"EXCEPTION: {e}\n")
        return False, 0.0


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    skip_existing = "--skip-existing" in sys.argv
    src_dir = args[0] if len(args) > 0 else "app/articles"
    worker = int(args[1]) if len(args) > 1 else 4
    subset = set(args[2].split(",")) if len(args) > 2 else None

    out_dir = os.path.join(src_dir, "pdf_export")
    log_dir = os.path.join(out_dir, "_logs")
    work_dir = os.path.join(out_dir, "_work")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)
    header = os.path.join(log_dir, "_header.tex")
    with open(header, "w", encoding="utf-8") as f:
        f.write(HEADER_TEX)

    md_files = sorted(f for f in os.listdir(src_dir) if f.lower().endswith(".md"))
    if not md_files:
        print("no md files found in", src_dir)
        return

    args = [a for a in sys.argv[1:] if not a.startswith("--")]

    results = {"ok": [], "fail": [], "skip": []}
    jobs = []
    with ThreadPoolExecutor(max_workers=worker) as pool:
        for fname in md_files:
            if subset and fname not in subset:
                continue
            out_pdf = os.path.join(out_dir, fname[:-3] + ".pdf")
            log_path = os.path.join(log_dir, fname[:-3] + ".log")
            if skip_existing and os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 1000:
                results["skip"].append(fname)
                continue
            jobs.append((fname, pool.submit(convert_one,
                        os.path.join(src_dir, fname), out_pdf, log_path, header, work_dir)))
        for fname, fut in jobs:
            ok, dur = fut.result()
            if ok:
                results["ok"].append(fname)
                print(f"  OK   {dur:5.1f}s  {fname}")
            else:
                results["fail"].append(fname)
                print(f"  FAIL {dur:5.1f}s  {fname}")

    print("\n==== SUMMARY ====")
    print(f"total: {len(jobs)}  ok: {len(results['ok'])}  fail: {len(results['fail'])}  skipped(existing): {len(results['skip'])}")
    if results["fail"]:
        print("failed files:")
        for f in results["fail"]:
            print("  -", f)
    print(f"output dir: {out_dir}")


if __name__ == "__main__":
    main()
