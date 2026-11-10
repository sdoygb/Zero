#!/usr/bin/env python3
"""Kaggle 竞赛扫描器（匿名，无需 API key）。

背景：Kaggle 官方 kaggle-api 的 /api/v1/competitions/list 需要认证（401），
但站点前端用的内部 API 对公开列表是匿名可读的。本脚本就调它。

用法：
    python3 tools/kaggle_scan.py                 # 当前 ACTIVE 竞赛（含奖金/队数/截止）
    python3 tools/kaggle_scan.py --json          # 同上，输出 JSON
    python3 tools/kaggle_scan.py --lb <id> [--top 12]   # 某赛榜单分布

注：数据是抓取时刻的快照；榜单分数含义以各赛规则页为准。
"""
import argparse
import json
import statistics
import subprocess

BASE = "https://www.kaggle.com/api/i/"
HEADERS = {"User-Agent": "Mozilla/5.0"}


def rpc(method, payload):
    """优先 requests；没装就退回 curl。

    注意：标准库 urllib 打这个端点会被边缘节点回 404（同样的头，curl 却是 200），
    所以这里不用 urllib。
    """
    url = BASE + method
    body = json.dumps(payload)
    try:
        import requests
        r = requests.post(url, data=body, headers={**HEADERS, "Content-Type": "application/json"}, timeout=30)
        r.raise_for_status()
        return r.json()
    except ImportError:
        out = subprocess.run(
            ["curl", "-s", "-X", "POST", url, "-H", "Content-Type: application/json",
             "-H", f"User-Agent: {HEADERS['User-Agent']}", "--data", body],
            capture_output=True, text=True, check=True).stdout
        return json.loads(out)


def active():
    d = rpc("competitions.CompetitionService/ListCompetitions",
            {"selector": {"listOption": "ACTIVE"}, "pageSize": 50})
    return d.get("competitions", [])


def leaderboard(cid, page_size=50):
    d = rpc("competitions.LeaderboardService/GetLeaderboard",
            {"competitionId": cid, "pageSize": page_size})
    rows = d.get("publicLeaderboard", [])
    vals = [float(r["displayScore"]) for r in rows if r.get("displayScore") not in (None, "")]
    return rows, vals


def fmt(c):
    r = c.get("reward") or {}
    prize = f"${r['quantity']:,}" if r.get("quantity") else "-"
    return (f"{c['deadline'][:10]} | {c['title'][:46]:46s} | {prize:>9s} | "
            f"teams={c.get('totalTeams', 0):>5} | id={c['id']} | {c['competitionName']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--lb", type=int, help="competition id")
    ap.add_argument("--top", type=int, default=12)
    a = ap.parse_args()

    if a.lb:
        rows, vals = leaderboard(a.lb)
        print(f"competition {a.lb}: scored={len(vals)}")
        for r in rows[:a.top]:
            print(f"  rank {r['rank']:>4}  {r['displayScore']}")
        if vals:
            s = sorted(vals)
            n = len(s)
            print(f"  min={s[0]:.6g} p1={s[n//100]:.6g} p10={s[n//10]:.6g} "
                  f"median={statistics.median(s):.6g} max={s[-1]:.6g}")
        return

    comps = active()
    if a.json:
        print(json.dumps(comps, ensure_ascii=False, indent=2)[:200000])
        return
    print(f"ACTIVE competitions: {len(comps)}  (snapshot)")
    for c in comps:
        print(fmt(c))


if __name__ == "__main__":
    main()
