#!/usr/bin/env python3
import argparse
import subprocess
import shutil
import re
from typing import List, Tuple

BAD_PATTERNS = (
    "Usage:",
    "Command Line Parsing Error",
    "Do not know what to do with argument",
    "runWolfPsortSummary",
)

# 解析形如：nucl 12.5, cyto_nucl 7.5
PAIR_RE = re.compile(r"([A-Za-z_]+)\s+([0-9]+(?:\.[0-9]+)?)")

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pep", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cmd", default="wolfpsort", help="wolfpsort command in PATH")
    ap.add_argument("--organism", default="plant", help="plant|animal|fungi")
    return ap.parse_args()

def looks_like_real_output(stdout: str) -> bool:
    if not stdout or not stdout.strip():
        return False
    for bad in BAD_PATTERNS:
        if bad in stdout:
            return False
    # 至少应当存在一行能提取到 loc+score
    for line in stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # 要有 seq_id + 至少一个 (loc score)
        parts = line.split(None, 1)
        if len(parts) < 2:
            continue
        if PAIR_RE.search(parts[1]):
            return True
    return False

def run_wolfpsort(cmd: str, organism: str, pep: str) -> str:
    if shutil.which(cmd) is None:
        raise SystemExit(f"ERROR: cannot find wolfpsort command in PATH: {cmd}")

    # 兼容两种常见调用：
    #   wolfpsort plant input.fa
    #   wolfpsort -organism plant input.fa   (有些实现才支持)
    candidates = [
        [cmd, organism, pep],
        [cmd, "-organism", organism, pep],
    ]

    tried = []
    for c in candidates:
        tried.append(" ".join(c))
        p = subprocess.run(c, capture_output=True, text=True)
        # 有些程序把信息打到 stdout，有些打到 stderr；我们都看
        stdout = (p.stdout or "").strip()
        stderr = (p.stderr or "").strip()
        merged = "\n".join([x for x in [stdout, stderr] if x])

        if p.returncode == 0 and looks_like_real_output(merged):
            return merged

    raise SystemExit("ERROR: wolfpsort failed or output not recognized. Tried:\n" + "\n".join(tried))

def parse_line(rest: str) -> List[Tuple[str, float]]:
    # 用正则提取所有 loc score，对逗号/分号等天然鲁棒
    pairs = []
    for loc, sc in PAIR_RE.findall(rest):
        try:
            pairs.append((loc, float(sc)))
        except ValueError:
            continue
    return pairs

def main():
    args = parse_args()
    out_text = run_wolfpsort(args.cmd, args.organism, args.pep)

    with open(args.out, "w") as o:
        o.write("seq_id\twolf_loc\twolf_score\twolf_scores\n")
        for line in out_text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # seq_id + rest
            parts = line.split(None, 1)
            if len(parts) < 2:
                continue
            sid, rest = parts[0], parts[1].strip()

            pairs = parse_line(rest)
            if not pairs:
                # 保留原始 rest 以便排查
                o.write(f"{sid}\tNA\tNA\t{rest}\n")
                continue

            best_loc, best_sc = max(pairs, key=lambda x: x[1])
            score_str = ";".join([f"{loc}={sc:g}" for loc, sc in pairs])
            o.write(f"{sid}\t{best_loc}\t{best_sc:g}\t{score_str}\n")

if __name__ == "__main__":
    main()
