#!/usr/bin/env python3
import argparse, glob, os, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_one(axt, kaks, method):
    out = axt + ".kaks"
    subprocess.check_call(
        [kaks, "-i", axt, "-o", out, "-m", method],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    return out

def parse_kaks_file(fp):
    """
    KaKs_Calculator 输出常见格式（示例）：
    Sequence  Method  Ka  Ks  Ka/Ks
    geneA-geneB  YN    0.12  0.34  0.35
    或者可能没有表头/以 # 开头。
    """
    rows = []
    with open(fp) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#") or line.lower().startswith("sequence"):
                continue
            parts = line.split()
            # 期望至少 5 列：seq method ka ks kaks
            if len(parts) >= 5:
                seq, method, ka, ks, wk = parts[0], parts[1], parts[2], parts[3], parts[4]
                rows.append((seq, method, ka, ks, wk))
            # 兜底：如果只有 4 列（seq ka ks wk），method 用 NA
            elif len(parts) == 4:
                seq, ka, ks, wk = parts[0], parts[1], parts[2], parts[3]
                rows.append((seq, "NA", ka, ks, wk))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axt_dir", required=True)
    ap.add_argument("--kaks", required=True)
    ap.add_argument("--method", default="YN")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    axts = sorted(glob.glob(os.path.join(args.axt_dir, "*.axt")))
    if not axts:
        raise SystemExit(f"[ERROR] no .axt found in {args.axt_dir}")

    outfiles = []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = [ex.submit(run_one, axt, args.kaks, args.method) for axt in axts]
        for fu in as_completed(futs):
            outfiles.append(fu.result())

    with open(args.out, "w") as w:
        w.write("pair\tmethod\tKa\tKs\tKaKs\n")
        for fp in sorted(outfiles):
            pair = os.path.basename(fp).replace(".axt.kaks", "")
            for _seq, method, ka, ks, wk in parse_kaks_file(fp):
                w.write(f"{pair}\t{method}\t{ka}\t{ks}\t{wk}\n")

if __name__ == "__main__":
    main()
