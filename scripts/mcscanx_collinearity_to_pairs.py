#!/usr/bin/env python3
import argparse, random

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--col", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max_pairs", type=int, default=0, help="0 means no limit")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    pairs = []
    seen = set()

    with open(args.col) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # 典型格式：idx geneA geneB ...
            if len(parts) >= 3 and parts[0].isdigit():
                a, b = parts[1], parts[2]
                if a == b:
                    continue
                key = (a, b) if a < b else (b, a)
                if key in seen:
                    continue
                seen.add(key)
                pairs.append((a, b, "syntenic"))

    if args.max_pairs and len(pairs) > args.max_pairs:
        random.seed(args.seed)
        pairs = random.sample(pairs, args.max_pairs)

    with open(args.out, "w") as w:
        w.write("geneA\tgeneB\ttype\n")
        for a, b, t in pairs:
            w.write(f"{a}\t{b}\t{t}\n")

if __name__ == "__main__":
    main()
