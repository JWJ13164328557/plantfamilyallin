#!/usr/bin/env python3
import argparse
import re


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--collinearity", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--min_block_hits",
        type=int,
        default=5,
        help="Only keep blocks with >= this many hits/anchors (rough filter).",
    )
    args = ap.parse_args()

    # MCScanX .collinearity typical structure:
    # ## Alignment 1: score=... e_value=... N=23
    #  0-  0:   geneA   geneB   evalue
    #  0-  1:   geneA   geneB   evalue
    # (blank line)
    aln_re = re.compile(r"^##\s*Alignment\s+\d+.*\bN=(\d+)\b", re.I)

    pairs = []
    keep_block = False

    with open(args.collinearity) as f:
        for line in f:
            line = line.strip()
            if not line:
                keep_block = False
                continue

            # comment lines
            if line.startswith("#"):
                m = aln_re.match(line)
                if m:
                    n = int(m.group(1))
                    keep_block = n >= args.min_block_hits
                continue

            if not keep_block:
                continue

            parts = re.split(r"\s+", line)
            if len(parts) < 2:
                continue

            # Case 1 (your file): "0-  0:  geneA  geneB  evalue"
            # Then geneA=parts[2], geneB=parts[3]
            if len(parts) >= 4 and re.fullmatch(r"\d+-", parts[0]) and re.fullmatch(r"\d+:", parts[1]):
                g1, g2 = parts[2], parts[3]
            else:
                # Case 2 (some MCScanX variants): "geneA geneB evalue ..."
                g1, g2 = parts[0], parts[1]

            if not g1 or not g2:
                continue
            if g1 == g2:
                continue

            pairs.append((g1, g2))

    # de-duplicate: same pair can appear multiple times
    seen = set()
    uniq = []
    for a, b in pairs:
        key = (a, b) if a < b else (b, a)
        if key in seen:
            continue
        seen.add(key)
        uniq.append((a, b))

    with open(args.out, "w") as w:
        w.write("geneA\tgeneB\ttype\n")
        for a, b in uniq:
            w.write(f"{a}\t{b}\tsyntenic\n")


if __name__ == "__main__":
    main()
