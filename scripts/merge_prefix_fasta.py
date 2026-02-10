#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target_fa", required=True)
    ap.add_argument("--model_fa", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--target_prefix", default="Target|")
    ap.add_argument("--model_prefix", default="Model|")
    args = ap.parse_args()

    recs = []
    seen = set()

    def add(fa, prefix):
        nonlocal recs, seen
        for r in SeqIO.parse(fa, "fasta"):
            rid = r.id.split()[0]
            new_id = f"{prefix}{rid}"
            if new_id in seen:
                i = 2
                while f"{new_id}__{i}" in seen:
                    i += 1
                new_id = f"{new_id}__{i}"
            seen.add(new_id)
            r.id = new_id
            r.name = new_id
            r.description = ""
            recs.append(r)

    add(args.target_fa, args.target_prefix)
    add(args.model_fa, args.model_prefix)

    SeqIO.write(recs, args.out, "fasta")

if __name__ == "__main__":
    main()
