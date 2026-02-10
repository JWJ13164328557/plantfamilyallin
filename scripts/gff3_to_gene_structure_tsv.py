#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="GFF3 (clean.filtered.gff3)")
    ap.add_argument("--ids", required=True, help="final_family_members.list (transcript IDs)")
    ap.add_argument("--out", required=True, help="gene_structure.tsv: seq_id feature start end (bp)")
    ap.add_argument("--features", default="exon,CDS,five_prime_UTR,three_prime_UTR,UTR",
                    help="comma-separated features to export")
    return ap.parse_args()

def load_ids(path):
    s = set()
    with open(path) as f:
        for line in f:
            x = line.strip()
            if x:
                s.add(x)
    return s

def attr_get(attr, key):
    # key=...; pattern
    m = re.search(rf"(?:^|;){re.escape(key)}=([^;]+)", attr)
    return m.group(1) if m else None

def main():
    args = parse_args()
    keep = load_ids(args.ids)
    feats = set([x.strip() for x in args.features.split(",") if x.strip()])

    rows = []
    with open(args.gff, "r", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr = parts
            if ftype not in feats:
                continue

            # GFF3 子特征通常用 Parent 指向 transcript/mRNA
            parent = attr_get(attr, "Parent")
            if not parent:
                continue
            # Parent 可能是逗号分隔多个
            parents = parent.split(",")
            # 只保留属于家族 transcript 的行
            ok_parents = [p for p in parents if p in keep]
            if not ok_parents:
                continue

            try:
                s = int(start)
                e = int(end)
            except Exception:
                continue
            if e < s:
                s, e = e, s

            # 统一 UTR 命名（可选）
            feat_out = ftype
            if ftype in ("five_prime_UTR", "three_prime_UTR"):
                feat_out = "UTR"

            for pid in ok_parents:
                rows.append((pid, feat_out, s, e))

    rows = sorted(set(rows), key=lambda x: (x[0], x[1], x[2], x[3]))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as o:
        o.write("seq_id\tfeature\tstart\tend\n")
        for r in rows:
            o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")

if __name__ == "__main__":
    main()
