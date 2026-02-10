#!/usr/bin/env python3
import argparse, re

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domtbl", required=True, help="hmmscan --domtblout output")
    ap.add_argument("--out", required=True, help="output tsv: seq_id domain start end")
    ap.add_argument("--only_pfam", default="", help="comma-separated PFAM IDs, e.g. PF00010,PF00249")
    ap.add_argument("--min_iE", type=float, default=1e-3, help="keep if i-Evalue <= this")
    return ap.parse_args()

def norm_pfam(acc: str) -> str:
    # PF00010.32 -> PF00010
    if acc is None:
        return ""
    return acc.split(".")[0]

def main():
    args = parse_args()

    keep = set()
    if args.only_pfam.strip():
        keep = {x.strip() for x in args.only_pfam.split(",") if x.strip()}
        # 也归一化一下，防止用户写 PF00010.32
        keep = {norm_pfam(x) for x in keep}

    out_rows = []
    with open(args.domtbl, "r", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            # domtblout 至少应有 23 列（description 之后更多也没事）
            if len(parts) < 23:
                continue

            # hmmscan domtblout:
            # 0 target_name
            # 1 target_accession (e.g. PF00010.32)
            # 3 query_name (sequence id)
            # 12 i-Evalue
            # 17 ali_from
            # 18 ali_to
            target_name = parts[0]
            target_acc  = parts[1]
            query_name  = parts[3]

            try:
                i_evalue = float(parts[12])
            except Exception:
                continue

            try:
                ali_from = int(parts[17])
                ali_to   = int(parts[18])
            except Exception:
                continue

            if i_evalue > args.min_iE:
                continue

            pf = norm_pfam(target_acc)
            if keep and pf not in keep:
                continue

            # domain 列你可以用 pf 或 target_name，看你画图想显示什么
            domain = pf if pf else target_name
            start = min(ali_from, ali_to)
            end   = max(ali_from, ali_to)
            out_rows.append((query_name, domain, start, end))

    # 写出（即使为空也写表头）
    with open(args.out, "w") as o:
        o.write("seq_id\tdomain\tstart\tend\n")
        for r in out_rows:
            o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")

if __name__ == "__main__":
    main()
