#!/usr/bin/env python3
import argparse
from collections import defaultdict
from Bio import SeqIO

def get_prefix(rid: str) -> str:
    rid = rid.split()[0]
    return rid.rsplit(".", 1)[0] if "." in rid else rid

def load_best_by_prefix(fa_path: str):
    """
    return:
      by_prefix: dict[prefix] -> dict[iso_id] -> SeqRecord
    """
    by_prefix = defaultdict(dict)
    for r in SeqIO.parse(fa_path, "fasta"):
        rid = r.id.split()[0]
        prefix = get_prefix(rid)
        by_prefix[prefix][rid] = r
    return by_prefix

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cds", required=True)
    ap.add_argument("--pep", required=True)
    ap.add_argument("--out_cds", required=True)
    ap.add_argument("--out_pep", required=True)
    ap.add_argument("--out_map", required=True)
    ap.add_argument("--prefer", choices=["cds", "pep"], default="cds",
                    help="choose longest isoform by cds length (default) or pep length")
    args = ap.parse_args()

    cds = load_best_by_prefix(args.cds)
    pep = load_best_by_prefix(args.pep)

    common_prefixes = sorted(set(cds.keys()) & set(pep.keys()))

    chosen = {}  # prefix -> (cds_rec, pep_rec)
    for prefix in common_prefixes:
        cds_ids = set(cds[prefix].keys())
        pep_ids = set(pep[prefix].keys())
        common_ids = sorted(cds_ids & pep_ids)
        if not common_ids:
            # prefix exists in both, but ids don't match exactly
            # fallback: pick longest independently (original behavior)
            best_c = max(cds[prefix].values(), key=lambda r: len(r.seq))
            best_p = max(pep[prefix].values(), key=lambda r: len(r.seq))
            chosen[prefix] = (best_c, best_p)
            continue

        # choose one isoform id shared by both cds & pep
        if args.prefer == "cds":
            best_id = max(common_ids, key=lambda rid: (len(cds[prefix][rid].seq), len(pep[prefix][rid].seq)))
        else:
            best_id = max(common_ids, key=lambda rid: (len(pep[prefix][rid].seq), len(cds[prefix][rid].seq)))

        chosen[prefix] = (cds[prefix][best_id], pep[prefix][best_id])

    # write map
    with open(args.out_map, "w") as f:
        f.write("gene_prefix\tcds_id\tpep_id\tcds_len\tpep_len\n")
        for prefix in sorted(chosen.keys()):
            c, p = chosen[prefix]
            f.write(f"{prefix}\t{c.id}\t{p.id}\t{len(c.seq)}\t{len(p.seq)}\n")

    # always create output files (even empty)
    SeqIO.write([chosen[p][0] for p in sorted(chosen.keys())], args.out_cds, "fasta")
    SeqIO.write([chosen[p][1] for p in sorted(chosen.keys())], args.out_pep, "fasta")

if __name__ == "__main__":
    main()
