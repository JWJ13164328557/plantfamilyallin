#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def read_list(p):
    s=set()
    with open(p) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            s.add(line.split()[0])
    return s

def write_list(s, out):
    with open(out,"w") as f:
        for x in sorted(s):
            f.write(x+"\n")

def blast_candidates(blast_tsv, out):
    s=set()
    with open(blast_tsv) as f:
        for line in f:
            a=line.rstrip("\n").split("\t")
            if len(a) < 12:
                continue
            s.add(a[1])  # sseqid
    write_list(s, out)

def final_members(blast, pfam, hmm, strategy, out_list, out_venn_tsv):
    A=read_list(blast)
    B=read_list(pfam)
    C=read_list(hmm)

    if strategy == "intersection":
        final = A & B & C
    elif strategy == "union":
        final = A | B | C
    elif strategy == "blast_and_domain":
        final = A & (B | C)
    else:
        raise SystemExit(f"Unknown strategy: {strategy}")

    all_ids = sorted(A | B | C)
    with open(out_venn_tsv,"w") as f:
        f.write("gene\tBLAST\tPFAM\tHMM\n")
        for g in all_ids:
            f.write(f"{g}\t{int(g in A)}\t{int(g in B)}\t{int(g in C)}\n")

    write_list(final, out_list)

def extract_fasta(fasta, ids, out):
    ids_set = read_list(ids)
    recs=[]
    for r in SeqIO.parse(fasta, "fasta"):
        rid = r.id.split()[0]
        if rid in ids_set or (rid.rsplit(".",1)[0] in ids_set):
            recs.append(r)
    SeqIO.write(recs, out, "fasta")

def main():
    ap=argparse.ArgumentParser()
    sub=ap.add_subparsers(dest="cmd", required=True)

    p1=sub.add_parser("blast_candidates")
    p1.add_argument("--blast_tsv", required=True)
    p1.add_argument("--out", required=True)

    p2=sub.add_parser("final_members")
    p2.add_argument("--blast", required=True)
    p2.add_argument("--pfam", required=True)
    p2.add_argument("--hmm", required=True)
    p2.add_argument("--strategy", required=True)
    p2.add_argument("--out_list", required=True)
    p2.add_argument("--out_venn_tsv", required=True)

    p3=sub.add_parser("extract_fasta")
    p3.add_argument("--fasta", required=True)
    p3.add_argument("--ids", required=True)
    p3.add_argument("--out", required=True)

    args=ap.parse_args()
    if args.cmd=="blast_candidates":
        blast_candidates(args.blast_tsv, args.out)
    elif args.cmd=="final_members":
        final_members(args.blast,args.pfam,args.hmm,args.strategy,args.out_list,args.out_venn_tsv)
    elif args.cmd=="extract_fasta":
        extract_fasta(args.fasta,args.ids,args.out)

if __name__=="__main__":
    main()
