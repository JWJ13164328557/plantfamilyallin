#!/usr/bin/env python3
import argparse

def read_set(p):
    s=set()
    with open(p) as f:
        for line in f:
            line=line.strip()
            if line:
                s.add(line.split()[0])
    return s

def parse_attr(attr):
    d={}
    for kv in attr.split(";"):
        if "=" in kv:
            k,v=kv.split("=",1)
            d[k]=v
    return d

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    genes=read_set(args.genes)
    out=[]

    with open(args.gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a=line.rstrip("\n").split("\t")
            if len(a) != 9:
                continue
            chrom, src, ftype, start, end, score, strand, phase, attr = a
            if ftype != "gene":
                continue
            info=parse_attr(attr)
            gid=info.get("ID") or info.get("Name")
            if not gid:
                continue
            gprefix = gid.rsplit(".",1)[0] if "." in gid else gid
            if gid in genes or gprefix in genes:
                out.append((chrom, int(start)-1, int(end), gid, "0", strand))

    out.sort(key=lambda x:(x[0], x[1], x[2]))
    with open(args.out,"w") as o:
        for r in out:
            o.write("\t".join(map(str,r))+"\n")

if __name__=="__main__":
    main()
