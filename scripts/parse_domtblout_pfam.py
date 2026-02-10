#!/usr/bin/env python3
import argparse

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--pfam_ids", required=True, help="comma separated PFxxxxx list")
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    want=set([x.strip() for x in args.pfam_ids.split(",") if x.strip()])
    hits=set()

    with open(args.domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a=line.strip().split()
            if len(a) < 23:
                continue
            target=a[0]
            query=a[3]
            pfid=query.split(".")[0]
            if pfid in want:
                hits.add(target)

    with open(args.out,"w") as o:
        for x in sorted(hits):
            o.write(x+"\n")

if __name__=="__main__":
    main()
