#!/usr/bin/env python3
import argparse

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    hits=set()
    with open(args.domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a=line.strip().split()
            if len(a) < 23:
                continue
            hits.add(a[0])

    with open(args.out,"w") as o:
        for x in sorted(hits):
            o.write(x+"\n")

if __name__=="__main__":
    main()
