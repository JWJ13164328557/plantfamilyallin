#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fimo", required=True)
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    df = pd.read_csv(args.fimo, sep="\t", comment="#")
    if df.shape[0] == 0:
        pd.DataFrame({"motif_id":[], "count":[]}).to_csv(args.out, sep="\t", index=False)
        return
    s = df.groupby("motif_id").size().reset_index(name="count").sort_values("count", ascending=False)
    s.to_csv(args.out, sep="\t", index=False)

if __name__=="__main__":
    main()
