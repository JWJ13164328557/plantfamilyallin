#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kaks_raw", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--min_ks", type=float, required=True)
    ap.add_argument("--max_ks", type=float, required=True)
    ap.add_argument("--max_w", type=float, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # 读 raw
    dt = pd.read_csv(args.kaks_raw, sep="\t", dtype=str)

    # 自动兼容 4列/5列：如果有 method，就丢掉它
    # 目标统一为: pair, Ka, Ks, KaKs
    if "method" in dt.columns:
        dt = dt[["pair", "Ka", "Ks", "KaKs"]]
    else:
        # 兼容旧格式：pair Ka Ks KaKs
        dt = dt.iloc[:, :4]
        dt.columns = ["pair", "Ka", "Ks", "KaKs"]

    # 转数值并过滤 NA/Inf
    for c in ["Ka", "Ks", "KaKs"]:
        dt[c] = pd.to_numeric(dt[c], errors="coerce")

    dt = dt.dropna(subset=["Ka", "Ks", "KaKs"])

    # 过滤阈值
    dt = dt[(dt["Ks"] >= args.min_ks) & (dt["Ks"] <= args.max_ks) & (dt["KaKs"] <= args.max_w)]

    # 加 type（如果你需要 pairs.tsv 的类型）
    # pairs.tsv: geneA geneB type
    pairs = pd.read_csv(args.pairs, sep="\t", dtype=str)
    if set(["geneA","geneB","type"]).issubset(pairs.columns):
        pairs["pair"] = pairs["geneA"] + "__" + pairs["geneB"]
        dt = dt.merge(pairs[["pair","type"]], on="pair", how="left")
    else:
        dt["type"] = "NA"

    dt.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
