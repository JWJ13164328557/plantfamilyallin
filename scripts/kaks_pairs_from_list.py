#!/usr/bin/env python3
import argparse
import itertools

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene_list", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    genes = []
    with open(args.gene_list) as f:
        for line in f:
            g = line.strip().split()[0]
            if g:
                genes.append(g)

    genes = sorted(set(genes))
    with open(args.out, "w") as w:
        w.write("geneA\tgeneB\ttype\n")
        for a, b in itertools.combinations(genes, 2):
            w.write(f"{a}\t{b}\tparalog_family\n")

if __name__ == "__main__":
    main()
