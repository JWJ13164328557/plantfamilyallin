#!/usr/bin/env python3
import argparse
import re
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA20)

# Kyte-Doolittle hydropathy (for GRAVY, if you ever need custom)
# but we still use ana.gravy() if available.

def aliphatic_index(seq: str) -> float:
    """
    Aliphatic index (Ikai 1980):
      AI = X(Ala) + 2.9*X(Val) + 3.9*(X(Ile)+X(Leu))
    where X(AA) = mole percent of that AA.
    """
    n = len(seq)
    if n == 0:
        return 0.0
    seq = seq.upper()
    xA = 100.0 * seq.count("A") / n
    xV = 100.0 * seq.count("V") / n
    xI = 100.0 * seq.count("I") / n
    xL = 100.0 * seq.count("L") / n
    return xA + 2.9 * xV + 3.9 * (xI + xL)

def sanitize_seq(raw: str):
    """
    Keep standard AA + X.
    Track invalid characters (e.g. '.', '-', etc).
    """
    raw = raw.strip().upper().replace("*", "")
    if not raw:
        return "", 0, 0, 0.0
    invalid = sum(1 for c in raw if (c not in AA_SET and c != "X"))
    x_count = raw.count("X")
    # remove invalid characters but keep X (ProteinAnalysis tolerates X poorly in some versions)
    # To be safest, remove X for property calculation, but keep X stats.
    cleaned = "".join([c for c in raw if c in AA_SET])
    invalid_frac = invalid / len(raw) if len(raw) else 0.0
    return cleaned, invalid, x_count, invalid_frac

def load_wolf_tsv(path: str):
    """
    wolfpsort_predict.py output:
      seq_id  wolf_loc  wolf_score  wolf_scores
    """
    d = {}
    if not path:
        return d
    with open(path) as f:
        header = f.readline()
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            seq_id = parts[0]
            loc = parts[1]
            score = parts[2]
            scores = parts[3] if len(parts) >= 4 else ""
            d[seq_id] = (loc, score, scores)
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pep", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--wolf_tsv", default="", help="Optional WoLF PSORT TSV to merge")
    args = ap.parse_args()

    wolf = load_wolf_tsv(args.wolf_tsv)

    # header
    cols = [
        "gene_id",
        "length_aa",
        "mw",
        "pi",
        "instability_index",
        "aliphatic_index",
        "gravy",
        "aromaticity",
        "charge_pH7",
        "ext_coeff_reduced",
        "ext_coeff_cystines",
        "A280_1mgml_reduced",
        "A280_1mgml_cystines",
        "x_count",
        "invalid_char_count",
        "invalid_char_frac",
    ]
    # 20 AA composition percent
    cols += [f"pct_{aa}" for aa in AA20]

    # wolfpsort columns
    cols += ["wolf_loc", "wolf_score", "wolf_scores"]

    with open(args.out, "w") as o:
        o.write(",".join(cols) + "\n")

        for r in SeqIO.parse(args.pep, "fasta"):
            raw = str(r.seq)
            seq, invalid_n, x_n, invalid_frac = sanitize_seq(raw)
            if not seq:
                continue

            ana = ProteinAnalysis(seq)

            length = len(seq)
            mw = ana.molecular_weight()
            pi = ana.isoelectric_point()
            instab = ana.instability_index()

            # Biopython may or may not have gravy in older versions; most have it.
            try:
                gravy = ana.gravy()
            except Exception:
                gravy = float("nan")

            arom = ana.aromaticity()

            # charge at pH 7
            try:
                charge7 = ana.charge_at_pH(7.0)
            except Exception:
                charge7 = float("nan")

            # extinction coefficient
            # returns (reduced, cystines)
            try:
                ext_red, ext_cys = ana.molar_extinction_coefficient()
            except Exception:
                ext_red, ext_cys = (float("nan"), float("nan"))

            # A280 for 1 mg/mL (1 g/L), pathlength 1 cm:
            # A = ext_coeff * (1 g/L)/MW = ext_coeff / MW
            a280_red = ext_red / mw if mw and ext_red == ext_red else float("nan")
            a280_cys = ext_cys / mw if mw and ext_cys == ext_cys else float("nan")

            ai = aliphatic_index(seq)

            # AA composition (percent)
            aa_pct = ana.get_amino_acids_percent()  # dict {AA: fraction}
            pct_list = [100.0 * aa_pct.get(aa, 0.0) for aa in AA20]

            wolf_loc, wolf_score, wolf_scores = ("NA", "NA", "NA")
            if r.id in wolf:
                wolf_loc, wolf_score, wolf_scores = wolf[r.id]

            row = [
                r.id,
                str(length),
                f"{mw:.3f}",
                f"{pi:.3f}",
                f"{instab:.3f}",
                f"{ai:.3f}",
                f"{gravy:.3f}" if gravy == gravy else "NA",
                f"{arom:.5f}",
                f"{charge7:.3f}" if charge7 == charge7 else "NA",
                str(ext_red),
                str(ext_cys),
                f"{a280_red:.6f}" if a280_red == a280_red else "NA",
                f"{a280_cys:.6f}" if a280_cys == a280_cys else "NA",
                str(x_n),
                str(invalid_n),
                f"{invalid_frac:.5f}",
            ]
            row += [f"{v:.3f}" for v in pct_list]
            row += [wolf_loc, wolf_score, wolf_scores]

            o.write(",".join(row) + "\n")

if __name__ == "__main__":
    main()
