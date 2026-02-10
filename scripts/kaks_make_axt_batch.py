#!/usr/bin/env python3
import argparse, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def read_fasta(fp):
    seqs = {}
    name = None
    buf = []
    with open(fp) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.replace(" ", "").replace("\r", ""))
        if name:
            seqs[name] = "".join(buf).upper()
    return seqs

def translate(cds):
    aa = []
    cds = cds.upper().replace("U", "T")
    # 尽量按三联体翻译；不够 3 的尾巴直接丢弃
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3]
        if any(c not in "ACGT" for c in codon):
            aa.append("X")
        else:
            aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)

def parse_pairs(fp):
    pairs = []
    with open(fp) as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue
            a, b, t = line.rstrip("\n").split("\t")[:3]
            pairs.append((a, b, t))
    return pairs

def fasta2kaks_axt(two_seq_fa: str, out_axt: str, nameA: str, nameB: str):
    """
    将 pal2nal 输出的 codon alignment fasta(2条序列) 写成 KaKs_Calculator 常用 AXT：
      >A-B
      SEQ_A
      SEQ_B
    """
    seqs = read_fasta(two_seq_fa)
    if len(seqs) != 2:
        raise ValueError(f"codon alignment not 2 seqs: {two_seq_fa} (got {len(seqs)})")

    # pal2nal 输出的 header 通常就是 gene id
    # 我们按 nameA/nameB 顺序写，避免 dict 无序导致对换
    if nameA not in seqs or nameB not in seqs:
        # 兜底：取任意两个
        keys = list(seqs.keys())
        s1, s2 = seqs[keys[0]], seqs[keys[1]]
    else:
        s1, s2 = seqs[nameA], seqs[nameB]

    if len(s1) != len(s2):
        raise ValueError(f"alignment length mismatch: {two_seq_fa}")

    # 密码子对齐一般应为 3 的倍数（含 gap 也占位）
    if len(s1) % 3 != 0 or len(s2) % 3 != 0:
        # 不强制退出也行，但 KaKs 可能会报错；这里直接报错更干净
        raise ValueError(f"codon aln length not multiple of 3: {two_seq_fa} len={len(s1)}")

    # 只保留 ACGTN-（pal2nal 可能会输出 N、-）
    def sanitize(x):
        x = x.upper()
        ok = set("ACGTN-")
        return "".join((c if c in ok else "N") for c in x)

    s1 = sanitize(s1)
    s2 = sanitize(s2)

    with open(out_axt, "w") as w:
        w.write(f">{nameA}-{nameB}\n")
        w.write(s1 + "\n")
        w.write(s2 + "\n")

def build_one(pair, cds_map, outdir, mafft, pal2nal):
    geneA, geneB, _type = pair
    if geneA not in cds_map or geneB not in cds_map:
        return (None, f"missing CDS for {geneA} or {geneB}")

    cdsA = cds_map[geneA]
    cdsB = cds_map[geneB]

    base = f"{geneA}__{geneB}".replace("|","_")
    tmp_dir = os.path.join(outdir, "_tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    cds_fp = os.path.join(tmp_dir, base + ".cds.fa")
    with open(cds_fp, "w") as w:
        w.write(f">{geneA}\n{cdsA}\n>{geneB}\n{cdsB}\n")

    pep_fp = os.path.join(tmp_dir, base + ".pep.fa")
    with open(pep_fp, "w") as w:
        w.write(f">{geneA}\n{translate(cdsA)}\n>{geneB}\n{translate(cdsB)}\n")

    pep_aln = os.path.join(tmp_dir, base + ".pep.aln.fa")
    with open(pep_aln, "w") as out:
        subprocess.check_call([mafft, "--auto", pep_fp], stdout=out, stderr=subprocess.DEVNULL)

    codon_aln = os.path.join(tmp_dir, base + ".codon.aln.fa")
    with open(codon_aln, "w") as out:
        subprocess.check_call(
            ["perl", pal2nal, pep_aln, cds_fp, "-output", "fasta"],
            stdout=out, stderr=subprocess.DEVNULL
        )

    axt_fp = os.path.join(outdir, base + ".axt")
    fasta2kaks_axt(codon_aln, axt_fp, geneA, geneB)

    return (axt_fp, None)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--cds_fa", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mafft", default="mafft")
    ap.add_argument("--pal2nal", required=True)
    # 兼容旧参数：不再需要，但保留不报错
    ap.add_argument("--axtconvertor", required=False, default="")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cds_map = read_fasta(args.cds_fa)
    pairs = parse_pairs(args.pairs)

    ok_fp = os.path.join(args.outdir, "ok.tsv")
    fail_fp = os.path.join(args.outdir, "failed.tsv")
    ok = 0
    fail = 0

    with open(ok_fp, "w") as okw, open(fail_fp, "w") as fw:
        okw.write("geneA\tgeneB\taxt\n")
        fw.write("geneA\tgeneB\treason\n")

        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = {ex.submit(build_one, p, cds_map, args.outdir, args.mafft, args.pal2nal): p for p in pairs}
            for fu in as_completed(futs):
                geneA, geneB, _t = futs[fu]
                try:
                    axt, err = fu.result()
                    if axt:
                        ok += 1
                        okw.write(f"{geneA}\t{geneB}\t{axt}\n")
                    else:
                        fail += 1
                        fw.write(f"{geneA}\t{geneB}\t{err}\n")
                except Exception as e:
                    fail += 1
                    fw.write(f"{geneA}\t{geneB}\tEXCEPTION: {e}\n")

    if ok == 0:
        print("[ERROR] No AXT generated. Check gene IDs between pairs and CDS fasta.", file=sys.stderr)
        sys.exit(2)

    print(f"[INFO] AXT generated: {ok}, failed: {fail}", file=sys.stderr)

if __name__ == "__main__":
    main()
