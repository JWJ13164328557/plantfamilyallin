#!/usr/bin/env python3
import argparse, re
from Bio import SeqIO

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pep", required=True)
    ap.add_argument("--meme_txt", required=True)
    ap.add_argument("--out_hits", required=True)
    ap.add_argument("--out_lens", required=True)
    return ap.parse_args()

def read_lengths(pep_fa):
    lens = {}
    for r in SeqIO.parse(pep_fa, "fasta"):
        seq = str(r.seq).replace("*", "")
        rid = r.id.split()[0]
        if seq:
            lens[rid] = len(seq)
    return lens

# MOTIF <consensus> MEME-1 width = 29 ...
RE_MOTIF = re.compile(r"^MOTIF\s+\S+.*\b(MEME-\d+)\b.*?\bwidth\s*=\s*(\d+)\b", re.I)

# Motif <consensus> MEME-1 sites sorted by position
RE_SITES_HEADER = re.compile(r"^Motif\s+\S+.*\b(MEME-\d+)\b\s+sites\s+sorted\s+by\s+position", re.I)

# 表头里常见字段（不同 MEME 版本略有差异）
# 常见：Sequence name | Start | P-value | ...
RE_TABLE_HEADER = re.compile(r"^Sequence\s+name", re.I)

def meme_id_to_num(tag: str) -> str:
    m = re.search(r"MEME-(\d+)$", tag)
    return m.group(1) if m else tag

def parse_int_any(s: str):
    """
    从字符串里抓第一个整数（允许 +/- 以及括号等）
    """
    m = re.search(r"[-+]?\d+", s)
    return int(m.group(0)) if m else None

def main():
    args = parse_args()
    lens = read_lengths(args.pep)

    # 输出蛋白长度
    with open(args.out_lens, "w") as o:
        o.write("seq_id\tlength\n")
        for k in sorted(lens.keys()):
            o.write(f"{k}\t{lens[k]}\n")

    motif_width = {}   # motif_id(num str) -> width int
    hits = set()       # 用 set 去重 (seq_id, motif, start, end)

    in_table = False
    curr_motif = None
    start_col = None   # start 列索引（从 0 开始）

    with open(args.meme_txt, "r", errors="ignore") as f:
        for line in f:
            s = line.rstrip("\n")

            # 1) motif 宽度
            m = RE_MOTIF.match(s.strip())
            if m:
                mid = meme_id_to_num(m.group(1))  # MEME-1 -> "1"
                w = int(m.group(2))
                motif_width[mid] = w
                continue

            # 2) sites 表开始
            mh = RE_SITES_HEADER.match(s.strip())
            if mh:
                curr_motif = meme_id_to_num(mh.group(1))
                in_table = True
                start_col = None
                continue

            if not in_table:
                continue

            st = s.strip()
            low = st.lower()

            # 3) 表结束判断（进入其他段落）
            if (st.startswith("Combined") or st.startswith("SUMMARY")
                or "position-specific" in low
                or (low.startswith("motif ") and "sites sorted by position" not in low)):
                in_table = False
                curr_motif = None
                start_col = None
                continue

            # 4) 跳过空行/分隔线
            if not st or st.startswith("-"):
                continue

            # 5) 捕获表头行：确定 start 列在哪
            # 典型表头：Sequence name  Start  P-value  ...
            if RE_TABLE_HEADER.match(st):
                cols = re.split(r"\s{2,}|\t+", st.strip())
                # 找 “Start” 的列
                for i, c in enumerate(cols):
                    if c.lower() == "start":
                        start_col = i
                        break
                continue

            # 6) 跳过子表头
            if low.startswith("sequence") or low.startswith("-------"):
                continue

            # 7) 解析数据行
            if curr_motif is None:
                continue

            # MEME 行通常用多个空格对齐，这里用“>=2 空格/Tab”切割更稳
            parts = re.split(r"\s{2,}|\t+", st)
            if len(parts) < 2:
                continue

            seq_id = parts[0].lstrip(">").strip().rstrip("*")
            if seq_id not in lens:
                continue

            # 取 start
            start = None
            if start_col is not None and start_col < len(parts):
                start = parse_int_any(parts[start_col])
            if start is None:
                # fallback：从剩余字段里找第一个整数
                for tok in parts[1:]:
                    start = parse_int_any(tok)
                    if start is not None:
                        break
            if start is None:
                continue

            w = motif_width.get(curr_motif)
            if w is None:
                continue

            end = start + w - 1

            # 越界过滤：保证在蛋白长度内
            L = lens[seq_id]
            if start < 1 or end < start or end > L:
                continue

            motif_name = f"Motif{curr_motif}"
            hits.add((seq_id, motif_name, int(curr_motif), start, end))

    # 输出 hits（按 seq_id 再按 motif_id 再按 start 排序）
    with open(args.out_hits, "w") as o:
        o.write("seq_id\tmotif\tmotif_id\tstart\tend\n")
        for seq_id, motif_name, motif_id, start, end in sorted(hits, key=lambda x: (x[0], x[2], x[3], x[4])):
            o.write(f"{seq_id}\t{motif_name}\t{motif_id}\t{start}\t{end}\n")

if __name__ == "__main__":
    main()
