from typing import List, Dict
from .io import middle_window

def search_one_read(
    read_id: str,
    read_seq: str,
    idx_fwd,
    idx_rev,
    genome_name: str,
) -> List[Dict]:
    """
    Search order:
      1) middle-50 window exact
      2) right-trim (keep left end) 50→20
      3) left-trim (keep right end) 50→20
      4) repeat on reverse-complement genome index

    Multimapping rule:
      - Determine the longest exact match length L (>=20) for the read.
      - Report ALL matches of that L-nt sequence (both strands).
      - If that L-nt sequence occurs >3 times, discard the read entirely.
    """
    hits: List[Dict] = []

    def try_windows(seq: str, strand: str):
        nonlocal hits
        s, e = middle_window(seq, 50)
        L = e - s

        # Step 1: Full window
        _try_subseq(read_id, seq, s, e, strand, idx_fwd if strand == '+' else idx_rev, genome_name, hits)

        # Step 2: right-trim keeping left end fixed
        if len(hits) == 0:
            for new_len in range(L - 1, 19, -1):  # 49..20
                _try_subseq(read_id, seq, s, s + new_len, strand, idx_fwd if strand == '+' else idx_rev, genome_name, hits)
                if hits:
                    break

        # Step 3: left-trim keeping right end fixed
        if len(hits) == 0:
            for shift in range(1, 50):  # safe upper bound
                if e - (s + shift) < 20:
                    break
                _try_subseq(read_id, seq, s + shift, e, strand, idx_fwd if strand == '+' else idx_rev, genome_name, hits)
                if hits:
                    break

    # Forward genome first; stop at first success so we know the current longest
    try_windows(read_seq, '+')

    # If none on forward, try reverse-complement genome
    if not hits:
        try_windows(read_seq, '-')

    # If still no hits, return empty
    if not hits:
        return []

    # --- Multimapping rule implementation ---
    # Determine the longest match length L among found hits
    max_len = max(h["match_len"] for h in hits)

    # Now we must gather ALL matches of that exact L-nt sequence on BOTH strands.
    # Identify the L-nt sequence(s) we matched (there can be >1 seq if different subseqs tie at L)
    longest_seqs = {h["sequence"] for h in hits if h["match_len"] == max_len}

    # Query both indexes for every such L-nt sequence to enumerate ALL genomic sites
    # (This ensures we don't accidentally limit to the first bucket we checked.)
    all_longest_hits: List[Dict] = []
    for lseq in longest_seqs:
        for strand, index in (('+', idx_fwd), ('-', idx_rev)):
            locs = index.lookup(lseq)
            for pos in locs:
                all_longest_hits.append({
                    "read_id": read_id,
                    "strand": strand,
                    "genome": genome_name,
                    "genome_start": pos,                 # 0-based
                    "genome_end": pos + len(lseq),       # exclusive
                    "read_start": None,                   # optional; not critical for rule
                    "read_end": None,                     # optional; not critical for rule
                    "sequence": lseq,
                    "match_len": len(lseq)
                })

    # Deduplicate by strand+coords+sequence
    dedup = []
    seen = set()
    for h in all_longest_hits:
        key = (h["strand"], h["genome_start"], h["genome_end"], h["sequence"])
        if key not in seen:
            seen.add(key)
            dedup.append(h)

    # Count distinct sites
    unique_sites = {(h["strand"], h["genome_start"], h["genome_end"]) for h in dedup}
    if len(unique_sites) > 3:
        # Highly repetitive → discard the read entirely
        return []

    # Otherwise, return all longest sites (1–3)
    # Sort for deterministic output: by strand, then genome_start
    dedup.sort(key=lambda x: (x["strand"], x["genome_start"]))
    return dedup

def _try_subseq(read_id, read_seq, rs, re, strand, index, genome_name, hits_out):
    subseq = read_seq[rs:re]
    if len(subseq) < 20:
        return
    locs = index.lookup(subseq)
    for pos in locs:
        hits_out.append({
            "read_id": read_id,
            "strand": strand,
            "genome": genome_name,
            "genome_start": pos,             # 0-based
            "genome_end": pos + len(subseq), # exclusive
            "read_start": rs,
            "read_end": re,
            "sequence": subseq,
            "match_len": len(subseq)
        })
