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
    # --- find all longest L-mers within the read ---
    max_len = max(h["match_len"] for h in hits)

    def all_lmers(seq: str, L: int):
        """Return every distinct L-mer substring of seq."""
        s = set()
        for i in range(0, len(seq) - L + 1):
            frag = seq[i:i+L]
            if "N" not in frag:
                s.add(frag)
        return s

    longest_seqs = all_lmers(read_seq, max_len)
    if not longest_seqs:
        return []

    all_hits = []
    for frag in longest_seqs:
        # lookup both strands
        fwd_sites = idx_fwd.lookup(frag)
        rev_sites = idx_rev.lookup(frag)

        # if this fragment maps >3 total times → discard read entirely
        if len(fwd_sites) + len(rev_sites) > 3:
            return []

        for pos in fwd_sites:
            all_hits.append({
                "read_id": read_id,
                "strand": "+",
                "genome": genome_name,
                "genome_start": pos,
                "genome_end": pos + max_len,
                "sequence": frag,
                "match_len": max_len
            })
        for pos in rev_sites:
            all_hits.append({
                "read_id": read_id,
                "strand": "-",
                "genome": genome_name,
                "genome_start": pos,
                "genome_end": pos + max_len,
                "sequence": frag,
                "match_len": max_len
            })

    # deduplicate and sort
    seen, dedup = set(), []
    for h in all_hits:
        key = (h["strand"], h["genome_start"], h["sequence"])
        if key not in seen:
            seen.add(key)
            dedup.append(h)
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
