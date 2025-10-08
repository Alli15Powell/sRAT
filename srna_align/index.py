from collections import defaultdict
from typing import Dict, List

# Index type:
# idx[k][prefix8][kmer] -> list of positions (0-based, inclusive start)
IndexType = Dict[int, Dict[str, Dict[str, List[int]]]]

def build_index(genome: str, min_k: int = 20, max_k: int = 50, prefix: int = 8) -> IndexType:
    """
    Build a multi-k exact-match index over the masked genome (A/C/G/T/N).
    Any k-mer containing 'N' is skipped.
    """
    idx: IndexType = {k: defaultdict(lambda: defaultdict(list)) for k in range(min_k, max_k + 1)}
    G = len(genome)
    for k in range(min_k, max_k + 1):
        bucket_k = idx[k]
        for i in range(0, G - k + 1):
            kmer = genome[i:i+k]
            if 'N' in kmer:
                continue
            p = kmer[:prefix]
            bucket_k[p][kmer].append(i)
    return idx

def lookup(idx: IndexType, subseq: str) -> List[int]:
    """Lookup positions (0-based starts) for 'subseq' in the index; empty if none."""
    k = len(subseq)
    table_k = idx.get(k)
    if table_k is None:
        return []
    bucket = table_k.get(subseq[:8], {})
    return bucket.get(subseq, [])
