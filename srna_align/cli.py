import argparse, os, time
from tqdm import tqdm
from .io import load_genome_fasta, mask_homopolymers, revcomp, iter_reads
from .index import build_index, lookup as index_lookup
from .search import search_one_read
from .export import write_excel

class IndexWrapper:
    """Simple wrapper to provide .lookup(kmer) API over the dict index."""
    def __init__(self, idx_dict): self.idx = idx_dict
    def lookup(self, kmer: str):
        return index_lookup(self.idx, kmer)

def main():
    p = argparse.ArgumentParser(description="sRAT – small RNA exact-match aligner")
    p.add_argument("--genome", required=True, help="Genome FASTA (single sequence)")
    p.add_argument("--reads", required=True, help="Reads FASTQ (or FASTA)")
    p.add_argument("--out", default="out", help="Output folder")
    args = p.parse_args()

    t0 = time.time()
    os.makedirs(args.out, exist_ok=True)

    # Load and mask genome
    gname, genome_raw = load_genome_fasta(args.genome)
    genome_fwd = mask_homopolymers(genome_raw)
    genome_rev = mask_homopolymers(revcomp(genome_raw))

    # Build indexes for 20..50 on both strands
    idx_fwd = IndexWrapper(build_index(genome_fwd, 20, 50))
    idx_rev = IndexWrapper(build_index(genome_rev, 20, 50))

    # Process reads
    all_hits = []
    for rid, rseq in tqdm(iter_reads(args.reads), desc="Aligning reads"):
        hits = search_one_read(rid, rseq, idx_fwd, idx_rev, gname)
        all_hits.extend(hits)

    # Output Excel
    reads_base = os.path.splitext(os.path.basename(args.reads))[0]
    genome_base = os.path.splitext(os.path.basename(args.genome))[0]
    out_path = os.path.join(args.out, f"{reads_base}_{genome_base}.xlsx")
    write_excel(all_hits, out_path)

    elapsed = time.time() - t0
    print(f"Done. Wrote {out_path}. Elapsed: {elapsed:.2f}s")

if __name__ == "__main__":
    main()
