from Bio import SeqIO
from Bio.Seq import Seq
import re

# 10+ identical bases (A/C/G/T) → mask with 'N'
HOMOPOLYMER_RE = re.compile(r'([ACGT])\1{9,}')

def load_genome_fasta(path: str):
    """Load a single-sequence genome FASTA. Returns (name, sequence_upper_T)."""
    recs = list(SeqIO.parse(path, "fasta"))
    if len(recs) != 1:
        raise ValueError("Genome FASTA should contain exactly one sequence")
    name = recs[0].id
    seq = str(recs[0].seq).upper().replace('U','T')
    return name, seq

def mask_homopolymers(seq: str) -> str:
    """Replace any 10+ homopolymer stretch with 'N' of the same length."""
    def repl(m):
        return "N" * (m.end() - m.start())
    return HOMOPOLYMER_RE.sub(repl, seq)

def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def iter_reads(path: str):
    """Yield (read_id, sequence_upper_T) from FASTQ or FASTA."""
    fmt = "fastq" if path.lower().endswith(("fastq", "fq")) else "fasta"
    for rec in SeqIO.parse(path, fmt):
        yield rec.id, str(rec.seq).upper().replace('U','T')

def middle_window(seq: str, win: int = 50) -> tuple[int, int]:
    """Return (start, end) 0-based of the middle window of length 'win'.
       If seq shorter than win, return entire seq."""
    n = len(seq)
    if n < win:
        return 0, n
    start = (n - win) // 2
    return start, start + win
