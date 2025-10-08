# sRAT — small RNA Alignment Tool

Exact-match small RNA aligner with:
- Middle-50 window search
- Right-trim (keep left end) → 50..20
- Left-trim (keep right end) → 50..20
- Reverse-complement genome pass
- Mask 10+ homopolymers to N (ignore hits there)
- **Multimapping rule:** Determine the longest exact match length **L ≥ 20** for each read.  
  - Report **all** matches of that L-nt sequence across forward and reverse-complement genomes.  
  - If that L-nt sequence occurs **> 3** times, **discard the entire read** (no output rows for it).
- Single Excel output: `read_id, strand, genome_coords, read_coords, sequence, match_len`

## Install
```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
source .venv/bin/activate
pip install -r requirements.txt
