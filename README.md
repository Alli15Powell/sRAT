# sRAT
### Multimapping logic
- Determine the **longest perfect match length L** (≥ 20 nt) found for the read.
- Collect **all** alignments of that exact L-nt sequence across both forward and reverse-complement genomes.
- If this L-nt sequence occurs **> 3 times**, the read is considered a *highly repetitive multimapper* and is **excluded entirely** from the results.
- Otherwise, all of its L-nt matches (1–3 total) are written to the Excel output.
