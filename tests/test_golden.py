import os
import pandas as pd
from srna_align.export import write_excel
from srna_align.search import search_one_read

def test_export_creates_excel(tmp_path):
    hits = [{
        "read_id": "r1",
        "strand": "+",
        "genome": "g",
        "genome_start": 100,
        "genome_end": 120,
        "read_start": 0,
        "read_end": 20,
        "sequence": "A"*20,
        "match_len": 20
    }]
    out = tmp_path / "test.xlsx"
    write_excel(hits, str(out))
    assert out.exists()
    df = pd.read_excel(out)
    assert list(df.columns) == ["read_id","strand","genome_coords","read_coords","sequence","match_len"]
    assert df.iloc[0]["genome_coords"] == "101-120"
    assert df.iloc[0]["read_coords"] == "1-20"

def test_multimap_discard_rule():
    # Fake index that returns specific positions for a given kmer; we bypass indexes by stubbing .lookup
    class FakeIdx:
        def __init__(self, positions): self.positions = positions
        def lookup(self, kmer): return list(self.positions)

    # read: any sequence; we just need search_one_read to call .lookup
    rid = "rX"
    rseq = "A"*60

    # Case 1: 4 sites for the longest match (≥20) -> discard entire read
    idx4 = FakeIdx([1,2,3,4])
    hits = search_one_read(rid, rseq, idx4, idx4, "g")
    assert hits == [], "Read with >3 identical longest-length hits must be discarded"

    # Case 2: 2 sites -> include both
    idx2 = FakeIdx([10,20])
    hits2 = search_one_read(rid, rseq, idx2, idx2, "g")
    assert len(hits2) == 2
    assert all(h["match_len"] >= 20 for h in hits2)
