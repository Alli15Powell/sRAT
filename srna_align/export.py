import pandas as pd

def write_excel(hits, out_path: str):
    """
    Write required columns to a single Excel:
      read_id, strand, genome_coords(1-based), read_coords(1-based), sequence, match_len
    If hits is empty, still create a sheet with headers.
    """
    cols = ["read_id","strand","genome_coords","read_coords","sequence","match_len", "reference_name"]
    if not hits:
        pd.DataFrame(columns=cols).to_excel(out_path, index=False)
        return
    
    rows = []
    for h in hits:
        g_start_1 = h["genome_start"] + 1
        g_end_1   = h["genome_end"]     # end already exclusive in 0-based → 1-based end is same value
        # read_start/read_end may be None if we re-enumerated longest hits globally
        if h.get("read_start") is None or h.get("read_end") is None:
            read_coords = ""  # not strictly required for spec after global re-enumeration
        else:
            r_start_1 = h["read_start"] + 1
            r_end_1   = h["read_end"]
            read_coords = f"{r_start_1}-{r_end_1}"

        rows.append({
            "read_id": h["read_id"],
            "strand": h["strand"],
            "genome_coords": f"{g_start_1}-{g_end_1}",
            "read_coords":   read_coords,
            "sequence": h["sequence"],
            "match_len": h["match_len"],
            "reference_name": h["genome"],
        })
    pd.DataFrame(rows, columns=cols).to_excel(out_path, index=False)
