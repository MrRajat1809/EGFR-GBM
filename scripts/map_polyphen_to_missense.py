import csv

# -------- FILES --------
BASE_FILE = "/mnt/d/EGFR_pipeline/reference/EGFR_All_Missense_nsSNPs_GRCh37.csv"
POLYPHEN_FILE = "/mnt/d/EGFR_pipeline/software/polyphen/EGFR_PolyPhenprediction.txt"
OUTPUT_FILE = "/mnt/d/EGFR_pipeline/results/EGFR_Missense_nsSNPs_with_PolyPhen.csv"

# -------- FIXED-WIDTH PARSER (robust for "possibly damaging") --------
def build_fixed_width_schema(header_line: str):
    """
    Build (name, start, end) slices for fixed-width columns based on the header line.
    """
    tokens = []
    i = 0
    n = len(header_line)

    while i < n:
        if header_line[i].isspace():
            i += 1
            continue
        start = i
        while i < n and not header_line[i].isspace():
            i += 1
        name = header_line[start:i]
        tokens.append((name, start))

    # Convert starts to (start,end)
    schema = []
    for idx, (name, start) in enumerate(tokens):
        end = tokens[idx + 1][1] if idx + 1 < len(tokens) else n
        schema.append((name, start, end))
    return schema

def parse_fixed_width_line(line: str, schema):
    row = {}
    for name, start, end in schema:
        row[name] = line[start:end].strip()
    return row

# -------- 1) READ POLYPHEN OUTPUT INTO LOOKUP DICT --------
# Key on (protein_position, refAA, altAA) using the submitted fields: o_pos, o_aa1, o_aa2
polyphen_map = {}
# ... (Files and header cleaning) ...

polyphen_map = {}

with open(POLYPHEN_FILE, "r", encoding="utf-8", errors="replace") as f:
    # 1. Find and clean the header line
    header_line = ""
    for line in f:
        if line.startswith("#o_acc"):
            header_line = line.lstrip("#").strip()
            break
    
    if not header_line:
        raise RuntimeError("Could not find PolyPhen header line starting with #o_acc")

    # 2. Use DictReader with tab delimiter
    # We split the header by tabs to get the exact column names
    fieldnames = [col.strip() for col in header_line.split('\t') if col.strip()]
    
    # Note: PolyPhen files sometimes have slightly inconsistent tab counts in headers.
    # If the standard DictReader fails, we use a simple split approach:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
            
        # Split by tab and strip whitespace from each value
        values = [val.strip() for val in line.split('\t')]
        
        # Create a dictionary mapping header names to values
        d = dict(zip(fieldnames, values))

        pos = d.get("o_pos", "")
        ref = d.get("o_aa1", "")
        alt = d.get("o_aa2", "")

        if pos and ref and alt:
            key = (pos, ref, alt)
            polyphen_map[key] = {
                "prediction": d.get("prediction", "NA"),
                "pph2_prob": d.get("pph2_prob", "NA"),
            }

# -------- 2) JOIN ONTO BASE CSV --------
with open(BASE_FILE, "r", newline="") as inp, open(OUTPUT_FILE, "w", newline="") as out:
    reader = csv.DictReader(inp)

    new_cols = ["PolyPhen_prediction", "PolyPhen_score"]
    fieldnames = reader.fieldnames + new_cols

    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        pos_str = row.get("Protein_position", "")
        aa_change = row.get("AA_change", "")

        # Default NA
        row["PolyPhen_prediction"] = "NA"
        row["PolyPhen_score"] = "NA"

        # Skip ambiguous ranges like "668-669"
        if not pos_str.isdigit():
            writer.writerow(row)
            continue

        try:
            ref_aa, alt_aa = aa_change.split("/")
        except ValueError:
            writer.writerow(row)
            continue

        key = (pos_str, ref_aa, alt_aa)
        hit = polyphen_map.get(key)

        if hit:
            pred = hit.get("prediction", "").strip()
            score = hit.get("pph2_prob", "").strip()

            row["PolyPhen_prediction"] = pred if pred else "NA"
            row["PolyPhen_score"] = score if score else "NA"

        writer.writerow(row)

print(f"Done. Wrote: {OUTPUT_FILE}")
