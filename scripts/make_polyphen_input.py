import csv

INPUT = "/mnt/d/EGFR_pipeline/reference/EGFR_All_Missense_nsSNPs_GRCh37.csv"
OUTPUT = "/mnt/d/EGFR_pipeline/software/polyphen/EGFR_PolyPhen_input.txt"

with open(INPUT) as f, open(OUTPUT, "w") as out:
    reader = csv.DictReader(f)
    for row in reader:
        pos = row["Protein_position"]
        aa_change = row["AA_change"]

        # Skip ambiguous ranges like 668-669
        if not pos.isdigit():
            continue

        try:
            ref_aa, alt_aa = aa_change.split("/")
        except ValueError:
            continue

        out.write(f"P00533 {pos} {ref_aa} {alt_aa}\n")
