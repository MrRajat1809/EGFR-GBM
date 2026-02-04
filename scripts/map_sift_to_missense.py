import csv

# -------- FILES --------
MISSENSE_FILE = "/mnt/d/EGFR_pipeline/reference/EGFR_All_Missense_nsSNPs_GRCh37.csv"
SIFT_FILE = "/mnt/d/EGFR_pipeline/software/sift4g/EGFR_SIFTprediction.txt"
OUTPUT_FILE = "/mnt/d/EGFR_pipeline/results/EGFR_Missense_nsSNPs_with_SIFT.csv"

# -------- SIFT MATRIX PARSING --------
# Order of amino acids in SIFT matrix header
AA_ORDER = [
    "A","B","C","D","E","F","G","H","I","K","L","M",
    "N","P","Q","R","S","T","V","W","X","Y","Z","*","-"
]

# position (1-based) -> {AA: score}
sift_matrix = {}

with open(SIFT_FILE) as f:
    position = 0
    for line in f:
        line = line.strip()

        # Skip headers / metadata lines
        if not line or line[0].isalpha():
            continue

        scores = line.split()
        if len(scores) != len(AA_ORDER):
            continue

        position += 1
        sift_matrix[position] = {
            aa: float(score) for aa, score in zip(AA_ORDER, scores)
        }

# -------- MAP TO MISSENSE SNPs --------
with open(MISSENSE_FILE) as inp, open(OUTPUT_FILE, "w", newline="") as out:
    reader = csv.DictReader(inp)
    fieldnames = reader.fieldnames + ["SIFT_score", "SIFT_prediction"]
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        pos_str = row["Protein_position"]

        # Handle ranges like "668-669" or any non-integer
        if not pos_str.isdigit():
            row["SIFT_score"] = "NA"
            row["SIFT_prediction"] = "NA"
            writer.writerow(row)
            continue

        pos = int(pos_str)

        aa_change = row["AA_change"]

        # Expected format: L/R
        try:
            ref_aa, alt_aa = aa_change.split("/")
        except ValueError:
            row["SIFT_score"] = "NA"
            row["SIFT_prediction"] = "NA"
            writer.writerow(row)
            continue

        score = None
        if pos in sift_matrix and alt_aa in sift_matrix[pos]:
            score = sift_matrix[pos][alt_aa]

        if score is None:
            row["SIFT_score"] = "NA"
            row["SIFT_prediction"] = "NA"
        else:
            row["SIFT_score"] = f"{score:.4f}"
            row["SIFT_prediction"] = (
                "Deleterious" if score <= 0.05 else "Tolerated"
            )

        writer.writerow(row)
