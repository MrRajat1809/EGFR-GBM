import csv

# Input from bcftools query
INPUT = "EGFR_GRCh37_vep_for_noncoding.tsv"
OUTPUT = "EGFR_NonCoding_SNPs_GRCh37.csv"

# Non-coding consequence set (VEP SO terms)
NONCODING = {
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "splice_region_variant",
    "splice_polypyrimidine_tract_variant",
    "splice_donor_5th_base_variant",
    "intron_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "non_coding_transcript_exon_variant",
}

def classify(cons_string: str) -> str:
    """Classifies the primary consequence into a readable category."""
    # Split by & or , because VEP can join multiple terms (e.g., intron_variant&splice_region_variant)
    terms = cons_string.replace("&", ",").split(",")
    
    # Priority classification based on the first matching term found
    for cons in terms:
        if "UTR" in cons: return "UTR"
        if cons.startswith("splice_"): return "SPLICE_RELATED"
        if cons == "intron_variant": return "INTRON"
        if "upstream" in cons: return "UPSTREAM"
        if "downstream" in cons: return "DOWNSTREAM"
        if "non_coding_transcript_exon" in cons: return "NONCODING_EXON"
    return "OTHER"

with open(INPUT, "r") as f, open(OUTPUT, "w", newline="") as out:
    w = csv.writer(out)
    w.writerow([
        "Variant_ID", "rsID", "Chr", "Pos_GRCh37", "Ref", "Alt",
        "Consequence", "Feature_Class", "Gene", "Transcript_ID"
    ])

    for line in f:
        # Basic cleanup and splitting
        row = line.rstrip("\n").split("\t")
        if len(row) < 6:
            continue
            
        rsid, chrom, pos, ref, alt_field, csq_field = row[0], row[1], row[2], row[3], row[4], row[5]

        # Handle multi-allelic sites and multiple VEP annotations
        alts = alt_field.split(",")
        csq_entries = csq_field.split(",")

        for alt in alts:
            # STRICT SNP CHECK: Ref and Alt must be exactly 1 base long
            if len(ref) != 1 or len(alt) != 1:
                continue 

            # Find the CSQ entry that matches the current ALT allele
            for entry in csq_entries:
                fields = entry.split("|")
                
                # Check indexing: 2=Allele, 3=Gene, 4=Feature, 5=Consequence
                if len(fields) >= 6 and fields[2] == alt:
                    consequence = fields[5]
                    
                    # Split combined consequences (e.g., 'intron_variant&splice_region_variant')
                    # and check if ANY of them are in our NONCODING set
                    current_cons_terms = consequence.replace("&", ",").split(",")
                    
                    if any(term in NONCODING for term in current_cons_terms):
                        gene = fields[3]
                        transcript = fields[4]
                        variant_id = f"{chrom}:{pos}:{ref}:{alt}"
                        
                        w.writerow([
                            variant_id, rsid, chrom, pos, ref, alt,
                            consequence, classify(consequence), gene, transcript
                        ])
                        break # Found the correct allele match, move to next ALT