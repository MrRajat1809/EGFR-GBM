STEP 2A — Functional Impact Prediction Using SIFT (Missense nsSNPs)
Objective

To evaluate the functional impact of EGFR missense nsSNPs using SIFT, an evolutionary conservation–based predictor, and integrate SIFT scores into the curated missense variant dataset.

Software and resources

Tool: SIFT4G

Protein: EGFR canonical sequence (UniProt P00533)

Genome build context: GRCh37

Variant set: 2,578 EGFR missense nsSNPs (validated against VEP summary)

Procedure
1. Installation of SIFT4G

SIFT4G was installed locally by cloning the official GitHub repository and compiling the source using make.

cd /mnt/d/EGFR_pipeline/software
git clone --recursive https://github.com/rvaser/sift4g.git
cd sift4g
make

Successful installation was confirmed by the presence of the executable bin/sift4g.

2. Preparation of protein inputs

Canonical EGFR protein sequence was downloaded from UniProt:

P00533.fasta

cd /mnt/d/EGFR_pipeline/data/uniprot
wget https://rest.uniprot.org/uniprotkb/P00533.fasta

Protein database for conservation analysis was prepared using:

cd /mnt/d/EGFR_pipeline/data/swissprot
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

UniProt Swiss-Prot FASTA (uniprot_sprot.fasta)

A full GRCh37 reference FASTA (~3 GB) was also downloaded earlier to ensure consistency with upstream genome-based annotations.

3. Execution of SIFT4G

SIFT4G was run in protein-centric mode, using:

 cd /mnt/d/EGFR_pipeline/software/sift4g

./bin/sift4g \
  -q /mnt/d/EGFR_pipeline/data/uniprot/P00533.fasta \
  -d /mnt/d/EGFR_pipeline/data/swissprot/uniprot_sprot.fasta \
  > EGFR_sift4g_raw.out

P00533.fasta as the query

uniprot_sprot.fasta as the reference database

Output was generated as a position-wise amino-acid tolerance matrix (.SIFTprediction file), covering 1,212 EGFR protein positions.

4. Post-processing and variant mapping

A custom Python script(map_sift_to_missense.py) was written to:

Parse the SIFT4G amino-acid substitution matrix

Map SIFT scores to each missense nsSNP using protein position and alternate amino acid

Assign functional labels using the standard threshold:

Score ≤ 0.05 → Deleterious

Score > 0.05 → Tolerated

Ambiguous protein position ranges (e.g. 668–669) were conservatively assigned NA.

Results

Total missense nsSNPs annotated: 2,578

SIFT classification:

Deleterious: 1,128

Tolerated: 1,442

NA (ambiguous positions): 8

Final output file:

EGFR_Missense_nsSNPs_with_SIFT.csv

Summary:
• EGFR locus (GRCh37) corrected to chr7:55,086,710–55,279,321 to ensure full 3′ UTR coverage.
• Re-ran VEP annotation and consequence extraction after locus correction.
• Final missense nsSNP count: 2,578 (matches VEP “Consequences (all)” summary).
• SIFT4G run on canonical EGFR protein (UniProt P00533).
• SIFT results:
    – Deleterious: 1,128
    – Tolerated: 1,442
    – NA (ambiguous protein positions): 8
• NA values retained to avoid incorrect positional mapping.
• Final output: EGFR_Missense_nsSNPs_with_SIFT.csv
