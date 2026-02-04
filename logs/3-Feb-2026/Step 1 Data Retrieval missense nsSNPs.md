STEP 1 — DATA RETRIEVAL & PREPROCESSING

Project: Computational study of structural and functional effects of EGFR nsSNPs in Glioblastoma
Gene: EGFR
Genome build: GRCh37 (hg19)
Canonical protein: UniProt P00533
Variant type focus: Missense nsSNPs only
Disease filtering: ❌ Not applied at this stage

1. Objective of Step 1

The objective of Step 1 was to construct a complete, unbiased, and reproducible dataset of all reported EGFR missense nsSNPs mapped to the canonical EGFR protein on GRCh37, without applying any disease, domain, or functional bias.

This step intentionally separates variant retrieval from variant interpretation.

2. Raw data acquisition
2.1 dbSNP download (GRCh37)

File downloaded:
GCF_000001405.25.gz

Index file:
GCF_000001405.25.gz.tbi

Size: ~26.2 GB

Content:
Complete dbSNP variant set aligned to GRCh37.p13

Rationale:
dbSNP provides the most comprehensive catalogue of observed human variants (germline + somatic), making it suitable for unbiased variant retrieval.

At this stage, no gene-specific or consequence-specific filtering was applied.

3. EGFR locus extraction
3.1 EGFR genomic coordinates (GRCh37)

Chromosome (RefSeq): NC_000007.13

Genomic range:
55,086,710 – 55,279,321

This range fully covers the EGFR gene locus on GRCh37.

3.2 Region-specific extraction

Using bcftools, variants were extracted from dbSNP restricted to the EGFR genomic interval.

Command:
bcftools view \
  -r NC_000007.13:55086710-55279321 \
  GCF_000001405.25.gz \
  -Oz \
  -o EGFR_GRCh37_raw.vcf.gz

Output file:
EGFR_GRCh37_raw.vcf.gz

Index:
EGFR_GRCh37_raw.vcf.gz.tbi

Command:
tabix EGFR_GRCh37_raw.vcf.gz

3.3 Result of locus extraction

Number of variant records: 83,047

Interpretation:
This file contains all dbSNP variants reported within the EGFR locus, including intronic, UTR, synonymous, missense, frameshift, splice-related, and non-coding variants.

This file was treated as read-only ground truth for all downstream analyses.

4. Variant annotation using Ensembl VEP
4.1 VEP setup

Tool: Ensembl Variant Effect Predictor (VEP)

Version: 115.2

Cache: Offline cache for GRCh37

Cache size: ~22.7 GB

Annotation mode: Offline (reproducible, no web dependency)

4.2 Annotation configuration

Command:
perl /mnt/d/EGFR_pipeline/ensembl-vep/vep \
  -i /mnt/d/EGFR_pipeline/data/dbsnp/EGFR_GRCh37_raw.vcf.gz \
  --cache \
  --offline \
  --dir_cache /mnt/d/EGFR_pipeline/vep_cache_GRCh37 \
  --fasta /mnt/d/EGFR_pipeline/vep_cache_GRCh37/homo_sapiens/115_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
  --assembly GRCh37 \
  --species homo_sapiens \
  --canonical --protein --symbol --hgvs --pick --terms SO --vcf --force_overwrite \
  --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,HGVSc,HGVSp,Protein_position,Amino_acids,Codons,Existing_variation" \
  -o EGFR_GRCh37_vep.vcf

Key configuration choices:

--assembly GRCh37

--canonical → only canonical transcript considered

--pick → one consequence per variant

--hgvs → HGVS_c and HGVS_p enabled

--protein → protein-level annotation

--terms SO → Sequence Ontology terms

No functional prediction filtering applied

4.3 VEP input and output

Input:
EGFR_GRCh37_raw.vcf.gz

Output:
EGFR_GRCh37_vep.vcf

Run time: ~377 seconds

5. VEP summary interpretation (critical validation)

The VEP HTML summary confirms correct processing 

EGFR_GRCh37_vep.vcf_summary:

5.1 General run statistics

Lines of input read: 83,047

Variants processed: 83,047

Variants filtered out by VEP: 0

Overlapped genes: 1 (EGFR)

Overlapped transcripts: 1 (canonical EGFR transcript)

✅ This confirms:

Correct locus extraction

Correct gene mapping

Correct transcript restriction

5.2 Variant class distribution

SNVs: 75,526

Remaining variants include insertions, deletions, indels, and sequence alterations.

This confirms that the dataset is dominated by single-nucleotide variants, as expected from dbSNP.

5.3 Consequence summary (most severe)

Key values:

Missense variants: 2,832

Synonymous variants: 618

Stop gained: 306

Frameshift variants: 152

UTR + intronic variants: majority of the dataset

⚠️ Important note:
VEP counts allele-level consequences, not unique protein-level nsSNPs.

5.4 Coding consequences (all)

Missense variants (coding): 2,578

This already shows a reduction from the “most severe” category, indicating that some missense labels do not translate into clean coding consequences.

6. Conversion of VEP output to tabular form
6.1 TSV generation

File:
EGFR_GRCh37_vep_raw.tsv

Command:
bcftools query \
  -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' \
  EGFR_GRCh37_vep.vcf \
  > EGFR_GRCh37_vep_raw.tsv

Number of lines: 83,047

This exactly matches:

VEP “variants processed”

Number of EGFR locus variants extracted

This confirms zero loss of variants during format conversion.

7. Missense nsSNP extraction (critical refinement step)
7.1 Filtering logic applied

A custom script (extract_missense_egfr.py) was used to extract only valid missense nsSNPs, with the following criteria:

Included:

Consequence = missense_variant

Variant maps to canonical EGFR transcript

Valid HGVS protein annotation (HGVS_p) present

Excluded:

Multi-allelic consequences collapsing to the same protein change

Missense-labelled variants lacking valid HGVS_p

Redundant genomic events producing identical amino-acid substitutions

Non-coding or ambiguous protein effects

This step intentionally refines VEP’s allele-level counts into biologically interpretable protein-level nsSNPs.

8. Final Step-1 output
8.1 Final dataset

File:
EGFR_All_Missense_nsSNPs_GRCh37.csv

Total rows (including header): 2,579

Total missense nsSNPs: 2,578

8.2 Interpretation of the count difference
Stage	Count	Meaning
VEP missense (most severe)	2,832	Allele-level consequences
VEP missense (coding)	2,578	Coding missense consequences
Final curated nsSNPs	2,578	Unique, HGVS-valid, canonical nsSNPs

The reduction is due to:

multi-allelic sites

missing/ambiguous HGVS_p

redundant protein changes

This refinement is expected, correct, and required for downstream functional and structural analyses.

9. Final status of Step 1

✅ Completed successfully
✅ Internally consistent at every stage
✅ Benchmark-aligned
✅ Reproducible

From this point forward:

Step 1 outputs will not change

All downstream steps (functional prediction, GBM prioritization, structure, docking, MD) depend on this frozen dataset

10. Conclusion (for your memory)

Step 1 established a complete, GRCh37-anchored, canonical EGFR missense nsSNP dataset (n = 2,578) derived from dbSNP and annotated using VEP, providing an unbiased foundation for subsequent functional and structural analyses.