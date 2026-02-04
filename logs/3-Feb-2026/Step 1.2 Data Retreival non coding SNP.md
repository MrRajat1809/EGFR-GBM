Objective

To retrieve and curate all single-nucleotide non-coding (regulatory) variants within the EGFR locus (GRCh37) for downstream regulatory impact analysis, in parallel with the coding nsSNP arm.

Data source and inputs

Primary source: dbSNP (GRCh37.p13)

Annotation tool: Ensembl Variant Effect Predictor (VEP v115.2, offline cache)

Input file:

EGFR_GRCh37_vep.vcf


(Previously generated in Step 1.1 from EGFR locus extraction)

EGFR genomic locus (GRCh37):

Contig: NC_000007.13

Coordinates: 55,086,724 – 55,279,321

Transcript policy: Canonical EGFR transcript only

Variant type policy: SNPs only (single-nucleotide substitutions)

Extraction strategy

Non-coding SNPs were extracted by re-filtering the existing VEP-annotated EGFR variant set, without downloading new data.

Inclusion criteria (VEP consequence terms)

Variants annotated with any of the following non-coding regulatory consequences were retained:

intron_variant

5_prime_UTR_variant

3_prime_UTR_variant

upstream_gene_variant

downstream_gene_variant

splice_region_variant

splice_polypyrimidine_tract_variant

non_coding_transcript_exon_variant

Exclusion criteria

Variants were excluded if they were:

Protein-altering (e.g. missense_variant, synonymous_variant, frameshift_variant, stop_gained)

Non-SNVs (insertions, deletions, indels)

Ambiguous or non-canonical transcript mappings

This ensured a pure regulatory SNP dataset, fully separated from the coding nsSNP arm.

script: /mnt/d/EGFR_pipeline/scripts/extract_noncoding_snps_egfr.py

Output file

Final file generated:

EGFR_NonCoding_SNPs_GRCh37.csv


Each row represents:

One single-nucleotide regulatory variant within the EGFR locus, mapped to the canonical transcript and annotated with a non-coding functional category.

Final counts and distribution

Total non-coding SNPs: 72,179

Breakdown by consequence type:

Intronic variants: 69,161

3′ UTR variants: 2124

5′ UTR variants: 167

Upstream gene variants: 61

Interpretation and validation

The dataset is intron-dominated, which is biologically expected given the large intronic span of EGFR.

UTR and upstream variants form a small but high-priority regulatory subset.

Counts are lower than raw VEP summary statistics because:

only SNPs (not indels) were retained

multi-allelic and redundant annotations were collapsed

canonical transcript filtering was enforced

The relative distribution of consequence classes matches the VEP summary, confirming that biological signal was preserved while technical noise was removed.

Step 1.2 status

✅ Completed successfully
✅ Consistent with Step 1.1 coding retrieval
✅ Benchmark-aligned regulatory dataset ready for tool-based prioritization