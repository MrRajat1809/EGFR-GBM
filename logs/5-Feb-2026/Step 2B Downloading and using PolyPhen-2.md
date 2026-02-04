STEP 2B — PolyPhen-2 Functional Impact Prediction (Missense nsSNPs)
Objective

To independently assess the structural and functional impact of EGFR missense nsSNPs using PolyPhen-2, and integrate PolyPhen evidence onto the base missense variant dataset without cross-contaminating other tool outputs.

Input dataset

Base file:

EGFR_All_Missense_nsSNPs_GRCh37.csv


Variants: 2,578 EGFR missense nsSNPs

Genome build: GRCh37

Protein: EGFR canonical sequence (UniProt P00533)

Procedure
1. Batch input preparation

A custom Python script(make_polyphen_input.py) was used to generate PolyPhen-2 batch input from the base missense dataset.

Input format:

P00533 <Protein_position> <RefAA> <AltAA>


Ambiguous protein positions (e.g. ranges such as 668–669) were excluded.

Total batch input variants: 2,570

2. PolyPhen-2 execution

PolyPhen-2 Batch Query mode was used via the official web server.

Parameters:

Classifier: HumDiv

Genome assembly: GRCh37/hg19

Transcript: Canonical

Variant type: Missense

Output file obtained:

EGFR_PolyPhenprediction.txt

3. Output parsing and mapping

PolyPhen-2 output was provided in a fixed-width, column-aligned format.

A custom Python parser(map_polyphen_to_missense.py) was implemented to:

Extract submitted mutation identifiers (o_pos, o_aa1, o_aa2)

Retrieve PolyPhen-2 predictions and probabilities (prediction, pph2_prob)

Map results back onto the base missense dataset using:

(Protein_position, RefAA, AltAA)


Variants without PolyPhen-2 scores were retained and labeled as NA.

Results

Total missense nsSNPs: 2,578

Successfully mapped: 2,570

NA (ambiguous positions): 8

PolyPhen-2 classification (mapped variants)
Class	Count
Benign	1,059
Possibly damaging	387
Probably damaging	1,124
Output file
EGFR_Missense_nsSNPs_with_PolyPhen.csv


Contains original missense variant information plus:

PolyPhen_prediction

PolyPhen_score

Step 2B status

✅ PolyPhen-2 executed successfully
✅ Results mapped cleanly to base dataset
✅ Counts consistent with VEP missense summary
✅ Independent evidence layer ready for consensus analysis