# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
use_working_directory('.')

import pandas as pd
import re
import os

# Paths
INPUT_CSV = project_path('outputs/variant_prioritization/EGFR_TCGA_Glioblastoma_Master_Consensus.csv')
OUTPUT_FASTA = project_path('outputs/structural_analysis/phenotypic_effect/EGFR_42_MutPred2_Input.fasta')

# 3-to-1 letter mapping
aa_map = {
    'Ala':'A', 'Arg':'R', 'Asn':'N', 'Asp':'D', 'Cys':'C',
    'Gln':'Q', 'Glu':'E', 'Gly':'G', 'His':'H', 'Ile':'I',
    'Leu':'L', 'Lys':'K', 'Met':'M', 'Phe':'F', 'Pro':'P',
    'Ser':'S', 'Thr':'T', 'Trp':'W', 'Tyr':'Y', 'Val':'V'
}

# Full Human EGFR Sequence (UniProt P00533)
EGFR_SEQ = (
    "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVER"
    "IPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPN"
    "GSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKK"
    "CPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQEL"
    "DILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKI"
    "ISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDG"
    "PHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRHIVRKRTLRRLLQER"
    "ELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTS"
    "TVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMA"
    "LESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLV"
    "IQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALT"
    "EDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDF"
    "FPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA"
)

def main():
    if not os.path.exists(INPUT_CSV):
        print(f"Error: Could not find {INPUT_CSV}")
        return

    df = pd.read_csv(INPUT_CSV)
    drivers = df[df['Final_Status'] == 'High-Confidence Driver'].copy()
    
    mut_list = []
    for _, row in drivers.iterrows():
        match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', str(row['HGVSp']))
        if match:
            wt_3, pos, mut_3 = match.groups()
            if wt_3 in aa_map and mut_3 in aa_map:
                mut_list.append(f"{aa_map[wt_3]}{pos}{aa_map[mut_3]}")

    # Create the single header line with mutations separated by spaces
    # Format: >ProteinID Mut1 Mut2 Mut3 ...
    header = f">EGFR_P00533 {' '.join(mut_list)}"

    with open(OUTPUT_FASTA, 'w') as f:
        f.write(header + "\n")
        f.write(EGFR_SEQ + "\n")

    print(f"\nSUCCESS! MutPred2 input file generated.")
    print(f"Header preview: {header[:100]}...")
    print(f"File saved to: {OUTPUT_FASTA}")

if __name__ == "__main__":
    main()