# Locate shared paths when this script is executed directly.
from pathlib import Path as _Path
import sys as _sys
_SRC = next(p for p in _Path(__file__).resolve().parents if p.name == "src")
_sys.path.insert(0, str(_SRC / "common"))
import importlib as _importlib
_paths = _importlib.import_module("01_project_paths")
project_path = _paths.project_path
use_working_directory = _paths.use_working_directory
use_working_directory('outputs/noncoding')

import torch
from pyfaidx import Fasta
import os
import gc
from enformer_pytorch import Enformer, seq_indices_to_one_hot

print("====================================================")
print("🚀 INITIALIZING ENFORMER IN LOCAL DOCKER (CPU MODE)")
print("====================================================")

device = torch.device('cpu')

print("🏗️ Building Enformer architecture manually...")
model = Enformer.from_hparams(
    dim = 1536,
    depth = 11,
    heads = 8,
    output_heads = dict(human = 5313, mouse = 1643),
    target_length = 896,
    use_tf_gamma = True 
)
model.all_tied_weights_keys = []

local_weights = project_path('data/reference/noncoding/pytorch_model.bin')
if not os.path.exists(local_weights):
    print("❌ ERROR: pytorch_model.bin not found!")
    exit()

state_dict = torch.load(local_weights, map_location='cpu', weights_only=True)
model.load_state_dict(state_dict)
model = model.to(device)
model.eval()

# --- THE GEOGRAPHY FIX ---
CHROMOSOME = 'chr7'
# Centered exactly between EGFR TSS (55086714) and Enhancer (55191823)
MIDPOINT = 55139268  
SEQUENCE_LENGTH = 196608
FASTA_PATH = project_path('data/reference/noncoding/hg38.fa')

# 0-based start position for sequence extraction
start_pos = MIDPOINT - (SEQUENCE_LENGTH // 2) 
end_pos = MIDPOINT + (SEQUENCE_LENGTH // 2)

print(f"\n🧬 Extracting 196kb window around {CHROMOSOME}:{start_pos}-{end_pos}...")
# Both the original interval and pyfaidx slicing use zero-based, end-exclusive coordinates.
with Fasta(FASTA_PATH, as_raw=True, sequence_always_upper=True) as fasta_file:
    wt_sequence = fasta_file[CHROMOSOME][start_pos:end_pos]

# --- EXACT 0-BASED STRING INDICES ---
idx_509 = 55191816 - start_pos - 1
idx_524 = 55191831 - start_pos - 1

mut_509_list = list(wt_sequence)
mut_509_list[idx_509] = 'G'
mut_509_seq = "".join(mut_509_list)

mut_524_list = list(wt_sequence)
mut_524_list[idx_524] = 'A'
mut_524_seq = "".join(mut_524_list)

print(f"✅ WT Allele at 509: {wt_sequence[idx_509]} | Mutant: {mut_509_seq[idx_509]}")
print(f"✅ WT Allele at 524: {wt_sequence[idx_524]} | Mutant: {mut_524_seq[idx_524]}\n")

def str_to_one_hot(seq_str):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    indices = [mapping.get(nuc, 4) for nuc in seq_str]
    tensor_indices = torch.tensor(indices).unsqueeze(0).to(device)
    return seq_indices_to_one_hot(tensor_indices)

# --- PROMOTER TARGETING ---
TSS_BIN_START = 36
TSS_BIN_END = 39

print("🧠 Simulating Wild-Type...")
with torch.no_grad():
    wt_out = model(str_to_one_hot(wt_sequence))['human'][0]
    # Extract the max expression across the promoter bins for ALL 5313 tracks
    wt_promoter = torch.max(wt_out[TSS_BIN_START:TSS_BIN_END, :], dim=0)[0] # Shape: (5313,)
    del wt_out
    gc.collect()
print("   -> Wild-Type tensor mapped and cleared from RAM.")

print("\n🧠 Simulating Mutant 55259509 (SOX2 Gain)...")
with torch.no_grad():
    mut_509_out = model(str_to_one_hot(mut_509_seq))['human'][0]
    mut_509_promoter = torch.max(mut_509_out[TSS_BIN_START:TSS_BIN_END, :], dim=0)[0]
    del mut_509_out
    gc.collect()
print("   -> Mutant 509 tensor mapped and cleared from RAM.")

print("\n🧠 Simulating Mutant 55259524 (OLIG2 Gain)...")
with torch.no_grad():
    mut_524_out = model(str_to_one_hot(mut_524_seq))['human'][0]
    mut_524_promoter = torch.max(mut_524_out[TSS_BIN_START:TSS_BIN_END, :], dim=0)[0]
    del mut_524_out
    gc.collect()
print("   -> Mutant 524 tensor mapped and cleared from RAM.")

# --- UNBIASED GLOBAL TRACK SCAN ---
print("\n🌍 Scanning all 5,313 human epigenetic tracks for EGFR activation...")

# Calculate Fold Changes across all tracks (Adding 1e-4 to prevent division by zero)
fc_509 = mut_509_promoter / (wt_promoter + 1e-4)
fc_524 = mut_524_promoter / (wt_promoter + 1e-4)

# Filter out noise: Only look at tracks where the mutant actually drives real expression (> 0.5)
valid_tracks_509 = torch.where(mut_509_promoter > 0.5)[0]
valid_tracks_524 = torch.where(mut_524_promoter > 0.5)[0]

if len(valid_tracks_509) > 0:
    best_track_509 = valid_tracks_509[torch.argmax(fc_509[valid_tracks_509])].item()
else:
    best_track_509 = torch.argmax(fc_509).item()

if len(valid_tracks_524) > 0:
    best_track_524 = valid_tracks_524[torch.argmax(fc_524[valid_tracks_524])].item()
else:
    best_track_524 = torch.argmax(fc_524).item()

# Get the final numbers for the best tracks
wt_baseline_for_509 = wt_promoter[best_track_509].item()
best_mut_509_expr = mut_509_promoter[best_track_509].item()
best_fc_509 = best_mut_509_expr / (wt_baseline_for_509 + 1e-4)

wt_baseline_for_524 = wt_promoter[best_track_524].item()
best_mut_524_expr = mut_524_promoter[best_track_524].item()
best_fc_524 = best_mut_524_expr / (wt_baseline_for_524 + 1e-4)

print("\n====================================================")
print("📊 UNBIASED GLOBAL IN SILICO PREDICTION")
print("====================================================")
print(f"🔥 MUTANT 509 (SOX2 Gain) - Biggest impact found on Track {best_track_509}")
print(f"   Wild-Type Baseline:  {wt_baseline_for_509:.4f}")
print(f"   Mutant Expression:   {best_mut_509_expr:.4f} --> {best_fc_509:.2f}x Fold Change")
print("----------------------------------------------------")
print(f"🔥 MUTANT 524 (OLIG2 Gain) - Biggest impact found on Track {best_track_524}")
print(f"   Wild-Type Baseline:  {wt_baseline_for_524:.4f}")
print(f"   Mutant Expression:   {best_mut_524_expr:.4f} --> {best_fc_524:.2f}x Fold Change")
print("====================================================")
