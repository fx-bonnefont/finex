"""Chemins et constantes partagés pour le pipeline finex."""

from pathlib import Path

# Répertoires
DATA_DIR = Path("data")
OUTPUT_DIR = Path("output")

# Fichiers récepteur
RECEPTOR_PDB = DATA_DIR / "receptor_xenoAMP.pdb"
RECEPTOR_PDBQT = DATA_DIR / "receptor_xenoAMP.pdbqt"

# Répertoires ligands
SDF_DROP_DIR = DATA_DIR / "sdf_drop"
LIGANDS_IN_DIR = DATA_DIR / "ligands_in"
LIGANDS_PDBQT_DIR = DATA_DIR / "ligands_pdbqt"

# Résultats docking
DOCKING_RESULTS_DIR = DATA_DIR / "docking_results"

# Sorties utilisateur
RESULTS_XLSX = OUTPUT_DIR / "results.xlsx"
VIZ_DIR = OUTPUT_DIR / "viz"

# Séquence cible
FASTA_SEQUENCE = "KSTNLVKNKCVNFNFNGLTGTGVLTESNKK"

# API OpenFold3
API_URL = "https://health.api.nvidia.com/v1/biology/openfold/openfold3/predict"
POLL_INTERVAL = 10
MAX_WAIT = 600

# Docking
DEFAULT_EXHAUSTIVENESS = 8

# Grid box (calibré sur résidus 546-550 du peptide xenoAMP)
GRID_CENTER_X = -1.0
GRID_CENTER_Y = 5.0
GRID_CENTER_Z = -2.5
GRID_SIZE_X = 20.0
GRID_SIZE_Y = 20.0
GRID_SIZE_Z = 20.0
