"""Génération des fichiers PDB de visualisation (récepteur + meilleure pose)."""

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def prepare_viz_for_ligand(
    ligand_name: str,
    results_dir: Path,
    receptor_pdb: Path,
    viz_dir: Path,
) -> Path:
    """Extrait la meilleure pose et crée un PDB complexe récepteur+ligand."""
    viz_dir.mkdir(parents=True, exist_ok=True)

    input_pdbqt = results_dir / f"{ligand_name}_out.pdbqt"
    if not input_pdbqt.exists():
        raise FileNotFoundError(f"Résultat de docking introuvable : {input_pdbqt}")

    ligand_pdb = viz_dir / f"{ligand_name}_pose1.pdb"
    complex_pdb = viz_dir / f"{ligand_name}_complex.pdb"

    # Extraire la meilleure pose avec obabel
    cmd = ["obabel", "-ipdbqt", str(input_pdbqt), "-opdb", "-O", str(ligand_pdb), "-f", "1", "-l", "1"]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    if result.returncode != 0:
        raise RuntimeError(f"Erreur obabel pour {ligand_name} : {result.stderr}")

    # Fusionner récepteur + ligand
    with open(complex_pdb, "w") as out_f:
        out_f.write(receptor_pdb.read_text())
        out_f.write(ligand_pdb.read_text())

    logger.info("Visualisation créée : %s", complex_pdb)
    return complex_pdb


def prepare_viz_top_n(
    scores: list[tuple[str, float]],
    top_n: int,
    results_dir: Path,
    receptor_pdb: Path,
    viz_dir: Path,
) -> list[Path]:
    """Génère les fichiers de visualisation pour les top N molécules."""
    top_ligands = scores[:top_n]
    paths = []

    for ligand_name, score in top_ligands:
        logger.info("Préparation visualisation : %s (score : %.1f kcal/mol)", ligand_name, score)
        path = prepare_viz_for_ligand(ligand_name, results_dir, receptor_pdb, viz_dir)
        paths.append(path)

    return paths
