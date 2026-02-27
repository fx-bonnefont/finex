"""Criblage virtuel avec AutoDock Vina."""

import logging
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from pathlib import Path

from tqdm import tqdm

from .config import (
    GRID_CENTER_X, GRID_CENTER_Y, GRID_CENTER_Z,
    GRID_SIZE_X, GRID_SIZE_Y, GRID_SIZE_Z,
)

logger = logging.getLogger(__name__)
VINA_EXECUTABLE = "vina"
DEFAULT_MAX_TORSIONS = 20


def _count_torsions(pdbqt_path: Path) -> int:
    """Lit le TORSDOF d'un fichier PDBQT."""
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("TORSDOF"):
            return int(line.split()[1])
    return 0


def _is_valid_result(pdbqt_path: Path) -> bool:
    """Vérifie qu'un fichier _out.pdbqt est complet (non tronqué par Ctrl+C)."""
    try:
        content = pdbqt_path.read_text()
        if not content or len(content) < 100:
            return False
        has_score = "REMARK VINA RESULT" in content
        has_end = "ENDMDL" in content
        return has_score and has_end
    except Exception:
        return False


def check_vina() -> None:
    """Vérifie que Vina est dans le PATH."""
    try:
        result = subprocess.run(
            [VINA_EXECUTABLE, "--version"], capture_output=True, text=True, timeout=10
        )
        logger.info("AutoDock Vina détecté : %s", result.stdout.strip() or result.stderr.strip())
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Exécutable '{VINA_EXECUTABLE}' introuvable. Installez AutoDock Vina :\n"
            "  https://vina.scripps.edu/downloads/\n"
            "  Conda : conda install -c conda-forge autodock-vina"
        )


def _dock_one_ligand(ligand_path_str: str, receptor_pdbqt: str, results_dir: str, exhaustiveness: int) -> dict:
    """Worker : docke un ligand contre le récepteur."""
    ligand_path = Path(ligand_path_str)
    stem = ligand_path.stem
    out_pdbqt = Path(results_dir) / f"{stem}_out.pdbqt"

    cmd = [
        VINA_EXECUTABLE,
        "--receptor", receptor_pdbqt,
        "--ligand", str(ligand_path),
        "--out", str(out_pdbqt),
        "--center_x", str(GRID_CENTER_X),
        "--center_y", str(GRID_CENTER_Y),
        "--center_z", str(GRID_CENTER_Z),
        "--size_x", str(GRID_SIZE_X),
        "--size_y", str(GRID_SIZE_Y),
        "--size_z", str(GRID_SIZE_Z),
        "--exhaustiveness", str(exhaustiveness),
    ]

    # Timeout proportionnel à l'exhaustiveness (base 300s pour e=8, linéaire)
    timeout = max(600, int(300 * exhaustiveness / 8))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            return {"ligand": stem, "status": "error", "msg": result.stderr.strip()[:500]}
        return {"ligand": stem, "status": "ok"}
    except subprocess.TimeoutExpired:
        return {"ligand": stem, "status": "timeout"}
    except Exception as e:
        return {"ligand": stem, "status": "error", "msg": str(e)}


def dock_all_ligands(
    receptor_pdbqt: Path,
    ligands_dir: Path,
    results_dir: Path,
    exhaustiveness: int,
    skip_existing: bool = True,
    only_ligands: set[str] | None = None,
    max_torsions: int = DEFAULT_MAX_TORSIONS,
) -> list[dict]:
    """Lance le docking de tous les ligands (ou un sous-ensemble). Retourne la liste des erreurs."""
    if not receptor_pdbqt.is_file():
        raise FileNotFoundError(f"Récepteur '{receptor_pdbqt}' introuvable.")

    ligand_files = sorted(ligands_dir.glob("*.pdbqt"))
    if not ligand_files:
        raise FileNotFoundError(f"Aucun fichier .pdbqt dans '{ligands_dir}'.")

    if only_ligands is not None:
        ligand_files = [f for f in ligand_files if f.stem in only_ligands]
        if not ligand_files:
            raise FileNotFoundError("Aucun des ligands demandés n'a de fichier .pdbqt.")

    # Filtrer les molécules trop flexibles
    filtered = []
    skipped_torsions = 0
    for lig in ligand_files:
        t = _count_torsions(lig)
        if t > max_torsions:
            skipped_torsions += 1
        else:
            filtered.append(lig)
    if skipped_torsions:
        logger.info("Skip de %d ligand(s) avec >%d torsions (trop flexibles pour Vina).", skipped_torsions, max_torsions)
    ligand_files = filtered

    results_dir.mkdir(parents=True, exist_ok=True)

    if skip_existing:
        to_dock = []
        skipped = 0
        corrupt = 0
        for lig in ligand_files:
            out_file = results_dir / f"{lig.stem}_out.pdbqt"
            if out_file.exists():
                if _is_valid_result(out_file):
                    skipped += 1
                else:
                    out_file.unlink()
                    corrupt += 1
                    to_dock.append(lig)
            else:
                to_dock.append(lig)
        if skipped:
            logger.info("Skip de %d ligand(s) déjà dockés.", skipped)
        if corrupt:
            logger.warning("Suppression de %d résultat(s) corrompus (re-docking).", corrupt)
    else:
        to_dock = ligand_files

    if not to_dock:
        logger.info("Tous les ligands sont déjà dockés.")
        return []

    num_workers = os.cpu_count()
    logger.info(
        "Docking : %d ligand(s), %d workers, exhaustiveness=%d",
        len(to_dock), num_workers, exhaustiveness,
    )
    logger.info(
        "Grid Box : center=(%.1f, %.1f, %.1f), size=(%.1f, %.1f, %.1f)",
        GRID_CENTER_X, GRID_CENTER_Y, GRID_CENTER_Z,
        GRID_SIZE_X, GRID_SIZE_Y, GRID_SIZE_Z,
    )

    worker = partial(
        _dock_one_ligand,
        receptor_pdbqt=str(receptor_pdbqt),
        results_dir=str(results_dir),
        exhaustiveness=exhaustiveness,
    )

    errors = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(worker, str(lig)): lig.stem for lig in to_dock}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Docking"):
            res = future.result()
            if res["status"] != "ok":
                errors.append(res)

    if errors:
        logger.warning("%d docking(s) en erreur :", len(errors))
        for err in errors:
            logger.warning("  - %s : %s", err["ligand"], err.get("msg", err["status"]))

    succeeded = len(to_dock) - len(errors)
    logger.info("Docking terminé : %d/%d réussis.", succeeded, len(to_dock))
    return errors
