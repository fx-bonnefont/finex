"""Analyse des résultats de docking et export Excel."""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def load_drug_names(ligands_dir: Path) -> dict[str, str]:
    """Lit la première ligne de chaque SDF pour récupérer le nom du médicament.

    Retourne un dict {chembl_id: drug_name}.
    """
    names = {}
    for sdf in ligands_dir.glob("*.sdf"):
        chembl_id = sdf.stem
        try:
            first_line = sdf.read_text().split("\n", 1)[0].strip()
            if first_line:
                names[chembl_id] = first_line
        except Exception:
            pass
    return names


def parse_all_scores(results_dir: Path) -> list[tuple[str, float]]:
    """Parse tous les fichiers _out.pdbqt et retourne les scores triés (meilleur en premier)."""
    scores = []

    for filepath in sorted(results_dir.glob("*_out.pdbqt")):
        ligand_name = filepath.stem.replace("_out", "")

        with open(filepath) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT"):
                    parts = line.split()
                    try:
                        score = float(parts[3])
                        scores.append((ligand_name, score))
                    except (IndexError, ValueError):
                        pass
                    break

    scores.sort(key=lambda x: x[1])
    return scores


def export_xlsx(
    scores: list[tuple[str, float]],
    output_path: Path,
    drug_names: dict[str, str] | None = None,
) -> Path:
    """Exporte les scores dans un fichier Excel avec noms de médicaments."""
    import pandas as pd

    df = pd.DataFrame(scores, columns=["chembl_id", "score_kcal_mol"])
    df.insert(0, "rank", range(1, len(df) + 1))

    if drug_names:
        df.insert(2, "drug_name", df["chembl_id"].map(drug_names).fillna(""))
    else:
        df.insert(2, "drug_name", "")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(output_path, index=False, engine="openpyxl", sheet_name="Docking Results")
    logger.info("Excel exporté : %s (%d molécules).", output_path, len(df))
    return output_path
