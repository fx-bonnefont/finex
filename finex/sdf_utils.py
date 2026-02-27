"""Utilitaires pour splitter des fichiers SDF multi-molécules."""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def split_sdf(input_sdf: Path, output_dir: Path) -> int:
    """Sépare un SDF multi-molécules en fichiers individuels. Retourne le nombre écrit."""
    from rdkit import Chem

    output_dir.mkdir(parents=True, exist_ok=True)
    supplier = Chem.SDMolSupplier(str(input_sdf))
    count = 0

    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name").strip() else ""
        if name:
            safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
        else:
            safe_name = f"ligand_{i}"

        out_path = output_dir / f"{safe_name}.sdf"
        # Éviter les collisions de noms
        if out_path.exists():
            out_path = output_dir / f"{safe_name}_{i}.sdf"

        writer = Chem.SDWriter(str(out_path))
        writer.write(mol)
        writer.close()
        count += 1

    logger.info("Split de '%s' : %d molécule(s) écrites dans '%s'.", input_sdf.name, count, output_dir)
    return count


def auto_split_drop_dir(drop_dir: Path, ligands_dir: Path) -> int:
    """Scanne drop_dir pour des SDF multi-molécules et les split dans ligands_dir."""
    if not drop_dir.exists():
        return 0

    total = 0
    for sdf_file in sorted(drop_dir.glob("*.sdf")):
        total += split_sdf(sdf_file, ligands_dir)

    return total
