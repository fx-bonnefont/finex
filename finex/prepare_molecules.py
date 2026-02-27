"""Conversion du récepteur et des ligands aux formats requis pour le docking."""

import logging
import subprocess
from multiprocessing import Pool, cpu_count
from pathlib import Path

from tqdm import tqdm

logger = logging.getLogger(__name__)


def check_obabel() -> None:
    """Vérifie que obabel est dans le PATH."""
    try:
        result = subprocess.run(["obabel", "-V"], capture_output=True, text=True, timeout=10)
        logger.info("OpenBabel détecté : %s", result.stdout.strip().split("\n")[0])
    except FileNotFoundError:
        raise FileNotFoundError(
            "obabel introuvable. Installez OpenBabel :\n"
            "  macOS  : brew install open-babel\n"
            "  Ubuntu : sudo apt-get install openbabel\n"
            "  Conda  : conda install -c conda-forge openbabel"
        )


def _generate_3d_sdf(args: tuple) -> dict:
    """Worker : génère les coordonnées 3D et écrit le SDF. Exécuté en parallèle."""
    chembl_id, name, smiles, output_dir = args
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    # Supprimer les warnings UFF/MMFF (atomes métalliques, charges exotiques)
    RDLogger.DisableLog("rdApp.*")

    try:
        # Garder uniquement le fragment principal (éliminer les contre-ions de type sel)
        if "." in smiles:
            smiles = max(smiles.split("."), key=len)

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"status": "skip"}

        # Filtrer les molécules contenant des métaux (inutiles pour le docking)
        organic_atoms = {1, 5, 6, 7, 8, 9, 15, 16, 17, 34, 35, 53}  # H,B,C,N,O,F,P,S,Cl,Se,Br,I
        if not all(a.GetAtomicNum() in organic_atoms for a in mol.GetAtoms()):
            return {"status": "skip"}

        mol = Chem.AddHs(mol)
        res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if res == -1:
            return {"status": "skip"}
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        mol.SetProp("_Name", name)
        mol.SetProp("chembl_id", chembl_id)

        sdf_path = Path(output_dir) / f"{chembl_id}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()
        return {"status": "ok"}
    except Exception:
        return {"status": "error"}


def fetch_fda_drugs(output_dir: Path) -> int:
    """Télécharge les médicaments FDA-approved depuis ChEMBL et génère les SDF 3D en parallèle."""
    from chembl_webresource_client.new_client import new_client
    from chembl_webresource_client.settings import Settings

    # Augmenter la taille des pages API pour réduire les requêtes réseau
    Settings.Instance().MAX_LIMIT = 1000

    logger.info("Phase 1/2 : Téléchargement des SMILES depuis ChEMBL (max_phase=4)...")
    molecule = new_client.molecule
    molecule.set_format("json")
    approved = molecule.filter(max_phase=4, molecule_type="Small molecule").only(
        ["molecule_chembl_id", "pref_name", "molecule_structures"]
    )

    # Collecter toutes les molécules d'abord (I/O réseau)
    entries = []
    for m in tqdm(approved, desc="Téléchargement ChEMBL"):
        structs = m.get("molecule_structures")
        if not structs or not structs.get("canonical_smiles"):
            continue
        entries.append((
            m["molecule_chembl_id"],
            m.get("pref_name") or m["molecule_chembl_id"],
            structs["canonical_smiles"],
        ))

    logger.info("%d molécules avec SMILES récupérées.", len(entries))

    # Phase 2 : génération 3D en parallèle (CPU-bound)
    output_dir.mkdir(parents=True, exist_ok=True)
    num_workers = cpu_count()
    logger.info("Phase 2/2 : Génération 3D avec %d workers...", num_workers)

    work_items = [(cid, name, smi, str(output_dir)) for cid, name, smi in entries]
    ok = 0
    errors = 0

    with Pool(processes=num_workers) as pool:
        for res in tqdm(
            pool.imap_unordered(_generate_3d_sdf, work_items, chunksize=16),
            total=len(work_items),
            desc="Génération 3D",
        ):
            if res["status"] == "ok":
                ok += 1
            elif res["status"] == "error":
                errors += 1

    logger.info("ChEMBL : %d molécules FDA écrites dans '%s' (%d erreurs ignorées).", ok, output_dir, errors)
    return ok


def convert_receptor(pdb_path: Path, pdbqt_path: Path) -> Path:
    """Convertit le PDB du récepteur en PDBQT."""
    if not pdb_path.is_file():
        raise FileNotFoundError(f"Fichier récepteur '{pdb_path}' introuvable.")

    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["obabel", str(pdb_path), "-O", str(pdbqt_path), "-xr", "--partialcharge", "gasteiger"]
    logger.info("Conversion du récepteur : %s -> %s", pdb_path, pdbqt_path)
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

    if result.returncode != 0:
        raise RuntimeError(f"Erreur conversion récepteur :\n{result.stderr}")

    logger.info("Récepteur converti avec succès.")
    return pdbqt_path


def _is_valid_pdbqt(path: Path) -> bool:
    """Vérifie qu'un PDBQT ligand est valide (une seule molécule, atomes + structure ROOT/TORSDOF)."""
    try:
        text = path.read_text()
        return (
            len(text) > 200
            and ("ATOM" in text or "HETATM" in text)
            and "ROOT" in text
            and "ENDROOT" in text
            and text.count("TORSDOF") == 1
        )
    except Exception:
        return False


def _convert_one_ligand(args: tuple) -> dict:
    """Worker : convertit un SDF en PDBQT. Reçoit (sdf_path, out_dir) en tuple."""
    import tempfile
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    sdf_path, out_dir = Path(args[0]), Path(args[1])
    out_path = out_dir / sdf_path.with_suffix(".pdbqt").name

    # Prétraitement : extraire le fragment le plus grand du SDF (contre-ions, sels)
    sdf_to_use = sdf_path
    tmp_path = None
    try:
        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=True)
        if mol is not None:
            frags = Chem.GetMolFrags(mol, asMols=True)
            if len(frags) > 1:
                largest = max(frags, key=lambda m: m.GetNumAtoms())
                fd, tmp_str = tempfile.mkstemp(suffix=".sdf", prefix="_tmp_")
                import os; os.close(fd)
                tmp_path = Path(tmp_str)
                writer = Chem.SDWriter(str(tmp_path))
                writer.write(largest)
                writer.close()
                sdf_to_use = tmp_path
    except Exception:
        pass

    cmd = [
        "obabel", str(sdf_to_use), "-O", str(out_path),
        "--partialcharge", "gasteiger", "-h",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if tmp_path and tmp_path.exists():
            tmp_path.unlink()
        if result.returncode != 0:
            return {"file": sdf_path.name, "status": "error", "msg": result.stderr.strip()}
        if not _is_valid_pdbqt(out_path):
            if out_path.exists():
                out_path.unlink()
            return {"file": sdf_path.name, "status": "error", "msg": "obabel: PDBQT invalide (0 atomes ou structure incorrecte)"}
        return {"file": sdf_path.name, "status": "ok"}
    except subprocess.TimeoutExpired:
        if tmp_path and tmp_path.exists():
            tmp_path.unlink()
        return {"file": sdf_path.name, "status": "timeout"}
    except Exception as e:
        if tmp_path and tmp_path.exists():
            tmp_path.unlink()
        return {"file": sdf_path.name, "status": "error", "msg": str(e)}


def convert_ligands(sdf_dir: Path, pdbqt_dir: Path, skip_existing: bool = True) -> list[dict]:
    """Convertit tous les SDF en PDBQT. Skip ceux qui existent déjà si demandé."""
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    sdf_files = sorted(sdf_dir.glob("*.sdf"))

    if not sdf_files:
        logger.warning("Aucun fichier .sdf trouvé dans '%s'.", sdf_dir)
        return []

    if skip_existing:
        to_convert = []
        skipped = 0
        corrupt = 0
        for sdf in sdf_files:
            target = pdbqt_dir / sdf.with_suffix(".pdbqt").name
            if target.exists():
                # Vérifier que le PDBQT n'est pas vide/tronqué/malformé
                if _is_valid_pdbqt(target):
                    skipped += 1
                else:
                    target.unlink()
                    corrupt += 1
                    to_convert.append(sdf)
            else:
                to_convert.append(sdf)
        if skipped:
            logger.info("Skip de %d ligand(s) déjà convertis.", skipped)
        if corrupt:
            logger.warning("Suppression de %d PDBQT corrompus (re-conversion).", corrupt)
    else:
        to_convert = sdf_files

    if not to_convert:
        logger.info("Tous les ligands sont déjà convertis.")
        return []

    num_workers = cpu_count()
    logger.info("Conversion de %d ligand(s) avec %d workers...", len(to_convert), num_workers)

    work_items = [(str(f), str(pdbqt_dir)) for f in to_convert]
    errors = []

    with Pool(processes=num_workers) as pool:
        results = pool.imap_unordered(_convert_one_ligand, work_items)
        for res in tqdm(results, total=len(to_convert), desc="Conversion ligands"):
            if res["status"] != "ok":
                errors.append(res)

    if errors:
        from collections import Counter
        counts = Counter(err.get("msg", err["status"]) for err in errors)
        for msg, count in counts.most_common():
            logger.warning("%d ligand(s) ignoré(s) : %s", count, msg)

    converted = len(to_convert) - len(errors)
    logger.info("Conversion terminée : %d/%d ligand(s) convertis.", converted, len(to_convert))
    return errors
