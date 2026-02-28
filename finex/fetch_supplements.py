"""Génération 3D des composés issus de CMAUP et FooDB (fichiers locaux).

Sources :
- CMAUP v2.0  : plantes médicinales → principes actifs (60k composés, 7.8k plantes)
- FooDB       : aliments → composés alimentaires (70k composés, ~1k aliments)

Fichiers attendus dans data/supplements_cache/ :
    cmaup/
        CMAUPv2.0_download_Plants.txt
        CMAUPv2.0_download_Ingredients_onlyActive.txt
        CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt
    foodb/
        Food.csv
        Compound.csv
        Content.csv
"""

import csv
import logging
from multiprocessing import Pool, Process, Queue, cpu_count
from pathlib import Path

from tqdm import tqdm

logger = logging.getLogger(__name__)

# --- Chemins vers les fichiers locaux ---
CACHE_DIR = Path("data") / "supplements_cache"

CMAUP_DIR = CACHE_DIR / "cmaup"
CMAUP_PLANTS_FILE = CMAUP_DIR / "CMAUPv2.0_download_Plants.txt"
CMAUP_INGREDIENTS_FILE = CMAUP_DIR / "CMAUPv2.0_download_Ingredients_onlyActive.txt"
CMAUP_ASSOC_FILE = CMAUP_DIR / "CMAUPv2.0_download_Plant_Ingredient_Associations_onlyActiveIngredients.txt"

FOODB_DIR = CACHE_DIR / "foodb"
FOODB_FOOD_FILE = FOODB_DIR / "Food.csv"
FOODB_COMPOUND_FILE = FOODB_DIR / "Compound.csv"
FOODB_CONTENT_FILE = FOODB_DIR / "Content.csv"

# Atomes organiques autorisés (même filtre que prepare_molecules)
_ORGANIC_ATOMS = {1, 5, 6, 7, 8, 9, 15, 16, 17, 34, 35, 53}


# ---------------------------------------------------------------------------
# Génération 3D (worker multiprocessing)
# ---------------------------------------------------------------------------

def _generate_3d_worker(args: tuple) -> dict:
    """Worker : SMILES → SDF 3D. Même logique que prepare_molecules._generate_3d_sdf."""
    compound_id, name, smiles, output_dir = args
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    try:
        if "." in smiles:
            smiles = max(smiles.split("."), key=len)

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"status": "skip", "id": compound_id}

        if not all(a.GetAtomicNum() in _ORGANIC_ATOMS for a in mol.GetAtoms()):
            return {"status": "skip", "id": compound_id}

        # Filtre torsions AVANT l'embedding (calcul 2D, instantané)
        # Cohérent avec le filtre --max-torsions du docking (défaut 20)
        from rdkit.Chem import Lipinski
        if Lipinski.NumRotatableBonds(mol) > 15:
            return {"status": "skip", "id": compound_id}

        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) == -1:
            return {"status": "skip", "id": compound_id}
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

        mol.SetProp("_Name", name[:80])
        mol.SetProp("compound_id", compound_id)

        sdf_path = Path(output_dir) / f"{compound_id}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()
        return {"status": "ok"}
    except Exception:
        return {"status": "skip", "id": compound_id}


def _run_3d_generation(entries: list[tuple], output_dir: Path) -> int:
    """Lance la génération 3D en parallèle. Retourne le nombre de succès.

    Les composés définitivement non-générables (métal, SMILES invalide, embedding
    impossible) reçoivent un marqueur .skip pour ne plus être retentés.
    """
    if not entries:
        return 0
    output_dir.mkdir(parents=True, exist_ok=True)

    # Petits lots : traitement séquentiel avec timeout par molécule
    # Grands lots : Pool multiprocessing
    use_pool = len(entries) >= 200

    ok = skipped = 0

    def _record(res: dict) -> None:
        nonlocal ok, skipped
        if res["status"] == "ok":
            ok += 1
        else:
            skipped += 1
            cid = res.get("id", "")
            if cid:
                (output_dir / f"{cid}.skip").touch()

    if use_pool:
        with Pool(processes=cpu_count()) as pool:
            for res in tqdm(
                pool.imap_unordered(_generate_3d_worker, entries, chunksize=16),
                total=len(entries),
                desc="Génération 3D",
            ):
                _record(res)
    else:
        # Timeout par molécule via Process individuel (évite les blocages RDKit)
        _TIMEOUT = 30  # secondes max par composé
        for entry in tqdm(entries, total=len(entries), desc="Génération 3D"):
            q: Queue = Queue()
            p = Process(target=_worker_enqueue, args=(entry, q))
            p.start()
            p.join(timeout=_TIMEOUT)
            if p.is_alive():
                p.kill()
                p.join()
                _record({"status": "skip", "id": entry[0]})
            else:
                _record(q.get() if not q.empty() else {"status": "skip", "id": entry[0]})

    if skipped:
        logger.info("%d composé(s) ignorés (métal, embedding impossible).", skipped)
    return ok


def _worker_enqueue(args: tuple, q: Queue) -> None:
    """Wrapper pour exécuter _generate_3d_worker et mettre le résultat dans une Queue."""
    q.put(_generate_3d_worker(args))


# ---------------------------------------------------------------------------
# CMAUP
# ---------------------------------------------------------------------------

def _check_cmaup_files() -> None:
    """Vérifie que les fichiers CMAUP locaux sont présents."""
    for f in (CMAUP_PLANTS_FILE, CMAUP_INGREDIENTS_FILE, CMAUP_ASSOC_FILE):
        if not f.exists():
            raise FileNotFoundError(
                f"Fichier CMAUP manquant : {f}\n"
                "Téléchargez les fichiers depuis https://bidd.group/CMAUP/download.html\n"
                f"et placez-les dans {CMAUP_DIR}/"
            )


def list_cmaup_plants(n: int = 0) -> list[str]:
    """Retourne la liste des plantes disponibles dans CMAUP.

    Args:
        n: si > 0, retourne les n premières. 0 = tout.
    """
    _check_cmaup_files()
    plants = set()
    with CMAUP_PLANTS_FILE.open(encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row.get("Plant_Name", "").strip()
            if name:
                plants.add(name)
    result = sorted(plants)
    return result[:n] if n > 0 else result


def _cmaup_filter_ids(plant_filter: set[str]) -> set[str]:
    """Plant_Name → Plant_ID → Ingredient_ID (= np_id).

    La correspondance est insensible à la casse et partielle (substring).
    """
    filter_lower = {p.lower() for p in plant_filter}

    # Étape 1 : Plant_Name/Species_Name → Plant_ID
    matched_plant_ids: set[str] = set()
    with CMAUP_PLANTS_FILE.open(encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for col in ("Plant_Name", "Species_Name", "Genus_Name", "Family_Name"):
                val = row.get(col, "").lower()
                if val and any(term in val for term in filter_lower):
                    pid = row.get("Plant_ID", "").strip()
                    if pid:
                        matched_plant_ids.add(pid)
                    break

    if not matched_plant_ids:
        return set()

    # Étape 2 : Plant_ID → Ingredient_ID (np_id)
    # Le fichier d'associations n'a PAS de ligne d'en-tête — deux colonnes tab-séparées
    ingredient_ids: set[str] = set()
    with CMAUP_ASSOC_FILE.open(encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=["Plant_ID", "Ingredient_ID"])
        for row in reader:
            if row["Plant_ID"].strip() in matched_plant_ids:
                ing_id = row["Ingredient_ID"].strip()
                if ing_id:
                    ingredient_ids.add(ing_id)

    return ingredient_ids


def fetch_cmaup(
    output_dir: Path,
    plant_filter: set[str] | None = None,
) -> int:
    """Génère les SDF 3D des ingrédients actifs CMAUP dans `output_dir`.

    Args:
        output_dir   : répertoire de sortie des SDF (ex. data/ligands_in)
        plant_filter : termes de recherche sur les noms de plantes (insensible
                       à la casse, correspondance partielle).
                       Ex. {"curcuma longa", "ginkgo"}.
                       Si None, toutes les plantes sont incluses.
    Returns:
        Nombre de SDF générés avec succès.
    """
    _check_cmaup_files()

    allowed_ids: set[str] | None = None
    if plant_filter is not None:
        allowed_ids = _cmaup_filter_ids(plant_filter)
        logger.info(
            "Filtre CMAUP : %d ingrédient(s) pour %d terme(s) de plante.",
            len(allowed_ids), len(plant_filter),
        )
        if not allowed_ids:
            logger.warning(
                "Aucune plante trouvée pour : %s\n"
                "Listez les plantes avec : python finex.py fetch-supplements --list-plants",
                ", ".join(plant_filter),
            )
            return 0

    entries = _cmaup_parse_ingredients(output_dir, allowed_ids)
    logger.info("CMAUP : %d composé(s) à générer (existants skippés).", len(entries))
    ok = _run_3d_generation(entries, output_dir)
    logger.info("CMAUP : %d SDF générés avec succès.", ok)
    return ok


def _cmaup_parse_ingredients(
    output_dir: Path,
    allowed_ids: set[str] | None,
) -> list[tuple]:
    """Parse Ingredients_onlyActive.txt (colonnes : np_id, pref_name, …, SMILES)."""
    entries = []
    with CMAUP_INGREDIENTS_FILE.open(encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            raw_id = row.get("np_id", "").strip()
            if not raw_id:
                continue

            if allowed_ids is not None and raw_id not in allowed_ids:
                continue

            compound_id = f"CMAUP_{raw_id}"
            if (output_dir / f"{compound_id}.sdf").exists() or (output_dir / f"{compound_id}.skip").exists():
                continue

            smiles = row.get("SMILES", "").strip()
            if not smiles or smiles in ("", "N/A", "None", "nan"):
                continue

            name = row.get("pref_name", "").strip() or row.get("iupac_name", "").strip() or compound_id

            entries.append((compound_id, name, smiles, str(output_dir)))
    return entries


# ---------------------------------------------------------------------------
# FooDB
# ---------------------------------------------------------------------------

def _check_foodb_files() -> None:
    """Vérifie que les fichiers FooDB locaux sont présents."""
    for f in (FOODB_FOOD_FILE, FOODB_COMPOUND_FILE, FOODB_CONTENT_FILE):
        if not f.exists():
            raise FileNotFoundError(
                f"Fichier FooDB manquant : {f}\n"
                "Téléchargez les CSV depuis https://foodb.ca/downloads\n"
                f"et placez-les dans {FOODB_DIR}/"
            )


def list_foodb_foods(n: int = 0) -> list[str]:
    """Retourne la liste des aliments disponibles dans FooDB.

    Args:
        n: si > 0, retourne les n premiers. 0 = tout.
    """
    _check_foodb_files()
    foods = []
    with FOODB_FOOD_FILE.open(encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("name", "").strip()
            if name:
                foods.append(name)
    result = sorted(set(foods))
    return result[:n] if n > 0 else result


def _foodb_filter_compound_ids(food_filter: set[str]) -> set[str]:
    """food_name → food_id → source_id (compound_id) via Content.csv."""
    filter_lower = {f.lower() for f in food_filter}

    # Food.csv : trouver les food_id correspondant aux aliments filtrés
    food_ids: set[str] = set()
    with FOODB_FOOD_FILE.open(encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("name", "").lower()
            sci = row.get("name_scientific", "").lower()
            if any(term in name or term in sci for term in filter_lower):
                fid = row.get("id", "").strip()
                if fid:
                    food_ids.add(fid)

    if not food_ids:
        return set()

    # Content.csv : food_id → source_id où source_type == "Compound"
    compound_ids: set[str] = set()
    with FOODB_CONTENT_FILE.open(encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("source_type", "") == "Compound" and row.get("food_id", "") in food_ids:
                cid = row.get("source_id", "").strip()
                if cid:
                    compound_ids.add(cid)
    return compound_ids


def fetch_foodb(
    output_dir: Path,
    food_filter: set[str] | None = None,
) -> int:
    """Génère les SDF 3D des composés FooDB dans `output_dir`.

    Args:
        output_dir  : répertoire de sortie des SDF
        food_filter : noms d'aliments (insensible à la casse, correspondance
                      partielle). Ex. {"turmeric", "green tea", "garlic"}.
                      Si None, tous les composés sont inclus.
    Returns:
        Nombre de SDF générés avec succès.
    """
    _check_foodb_files()

    allowed_compound_ids: set[str] | None = None
    if food_filter is not None:
        allowed_compound_ids = _foodb_filter_compound_ids(food_filter)
        logger.info(
            "Filtre FooDB : %d composé(s) pour %d terme(s) d'aliment.",
            len(allowed_compound_ids), len(food_filter),
        )
        if not allowed_compound_ids:
            logger.warning(
                "Aucun aliment trouvé pour : %s\n"
                "Listez les aliments avec : python finex.py fetch-supplements --list-foods",
                ", ".join(food_filter),
            )
            return 0

    entries = _foodb_parse_compounds(output_dir, allowed_compound_ids)
    logger.info("FooDB : %d composé(s) à générer (existants skippés).", len(entries))
    ok = _run_3d_generation(entries, output_dir)
    logger.info("FooDB : %d SDF générés avec succès.", ok)
    return ok


def _foodb_parse_compounds(
    output_dir: Path,
    allowed_ids: set[str] | None,
) -> list[tuple]:
    """Parse Compound.csv FooDB → liste de tuples (id, name, smiles, dir).

    Note : les colonnes du CSV FooDB sont décalées par rapport au header —
    le SMILES réel est à l'index 7 (colonne nommée 'cas_number' dans le header).
    On utilise un accès positionnel pour contourner ce désalignement.
    """
    entries = []
    with FOODB_COMPOUND_FILE.open(encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if len(row) < 8:
                continue

            raw_id = row[0].strip()
            if not raw_id:
                continue

            if allowed_ids is not None and raw_id not in allowed_ids:
                continue

            public_id = row[1].strip() or raw_id
            compound_id = f"FOODB_{public_id}"
            if (output_dir / f"{compound_id}.sdf").exists() or (output_dir / f"{compound_id}.skip").exists():
                continue

            smiles = row[7].strip()  # SMILES réel (décalage CSV FooDB)
            if not smiles or smiles in ("", "N/A", "None"):
                continue

            name = row[2].strip() or compound_id

            entries.append((compound_id, name, smiles, str(output_dir)))
    return entries
