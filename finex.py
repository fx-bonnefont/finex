#!/usr/bin/env python3
"""finex — Pipeline de criblage virtuel pour le repositionnement de médicaments."""

import argparse
import logging
import sys

from finex.config import (
    DEFAULT_EXHAUSTIVENESS,
    DOCKING_RESULTS_DIR,
    LIGANDS_IN_DIR,
    LIGANDS_PDBQT_DIR,
    RECEPTOR_PDB,
    RECEPTOR_PDBQT,
    RESULTS_XLSX,
    SDF_DROP_DIR,
    VIZ_DIR,
)



def setup_logging(verbose: bool = False):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def cmd_run(args):
    """Exécute le pipeline complet de docking."""
    from finex.analyze_results import export_xlsx, load_drug_names, parse_all_scores
    from finex.fetch_target import fetch_target
    from finex.prepare_molecules import (
        check_obabel,
        convert_ligands,
        convert_receptor,
        fetch_fda_drugs,
    )
    from finex.run_docking import check_vina, dock_all_ligands
    from finex.sdf_utils import auto_split_drop_dir

    logger = logging.getLogger("finex")

    # Étape 0 : vérifier les outils externes
    check_obabel()
    check_vina()

    # Étape 1 : récupérer la structure du récepteur
    if RECEPTOR_PDB.exists():
        logger.info("Récepteur PDB existant, skip du fetch.")
    else:
        logger.info("Récupération de la structure 3D via OpenFold3...")
        fetch_target(RECEPTOR_PDB)

    # Étape 2 : reset si --fresh-db demandé
    if args.fresh_db:
        import shutil
        logger.info("--fresh-db : suppression des données existantes...")
        for d in (LIGANDS_IN_DIR, LIGANDS_PDBQT_DIR, DOCKING_RESULTS_DIR):
            if d.exists():
                shutil.rmtree(d)
                logger.info("  Supprimé : %s", d)

    # Étape 3 : auto-split des SDF déposés dans sdf_drop/
    SDF_DROP_DIR.mkdir(parents=True, exist_ok=True)
    split_count = auto_split_drop_dir(SDF_DROP_DIR, LIGANDS_IN_DIR)
    if split_count > 0:
        logger.info("Auto-split : %d molécule(s) importées depuis sdf_drop/.", split_count)

    # Étape 4 : télécharger les médicaments FDA si aucun ligand disponible
    LIGANDS_IN_DIR.mkdir(parents=True, exist_ok=True)
    sdf_files = list(LIGANDS_IN_DIR.glob("*.sdf"))
    if not sdf_files:
        logger.info("Aucun SDF trouvé, téléchargement FDA depuis ChEMBL...")
        fetch_fda_drugs(LIGANDS_IN_DIR)
        sdf_files = list(LIGANDS_IN_DIR.glob("*.sdf"))

    if not sdf_files:
        logger.error("Aucun fichier SDF disponible. Abandon.")
        sys.exit(1)

    logger.info("%d ligand(s) SDF disponibles.", len(sdf_files))

    # Étape 4 : convertir le récepteur PDB → PDBQT
    if RECEPTOR_PDBQT.exists():
        logger.info("Récepteur PDBQT existant, skip de la conversion.")
    else:
        convert_receptor(RECEPTOR_PDB, RECEPTOR_PDBQT)

    # Étape 5 : convertir les ligands SDF → PDBQT
    convert_ligands(LIGANDS_IN_DIR, LIGANDS_PDBQT_DIR, skip_existing=True)

    # Étape 6 : docking
    dock_all_ligands(
        receptor_pdbqt=RECEPTOR_PDBQT,
        ligands_dir=LIGANDS_PDBQT_DIR,
        results_dir=DOCKING_RESULTS_DIR,
        exhaustiveness=args.exhaustiveness,
        skip_existing=True,
        max_torsions=args.max_torsions,
    )

    # Étape 7 : analyse et export Excel
    scores = parse_all_scores(DOCKING_RESULTS_DIR)
    if not scores:
        logger.error("Aucun score de docking trouvé.")
        sys.exit(1)

    drug_names = load_drug_names(LIGANDS_IN_DIR)
    xlsx_path = export_xlsx(scores, RESULTS_XLSX, drug_names=drug_names)
    logger.info("Résultats exportés : %s (%d molécules classées).", xlsx_path, len(scores))

    # Afficher le top 10
    print(f"\n{'Rang':<6} {'ChEMBL ID':<20} {'Médicament':<30} {'Score (kcal/mol)'}")
    print("-" * 80)
    for i, (chembl_id, score) in enumerate(scores[:10], 1):
        name = drug_names.get(chembl_id, "")
        print(f"{i:<6} {chembl_id:<20} {name:<30} {score:.1f}")
    print(f"\nClassement complet : {xlsx_path}")


def cmd_refine(args):
    """Relance le docking sur les N meilleurs candidats avec un exhaustiveness plus élevé."""
    from finex.analyze_results import export_xlsx, load_drug_names, parse_all_scores
    from finex.run_docking import check_vina, dock_all_ligands

    logger = logging.getLogger("finex")

    check_vina()

    # Lire les scores existants
    scores = parse_all_scores(DOCKING_RESULTS_DIR)
    if not scores:
        logger.error("Aucun résultat de docking existant. Lancez d'abord : python finex.py run")
        sys.exit(1)

    top_n = args.top
    top_ligands = [name for name, _ in scores[:top_n]]
    logger.info(
        "Refine : re-docking des %d meilleurs candidats avec exhaustiveness=%d",
        len(top_ligands), args.exhaustiveness,
    )

    # Sauvegarder les anciens résultats avant de relancer
    backup_dir = DOCKING_RESULTS_DIR / "backup_pre_refine"
    backup_dir.mkdir(parents=True, exist_ok=True)
    backed_up = 0
    for name in top_ligands:
        out_file = DOCKING_RESULTS_DIR / f"{name}_out.pdbqt"
        if out_file.exists():
            import shutil
            shutil.copy2(out_file, backup_dir / out_file.name)
            out_file.unlink()
            backed_up += 1
    logger.info("Backup de %d résultat(s) dans %s, lancement du re-docking...", backed_up, backup_dir)

    # Relancer le docking uniquement sur ces ligands
    errors = dock_all_ligands(
        receptor_pdbqt=RECEPTOR_PDBQT,
        ligands_dir=LIGANDS_PDBQT_DIR,
        results_dir=DOCKING_RESULTS_DIR,
        exhaustiveness=args.exhaustiveness,
        skip_existing=False,
        only_ligands=set(top_ligands),
    )

    # Restaurer les backups pour les ligands en échec (timeout, erreur)
    failed_names = {e["ligand"] for e in errors}
    restored = 0
    for name in failed_names:
        backup_file = backup_dir / f"{name}_out.pdbqt"
        target_file = DOCKING_RESULTS_DIR / f"{name}_out.pdbqt"
        if backup_file.exists() and not target_file.exists():
            import shutil
            shutil.copy2(backup_file, target_file)
            restored += 1
    if restored:
        logger.info("Restauration de %d résultat(s) depuis le backup (timeout/erreur).", restored)

    # Re-analyser et exporter
    scores = parse_all_scores(DOCKING_RESULTS_DIR)
    drug_names = load_drug_names(LIGANDS_IN_DIR)
    xlsx_path = export_xlsx(scores, RESULTS_XLSX, drug_names=drug_names)

    print(f"\n{'Rang':<6} {'ChEMBL ID':<20} {'Médicament':<30} {'Score (kcal/mol)'}")
    print("-" * 80)
    for i, (chembl_id, score) in enumerate(scores[:top_n], 1):
        name = drug_names.get(chembl_id, "")
        print(f"{i:<6} {chembl_id:<20} {name:<30} {score:.1f}")
    print(f"\nClassement mis à jour : {xlsx_path}")


def cmd_fetch_supplements(args):
    """Télécharge CMAUP et/ou FooDB et génère les SDF 3D."""
    from finex.fetch_supplements import (
        fetch_cmaup,
        fetch_foodb,
        list_cmaup_plants,
        list_foodb_foods,
    )

    logger = logging.getLogger("finex")

    # --- Mode listing ---
    if args.list_plants:
        plants = list_cmaup_plants()
        print(f"\n{len(plants)} plantes disponibles dans CMAUP :\n")
        for p in plants:
            print(f"  {p}")
        return

    if args.list_foods:
        foods = list_foodb_foods()
        print(f"\n{len(foods)} aliments disponibles dans FooDB :\n")
        for f in foods:
            print(f"  {f}")
        return

    # --- Parsing du filtre ---
    plant_filter = (
        {t.strip() for t in args.filter.split(",")} if args.filter and args.source != "foodb"
        else None
    )
    food_filter = (
        {t.strip() for t in args.filter.split(",")} if args.filter and args.source != "cmaup"
        else None
    )

    total = 0
    LIGANDS_IN_DIR.mkdir(parents=True, exist_ok=True)

    if args.source in ("cmaup", "all"):
        logger.info("=== Source : CMAUP v2.0 ===")
        n = fetch_cmaup(LIGANDS_IN_DIR, plant_filter=plant_filter)
        total += n

    if args.source in ("foodb", "all"):
        logger.info("=== Source : FooDB ===")
        n = fetch_foodb(LIGANDS_IN_DIR, food_filter=food_filter)
        total += n

    logger.info("Total : %d SDF générés dans '%s'.", total, LIGANDS_IN_DIR)
    if total > 0:
        print(f"\n{total} molécule(s) prêtes dans {LIGANDS_IN_DIR}/")
        print("Lancez le docking avec : python finex.py run")


def cmd_report(args):
    """Recompile le fichier Excel avec les résultats de docking disponibles."""
    from finex.analyze_results import export_xlsx, load_drug_names, parse_all_scores

    logger = logging.getLogger("finex")

    scores = parse_all_scores(DOCKING_RESULTS_DIR)
    if not scores:
        logger.error("Aucun résultat de docking trouvé. Lancez d'abord : python finex.py run")
        sys.exit(1)

    drug_names = load_drug_names(LIGANDS_IN_DIR)
    xlsx_path = export_xlsx(scores, RESULTS_XLSX, drug_names=drug_names)

    total_sdf = len(list(LIGANDS_IN_DIR.glob("*.sdf")))
    logger.info("Rapport : %d/%d molécules dockées.", len(scores), total_sdf)
    logger.info("Fichier mis à jour : %s", xlsx_path)

    top_n = min(args.top, len(scores))
    print(f"\n{'Rang':<6} {'ID':<22} {'Nom':<30} {'Score (kcal/mol)'}")
    print("-" * 82)
    for i, (chembl_id, score) in enumerate(scores[:top_n], 1):
        name = drug_names.get(chembl_id, "")
        print(f"{i:<6} {chembl_id:<22} {name:<30} {score:.1f}")
    print(f"\n{len(scores)}/{total_sdf} molécules dockées — classement complet : {xlsx_path}")


def cmd_viz(args):
    """Génère les fichiers de visualisation pour les top N candidats."""
    from finex.analyze_results import load_drug_names, parse_all_scores
    from finex.prepare_molecules import check_obabel
    from finex.prepare_viz import prepare_viz_top_n

    logger = logging.getLogger("finex")

    check_obabel()

    scores = parse_all_scores(DOCKING_RESULTS_DIR)
    if not scores:
        logger.error("Aucun score trouvé. Lancez d'abord : python finex.py run")
        sys.exit(1)

    drug_names = load_drug_names(LIGANDS_IN_DIR)

    paths = prepare_viz_top_n(
        scores=scores,
        top_n=args.top,
        results_dir=DOCKING_RESULTS_DIR,
        receptor_pdb=RECEPTOR_PDB,
        viz_dir=VIZ_DIR,
    )

    print(f"\n{len(paths)} fichier(s) de visualisation générés dans {VIZ_DIR}/ :")
    for path, (chembl_id, score) in zip(paths, scores[: args.top]):
        name = drug_names.get(chembl_id, "")
        print(f"  {path}  ({name}, {score:.1f} kcal/mol)")
    print("\nOuvrir avec PyMOL : pymol output/viz/*_complex.pdb")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="finex",
        description="Pipeline de criblage virtuel pour le repositionnement de médicaments (Long COVID / xenoAMP).",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Activer le mode debug")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- run ---
    run_parser = subparsers.add_parser("run", help="Lancer le pipeline complet de docking")
    run_parser.add_argument(
        "-e", "--exhaustiveness", type=int, default=DEFAULT_EXHAUSTIVENESS,
        help=f"Paramètre exhaustiveness de Vina (défaut : {DEFAULT_EXHAUSTIVENESS})",
    )
    run_parser.add_argument(
        "--max-torsions", type=int, default=15,
        help="Exclure les ligands avec plus de N torsions (défaut : 15)",
    )
    run_parser.add_argument(
        "--fresh-db", action="store_true",
        help="Retélécharger la base FDA complète depuis ChEMBL (supprime les anciens ligands/résultats)",
    )
    run_parser.set_defaults(func=cmd_run)

    # --- refine ---
    refine_parser = subparsers.add_parser("refine", help="Re-docker les top N candidats avec un exhaustiveness élevé")
    refine_parser.add_argument(
        "-e", "--exhaustiveness", type=int, default=128,
        help="Paramètre exhaustiveness de Vina (défaut : 128)",
    )
    refine_parser.add_argument(
        "--top", type=int, default=100,
        help="Nombre de meilleurs candidats à re-docker (défaut : 100)",
    )
    refine_parser.set_defaults(func=cmd_refine)

    # --- fetch-supplements ---
    sup_parser = subparsers.add_parser(
        "fetch-supplements",
        help="Télécharger CMAUP / FooDB et générer les SDF 3D",
    )
    sup_parser.add_argument(
        "--source",
        choices=["cmaup", "foodb", "all"],
        default="all",
        help="Source à télécharger : cmaup, foodb, ou all (défaut : all)",
    )
    sup_parser.add_argument(
        "--filter",
        metavar="TERMES",
        default=None,
        help=(
            "Filtrer par plante/aliment (virgule-séparé, insensible à la casse). "
            "Ex. : --filter 'curcuma longa,ginkgo' (CMAUP) "
            "ou --filter 'turmeric,garlic' (FooDB)"
        ),
    )
    sup_parser.add_argument(
        "--list-plants",
        action="store_true",
        help="Lister les plantes disponibles dans CMAUP et quitter",
    )
    sup_parser.add_argument(
        "--list-foods",
        action="store_true",
        help="Lister les aliments disponibles dans FooDB et quitter",
    )
    sup_parser.set_defaults(func=cmd_fetch_supplements)

    # --- report ---
    report_parser = subparsers.add_parser("report", help="Recompiler le fichier Excel avec les résultats disponibles")
    report_parser.add_argument(
        "--top", type=int, default=20,
        help="Nombre de meilleurs candidats à afficher (défaut : 20)",
    )
    report_parser.set_defaults(func=cmd_report)

    # --- viz ---
    viz_parser = subparsers.add_parser("viz", help="Générer les fichiers de visualisation")
    viz_parser.add_argument(
        "--top", type=int, default=5,
        help="Nombre de molécules à visualiser (défaut : 5)",
    )
    viz_parser.set_defaults(func=cmd_viz)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    setup_logging(getattr(args, "verbose", False))
    args.func(args)


if __name__ == "__main__":
    main()
