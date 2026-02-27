"""Appel à l'API OpenFold3 pour prédire la structure 3D du peptide cible."""

import json
import logging
import os
import time
from pathlib import Path

import requests
from dotenv import load_dotenv

from .config import API_URL, FASTA_SEQUENCE, MAX_WAIT, POLL_INTERVAL

load_dotenv()
logger = logging.getLogger(__name__)


def _get_api_key() -> str:
    key = os.environ.get("NVIDIA_API_KEY")
    if not key:
        raise EnvironmentError("Variable d'environnement NVIDIA_API_KEY non définie.")
    return key


def _build_payload() -> dict:
    msa_csv = f"key,sequence\n-1,{FASTA_SEQUENCE}"
    return {
        "request_id": "xenoAMP_529_558",
        "inputs": [
            {
                "input_id": "xenoAMP_529_558",
                "molecules": [
                    {
                        "type": "protein",
                        "id": "A",
                        "sequence": FASTA_SEQUENCE,
                        "msa": {
                            "main_db": {
                                "csv": {
                                    "alignment": msa_csv,
                                    "format": "csv",
                                }
                            }
                        },
                    }
                ],
                "output_format": "pdb",
            }
        ],
    }


def _submit_prediction(api_key: str) -> requests.Response:
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json",
        "NVCF-POLL-SECONDS": "300",
    }
    logger.info("Envoi de la séquence FASTA à l'API OpenFold3...")
    logger.info("Séquence : %s (%d AA)", FASTA_SEQUENCE, len(FASTA_SEQUENCE))

    response = requests.post(API_URL, headers=headers, json=_build_payload(), timeout=360)
    response.raise_for_status()
    return response


def _poll_for_result(status_url: str, api_key: str) -> dict:
    headers = {"Authorization": f"Bearer {api_key}", "Accept": "application/json"}
    elapsed = 0
    while elapsed < MAX_WAIT:
        resp = requests.get(status_url, headers=headers, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        status = data.get("status", "").lower()

        if status == "completed" or "pdb" in data or "output" in data:
            logger.info("Prédiction terminée.")
            return data
        if status == "failed":
            raise RuntimeError(f"Prédiction échouée côté serveur : {data}")

        logger.info("Statut : %s — nouvelle vérification dans %ds...", status, POLL_INTERVAL)
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL

    raise TimeoutError(f"La prédiction n'a pas abouti en {MAX_WAIT}s.")


def _extract_pdb(data: dict) -> str:
    pdb_keys = ("pdb", "pdb_string", "structure")

    for key in pdb_keys:
        if key in data:
            return data[key]

    outputs = data.get("outputs", [])
    if isinstance(outputs, list) and outputs:
        first_output = outputs[0]
        if isinstance(first_output, dict):
            for key in pdb_keys:
                if key in first_output:
                    return first_output[key]
            structs = first_output.get("structures_with_scores", [])
            if isinstance(structs, list) and structs:
                first_struct = structs[0]
                if isinstance(first_struct, dict):
                    for key in pdb_keys:
                        if key in first_struct:
                            return first_struct[key]

    raise RuntimeError(
        f"Impossible d'extraire le PDB. Clés : {list(data.keys())}. "
        f"Réponse (tronquée) : {json.dumps(data, indent=2)[:1000]}"
    )


def fetch_target(output_path: Path) -> Path:
    """Prédit la structure 3D via OpenFold3 et sauvegarde le PDB."""
    api_key = _get_api_key()
    response = _submit_prediction(api_key)
    data = response.json()

    if "status_url" in data or "reqId" in data:
        status_url = data.get("status_url", f"{API_URL}/status/{data.get('reqId')}")
        logger.info("Tâche asynchrone détectée. Polling sur : %s", status_url)
        data = _poll_for_result(status_url, api_key)

    pdb_content = _extract_pdb(data)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(pdb_content)
    logger.info("Structure 3D sauvegardée dans '%s' (%d octets).", output_path, len(pdb_content))
    return output_path
