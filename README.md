# CONTEXTE DU PROJET : High-Throughput Virtual Screening (HTVS) in silico pour le repositionnement de médicaments (Covid Long)

## 1. Contexte Biologique et Pathologique
Ce projet vise à trouver une molécule (médicament approuvé par la FDA) capable de neutraliser un mécanisme destructeur identifié dans les cas de Covid long sévère. Le mécanisme est détaillé dans l'étude PNAS jointe dans ce répertoire (Zhang et al., 2026).
* L'étude démontre que la digestion protéolytique de la protéine Spike du SARS-CoV-2 par l'hôte produit des fragments peptidiques toxiques appelés "xenoAMPs".
* Ces xenoAMPs ont la capacité physique d'induire une courbure de Gauss négative (NGC) sur les membranes lipidiques, ce qui leur permet de former des pores transmembranaires destructeurs.
* Ce mécanisme cible spécifiquement et décime les cellules immunitaires présentant des morphologies complexes riches en NGC, comme les cellules dendritiques plasmacytoïdes (pDC), les lymphocytes T CD4+ et les lymphocytes T CD8+.
* La destruction des CD8+ empêche le contrôle des virus latents (comme l'EBV) et provoque une hyper-inflammation chronique.

## 2. La Cible Moléculaire (Le Récepteur)
Nous ciblons un fragment spécifique identifié dans l'étude : le **xenoAMP(S)**.
* **Position initiale :** Acides aminés 529 à 558 de la protéine Spike du SARS-CoV-2.
* **Séquence FASTA (30 AA) :** `KSTNLVKNKCVNFNFNGLTGTGVLTESNKK`
* **Zone active (Binding Pocket cible) :** L'étude montre que ce peptide s'insère dans la bicouche lipidique et stabilise les pores grâce à des résidus hydrophiles spécifiques (THR547, GLY548, THR549, GLY550) et un résidu hydrophobe (LEU546). C'est cette région précise que notre ligand devra cibler et bloquer (encombrement stérique).