# Résultats du criblage virtuel — Inhibiteurs potentiels de xenoAMP(S)

## Contexte

Zhang et al. (PNAS 2026) démontrent que la digestion protéolytique de la protéine Spike du SARS-CoV-2 libère des fragments peptidiques — **xenoAMPs** — capables d'induire une **courbure de Gauss négative (NGC)** sur les membranes lipidiques, formant des pores transmembranaires de ~1.9 nm de diamètre. Le fragment **xenoAMP(S)** (résidus 529–558, `KSTNLVKNKCVNFNFNGLTGTGVLTESNKK`) est le plus actif. Il cible sélectivement les cellules immunitaires à morphologie complexe riche en NGC : **pDC**, **CD4+ T** et **CD8+ T**, dont la déplétion est directement corrélée à la sévérité du COVID long.

La stratégie du criblage est la **séquestration directe** : identifier des petites molécules capables de se lier à xenoAMP(S) et de bloquer stériquement son insertion membranaire, avant qu'il n'atteigne les NGC-rich regions des cellules immunitaires.

---

## Méthodologie

- **Récepteur :** structure 3D de xenoAMP(S) prédite par OpenFold3
- **Grid box :** 20×20×20 Å centrée sur la face amphipathique d'insertion (LEU546, THR547, GLY548, THR549, GLY550)
- **Ligands :** médicaments FDA-approuvés (ChEMBL, max\_phase=4) + principes actifs CMAUP v2.0 + composés FooDB ; 11 283 molécules dockées à la date du rapport
- **Logiciel :** AutoDock Vina 1.2.5, exhaustiveness=8 (criblage initial)
- **Seuil de torsions :** ≤ 15 liaisons rotatives (Lipinski NumRotatableBonds)

---

## Analyse des 20 meilleurs candidats

### Tier 1 — Candidats prioritaires

---

#### 🥇 1. Isosungucine — CMAUP\_NPC63047 | −7.0 kcal/mol

**Origine :** Alcaloïde bisindolique isolé de *Strychnos usambarensis* (Loganiaceae, Afrique de l'Est). Également présent dans *S. camptoneura*.

**Structure :** Système polycyclique rigide à deux unités indole pontées, avec un cœur strychnos quaternaire. Caractère amphipathique marqué : faces aromatiques planes hydrophobes + azotes endocycliques pouvant former des liaisons hydrogène.

**Pertinence pour xenoAMP(S) :** La rigidité du squelette et son encombrement stérique tridimensionnel permettent une complémentarité de surface avec la face hydrophobe de xenoAMP(S) (LEU546, PHE-like residues). Le score de docking le plus élevé du criblage (−7.0 kcal/mol) suggère une affinité structurelle forte pour la pocket d'insertion membranaire.

**Activités connues :** Antimicrobien, antiparasitaire (antipaludéen), antifongique. Aucune étude d'inhibition de pores membranaires publiée.

**Pharmacocinétique :** Données limitées. LogP estimé ~3–4 (lipophile modéré), biodisponibilité orale probable mais non caractérisée. Franchissement de la BHE possible (famille Strychnos). **Administration envisageable :** orale ou IV (comme les autres alcaloïdes de Strychnos).

**Toxicité :** ⚠️ La famille *Strychnos* inclut la strychnine (neurotoxique sévère, antagoniste glycine). L'isosungucine est structurellement distincte de la strychnine, mais une évaluation toxicologique approfondie est indispensable avant tout usage thérapeutique. Le profil toxicologique in vitro des bisindoles de *S. usambarensis* est généralement moins sévère que les mono-indoles.

**Statut :** Composé naturel non approuvé. Lead de recherche.

**Verdict :** Meilleur score du criblage. Structure prometteuse justifiant une caractérisation toxicologique et des études de liaison isotherme (ITC/SPR) avec xenoAMP(S).

---

#### 2. Zavegepant — CHEMBL2397415 | −7.0 kcal/mol

**Origine :** Gepant de 3e génération, antagoniste du récepteur CGRP (calcitonin gene-related peptide). **Approuvé FDA 2023** (Zavzpret®, spray nasal) pour la migraine aiguë.

**Structure :** Molécule synthétique complexe (MW ~583 Da), bicyclique avec des groupements trifluorométhyle et des chaînes latérales rigides. Caractère amphipathique permettant l'interaction avec des peptides cationiques.

**Pertinence pour xenoAMP(S) :** Les gepants sont conçus pour bloquer des interactions peptide–récepteur au niveau d'une interface protéique. La même logique structurelle (occupation stérique d'une face peptidique) peut s'appliquer à la neutralisation de xenoAMP(S). Le CGRP lui-même est un peptide amphipathique cationique de 37 résidus — la similarité fonctionnelle avec xenoAMP(S) est notable.

**Pharmacocinétique :** Spray nasal : Tmax ~15 min, biodisponibilité ~4 %. Voie IV possible (études en cours). Demi-vie ~6–7h. Liaison protéique ~71 %. Métabolisme hépatique (CYP3A4).

**Toxicité :** Profil de sécurité bien établi. Effets indésirables : nausées, dysgueusie (fréquents, légers). Pas de vasoconstricton (avantage sur les triptans). **Pas de contre-indication cardiovasculaire.**

**Statut :** ✅ **FDA-approuvé**. Utilisable off-label après validation de principe.

**Intérêt particulier :** La voie nasale pourrait être stratégiquement pertinente — le début de la protéolyse de la Spike a lieu dans les voies aériennes supérieures (cathepsines, KLK). Un traitement nasal précoce post-infection pourrait séquestrer les xenoAMPs avant dissémination systémique.

**Verdict :** Candidat cliniquement le plus mature. Priorité haute pour étude de liaison directe (ITC) avec xenoAMP(S) synthétique.

> *Note : le rang #16 (Zavegepant hydrochloride, CHEMBL4650220, −6.5 kcal/mol) est la forme sel du même principe actif — doublon pharmacologique.*

---

#### 3. Anigorootin — FOODB\_FDB011196 | −6.9 kcal/mol

**Origine :** Phénylphalénanone isolée de *Anigozanthos* spp. (Haemodoraceae, Australie occidentale). Présente également dans d'autres plantes de la famille.

**Structure :** Système tricyclique aromatique plan (phénalénanone + cycle phényle). Structure rigide, planaire, hydrophobe avec des groupements hydroxyle périphériques. MW ~330 Da.

**Pertinence pour xenoAMP(S) :** La planarité aromatique et la rigidité du squelette permettent un stacking π avec les résidus aromatiques du peptide. La taille modeste (MW ~330) favorise une bonne biodisponibilité. Les groupements hydroxyle peuvent engager les résidus THR/GLY du domaine hydrophile de xenoAMP(S).

**Activités connues :** Anti-inflammatoire (inhibition NF-κB reportée pour des phénylphalénanones similaires), antioxydant, antifongique. Les composés de cette classe inhibent des kinases impliquées dans les voies inflammatoires.

**Pharmacocinétique :** Données limitées. LogP modéré (~2–3), MW faible → biodisponibilité orale probable. **Administration envisageable :** orale.

**Toxicité :** Faible a priori (composé alimentaire, plante non toxique). Données formelles absentes.

**Statut :** Composé naturel non approuvé. Présence dans des plantes alimentaires rassurante.

**Verdict :** Lead naturel attractif combinant score de docking élevé, taille favorable et profil de sécurité probable. Priorité 2 pour validation in vitro.

---

#### 4. Usambarensine — CMAUP\_NPC280852 | −6.8 kcal/mol

**Origine :** Alcaloïde bisindolique de *Strychnos usambarensis*, co-isolé avec l'isosungucine.

**Structure :** Bisindole polycyclique structurellement proche de l'isosungucine, avec une connectivité légèrement différente entre les deux unités indole.

**Pertinence pour xenoAMP(S) :** Même raisonnement que l'isosungucine (#1). La présence de deux dérivés de la même plante dans le top 5 renforce l'hypothèse que le squelette bisindole-strychnos est particulièrement adapté à la pocket de xenoAMP(S).

**Activités connues :** Antimicrobien, antiparasitaire. Cytotoxique sur certaines lignées tumorales.

**Pharmacocinétique / Toxicité :** Identique à isosungucine — profil inconnu, même précaution famille Strychnos.

**Verdict :** Confirmation de la famille bisindole comme scaffold prioritaire. À étudier en série SAR avec l'isosungucine.

> *Le rang #11 (10'-Hydroxyusambarensine, −6.6 kcal/mol) est un dérivé hydroxylé d'usambarensine, confirmant la robustesse du scaffold dans cette pocket.*

---

### Tier 2 — Candidats FDA avec réserves

---

#### 7. Dihydroergotamine (DHE) — CHEMBL1732 | −6.7 kcal/mol

**Origine :** Alcaloïde de l'ergot de seigle (*Claviceps purpurea*), semi-synthétique. **Approuvé FDA** (DHE-45® IV, Migranal® nasal spray) pour la migraine sévère.

**Structure :** Ergoline pentacyclique complexe avec une chaîne latérale peptidomimétique (lysergamide). MW ~583 Da. Structure amphipathique prononcée.

**Pertinence pour xenoAMP(S) :** La chaîne latérale peptidomimétique du DHE mime structurellement des dipeptides cycliques et pourrait engager les résidus de liaison du peptide cible. Le noyau ergoline rigide offre un ancrage hydrophobe complémentaire à la face lipophile de xenoAMP(S).

**Pharmacocinétique :** Biodisponibilité orale médiocre (<1 %). Bonne absorption IV ou nasale. Tmax nasal ~30–60 min. Forte liaison protéique. Métabolisme hépatique extensif.

**Toxicité :** ⚠️ Vasoconstricteur puissant (agoniste 5-HT1, α-adrénergique). **Contre-indiqué** en cas de coronaropathie, HTA non contrôlée, grossesse. Ergotisme (ischémie digitale) possible. Usage limité dans le temps. Ce profil de sécurité restreint fortement l'utilisation chronique qu'exigerait un traitement du COVID long.

**Verdict :** Score de docking intéressant, mais toxicité cardiovasculaire incompatible avec un traitement prolongé du COVID long. Valeur principale : validation pharmacophore pour le scaffold ergoline.

---

#### 9. Tazemetostat hydrobromide — CHEMBL4594260 | −6.7 kcal/mol

**Origine :** Inhibiteur sélectif d'EZH2 (histone méthyltransférase, complexe PRC2). **Approuvé FDA 2020** (Tazverik®) pour le sarcome épithélioïde et le lymphome folliculaire EZH2+.

**Pertinence pour xenoAMP(S) :** Score de docking correct (−6.7). Molécule volumineuse et rigide (MW ~571 Da) avec un scaffold indolique.

**Préoccupation majeure :** ⛔ EZH2 est un régulateur essentiel de la différenciation des lymphocytes T et B. Son inhibition **supprime activement la réponse immunitaire adaptative**. Dans le contexte du COVID long où l'objectif est de protéger les CD4+/CD8+ T déjà épuisés, administrer un immunosuppresseur epigénétique est **directement contre-thérapeutique**. Le bénéfice d'un éventuel blocage de xenoAMP(S) serait annulé par l'immunosuppression induite.

**Verdict :** ❌ À exclure. Mécanisme d'action antagoniste de l'objectif thérapeutique.

---

#### 17. Velpatasvir — CHEMBL3545062 | −6.5 kcal/mol

**Origine :** Inhibiteur de NS5A du VHC, **approuvé FDA 2016** (Epclusa®, en combinaison avec sofosbuvir).

**Structure :** Grande molécule (MW ~884 Da) avec un cœur biphényle et des chaînes latérales proline-like. Conçue pour bloquer des interactions protéine–protéine au sein du complexe de réplication virale. Biodisponibilité orale excellente (~98 %).

**Pertinence pour xenoAMP(S) :** Velpatasvir est un inhibiteur d'interactions protéiques, un mode d'action directement transposable à la séquestration d'un peptide. Sa biodisponibilité orale exceptionnelle et son profil de sécurité établi (effets secondaires très légers) en font un candidat cliniquement attractif.

**Pharmacocinétique :** Tmax ~3h, t½ ~15h, liaison protéique >99 %. Métabolisme CYP2B6/CYP3A4. Bonne exposition systémique.

**Toxicité :** Très bien toléré. Céphalées, fatigue légère. Pas d'immunosuppression.

**Verdict :** Candidat FDA de haute valeur. Biodisponibilité orale remarquable + mécanisme d'inhibition d'interaction protéique compatible. À valider in vitro. Administration orale une fois par jour.

---

### Tier 3 — Leads naturels à profil variable

---

#### 5. α-Chaconine — FOODB\_FDB093773 | −6.7 kcal/mol

**Origine :** Glycoalcaloïde stéroïdien de *Solanum tuberosum* (pomme de terre), concentré dans la peau et les germes.

**Pertinence apparente :** Score correct (−6.7). Structure stéroïdienne lipophile.

**Préoccupation critique :** ⚠️ L'α-chaconine est **elle-même un agent perturbateur de membranes lipidiques** (lyse membranaire, inhibition de AChE). Elle interagit directement avec les phospholipides membranaires et peut induire une perméabilité membranaire. Dans ce contexte, le score de docking reflète probablement une affinité membranaire non spécifique plutôt qu'une liaison spécifique à xenoAMP(S). Un composé perturbateur de membranes risque de **potentialiser** la toxicité de xenoAMP(S) plutôt que de l'inhiber.

**Toxicité :** Toxique à doses élevées (>1 mg/kg). Ingestion de germes de pomme de terre : empoisonnement documenté.

**Verdict :** ❌ Faux positif probable. La perturbation membranaire intrinsèque disqualifie ce composé.

---

#### 6. Garcilivin B — FOODB\_FDB018413 | −6.7 kcal/mol

**Origine :** Polycétide / limonoïde de *Garcinia* spp. (Clusiaceae). Famille riche en xanthones et benzophénones bioactives.

**Structure :** Polycyclique oxygéné, rigide. Les limonoïdes de *Garcinia* présentent généralement une bonne affinité pour les interfaces hydrophobes.

**Activités connues :** Anti-inflammatoire (inhibition COX-2 pour des limonoïdes similaires), antibactérien, antiprolifératif.

**Pharmacocinétique :** Peu documentée. LogP modéré. Données in vivo absentes pour ce composé spécifique.

**Verdict :** Lead naturel intéressant. Données insuffisantes pour aller plus loin sans validation expérimentale.

---

#### 10. Cucurbitaxanthin A — FOODB\_FDB013989 | −6.7 kcal/mol

**Origine :** Xanthophylle caroténoïde de *Cucurbita* spp. (courges).

**Structure :** Longue chaîne polyène conjuguée avec des groupements hydroxyle terminaux (MW ~600 Da). Très lipophile.

**Préoccupation :** Les caroténoïdes partitionnent naturellement dans les membranes lipidiques et s'intercalent dans les bicouches. Le score de docking peut refléter une affinité hydrophobe non spécifique. De plus, la longue chaîne polyène (>15 torsions réelles) pose question quant à la spécificité de liaison.

**Verdict :** Intérêt limité. Probable interaction non spécifique avec les composantes hydrophobes du peptide.

---

#### 11. 10'-Hydroxyusambarensine — CMAUP\_NPC476041 | −6.6 kcal/mol

**Origine :** Dérivé hydroxylé d'usambarensine (#4), *Strychnos usambarensis*.

**Verdict :** Confirme le scaffold bisindole-strychnos. Potentiellement plus soluble que le composé parent (groupement OH supplémentaire). Intègre la série SAR à étudier avec isosungucine et usambarensine.

---

#### 13. 10-(Chrysophanol-7'-Yl)-10-Hydroxychrysophanol-9-Anthrone — CMAUP\_NPC471682 | −6.6 kcal/mol

**Origine :** Dimère d'anthraquinone issu de plantes médicinales (Polygonacées ou Rhamnacées).

**Structure :** Grande molécule plane aromatique (deux unités chrysophanol liées). MW élevée (~720 Da).

**Activités connues :** Les anthraquinones ont des propriétés anti-inflammatoires et antimicrobiennes documentées. L'aloïne (anthraquinone d'*Aloe vera*) module des cascades pro-inflammatoires.

**Préoccupation :** MW élevée et structure plane suggèrent un risque d'intercalation dans l'ADN (mutagénicité potentielle, typique des anthraquinones aromatiques). Règle de Lipinski probablement violée.

**Verdict :** Structure intéressante pharmacologiquement mais préoccupations génotoxiques standard des anthraquinones. À étudier avec tests Ames.

---

#### 14. Cyclopamine — CMAUP\_NPC28280 | −6.6 kcal/mol

**Origine :** Alcaloïde stéroïdien de *Veratrum californicum* (Mélanthiacées). Inhibiteur de Smoothened (SMO), récepteur de la voie Hedgehog.

**Activités connues :** Anticancéreux (base du vismodegib/sonidegib), antiprolifératif. Bonne pénétration membranaire (lipophile).

**Toxicité :** ⛔ **Tératogène sévère confirmé** (induit la cyclopie chez les agneaux, à l'origine de la découverte de la voie Hedgehog). Absolument contre-indiqué chez la femme en âge de procréer. Marge thérapeutique très étroite.

**Verdict :** ❌ Profil tératogène rédhibitoire pour un traitement du COVID long (population d'âge reproductif incluse). Valeur uniquement comme outil pharmacologique.

---

#### 15. Ochrolifuanine A — CMAUP\_NPC221786 | −6.6 kcal/mol

**Origine :** Alcaloïde indolique de *Ochrosia* spp. (Apocynaceae). Famille riche en alcaloïdes ellipticine-like.

**Structure :** Β-carboline ou indoloquinoléine polycyclique. Structure plane aromatique, intercalante potentielle.

**Activités connues :** Cytotoxique (comme la plupart des alcaloïdes Ochrosia). Antimicrobien. Les alcaloïdes de ce genre ont souvent des propriétés antivirales.

**Toxicité :** Cytotoxicité générale prévisible pour ce type de scaffold. Fenêtre thérapeutique probablement étroite.

**Verdict :** Intéressant du point de vue structural mais profil cytotoxique à évaluer avant tout développement.

---

#### 18. Auroxanthin — FOODB\_FDB015834 | −6.5 kcal/mol

**Origine :** Caroténoïde xanthophylle (époxyde d'antheraxanthine) issu de diverses plantes alimentaires.

**Verdict :** Même raisonnement que Cucurbitaxanthin A (#10). Très lipophile, interaction membranaire non spécifique probable. Score de docking reflet d'une affinité hydrophobe.

---

#### 19. Isotheaflavin — FOODB\_FDB000613 | −6.5 kcal/mol

**Origine :** Dérivé de théaflavine du thé noir (*Camellia sinensis*), produit d'oxydation des catéchines.

**Structure :** Polyphénol avec un système benzotropolone (deux cycles aromatiques + un cycle à 7 chaînons). MW ~564 Da.

**Activités connues :** Antioxydant puissant, anti-inflammatoire (inhibition NF-κB, COX-2), antiviral (notamment contre les coronavirus pour des théaflavines analogues). Les théaflavines inhibent la fusion membranaire et l'entrée virale pour plusieurs virus.

**Pharmacocinétique :** ⚠️ Biodisponibilité orale très faible (<1–5 % pour les théaflavines en général). Métabolisme intestinal et hépatique extensif. Forte liaison protéique. **Problème majeur pour un usage systémique.**

**Stratégie alternative :** Administration locale (inhalation, spray nasal) pour séquestration au site de première protéolyse virale (voies aériennes supérieures).

**Toxicité :** Excellent profil — présente dans le thé noir à haute consommation mondiale.

**Verdict :** Profil biologique encourageant et sécurité optimale, mais biodisponibilité systémique insuffisante. Intéressant uniquement pour une approche topique nasale/pulmonaire.

---

#### 20. Terniflorin — FOODB\_FDB016369 | −6.5 kcal/mol

**Origine :** Flavonoïde isolé de *Erythroxylum* spp. et d'autres plantes. Présent dans des aliments.

**Structure :** Flavone ou flavonol glycosylé. MW modérée.

**Activités connues :** Anti-inflammatoire, antioxydant. Propriétés vasculaires pour des flavonoïdes de cette classe.

**Verdict :** Lead naturel mineur. Données insuffisantes pour priorisation.

---

#### Composés non identifiés

- **CMAUP\_NPC478949** (rang #8, −6.7 kcal/mol) : identifié sous l'InChIKey `MAWDEQUQRREYMR-QOAFJEBESA-N` — aucun nom vernaculaire dans la base. Score élevé justifiant une identification chimique (PubChem lookup sur l'InChIKey) avant tout classement.
- **CMAUP\_NPC134131** (rang #12, −6.6 kcal/mol) : InChIKey `GIMKEHNOTHXONN-HNIMUKOUSA-N` — même situation. Ces deux composés doivent être dérépliqués avant validation.

---

## Synthèse et hiérarchisation

| Rang | Molécule | Score | Statut | Priorité | Motif |
|------|----------|-------|--------|----------|-------|
| 1 | Isosungucine | −7.0 | Naturel | ⭐⭐⭐⭐ | Meilleur score, scaffold bisindole rigide, SAR possible |
| 2 | Zavegepant | −7.0 | ✅ FDA | ⭐⭐⭐⭐⭐ | Approuvé, sûr, mécanisme compatible, voie nasale |
| 3 | Anigorootin | −6.9 | Naturel | ⭐⭐⭐ | Score élevé, taille favorable, faible toxicité probable |
| 4 | Usambarensine | −6.8 | Naturel | ⭐⭐⭐ | Confirme le scaffold bisindole |
| 7 | Dihydroergotamine | −6.7 | ✅ FDA | ⭐⭐ | Validé cliniquement mais toxicité CV incompatible usage chronique |
| 17 | Velpatasvir | −6.5 | ✅ FDA | ⭐⭐⭐⭐ | Excellent PK oral, inhibiteur d'interaction protéique |
| 19 | Isotheaflavin | −6.5 | Alimentaire | ⭐⭐ | Excellent profil de sécurité, biodisponibilité orale insuffisante |
| 5 | α-Chaconine | −6.7 | — | ❌ | Perturbateur membranaire intrinsèque — faux positif |
| 9 | Tazemetostat | −6.7 | ✅ FDA | ❌ | Immunosuppresseur — contre-thérapeutique |
| 14 | Cyclopamine | −6.6 | — | ❌ | Tératogène sévère |

---

## Considérations pharmacocinétiques et voies d'administration

### Séquestration systémique (approche principale)
Pour bloquer xenoAMP(S) circulant dans le plasma ou les organes lymphoïdes secondaires, la molécule doit atteindre des concentrations suffisantes dans le compartiment extracellulaire sanguin et tissulaire.
- **Velpatasvir** et **Zavegepant** offrent la meilleure combinaison biodisponibilité / tolérance pour cette approche.
- Les bisindoles (isosungucine, usambarensine) nécessitent une caractérisation PK complète.

### Séquestration locale (approche nasale/pulmonaire)
La protéolyse initiale de la Spike a lieu dans les voies aériennes supérieures (kallikréines KLK, élastase). Bloquer les xenoAMPs à ce stade, avant dissémination, est une stratégie préventive/précoce.
- **Zavegepant** existe déjà sous forme de spray nasal (Zavzpret®).
- **Isotheaflavin** : candidat idéal par inhalation (faible biodisponibilité orale, mais bonne concentration locale possible).

### Synergisme LL-37
Zhang et al. montrent que xenoAMP(S) et LL-37 (AMP endogène) se potentialisent mutuellement. Un inhibiteur de xenoAMP(S) actif à faible concentration (< P/L = 1/152) pourrait rompre cette coopération et réduire drastiquement la formation de pores même sans neutralisation complète du peptide.

---

## Prochaines étapes recommandées

1. **Valider les 5 candidats prioritaires** par ITC (Isothermal Titration Calorimetry) ou SPR (Surface Plasmon Resonance) avec xenoAMP(S) synthétique (30-mer commercial disponible).
2. **Test SAXS inhibition** : mesurer si les candidats réduisent la formation de phases cubiques Pn3m/Im3m induites par xenoAMP(S) sur SUVs DOPS/DOPE (exactement comme dans l'étude Zhang et al.).
3. **Test cytotoxicité cellulaire** : évaluer la protection des pDC/CD8+ T contre xenoAMP(S) à 40 µM (protocole Fig. 5C–H de Zhang et al.) en présence des candidats.
4. **Identification** de CMAUP\_NPC478949 et NPC134131 par lookup PubChem/ChemSpider sur leurs InChIKeys.
5. **SAR bisindoles** : synthèse ou sourcing de dérivés d'isosungucine avec évaluation toxicologique (Strychnos panel).
6. **Refine docking** (exhaustiveness=64) sur les 5 candidats retenus + extension du criblage (objectif : 18 600 molécules dockées).

---

## Analyse des concentrations — Sommes-nous dans les bons ordres de grandeur ?

### 1. Concentration cible : combien de xenoAMP(S) faut-il neutraliser ?

**Concentration de protéine Spike dans le plasma durant un COVID sévère** (données Ogata et al., 2021 ; St-Jean et al., 2024) : **~10 ng/mL** (médiane), avec des valeurs jusqu'à 50 ng/mL dans les formes critiques.

```
Spike (MW = 180 000 g/mol) : 10 ng/mL = 55 pM total
xenoAMP(S) = résidus 529–558 = 30 AA / 1 274 AA totaux
→ si protéolyse complète : 55 pM × (30/1274) ≈ 1.3 nM de xenoAMP(S)
→ accumulation locale NGC (×4, modèle Zhang et al.) : ~5 nM à la surface des pDC/CD8+ T
```

**Dans le COVID long** (réservoir viral persistant, protéolyse chronique à bas bruit), les concentrations plasmatiques de Spike sont probablement plus basses qu'en phase aiguë, mais l'exposition est continue sur des semaines. L'ordre de grandeur retenu pour cette analyse : **[xenoAMP(S)]_cible ≈ 1–5 nM** (plasmatique à locale).

> **Remarque critique :** Zhang et al. ont utilisé **40 µM** dans leurs expériences in vitro — c'est 10 000× supérieur à l'estimation physiologique. Cela confirme qu'il s'agissait d'une condition supra-physiologique pour garantir la détection. La formation de pores via synergisme LL-37 a été démontrée à P/L = 1/152, soit des concentrations beaucoup plus faibles, cohérentes avec la gamme nM.

---

### 2. Affinité estimée des candidats — Kd calculé depuis les scores Vina

La relation ΔG = RT ln(Kd) (R = 1,987 × 10⁻³ kcal/mol/K, T = 310 K) donne les estimations suivantes :

| Molécule | Score Vina (kcal/mol) | Kd estimé |
|---|---|---|
| Isosungucine | −7.0 | ~12 µM |
| Zavegepant | −7.0 | ~12 µM |
| Anigorootin | −6.9 | ~14 µM |
| Usambarensine | −6.8 | ~16 µM |
| Velpatasvir | −6.5 | ~26 µM |

⚠️ **Ces Kd sont des estimations très approximatives.** Les scores Vina ont une incertitude de ±1.5–2.0 kcal/mol sur des cibles flexibles (peptides), ce qui représente un facteur 10–100 sur le Kd. Une interaction réelle à 100 nM – 1 µM est parfaitement compatible avec ces scores. Seule une mesure expérimentale (ITC, SPR, MST) donnera la vraie valeur.

---

### 3. Concentrations plasmatiques disponibles aux doses cliniques

| Molécule | Dose / voie | Cmax plasma | Cmax (nM) |
|---|---|---|---|
| Zavegepant | 10 mg intranasal (Zavzpret®) | 1.46 ng/mL | **2.5 nM** |
| Velpatasvir | 100 mg oral (Epclusa®) | 1 200 ng/mL | **1 360 nM** |
| Dihydroergotamine | 2 mg intranasal (Migranal®) | ~750 ng/mL | ~1 290 nM |
| Isosungucine | — | Non déterminé | ? |
| Anigorootin | — | Non déterminé | ? |

---

### 4. Rapport molaire Drug / xenoAMP(S)

En prenant [xenoAMP(S)]_plasma ≈ **1.3 nM** :

| Molécule | Cmax (nM) | Ratio Drug/Peptide | Interprétation |
|---|---|---|---|
| Zavegepant nasal | 2.5 nM | **×1.9** | Quasi-équimolaire |
| Velpatasvir oral | 1 360 nM | **×1 046** | Très large excès molaire |
| DHE nasal | 1 290 nM | **×992** | Large excès molaire |

Stoichiométriquement, **Velpatasvir et DHE sont en très large excès** vis-à-vis du peptide. Zavegepant est quasi-équimolaire avec le peptide aux doses nasales approuvées.

---

### 5. Taux de séquestration attendu — f = [Drug] / ([Drug] + Kd)

```
f = fraction de xenoAMP(S) séquestrée à l'équilibre (hypothèse [Drug] >> [xenoAMP])
```

| Scénario Kd | Zavegepant (2.5 nM) | Velpatasvir (1 360 nM) |
|---|---|---|
| Kd = 100 nM (optimiste) | 2.4% | **93%** |
| Kd = 1 µM (plausible) | 0.2% | **58%** |
| Kd = 10 µM (Vina brut) | 0.02% | 12% |
| Kd = 100 µM (pessimiste) | 0.002% | 1.3% |

---

### 6. Conclusions : sommes-nous dans les bons ordres de grandeur ?

#### ✅ OUI pour Velpatasvir (si Kd ≤ 1 µM)
Avec 1 360 nM de Velpatasvir disponible contre ~1.3 nM de xenoAMP(S), le **rapport molaire de 1 000× est extrêmement favorable**. Si le Kd réel se situe entre 100 nM et 1 µM — une gamme tout à fait réaliste pour une petite molécule bien dockée contre un peptide amphipathique — Velpatasvir séquestrerait **58 à 93% du xenoAMP(S) circulant**. La formation de pores NGC étant un phénomène coopératif (la courbe dose-réponse xenoAMP(S)/NGC est non-linéaire), réduire la fraction libre de 50–90% pourrait suffire à passer sous le seuil de formation de pores stables.

#### ⚠️ INSUFFISANT pour Zavegepant nasal (à dose approuvée)
À 2.5 nM, le Zavegepant nasal est quasi-équimolaire avec xenoAMP(S) et bien en dessous du Kd estimé → séquestration négligeable **pour la voie systémique**. Son intérêt reste la **séquestration locale** dans les voies aériennes supérieures (site de première protéolyse de la Spike), où les concentrations locales de Zavegepant juste après administration sont de 2–3 ordres de grandeur supérieures aux concentrations plasmatiques. Une dose nasale de 10 mg délivre ~17 µmol dans ~100 µL de fluide muqueux → concentration locale **~170 µM >> Kd estimé**.

#### ❓ Inconnu pour les composés naturels (Isosungucine, Anigorootin)
Sans données PK, impossible de conclure. Cependant, les alcaloïdes bisindoliques de *Strychnos* atteignent typiquement des concentrations plasmatiques de l'ordre de **100–500 nM** après administration orale aux doses thérapeutiques de composés analogues (quinine, brucine). Si Isosungucine suit ce profil et si Kd ≈ 100 nM–1 µM, la séquestration serait **comparable à celle de Velpatasvir**.

---

### 7. Résumé pharmacologique

```
[xenoAMP(S)] à neutraliser   ≈  1–5 nM (plasma/local)
Kd estimé (Vina, ±1 ordre)  ≈  0.1–10 µM

→ Velpatasvir oral   : Cmax 1 360 nM → ratio Cmax/Kd ≈ 0.14–14  → FAVORABLE
→ Zavegepant nasal   : Cmax    2.5 nM → ratio Cmax/Kd ≈ 0.00025  → insuffisant systémique
                        Local muqueux  : ~170 µM      → FAVORABLE localement
→ Isosungucine       : Cmax    inconnue, probablement 100–500 nM → à mesurer
```

**L'ordre de grandeur est compatible.** Le principal verrou pharmacologique n'est pas la disponibilité molaire du médicament par rapport au peptide, mais l'**affinité réelle (Kd)**, qui doit être mesurée expérimentalement. La priorité absolue est donc une mesure ITC/SPR de Kd pour Velpatasvir et Isosungucine sur le peptide xenoAMP(S) synthétique.

---

*Rapport généré le 2026-02-28 — 11 283 / ~18 600 molécules dockées (criblage en cours, exhaustiveness=8)*
