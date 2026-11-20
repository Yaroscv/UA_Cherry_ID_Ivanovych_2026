# UA_Cherry_ID_Ivanovych_2026

Data and analysis code for:

> **Ivanovych Y** (2026) Germplasm characterization and genetic structure of Ukrainian sweet cherry (*Prunus avium* L.) assessed by microsatellite and *S*-locus markers. *Genetic Resources and Crop Evolution*. https://doi.org/10.1007/s10722-026-02994-3

If you use these data or scripts, please cite the article (see [Citation](#citation)).

## Abstract

Ukrainian sweet cherry (*Prunus avium* L.) germplasm is among the larger national collections in Europe, yet its cultivar identities and recorded pedigrees have not been tested against molecular evidence, and the collections holding it are at risk. Ninety-four accessions were genotyped at 18 microsatellite (SSR) loci and three self-incompatibility (*S*) locus markers, with allele sizes at six shared loci calibrated to the European reference scale to permit direct cross-country comparison — the first such calibration applied to this germplasm. Identity testing resolved five duplicate cultivar names into three homonym pairs, two clonal pairs, and identified two additional accessions as mislabeled propagations, leaving 90 unique genotypes. Comparison with European reference collections revealed a reciprocal exchange of identity between the founder cultivars 'Drogans gelbe Knorpelkirsche' and 'Dönissens gelbe Knorpelkirsche' in both Ukrainian and European holdings. *S*-genotypes were determined for the 94 and 17 duplicate accessions, yielding 19 *S*-alleles, resolving six of eleven conflicts in the published record, and providing first records for 33 Ukrainian cultivars and selections. Parentage testing confirmed three pedigrees in full and one parent in 25 further pedigrees, rejected 13, and proposed four new assignments; 'Drogans gelbe Knorpel.' was confirmed as the parent of 17 pedigrees, and two distinct genotypes were shown to have been propagated under the name 'Valerii Chkalov'. Genetic structure reflected a narrow founder base. Cultivar labels cannot be assumed correct without verification, and conservation necessitates verification through convergent SSR, *S*-locus, and parentage evidence, with the individual tree code serving as the primary unit of record.

**Keywords:** cultivar identification · genetic diversity · germplasm management · microsatellite markers · parentage analysis · *Prunus avium* · self-incompatibility · Ukraine

## Repository structure

```
UA_Cherry_ID_Ivanovych_2026/
├── data/
│   ├── UKR_Sw_Cherry_ID_18_SSR_3_S-locus_data.txt   # 94 accessions × 18 SSR + 3 S-locus markers (diploid, two columns per locus)
│   ├── metadata.csv / metadata.txt                  # sample codes, original and corrected names, clone status, origin, S-genotype, CIG
│   ├── pedigree_v2.csv                              # recorded pedigrees with parent confirmation status (input for Fig. 1)
│   ├── BIC_table.csv                                # DAPC BIC values, K = 1–10
│   ├── DAPC_K3_Qscores_ed.csv                       # DAPC membership probabilities, K = 3 (n = 90)
│   ├── structure_results_K1-11.zip                  # STRUCTURE 2.3.4 output, K = 1–11 × 20 replicates
│   ├── CLUMPP/                                      # CLUMPP 1.1.2 consensus Q-matrices (K = 2, 4, 7) and paramfile
│   └── DARwin_PCoA/                                 # DARwin 6 dissimilarity / factorial-analysis output (internal use; not reported)
├── R_UA_Cherry_ID/
│   ├── 01.UA_Cherry_ID_DAPC.qmd                     # v3.1.0 — DAPC (adegenet): identity curation (n = 90), n.pca, replicated BIC, K = 3 / 4, xvalDapc
│   ├── 02.UA_Cherry_ID_Paternity.qmd                # v6.1.0 — direct-match and parentage search vs the ECPGR reference (Familias-style LR)
│   ├── 03.UA_Cherry_ID_Harmonise_ECPGR.py           # v1.0.1 — calibration of six SSR loci to the ECPGR scale (Supplementary Table S3)
│   └── data/                                        # inputs for 01–03: UKR genotypes and metadata, harmonised UKR query (90 × 6 loci), ECPGR reference files
├── UA_Cherry_ID_Ivanovych_2026.ipynb                # Python: structureHarvester / Evanno ΔK, CLUMPP indfile, bar plots, BIC plot, pedigree figure
├── LICENSE
└── README.md
```

## Data notes

- **Sample codes** (e.g. `DRO/1`, `VAL/3`, `#110`) identify individual trees and are the primary unit of record. Cultivar names in the data files are the labels under which material was received; corrected identities are given in `metadata.csv` (`posterior_name`, `harmonised_name`).
- **Missing data** are coded `*` in the genotype file and `0` in the harmonised query file.
- **Clone correction:** four duplicate accessions (`OTR/2`, `TAY/1`, `ZOR/1`, `OLE/2`) were removed before structure analyses, leaving n = 90.
- *S*-locus markers (PaConsI, PaConsII, F-box) were excluded from DAPC and STRUCTURE because they violate neutrality and linkage-equilibrium assumptions.
- Allele sizes at six SSR loci shared with the ECPGR panel were calibrated to the ECPGR reference scale (Supplementary Table S3 of the article); the harmonised profiles are in `R_UA_Cherry_ID/data/UKR_Cherry_90_harmonised_query.txt`. `03.UA_Cherry_ID_Harmonise_ECPGR.py` regenerates Table S3 and the 94 × 6 harmonised matrix from the raw genotypes; it checks its own counts against the published Table S3 (10 regression checks) and stops with an error if any differs.

## Code and versioning

Each script carries its version inside the file (`script_version` in the `.qmd` files, `VERSION` in the `.py` file), prints it when run, and writes it to a run record next to its outputs (`DAPC_run_info.txt`, `<out_prefix>_run_info.txt`, `manifest.json` with input and output checksums). A CHANGELOG at the top of each file lists every version. Version numbers: MAJOR = method change (results may differ), MINOR = new options or documentation (results unchanged), PATCH = fixes. File names carry no version; earlier versions are in the git history.

| Script | Version | Run | Outputs |
|---|---|---|---|
| `01.UA_Cherry_ID_DAPC.qmd` | 3.1.0 | RStudio / Quarto, *Run All*; set `proj_root` | `results/cherry_SSR_cc/` |
| `02.UA_Cherry_ID_Paternity.qmd` | 6.1.0 | RStudio / Quarto; set `proj_root` and `run_profile` (`ECPGR_vs_ECPGR` or `UKR_vs_ECPGR`) | `results/cherry_ID/` |
| `03.UA_Cherry_ID_Harmonise_ECPGR.py` | 1.0.1 | `python 03.UA_Cherry_ID_Harmonise_ECPGR.py` (Python ≥ 3.8, openpyxl) from `R_UA_Cherry_ID/` | `results/harmonise/` |

The Jupyter notebook (`UA_Cherry_ID_Ivanovych_2026.ipynb`) is kept as run for the figures.

## Software

| Analysis | Software |
|---|---|
| Parentage, Ukrainian panel (94 × 18 SSR) | Familias 3.4.1 |
| Allele-size calibration to the ECPGR scale (Table S3) | Python 3, openpyxl (`03.UA_Cherry_ID_Harmonise_ECPGR.py`) |
| Direct-match and parentage vs. ECPGR reference | R / Quarto (`02.UA_Cherry_ID_Paternity.qmd`) |
| Bayesian clustering | STRUCTURE 2.3.4; structureHarvester (Evanno ΔK); CLUMPP 1.1.2 |
| DAPC | R, adegenet (`01.UA_Cherry_ID_DAPC.qmd`) |
| NJ dendrogram, factorial analysis | DARwin 6 |
| Figures | Python (matplotlib, Graphviz), Google Colab |

## Laboratory work

SSR and *S*-locus genotyping was carried out at INRAE UMR 1332 BFP (Bordeaux, France) in 2016 during a Short-Term Scientific Mission of COST Action FA1104. Reference DNA was provided by the Centre de Ressources Génétiques (CRG) *Prunus* and INRAE Bordeaux.

## Third-party data

The ECPGR sweet cherry reference files in `R_UA_Cherry_ID/data/` are derived from the supplementary data of the two publications below and are included only so the analyses can be reproduced. They remain the work of their original authors; please cite these publications if you use them.

- Barreneche T, Cárcamo de la Concepción M, Blouin-Delmas M, et al. (2021) SSR-based analysis of genetic diversity and structure of sweet cherry (*Prunus avium* L.) from 19 countries in Europe. *Plants* 10:1983. https://doi.org/10.3390/plants10101983
- Ordidge M, Litthauer S, Venison E, Blouin-Delmas M, Fernandez-Fernandez F, Höfer M, Kägi C, Kellerhals M, Marchese A, Mariette S, Nybom H, Giovannini D (2021) Towards a joint international database: alignment of SSR marker data for European collections of cherry germplasm. *Plants* 10:1243. https://doi.org/10.3390/plants10061243

## Citation

> Ivanovych Y (2026) Germplasm characterization and genetic structure of Ukrainian sweet cherry (*Prunus avium* L.) assessed by microsatellite and *S*-locus markers. Genet Resour Crop Evol. https://doi.org/10.1007/s10722-026-02994-3

## Licence

© 2026 Yaroslav Ivanovych. The data and code in this repository are licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE). You may share and adapt them for any purpose, provided you credit the author by citing the article above. Third-party data are credited to their original authors (see above).

## Contact

Yaroslav Ivanovych — Department of Forest Mycology and Plant Pathology, Swedish University of Agricultural Sciences (SLU), Uppsala, Sweden; Institute of Horticulture, NAAS of Ukraine, Kyiv, Ukraine
yaroslav.ivanovych@slu.se · ORCID [0000-0002-4066-4106](https://orcid.org/0000-0002-4066-4106)
