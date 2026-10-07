# bile_acid_updated_library

Code for **"The evolving bileome: a pan-repository landscape of expanded bile acid and related steroidal signatures"**
(Mannochio-Russo *et al.*).

The six spectral libraries (GNPS-C23/C24/C24-ALKAMINE/C26/C27-BILE-ACIDS-MODIFICATIONS and
GNPS-C28-BILE-ALCOHOL-MODIFICATIONS) are available in GNPS2 and in Zenodo
([10.5281/zenodo.21941572](https://doi.org/10.5281/zenodo.21941572)), together with the annotation and provenance
tables that these scripts read and write.

The first part of the pipeline (merging the MassQL results, removing duplicates, precursor checks, ReDU-based
filters, delta-mass and precursor m/z filters) is in
[library_generation_from_MassQL](https://github.com/helenamrusso/library_generation_from_MassQL). This repository
covers scripts after that step, as well as the figures.

## Figures

| Folder | Manuscript item | Main file(s) |
|---|---|---|
| `Figures/Figure1B_SuppFig2_SuppFig3_body_part_distribution` | Fig. 1B, Supp. Figs 2–3 | `body_part_distribution_BA_update_remove_artifacts.ipynb` |
| `Figures/Figure1C_SuppFig4_taxonomic_trees` | Fig. 1C, Supp. Fig. 4 | `NCBI_Taxonomy_Megan.ipynb` (count tables); R tree scripts |
| `Figures/Figure2A_vertebrate_bile_use_case` | Fig. 2A (MSV000085120) | `bile_samples_QTOF_final2.ipynb` |
| `Figures/Figure2B_zoo_feces_use_case` | Fig. 2B (MSV000086131) | `animal_zoo_diet_final.ipynb` |
| `Figures/SuppFig5_spectral_entropy` | Supp. Fig. 5 | `export_mgf_entropy_tsv.py` + `spectral_entropy.py` (entropy), `calculate_number_fragments_mgf.ipynb`, `density_entropy_all_reorder.ipynb` |
| `Figures/SuppFig7_SuppFig8_alternative_ion_forms` | Supp. Figs 7–8; ion-form flags in the libraries | `ISF_plot.ipynb`, `artifact_screen.py`, `ISF_C27_tetra.ipynb` |

## Library generation after clustering and validation (`Other_codes`)

Folders are numbered in the order they were run.

| Folder | Manuscript section | Content |
|---|---|---|
| `01_query_specificity_FDR` | Technical validation: specificity of the queries | Run the MassQL queries against the GNPS public libraries |
| `02_mgf_utilities_pre_post_falcon` | Spectral clustering | Add `CHARGE`, renumber/sort/split large MGFs for Falcon, collect Falcon outputs, add scan numbers, recombine split outputs |
| `03_library_match_classification` | Confidence-aware molecular networking | GNPS2 library matches of the clustered spectra + manual classification (steroidal yes/no) |
| `04_entropy_MN_filters_and_deltas` | Spectral entropy; network-based filtering | Synthesis-dataset filter, entropy > 7.5 filter, SpecReBoot network filter (1/3 rule), final MGFs, delta counts |
| `05_in_silico_formula_annotation` | In silico molecular formula annotation | BUDDY download, SIRIUS rerun helper, formula/delta-formula annotation |
| `06_final_library_formatting` | Data records | Annotation tables, GNPS2 library templates and final MGFs, acquisition metadata (collision energy) |


## Running the notebooks

Every notebook has a short header (what it produces, its inputs) and, where it reads project files. The use-case notebooks (Fig. 2) download
the GNPS2 feature tables themselves and cache them in `data/`.

Analyses were run with Python 3.10.18 (see `requirements.txt`). The taxonomic trees were drawn in R 4.5.1 with
taxize, ape, readr, dplyr, tidyr, purrr, RColorBrewer and fields.
