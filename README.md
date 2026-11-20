# Cellular and Molecular Determinants of Human Kidney Disease

**Single-nuclei, ATAC-seq and Xenium based high-resolution spatial mapping of human kidney diseases**

Analysis and figure code of the manuscript (revised version, 2026): single-nucleus RNA-seq, single-nucleus ATAC-seq and
Xenium spatial transcriptomics of human kidney biopsies, and plasma proteomics of the UK Biobank.

## Authors

Saritas T\*, Boys C\*, Frommer L\*, Kuppe C\*, Gómez-Ochoa SA, Maryam S, van den Bosch TPP, Goepp V, van der Velde M, Buse M, Ulrich J, Menzel S, Schäfer G, Nlandu Q, Spaliviero M, Bleckwehl T, Shao X, Babler A, Krauß D, Li T, Malik MS, Mohanta D, Ijaz S, Schreibing T, Knoppert SN, Goldschmeding R, Nguyen TQ, Hoeft K, Doğan Dar E, Ziegler S, Perales-Patón J, Schneider KM, Schreibing F, Boor P, Ramirez Flores R, Reinders MEJ, Baan CC, Kranz J, Luedde T, Geisler F, Kollar M, Tesař V, Hrušková Z, Shaigan M, Costa I, Levinson RT, Smeets B, Chen Y, Schneider CV, Minnee RC, Schneider RK, Jansen J, Clahsen-van Groningen MC, Hayat S#, Saez-Rodriguez J#, Kramann R#

<img width="1262" height="583" alt="Study overview: data modalities, analysis steps and clinical endpoints" src="https://github.com/user-attachments/assets/3ce98ab0-a98f-4768-88de-49d363219c6a" />

*Study overview, as drawn for the first submission. The revised manuscript describes the Xenium panel with its 468 genes
and uses Factor 1 of the recomputed factor analysis in place of Factor 2.*

## Where to start

| you are looking for | go to |
|---|---|
| the code of a figure panel | [`Code/revision_2026/FIGURE_PANELS.md`](Code/revision_2026/FIGURE_PANELS.md): every panel of the main, Extended Data and Supplementary figures, who made it and the script that draws it |
| how an analysis was run | the README of its folder under [`Code/revision_2026/analyses/`](Code/revision_2026/analyses): purpose, run order, conda environment, inputs, outputs |
| the multicellular factor analysis (MCFA) | its own repository, [saezlab/kidney_biopsies_MCFA](https://github.com/saezlab/kidney_biopsies_MCFA) |
| the design of the Xenium gene panel | [`Code/Xenium_Panel_Design/`](Code/Xenium_Panel_Design) |

## Repository contents

| folder | contents |
|---|---|
| [`Code/revision_2026/`](Code/revision_2026) | Analyses and per-panel figure scripts of the revised manuscript: Xenium cell typing, factor 1 on Xenium, spatial neighbourhoods and cell-cell interactions, NicheCompass niches, glomerular niche 20, PEC and immune subclusters with the PEC trajectory, PT→aPT gene-regulatory network, snRNA-seq overview, Xenium pseudobulk tables for the MCFA, cell-cell interactions within Xenium niches (LIANA, CrossTalkeR), gene expression in PT and aPT cells per core, Xenium quality metrics per core, bulk RNA-seq (DESeq2). `analyses/` holds the analyses in pipeline order, `figures/` one entry script per panel, `config/` all data paths, `envs/` the conda environments. |
| [`Code/Xenium_Panel_Design/`](Code/Xenium_Panel_Design) | Gene selection for the custom 468-gene Xenium panel (GT3ZEN) and all output files of the 10x Genomics Xenium Panel Designer. |
| [`Code/Preprocessing_Xenium/`](Code/Preprocessing_Xenium) | Xenium preprocessing: assignment of cells to tissue cores, quality metrics, marker genes, first cell-type annotation. The quality-metric plots of the revised figures are in `Code/revision_2026/analyses/14_xenium_qc_per_core/`. |
| [`Code/CCI_Niches/`](Code/CCI_Niches) | Cell-cell interactions per spatial niche (LIANA, CrossTalkeR), first-submission version. The version of the revised figures is in `Code/revision_2026/analyses/12_xenium_niche_cci_liana_crosstalker/`. |
| [`Code/ukbb_processing/`](Code/ukbb_processing) | UK Biobank plasma proteomics: cohort, incident chronic kidney disease, specificity of the protein signature. |
| [`Code/cell_boundaries/`](Code/cell_boundaries), [`Code/plot_egfr_f2_on_spatial/`](Code/plot_egfr_f2_on_spatial) | Two scripts of the first submission that are still needed: `1_preprocess_xenium.py` and `1_preproc_cell_boundaries.py` write the cell polygons of the 2025 three-region Xenium object that Fig. 6g (right) reads. |

## Running the code

The scripts in `Code/revision_2026` read every data location from environment variables:

```bash
cd Code/revision_2026
source config/env.sh             # sets KB_REPO, PYTHONPATH and the data locations
export KB_DATA=/path/to/data     # override any location
python figures/Fig3/Fig3c_*.py   # writes the panel to figure_output/Fig3/
```

The variables, the conda environments and what to expect when re-drawing a panel are described in
[`Code/revision_2026/README.md`](Code/revision_2026/README.md). Heavy steps ran as SLURM jobs on CPUs; no GPU is needed.
The other folders hold the scripts and notebooks as their authors ran them; input paths are set inside each file.

## Data

The code needs the processed snRNA-seq, snATAC-seq and Xenium objects and the clinical tables of the study; see the Data
availability statement of the paper. UK Biobank participant data are available only through UK Biobank and are not part
of this repository.

## Code of the first submission

Code that was re-implemented for the revision is no longer in this repository:

- NicheCompass and PEC trajectory (`Code/Niche_analysis/`, `Code/Trajectory/`): now
  `Code/revision_2026/analyses/01_xenium_preprocessing_niches/` and
  `Code/revision_2026/analyses/09_snrna_pec_immune_subclusters/`.
- Xenium neighbourhoods and cell-cell interactions (`Code/nh_analysis/`, `Code/spatial_cci/`, `Code/utils/` and the
  plotting scripts of `Code/cell_boundaries/` and `Code/plot_egfr_f2_on_spatial/`): now
  `Code/revision_2026/analyses/03_xenium_factor1_neighbourhoods_cci/`.
- Milo scripts (`Code/milo/`): removed because no figure of the revised manuscript uses them.

## License

This project is distributed under the GNU General Public License v3.0. See [`LICENSE`](LICENSE) for full details.
