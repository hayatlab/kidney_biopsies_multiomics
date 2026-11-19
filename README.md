# Cellular and Molecular Determinants of Human Kidney Disease
# Single-nuclei, ATAC-seq and Xenium based high-resolution spatial mapping of human kidney diseases

Saritas T1,2*, Boys C3*, Kuppe C1*, Frommer L1, Maryam S1, Goepp V1, SA Gómez-Ochoa4, Bleckwehl T1, Xiaohang Shao1, Malik MS1, Max van der Velde1, Mohanta D1,2, Ijaz S1, Schreibing T1, Knoppert SN5, Goldschmeding R5, Nguyen TQ5, Konrad Hoeft1, van den Bosch TPP6, Elif Doğan Dar3, Ziegler S1, Perales-Patón J3, Schneider KM7, Schreibing F1, Schneider RK8,9,10, Boor P11, Ramirez Flores R12, Levinson RT3,4, Schneider CV7, Minnee R13, Jansen J1, Clahsen-van Groningen MC1,14, Hayat S1#, Saez-Rodriguez J3,12#, Kramann R1,2#

<img width="1262" height="583" alt="image" src="https://github.com/user-attachments/assets/3ce98ab0-a98f-4768-88de-49d363219c6a" />

## Code in this repository

| folder | contents |
|---|---|
| [`Code/revision_2026/`](Code/revision_2026) | Analyses and per-panel figure scripts of the revised manuscript: Xenium cell typing, factor 1 on Xenium, spatial neighbourhoods and cell-cell interactions, NicheCompass niches, glomerular niche 20, PEC and immune subclusters with the PEC trajectory, PT→aPT gene-regulatory network. [`FIGURE_PANELS.md`](Code/revision_2026/FIGURE_PANELS.md) lists every figure panel, who made it and the script that draws it. |
| [`Code/Xenium_Panel_Design/`](Code/Xenium_Panel_Design) | Design of the custom 468-gene Xenium panel (GT3ZEN) and all output files of the 10x Genomics Xenium Panel Designer. |
| [`Code/Preprocessing_Xenium/`](Code/Preprocessing_Xenium) | Xenium preprocessing: core assignment, quality metrics, first cell-type annotation. |
| [`Code/CCI_Niches/`](Code/CCI_Niches) | Cell-cell interactions per spatial niche (LIANA, CrossTalkeR). |
| [`Code/ukbb_processing/`](Code/ukbb_processing) | UK Biobank plasma proteomics analyses. |
| [`Code/cell_boundaries/`](Code/cell_boundaries), [`Code/plot_egfr_f2_on_spatial/`](Code/plot_egfr_f2_on_spatial) | Two first-submission scripts that are still needed: `1_preprocess_xenium.py` and `1_preproc_cell_boundaries.py` write the cell polygons of the 2025 three-region Xenium object that Fig. 6g (right) reads. |

The first-submission NicheCompass and PEC trajectory code (`Code/Niche_analysis/`, `Code/Trajectory/`) was replaced by
`Code/revision_2026/analyses/01_xenium_preprocessing_niches/` and
`Code/revision_2026/analyses/09_snrna_pec_immune_subclusters/`.
The first-submission Xenium neighbourhood and cell-cell interaction scripts (`Code/nh_analysis/`, `Code/spatial_cci/`,
`Code/utils/` and the plotting scripts of `Code/cell_boundaries/` and `Code/plot_egfr_f2_on_spatial/`) were replaced by
`Code/revision_2026/analyses/03_xenium_factor1_neighbourhoods_cci/`. The Milo scripts (`Code/milo/`) were removed
because no figure of the revised manuscript uses them.

## License

This project is distributed under the GNU General Public License v3.0. See [`LICENSE`](LICENSE) for full details.
