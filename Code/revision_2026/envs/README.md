# Conda environments

`<name>.yml` = `conda env export --no-builds` of the environment the scripts were run with (12 Nov 2026, Aachen HPC,
CPU only). `<name>_R_packages.csv` = every R package installed in the R environments (incl. GitHub installs such as
scMEGA and ArchR). Create one with `conda env create -f envs/<name>.yml`; the SLURM wrappers activate them by name.

| environment | used for |
|---|---|
| agent_env | general Python analyses and figures (ATAC views 34/45-48, snRNA overview, biopsy-age QC, factor-1 list) |
| ViaEnv | PHATE / VIA / diffusion-map trajectories (ATAC joint pseudotime, PEC trajectory, gene trends) |
| nichecompass | NicheCompass models and downstream, Xenium factor scoring, NELL1, filtered DE tables |
| spatial_env / spatialdata | Xenium geometry, Voronoi neighbourhoods, cell-boundary maps, Xenium figures |
| cci_env | local cell-cell interaction scoring (CellPhoneDB v5 database) |
| kpmp_repurposing | scANVI label transfer (scvi-tools), PEC mapping onto Xenium |
| omicverse | joint PCA / Harmony / KNN label transfer when merging the five Xenium regions |
| pathfire | level2 subcluster sweep (snRNA) |
| decoupler | Venn of niche-20 genes vs factor 1 |
| de_env | MAST, fgsea, clusterProfiler (R) |
| covariate | limma-voom (R) |
| scmega2 | ArchR, Signac, Seurat, scMEGA eGRN (R) |
