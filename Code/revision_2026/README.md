# Code for the revised manuscript: analyses and figure panels (Leon Frommer)

Saritas, Boys, Frommer, Kuppe *et al.* — *Cellular and molecular determinants of human kidney disease: single-nuclei, ATAC-seq
and Xenium based high-resolution spatial mapping of human kidney diseases* (revision 2026).

This folder holds the code behind the figure panels produced by Leon Frommer (Hayat lab): the Xenium cell typing,
factor 1 on Xenium, spatial neighbourhoods and cell-cell interactions, NicheCompass niches, the glomerular niche 20 /
crescent analysis, the snRNA-seq PEC and immune subclusters, the PT→aPT snATAC/scMEGA gene-regulatory network, the
snRNA-seq overview, and the Xenium pseudobulk tables handed to the MCFA authors. [`FIGURE_PANELS.md`](FIGURE_PANELS.md) lists every panel of the paper, who made it, and the script
that draws it. Panels made by other authors (MCFA, atlas construction, UK Biobank, wet lab, a colleague's Xenium
notebooks) are listed there too, with a pointer to where their code lives.

## Layout

| folder | contents |
|---|---|
| `analyses/` | the analyses, in pipeline order (`01_` … `11_`). Each folder has a README: purpose, run order, conda environment, inputs, outputs, panels. Scripts keep the names they were run under (several import each other). |
| `figures/<figure>/` | one entry point per published panel, named after it (`Fig3/Fig3c_…`, `ED4/ED4d-e_…`, `S19/S19c-d_…`). Each reads saved analysis outputs and writes the panel to `$KB_FIG_OUT/<figure>/`; its header names the inputs, the plotting code, edits made afterwards in Illustrator, and how it was verified. |
| `config/` | `env.sh` (source this first), `kbpaths.py` / `kbpaths.R` (every data location in one place), `cores.py` (Xenium core ids). |
| `envs/` | conda environment exports and R package lists. |
| `FIGURE_PANELS.md` / `.csv` | panel → maker → entry script. |
| `MANIFEST.csv` | every file in `analyses/` with the workspace file it was copied from and checksums. |

## Running the code

```bash
source config/env.sh          # sets KB_REPO, PYTHONPATH and the data locations below
export KB_DATA=/path/to/data  # override any location before or after sourcing
python figures/Fig3/Fig3c_*.py   # writes figure_output/Fig3/...
```

All paths go through these variables (defaults: the Aachen HPC where the analyses ran):

| variable | meaning |
|---|---|
| `KB_ANALYSIS` | analysis workspaces; every analysis reads and writes under `$KB_ANALYSIS/<workspace>` (named in its README) |
| `KB_PROJECT` | project folder with the snRNA-seq objects and earlier NicheCompass runs |
| `KB_DATA`, `KB_RESULTS` | shared processed data and results (snRNA/snATAC objects, Xenium/NicheCompass objects) |
| `KB_XENIUM_RAW` | Xenium Onboard Analysis output bundles |
| `KB_FIG_OUT` | where figure entry points write (default `figure_output/`) |
| `KB_CORE_KEY`, `KB_EXCLUDED_CORES` | optional private files (see below) |

Heavy steps were run as SLURM jobs; the `.sh` wrappers next to the scripts carry the resources and the conda
environment (`#SBATCH` log paths are relative: create `logs/` in the directory you submit from). Everything ran on
CPUs (no GPU). Human GRCh38, HGNC gene symbols, random seed 42 unless a script says otherwise.

## Data

The code needs the processed objects of the study (snRNA-seq, snATAC-seq, Xenium) and clinical tables; see the Data
availability statement of the paper. The UK Biobank participant data used for Fig. 4e's protein table are available
only through UK Biobank and are not part of this repository.

## Patient identifiers

Internal Xenium core labels contain pathology numbers and are never written into this repository. Scripts refer to a
core by its public id (`xenium_core_id`, XEN###, as in the de-identified Xenium data); `config/cores.py` translates it
to the internal label only when a private key is supplied (`KB_CORE_KEY`). Cores excluded by quality control are read
from a private list (`KB_EXCLUDED_CORES`); none of them is part of the public Xenium data. Donors are referred to by
K-numbers.

## Reproducing the panels: what to expect

- Every entry script was re-run on the saved analysis outputs (12 Nov 2026). Of the 68 panel letters made with this code,
  46 came out pixel- or byte-identical to the files used in the published figures, 7 are identical except for the title,
  and the 4 panels of Supplementary Fig. 10 are identical when drawn from the earlier object they were made from (below);
  the `re-run check` column of `FIGURE_PANELS.md` gives the result and the reason for every other panel.
- Map panels (Fig. 3f, 6b, 6g; Extended Data Fig. 4d–e; Supplementary Fig. 13) now print the public core id in their
  title; everything below the title is identical.
- Some originals used unseeded randomness that cannot be replayed exactly: label placement (ggrepel, Supplementary
  Fig. 19b), jitter (Extended Data Fig. 2e–g), bootstrap confidence intervals (Supplementary Fig. 23h–i), and the VIA
  pseudotime of the PEC trajectory, whose fitted object was not saved (Fig. 6i and the gene trends of 6j are produced
  only by re-running the analysis step; their entry scripts say how).
- Fonts: several plotting scripts register Liberation Sans from `/usr/share/fonts/liberation-sans/` and otherwise fall
  back to DejaVu Sans, which shifts text and therefore the exact layout (seen for Fig. 6f and 6k on machines without the
  font).
- Supplementary Fig. 10 was drawn from the cell-type object before 182 aPT cells were re-assigned to PEC; its entry
  script reproduces it with `--xenium <that earlier object>` and draws the final object by default.

## How this folder was built

The scripts were copied from the analysis workspaces on 12 Nov 2026 (`analyses/11_xenium_nuclear_pseudobulks` on 19 Nov 2026). Hard-coded paths were replaced by the variables
above by an automatic rewrite whose equivalence was checked file by file (Python: identical syntax trees once the
variables take their default values; R and shell: exact reversal of every substitution); deliberate edits (removal of
identifiers, locating shipped instead of workspace copies) are listed per file in `MANIFEST.csv` (`patch_rules`). Every
figure entry point was re-run against the saved analysis outputs and compared with the file used in the published
figure; the result is in its header.
