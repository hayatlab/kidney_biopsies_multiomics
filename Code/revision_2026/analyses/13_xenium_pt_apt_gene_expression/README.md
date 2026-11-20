# 13 — Expression of single genes in PT and aPT cells per Xenium core

## Purpose
Compares the expression of a gene between PT and aPT cells of the Xenium data with the core as the unit (Fig. 3e:
RNF152 and NFKBIZ). Analysis by Sidrah Maryam (notebook `Fig3e.ipynb`); `gene_expression_per_core.py` holds the
functions of that notebook, unchanged except that the name of the core column is an argument.

Method: per cell log1p(count / library size × 10,000), library size = sum of the raw counts over the 468 panel genes;
mean per core and cell type (`level2_refined`, PT and aPT); one point per core and cell type, coloured by condition;
diamonds = mean of the core values per condition; bracket = two-sided Welch t-test PT vs aPT on the core values
(Benjamini-Hochberg over the tested pairs: one pair, P unchanged). Human, HGNC symbols.

## Analysis workspace
None: the figure script reads the Xenium object and writes only the panel.

## Run (source `config/env.sh` first)
| script | conda env | output |
|---|---|---|
| `python figures/Fig3/Fig3e_rnf152_nfkbiz_pt_vs_apt.py` (imports `gene_expression_per_core.py`; login node, about 2 min) | nichecompass | `$KB_FIG_OUT/Fig3/Fig3e_rnf152_nfkbiz_pt_vs_apt.pdf` |

## Input
| input | producer |
|---|---|
| `$KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad` (258,375 cells × 468 genes, `layers["counts"]`, `level2_refined`, `condition`) | `analyses/01_xenium_preprocessing_niches` |

## Result (20 Nov 2026, re-run of the figure script)
119 cores; 117 have PT cells and 119 have aPT cells (56,181 PT and 14,890 aPT cells).

| gene | mean of core values, PT | aPT | difference aPT − PT (95% CI) | Welch t-test |
|---|---|---|---|---|
| RNF152 | 3.97 | 1.48 | −2.49 (−2.69 to −2.28) | P = 3.3e-63 |
| NFKBIZ | 0.98 | 2.06 | +1.08 (+0.88 to +1.27) | P = 6.6e-23 |

Paired over the 117 cores with both cell types (not shown in the panel): RNF152 −2.47 (−2.60 to −2.34), NFKBIZ +1.06
(+0.91 to +1.21). The 24 condition means measured in the printed panel agree with the re-run within 0.004.

## Caveats (read before reusing)
- The test treats every core × cell type value as an independent observation. PT and aPT values come from the same
  cores, and several cores can come from one donor; the paired differences above respect the first of these, neither
  respects the second.
- **The gene list of the notebook as received (PROM1, NFKBIZ, RNF212B, GRAMD1B, SLC13A3, PDZD2) does not contain
  RNF152**; the published panel was drawn with the same functions for RNF152 and NFKBIZ.
- The notebook computes a one-way ANOVA P for the plot title, but its table holds empty core × condition combinations,
  so that P is not a number and no ANOVA P appears in the plot; nothing in the panel depends on it.
- Library size is the sum over the panel genes, so "expression" is relative to the 468-gene panel.
- The legend of the published panel says "CTRL_Xenium"; the object's label is `CTRL` (33 control cores).
