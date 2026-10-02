# Pathway-score mediation reproducibility code

R code and frozen reference results for:

> **When Does a Pathway Score Preserve a Mediated Effect? Aggregation Geometry, Minimax Adequacy, and Gap Inference for High-Dimensional Mediation**

Author: **Subir Hait**

## Scope of this repository

This repository reproduces the manuscript analyses **from the processed TCGA LUAD stage-1 object onward**. It contains the simulation code, three-gene application, genome-wide outcome-free audit, independent cross-screened audit, figure generation, and numerical verification against the frozen v8.4 manuscript results.

The required processed input is:

```text
tcga_lung/luad_stage1.rds
```

with list elements `Mval`, `ph`, `ann`, and `alpha`.

### Important provenance limitation

The exact historical raw-data preparation script that created `luad_stage1.rds` is not included in the recoverable analysis archive. This repository therefore **does not claim raw-data-to-results reproducibility**. It provides reproducibility from the frozen processed stage-1 object onward. The original preprocessing script should be added if it is recovered before public release.

The processed stage-1 object itself is **not tracked in Git**.

## Repository contents

```text
00_install_packages.R
01_validate_stage1_input.R
02_simulation_and_application.R
03_gap_inference_simulation.R
04_make_figures.R
05_genome_wide_audit.R
06_cross_screened_audit.R
07_verify_against_frozen_v8_4.R
99_session_info.R
run_all_after_stage1.R
run_full_verification.R
reference_results_v8_4/
tcga_lung/README.md
VERIFICATION.md
CITATION.cff
.gitignore
SHA256SUMS.txt
```

## Main analysis order

1. `00_install_packages.R` — install/check required packages.
2. `01_validate_stage1_input.R` — validate the stage-1 RDS object.
3. `02_simulation_and_application.R` — semi-synthetic study, estimator study, cautionary survival analysis, and three-gene application.
4. `03_gap_inference_simulation.R` — paired influence-function gap-inference simulation.
5. `05_genome_wide_audit.R` — genome-wide outcome-free audit.
6. `06_cross_screened_audit.R` — independent screen/estimate audit and pooled summaries.
7. `04_make_figures.R` — regenerate manuscript figures.
8. `99_session_info.R` — save the R environment.

Once `tcga_lung/luad_stage1.rds` is in place, the main workflow can be run with:

```r
source("run_all_after_stage1.R")
```

## Full verification run

For the frozen manuscript configuration:

```r
Sys.setenv(
  PATHWAY_SCORE_ROOT = getwd(),
  AUDIT_SPLITS = 25,
  AUDIT_INNER_SPLITS = 25,
  GENOME_AUDIT_GENE_BOOT = 5000,
  CROSS_AUDIT_GENE_BOOT = 2000,
  AUDIT_TOPK = "25,50,100,250",
  SIM_REPS = 800
)

source("run_full_verification.R")
```

The verifier compares manuscript-facing outputs against `reference_results_v8_4/` and checks representative RDS values and regenerated figures.

The final clean rerun on **2026-09-19** used R 4.6.0 on Windows 11 and returned:

```text
Checks: 48
PASS: 48
FAIL: 0
All verification checks passed.
```

See `VERIFICATION.md` for the tolerance note and verified configuration.

## Reproducibility settings

- Gap-inference simulation seed: `16092026`
- Gap-inference replicates per cell: `800`
- Genome-wide audit splits: `25`
- Inner audit splits: `25`
- Genome-wide descriptive gene bootstrap: `5000`
- Cross-screen gene bootstrap: `2000`
- Cross-screen top-K values: `25, 50, 100, 250`

## Verification tolerance

Deterministic quantities are checked under strict numerical tolerances. Only the two percentile-bootstrap endpoints `gene_boot_median_lo` and `gene_boot_median_hi` in `genome_wide_audit_benchmark.csv` use an absolute tolerance of `0.002`, reflecting small Monte Carlo/build-dependent bootstrap-quantile variation. Point estimates and all other outputs remain under the strict tolerance.

## Data source

The empirical application uses TCGA LUAD data obtained through UCSC Xena. Do not commit raw or processed patient-level data to this repository unless redistribution is explicitly permitted.

## Software environment

The successful verification run used R 4.6.0 on Windows 11. Key packages included `data.table`, `survival`, `pseudo`, `glmnet`, `dplyr`, `UCSCXenaTools`, and `limma`.

## License

The software and R code in this repository are released under the **MIT License**; see `LICENSE`. TCGA/UCSC Xena data and any local processed data objects are **not** relicensed by this repository and remain subject to their original terms of use. The manuscript text is not covered by the software MIT license.
