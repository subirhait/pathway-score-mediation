# v1.0.0 - Frozen reproducibility release

This is the frozen reproducibility release accompanying the manuscript:

**When Does a Pathway Score Preserve a Mediated Effect? Aggregation Geometry, Minimax Adequacy, and Gap Inference for High-Dimensional Mediation**

## Included

- R analysis scripts from the processed TCGA LUAD stage-1 object onward
- paired gap-inference simulation
- genome-wide outcome-free audit
- independent screen/estimate audit
- manuscript figure generation
- 16 frozen v8.4 reference-result tables
- numerical verification workflow

## Verification

Final clean verification on 2026-09-19 using R 4.6.0 on Windows 11:

```text
Checks: 48
PASS: 48
FAIL: 0
All verification checks passed.
```

The verifier retains strict tolerances for deterministic outputs. Only the two stochastic percentile-bootstrap endpoints in `genome_wide_audit_benchmark.csv` use an absolute tolerance of 0.002, as documented in `VERIFICATION.md`.

## Data

TCGA/UCSC Xena data are not redistributed. The repository begins from the local processed file `tcga_lung/luad_stage1.rds`. The historical raw-data preparation script that created that object is not present in the recoverable archive; this limitation is stated explicitly in the README.

## License

Software and R code: MIT License. Data remain subject to their original terms of use.
