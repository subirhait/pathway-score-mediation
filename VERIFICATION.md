# Verification record

## Final clean verification

Date: **2026-09-19**

Environment:

- R 4.6.0 (2026-04-24 ucrt)
- Windows 11 x64
- Genome-wide audit splits: 25
- Inner splits: 25
- Genome-wide descriptive gene bootstrap: 5,000
- Cross-screen gene bootstrap: 2,000
- Gap-inference simulation: 800 replicates per cell

Final automated result:

```text
Checks: 48
PASS: 48
FAIL: 0
All verification checks passed.
```

## Monte Carlo tolerance note

The verifier uses strict tolerances for deterministic quantities. The only exception is an absolute tolerance of **0.002** for the two stochastic percentile-bootstrap endpoints `gene_boot_median_lo` and `gene_boot_median_hi` in `genome_wide_audit_benchmark.csv`.

This exception was introduced after a complete clean rerun reproduced all point estimates and deterministic outputs, while the two bootstrap CI endpoints differed by less than 0.001 because of small Monte Carlo/build-dependent variation.

No manuscript point estimate was changed to obtain the 48/48 verification result.
