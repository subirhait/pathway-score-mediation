# TCGA LUAD processed input and generated outputs

Place the required processed stage-1 object here:

```text
luad_stage1.rds
```

The analysis expects the RDS object to contain `Mval`, `ph`, `ann`, and `alpha`.

This directory is also used by the scripts for generated genome-wide and cross-screen output files.

Patient-level/raw/processed TCGA data are intentionally not included in the public Git bundle. The exact historical script that originally created `luad_stage1.rds` is not present in the recoverable archive, so this repository begins at the processed stage-1 object.

## Licensing note

The MIT License in the repository applies to the software/code, not to TCGA/UCSC Xena data or any local processed data objects.
