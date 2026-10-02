# Public-release checklist

Before the GitHub/Zenodo release:

- [x] Full R rerun completed.
- [x] Automated numerical verification: 48/48 checks passed.
- [x] Frozen v8.4 reference-result tables included.
- [x] One-off reconciliation/debug scripts excluded from the public bundle.
- [x] Local Windows paths removed from release analysis scripts.
- [ ] Recover and add the original raw-data preparation script if available; otherwise retain the explicit stage-1 provenance limitation in README.
- [x] MIT software license added (`LICENSE`).
- [ ] Create GitHub repository and add its URL to README/CITATION metadata if desired.
- [ ] Create a versioned Zenodo release and DOI.
- [ ] Add the final GitHub/Zenodo links to the manuscript Data/Code Availability statement.
