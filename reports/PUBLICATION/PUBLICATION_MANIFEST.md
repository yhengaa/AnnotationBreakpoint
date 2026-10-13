# Publication snapshot

This repository is a compact distribution of AnnotationBreakpoint 0.2.0 source, examples, and selected numerical evidence. The downloadable public source package has SHA-256 `17e57844e6c8abad4147032bf87e1d3d0d0e03776f5226e5992be56eba0d04d3` and is available from the [v0.2.0 release](https://github.com/yhengaa/AnnotationBreakpoint/releases/tag/v0.2.0).

## Source provenance

The repository source was extracted from an earlier checked package (SHA-256 `733e005289dc9930b881f07fb1ebbb1c513424f746b0feda09024a5f7e33ce07`), whose source revision is `9f97b51c2b5c7bc8c15997022b513f7cb8c790e5`. That earlier artifact is provenance evidence rather than the public download. `package_source_parity.json` is its unchanged package-build receipt; it records 44 byte-identical package files and zero executable-source mismatches for the initial publication snapshot. Later public documentation revisions are recorded in ordinary Git history and do not change the executable source or numerical results.

Publication-specific differences from that checked package are limited to:

1. public-facing README, tutorial, license, notice, and evidence-index documentation;
2. omission of two private mandate-text files that are not executable package source;
3. publication `.gitignore` and `.Rbuildignore` files; and
4. selected compact evidence and publication metadata.

The package version is 0.2.0. The version-neutral `inst/NOTICE` retains the established dependency versions and attribution. `LICENSE.md` presents the full MIT terms using the established year and copyright holder; the R package `LICENSE` metadata remains unchanged.

## Evidence and exclusions

The selected evidence preserves the fixed 13-case comparator result, the FAM118A and PDXK upper-bound examples, the negative Hao5P demonstration, and the inconclusive scaling and public-utility conclusions. Numerical result tables and historical scientific result notes are unchanged.

Excluded material includes raw counts, caches, large result objects, job logs, failed partial saves, private coordination and runtime records, agent histories, full local Git history, and large archival bundles. These exclusions keep the public repository compact; they do not change the included scientific results.

`SHA256SUMS` at repository root hashes every included file except itself.
