# FAM118A walkthrough

## The question

The saved COMBAT example asks whether a positive FAM118A expression effect for
COVID_CRIT versus COVID_SEV in CD4.TCM cells remains above the registered 0.5
log2-effect threshold when selected cells may move between CD4.TCM and
CD4.NAIVE according to their registered allowed-label sets.

Under the observed labels, the fitted effect was `0.8106752561679119`. A
permitted 64-cell prefix produced `0.4731965904097605`, and an independent
replay confirmed the crossing. The cost of 64 is a verified upper bound: this
run established that 64 changes were sufficient, without certifying that no
smaller feasible assignment exists.

![FAM118A log2 effect decreased from 0.811 to 0.473 after 64 permitted label changes, crossing the registered threshold of 0.5.](assets/fam118a-sensitivity.png)

## Saved inputs

The registered analysis contains 28,121 cells from 16 donors, 20,615 RNA
features, and 8,898 cells with two permitted labels. It retains the original
integer counts, specimen records, allowed-label table, CD4.TCM/CD4.NAIVE
alternatives, global limit of 445 changes, and named donor and sample limits.
<details>
<summary>Semantic identities checked during reconstruction</summary>

| Component | SHA-256 |
| --- | --- |
| Counts | `3bde21f1e55309be2e8a262909dd20e0f39ef57e0609e85ec873e99c82893a0e` |
| Cells | `05ec140f67b8725d6cb322917f85bcea05317f1e20f156db7e3fe8846f0988fd` |
| Specimens | `ae1184622d10591c7a80be088788ff5a6da6cce06a0315813d345e668f924aab` |
| Allowed labels | `56606603ce0ebbe1604d09d4b78237d5f94c048204ea207650259f76dd4f1560` |
| Constraints | `80d387afb28734e7bcee816ad5276a62bc0e6a0b5aaccd6e332cd61e06762485` |

</details>

The source is the COMBAT Consortium CITE-seq dataset on
[Zenodo](https://doi.org/10.5281/zenodo.6120249), distributed under CC BY 4.0.

## Try the small bundled example

After installing AnnotationBreakpoint 0.3.0, run the bundled toy example:

```sh
AB_CLI="$(Rscript --vanilla -e 'cat(system.file("exec", "annotation-breakpoint", package="AnnotationBreakpoint"))')"
AB_EXAMPLE="$(Rscript --vanilla -e 'cat(system.file("extdata", "examples", "exact_toy", package="AnnotationBreakpoint"))')"

Rscript --vanilla "$AB_CLI" \
  --input-dir "$AB_EXAMPLE" \
  --layer-a \
  --output-dir ab-output/quickstart
```

This portable toy example demonstrates the input, search, replay, and reporting
flow with the `--layer-a` flag. It is separate from the edgeR public case above.

## Reconstruct and run the public example

The reconstruction requires the optional Bioconductor package `rhdf5`:

```sh
Rscript -e 'BiocManager::install("rhdf5")'
```

Download the registered public source file, which is about 6.4 GB:

```sh
curl --fail --location \
  --output COMBAT-CITESeq-DATA.h5ad \
  https://zenodo.org/api/records/6120249/files/COMBAT-CITESeq-DATA.h5ad/content
```

Build the exact registered input directories with the installed reconstruction
script:

```sh
AB_COMBAT_BUILDER="$(Rscript --vanilla -e 'cat(system.file("scripts", "build_combat_inputs.R", package="AnnotationBreakpoint"))')"

Rscript --vanilla "$AB_COMBAT_BUILDER" \
  --h5ad ./COMBAT-CITESeq-DATA.h5ad \
  --output-dir combat-inputs
```

Then run the default fixed-prefix search:

```sh
Rscript --vanilla "$AB_CLI" \
  --input-dir combat-inputs/inputs/FAM118A \
  --analysis-spec combat-inputs/analysis_specs/FAM118A.yaml \
  --output-dir combat-results/FAM118A
```

These commands are the shipped reconstruction and CLI workflow; they were not
rerun while preparing this walkthrough. The recorded 64-change result comes
from the accepted saved evidence in
[`results_13x3.tsv`](../reports/DECISIVE-COMPARATOR-001/results_13x3.tsv) and
[`method_b_prefixes.tsv`](../reports/DECISIVE-COMPARATOR-001/method_b_prefixes.tsv).

For full input details, see the
[COMBAT reconstruction guide](../inst/doc/COMBAT_REPRODUCTION.md). Cite this
software release with [doi:10.5281/zenodo.23353292](https://doi.org/10.5281/zenodo.23353292).
