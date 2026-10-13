# Reconstructing the saved COMBAT examples

AnnotationBreakpoint 0.3.0 includes a standalone builder for the public COMBAT
CITE-seq data on Zenodo. It reconstructs the exact 28,121-cell evaluation panel
used by the saved FAM118A and PDXK examples. The builder produces input files and
checks identities; it does not fit a model or run a search.

## Requirements and source

Use the package's registered R 4.5.x/Bioconductor 3.21 environment. The builder
also needs the optional Bioconductor package `rhdf5`:

```sh
Rscript -e 'BiocManager::install("rhdf5")'
```

Download `COMBAT-CITESeq-DATA.h5ad` from [Zenodo record
6120249](https://doi.org/10.5281/zenodo.6120249). The file is about 6.4 GB. Its
registered source identity is:

- bytes: `6409089483`
- MD5: `87c6b1a733ea1adc37c808d9a357de74`
- license: CC BY 4.0

```sh
curl --fail --location \
  --output COMBAT-CITESeq-DATA.h5ad \
  https://zenodo.org/api/records/6120249/files/COMBAT-CITESeq-DATA.h5ad/content
```

The builder checks the byte count and MD5 by default.

## Build the input directories

```sh
AB_COMBAT_BUILDER="$(Rscript --vanilla -e 'cat(system.file("scripts", "build_combat_inputs.R", package="AnnotationBreakpoint"))')"

Rscript --vanilla "$AB_COMBAT_BUILDER" \
  --h5ad /path/to/COMBAT-CITESeq-DATA.h5ad \
  --output-dir combat-inputs
```

The output contains separate `inputs/FAM118A` and `inputs/PDXK` directories,
matching analysis specifications under `analysis_specs`, semantic hashes, and a
file manifest. Each input directory has six files: sparse `counts.rds` plus
`cells.tsv`, `specimens.tsv`, `allowed_labels.tsv`, `claim.tsv`, and
`constraints.tsv`. Named donor and sample limits are JSON objects inside the
one-row constraints table.

`--skip-source-md5` is available only for a local file whose checksum has
already been independently verified. The builder still checks the byte count,
raw cell/feature mappings, ADT ambiguity rule, matrix structure, and all saved
semantic identities, and records that the source checksum step was skipped.

## Run an example

```sh
AB_CLI="$(Rscript --vanilla -e 'cat(system.file("exec", "annotation-breakpoint", package="AnnotationBreakpoint"))')"

Rscript --vanilla "$AB_CLI" \
  --input-dir combat-inputs/inputs/FAM118A \
  --analysis-spec combat-inputs/analysis_specs/FAM118A.yaml \
  --output-dir combat-results/FAM118A
```

The accepted saved results under the registered 445-change global limit and
named donor caps are replay-verified upper bounds:

| Gene | Original beta | Replayed beta | Changed labels |
| --- | ---: | ---: | ---: |
| FAM118A | 0.8106752562 | 0.4731965904 | 64 |
| PDXK | 0.7274641768 | 0.3373883736 | 32 |

The contrast is COVID_CRIT versus COVID_SEV in COMBAT CD4.TCM cells across 16
donors, with CD4.NAIVE as the permitted alternative and threshold 0.5. These
costs are verified upper bounds, not certified minima.

## Frozen reconstruction details

The package ships the accepted cell IDs and their raw H5AD row indices, exact
specimen metadata, ADT rule, analysis specification, named limits, and expected
semantic hashes. Raw RNA counts are rebuilt from `/layers/raw`; feature IDs and
types come from `/var`; CD45RA/CD45RO support is rebuilt from `/X`. The original
labels and the ambiguity flag remain separate fields in `cells.tsv`.

The frozen selection tables are derived from the COvid-19 Multi-omics Blood
ATlas (COMBAT) Consortium dataset *A blood atlas of COVID-19 defines hallmarks
of disease severity and specificity: Associated data* (DOI
10.5281/zenodo.6120249), distributed under CC BY 4.0. AnnotationBreakpoint
package code is distributed under the MIT License.
