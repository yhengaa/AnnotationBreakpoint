# Installed examples

`exact_toy/` is the frozen three-feature Layer-A example. `edger_witness/`
is the frozen 101-feature edgeR witness example and includes its authoritative
analysis specification. Each directory contains the six ordinary input tables;
neither contains a test-answer table. `provenance.tsv` binds every installed
file to its source path and SHA-256 digest.

The practical software default is `proxy_prefix`, labeled
**observed-centered proxy-ranked cumulative prefixes**. It fully refits the
available exact registered prefixes `1, 2, 4, 8, 16, 32, 64, 128`, including
those after an earlier crossing, and returns the earliest prefix whose witness
passes independent replay. An inadmissible prefix is rejected without
skip/fill, `outgoing_only` is rejected, and a no-hit remains unresolved. These tiny installed examples teach
the interface; the accepted FAM118A and PDXK saved examples remain in the
source repository and are not bundled here as large public inputs.
