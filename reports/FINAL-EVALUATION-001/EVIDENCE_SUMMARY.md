# Accepted evaluation evidence

This publication copy summarizes the accepted frozen evaluation without raw
counts, caches, job logs, or local runtime records. The numerical tables in
this directory are unchanged copies of the accepted project tables.

The scaling evaluation completed 984/984 registered tasks with all saved
verifications passing and no scientific failures. Greedy at 128 fits matched
the exact cost in all 20 repeated K=6/K=12 observations. The registered K=192
primary decision remained **inconclusive** because only 10 exact-attainable
parents were available, below the required denominator of 12. Within those ten
parents, all three methods found 2/10 replayed witnesses at 128 fits. At the
separate 1,024-fit coverage tier, greedy found 10/10, random 7/10, and pair
5/10. This is bounded evidence for the executed implementations and allowances,
not a general performance guarantee.

The public COMBAT evaluation completed 96/96 registered searches. Only three
claims originally met the fixed effect threshold. Four replayed witness
observations all concerned the near-threshold IL18R1 claim; the one-cell greedy
witness changed beta from 0.500445209221746 to 0.496022306993079. There were
zero strict added-information claims and one explanation/compression example.
With only three eligible claims against a registered minimum denominator of
ten, the public-utility conclusion was **inconclusive**. The observed evidence
did not support added public utility, and it is not a robustness certificate.

The later decisive comparator ended the methods claim. Its separate
[result note](../DECISIVE-COMPARATOR-001/RESULT_NOTE.md) and
[39-row table](../DECISIVE-COMPARATOR-001/results_13x3.tsv) are included in
this repository. No new fit or analysis was performed to create this
publication summary.

Included tables:

- `headline_summary.tsv`: accepted headline outcomes and interpretations.
- `public_eligible_claims.tsv`: the three originally eligible public claims.
- `public_witnesses.tsv`: all four replayed public witness observations.
- `scaling_k192_coverage.tsv`: primary and additional K=192 coverage tiers.
- `scaling_depth_denominators.tsv`: explicit depth denominators.
