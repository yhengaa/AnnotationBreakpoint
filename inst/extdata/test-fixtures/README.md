# Synthetic package test fixtures

These three tiny, deterministic Layer-A fixtures exercise the exact reference
implementation without external data or network access:

- `fragile`: one reassignment overturns the registered claim.
- `stable`: every feasible assignment preserves the registered claim.
- `coordinated`: the minimum overturn requires two coordinated moves.

The `expected_results.tsv` tables were retained unchanged from the accepted
historical test corpus. They are test expectations, not public benchmark
results or software output examples.
