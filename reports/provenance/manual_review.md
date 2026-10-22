# Manual Provenance Review

Date: 2026-06-01

This note interprets the regenerated mechanical audit in
`reports/provenance/audit.md` and
`reports/provenance/similarity_inventory.csv`. It is not a legal conclusion,
legal advice, or a provenance certificate.

## Count Reconciliation

The audit was regenerated with:

```bash
python3 scripts/generate_similarity_audit.py
```

Current mechanical counts:

| Risk | Count | Manual disposition |
|---|---:|---|
| High | 11 | All rows classified below |
| Medium | 58 | Keep in review queue; no exact code/comment/string shingles unless noted by the high table |
| Low-review | 19 | Track as weak mapped signals, especially on algorithmically important paths |
| Low | 10350 | Mechanical inventory only |
| Total pairs scanned | 10438 | Matches `audit.md` |

These counts match `reports/provenance/audit.md`.

## High-Risk Row Disposition

The `high` label means the mechanical tool found exact retained code-token
shingles, exact comment shingles, or shared strings. It does not by itself
prove source copying. Every high row remains listed here so release review does
not depend on memory.

| Rust file | ROOT file | Signal | Disposition | Rationale |
|---|---|---|---|---|
| `src/bin/ref_compare_runner.rs` | `test/testMinuit2.cxx` | code=22 | keep-validation | Reference comparison runner intentionally mirrors workload/test semantics; no comment or string shingles; not shipped library implementation logic. |
| `src/covariance_squeeze.rs` | `test/testMinimizer.cxx` | code=12 | keep-with-note | Overlap is with ROOT tests and generic covariance fixture logic; no comment or string shingles. Keep in provenance inventory. |
| `src/covariance_squeeze.rs` | `test/testADMinim.cxx` | code=12 | keep-with-note | Same covariance test-fixture overlap; not a mapped implementation source row; no copied comments or strings. |
| `src/covariance_squeeze.rs` | `test/testMinuit2.cxx` | code=6 | keep-with-note | Small exact token shingles against tests only; no copied comments or strings. |
| `src/bin/ref_compare_runner.rs` | `test/testUserFunc.cxx` | code=4 | keep-validation | Runner/test vocabulary overlap; no comment or string shingles. |
| `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad4F.h` | code=4 | keep-validation | Tutorial quadratic fixture overlap in validation code; no comment or string shingles. |
| `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad8F.h` | code=4 | keep-validation | Tutorial quadratic fixture overlap in validation code; no comment or string shingles. |
| `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad12F.h` | code=4 | keep-validation | Tutorial quadratic fixture overlap in validation code; no comment or string shingles. |
| `src/user_parameters.rs` | `src/MnUserTransformation.cxx` | comments=3 | document | Mapped high row is comment-shingle based. Current public wording uses ROOT as a numerical/API reference and avoids source-port language; keep source comments neutral. |
| `src/gradient/analytical.rs` | `test/testMinuit2.cxx` | code=3 | keep-validation | Shared analytical-gradient test vocabulary against a ROOT test file; no comment or string shingles. |
| `src/parameter.rs` | `test/testMinimizer.cxx` | code=3 | keep-api-vocabulary | Minuit parameter names are public API vocabulary; row is against tests, with no comment or string shingles. |

No high row currently points to the four implementation paths that previously
needed focused rewrite attention: `src/linesearch.rs`, `src/simplex/builder.rs`,
`src/transform/sin.rs`, and `src/migrad/builder.rs`.

## Rewritten Path Check

The four named paths were checked against the regenerated CSV. The important
post-burndown signal is that the relevant mapped/source rows have
`code_shingles = 0`, `comment_shingles = 0`, and `shared_string_literals = 0`.
Residual signal is generic structural similarity from the shared algorithm and
Minuit API vocabulary, so the disposition is `document` rather than `rewrite`.

| Rust path | ROOT row checked | Risk | Exact code/comment/string shingles | Disposition |
|---|---|---|---|---|
| `src/linesearch.rs` | `src/MnLineSearch.cxx` | medium | 0 / 0 / 0 | document residual generic line-search structure |
| `src/linesearch.rs` | `inc/Minuit2/MnLineSearch.h` | medium | 0 / 0 / 0 | document API mapping |
| `src/simplex/builder.rs` | `inc/Minuit2/SimplexBuilder.h` | medium | 0 / 0 / 0 | document API mapping |
| `src/simplex/builder.rs` | `src/SimplexBuilder.cxx` | low | 0 / 0 / 0 | implementation-source row no longer has exact retained shingles |
| `src/transform/sin.rs` | `inc/Minuit2/SinParameterTransformation.h` | medium | 0 / 0 / 0 | document mathematical bounded-transform formula |
| `src/transform/sin.rs` | `src/SinParameterTransformation.cxx` | low | 0 / 0 / 0 | implementation-source row no longer has exact retained shingles |
| `src/migrad/builder.rs` | `src/BFGSErrorUpdator.cxx` | medium | 0 / 0 / 0 | document residual DFP/BFGS algorithm structure |
| `src/migrad/builder.rs` | `inc/Minuit2/VariableMetricBuilder.h` | medium | 0 / 0 / 0 | document API mapping |
| `src/migrad/builder.rs` | `src/VariableMetricBuilder.cxx` | low | 0 / 0 / 0 | implementation-source row has generic structure only |

## Medium-Risk Policy

The 58 medium rows remain part of the release-review inventory. Most are mapped
ROOT/Rust API correspondences or shared mathematical algorithms. They should be
handled by source-neutral comments, tests, and reference-oracle wording, not by
claims of legal clearance. If a future review finds copied comments,
branch-for-branch structure, unusual copied local names, or non-mathematical
constants, create a targeted rewrite bead for that file.

## Claim Boundary

Public docs should continue to say:

- ROOT Minuit2 is used as a numerical/API reference and parity baseline.
- `minuit2-rs` is an independent Rust implementation of Minuit-style algorithms.
- The audit is provenance evidence for release review, not legal advice or legal
  certainty.

Current public-doc checks:

- `README.md` uses numerical-reference wording and says downstream users with
  strict provenance requirements should perform their own review.
- `reports/provenance/audit.md` says the mechanical report is not a legal
  conclusion.
- `reports/parity/dropin_compat.md` keeps compatibility claims scoped to
  measured surfaces.
