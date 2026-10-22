# Measured Compatibility Report

This report records the compatibility claim for the closeout. It is not an
unqualified drop-in replacement promise. The release-facing target is measured
compatibility on the surfaces covered by generated parity artifacts, Rust tests,
Python tests, and the Python differential harness.

Numerical-output parity remains tracked separately in `reports/verification/`.

## Claim Boundary

- C++/Rust: `minuit2-rs` provides Minuit2 capability parity through a Rust-native
  API. It does not provide a literal C++ source-shape facade for ROOT Minuit2.
- Python: the PyO3 binding provides a measured `iminuit.Minuit` subset. The
  subset is checked by `python/compat/diff_iminuit.py`; deferred APIs are listed
  below and must not be described as implemented.
- Provenance: ROOT Minuit2 is the numerical reference and API comparison point,
  not a claim of source derivation or legal clearance.

## Surface 1: Rust Crate API vs ROOT Minuit2 C++

Generated parity artifacts currently report:

| Metric | Count |
|---|---:|
| Total upstream symbols in scope | 415 |
| Implemented | 217 |
| Missing | 0 |
| Needs review | 170 |
| Intentionally skipped | 28 |

Source of truth: `reports/parity/gaps.md` and `reports/parity/functions.csv`,
generated from ROOT Minuit2 `v6-36-08`.

The `needs-review` bucket is not a known blocker list. It includes matcher
ambiguity and architectural redistribution where ROOT's class-shaped API maps to
Rust builders, state objects, and transformation types. Examples:

| ROOT C++ shape | minuit2-rs shape |
|---|---|
| `MnApplication` facade mutators/accessors | `MnMigrad`/`MnSimplex` builders plus `MnUserParameters` and `MnUserTransformation` |
| CamelCase methods such as `SetValue` | Rust snake_case methods such as `set_value` |
| In-place mutable application objects | Builder and ownership-oriented Rust flow |
| Constructor/destructor/operator symbols | Rust idioms or intentionally skipped entries |

If literal C++ source compatibility is required later, that should be a separate
facade-design bead because it would intentionally mirror more of ROOT's class
shape. It is outside this closeout target.

## Surface 2: Python Binding vs `iminuit.Minuit`

The Python closeout target is a tested `iminuit`-compatible subset, not complete
coverage of every `iminuit` convenience API.

### Implemented and Measured

| Area | Current status |
|---|---|
| Construction | `Minuit(fcn, x=..., y=...)`; positional starts from function signature; `name=[...]`; signature introspection |
| Constants and tunables | `Minuit.LEAST_SQUARES`, `Minuit.LIKELIHOOD`, `errordef`, `strategy`, `tol` |
| Value/error/fixed/limit views | `values`, `errors`, `fixed`, `limits` support int and str indexing, assignment, iteration, length, and `to_dict()` |
| Minimization chaining | `migrad()`, `simplex()`, `hesse()`, and `minos()` return `self` |
| Fit status | `fval`, `valid`, `accurate`, `nfcn`, `npar`, `nfit`, `parameters` |
| Result objects | `fmin`, `params`, `merrors`, covariance as nested Python lists |
| Error and profile tools | `minos(*parameters)`, `profile`, `mnprofile`, `mncontour`, `contour` |
| State helpers | `reset()` and `fixto(key, value)` |
| Extra surface | `global_cc` is exposed even though `iminuit` has no direct property with that name |

The differential harness runs identical checks against `iminuit.Minuit` and
`minuit2.Minuit`. The current checks are:

| Check | Covered behavior |
|---|---|
| `construct_kwargs` | keyword construction, `migrad`, values, `fval` |
| `construct_positional` | positional starts with names from function signature |
| `values_by_index` | `values[0]` and `values[1]` |
| `values_setitem_persists` | view assignment persists on the object |
| `limits_setitem_persists` | limit assignment and readback |
| `errordef_constant` | class constant plus property assignment |
| `errors_and_valid` | error view and validity |
| `fixed_view_assign` | fixed view assignment before minimization |
| `limits_view_assign` | one-sided limits through the limits view |
| `start_at_limit` | fit can move from an exact active bound |
| `fmin_fields` | `fmin.is_valid` and `fmin.edm` |
| `params_objects` | `params` object fields |
| `counts` | `npar`, `nfit`, and `parameters` |
| `hesse_covariance` | covariance after `hesse()` |
| `minos_merrors` | `minos()` populates `merrors` |
| `profile` | `profile(..., subtract_min=True)` arrays |
| `mncontour` | closed `N x 2` contour array |
| `contour_grid` | grid-shaped `contour` output |
| `reset` | state reset to initial values |
| `fixto` | fix a parameter to a supplied value |

### Deferred Python APIs

These are explicitly out of scope for the current closeout unless promoted by a
new bead:

| Deferred API | Current behavior |
|---|---|
| `scan(ncall=None)` brute-force global minimizer | Raises `NotImplementedError`; use `profile` or `mnprofile` for 1D scans |
| `mncontour(..., cl=...)` confidence-level scaling | Raises `NotImplementedError`; the contour level follows `errordef` |
| `grad`, `g2`, `hessian` constructor callbacks | Not accepted by the constructor |
| Single-array constructor form | Not documented or covered by the harness |
| `covariance` Matrix helper methods | Covariance is exposed as nested lists, without `.correlation()` helpers |
| `precision`, `print_level`, `throw_nan`, `ngrad`, `init_params`, `pos2var`, `var2pos` | Not exposed |
| Plotting, `draw_*`, `visualize`, `interactive`, `scipy`, `iminuit.cost` | Not part of this crate's current Python surface |

## Verification Commands

Regenerate parity artifacts:

```bash
python3 scripts/generate_parity_report.py
```

Run Python compatibility checks after building the extension:

```bash
python python/compat/diff_iminuit.py
pytest python/tests/test_smoke.py python/tests/test_phase2.py -q
```

Run the stale-claim search from the Beads acceptance criteria after editing.
