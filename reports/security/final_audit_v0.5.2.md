# Final security audit — minuit2 v0.5.2

Date: 2026-06-11. Performed at the v0.5.2 release (project entering
maintenance mode). Auditor: Claude (Fable 5), run locally against commit
`21f2abc` plus the maintenance-mode working tree.

## Scope

Dependency advisories, license/source policy, unsafe code, secrets in the
tree, GitHub Actions supply chain, and published-crate contents. Out of
scope: full git-history secret scan, fuzzing, and formal review of the
numerical code (covered instead by the ROOT differential harness and the
NIST certified-oracle tests).

## Results

| Check | Tool / method | Result |
|-------|---------------|--------|
| Vulnerability advisories | `cargo audit` (137 deps) | No vulnerabilities. One allowed warning: `paste` 1.0.15 unmaintained (RUSTSEC-2024-0436), transitive via `nalgebra → simba`; compile-time proc-macro only, allowlisted with rationale in `deny.toml`. |
| License / bans / sources | `cargo deny check` | advisories ok, bans ok, licenses ok, sources ok. One cosmetic warning: unused `Unlicense` allowance in `deny.toml`. |
| Unsafe code | grep over `src/` + `unsafe_code = "warn"` lint | Zero `unsafe` blocks in `src/`, matching the README claim. |
| Secrets | grep for keys/tokens/passwords/private keys across `src/`, `python/`, `scripts/`, `.github/`, manifests | No credentials; all matches are lexer/parser "token" identifiers. |
| Actions supply chain | grep `uses:` across `.github/workflows/*.yml` | Every third-party action pinned by 40-char commit SHA (hardening from v0.4.2 still intact). |
| Published crate contents | `cargo package --list` (69 files) | No internal planning files (`VISION.md`, `docs/goals*`, ledgers, `.beads/`, `.sc/`) in the package; the `include` allowlist holds. |

## Conclusion

No vulnerabilities, no policy violations, no secrets, no unpinned actions.
The single advisory is an unmaintained transitive proc-macro with no
exploit path, already documented and allowlisted. Status: **pass**.

## Standing notes for maintenance mode

- Re-run `cargo audit` / `cargo deny check` when bumping any dependency;
  both run in CI.
- If `nalgebra` adopts a maintained `paste` replacement, drop the
  RUSTSEC-2024-0436 allowlist entry.
- Remove the unused `Unlicense` entry from `deny.toml` on the next
  housekeeping pass (cosmetic).
