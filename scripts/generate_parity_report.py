#!/usr/bin/env python3
"""
Generate function-level parity artifacts for minuit2-rs.

Outputs:
  - reports/parity/functions.csv
  - reports/parity/gaps.md
  - reports/parity/missing.csv
  - reports/parity/needs_review.csv

Method:
  1. Parse in-scope entries from INVENTORY.md
  2. Fetch upstream files from a pinned Minuit2 commit
  3. Extract upstream function-like symbols (heuristic parser)
  4. Extract Rust symbols from src/ (excluding cfg(test) blocks)
  5. Map upstream symbol -> Rust symbol with confidence buckets
"""

from __future__ import annotations

import csv
import argparse
import re
import subprocess
import sys
import textwrap
import urllib.error
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
INVENTORY_PATH = REPO_ROOT / "INVENTORY.md"
REPORT_DIR = REPO_ROOT / "reports" / "parity"
CACHE_DIR = REPO_ROOT / ".cache" / "parity"

DEFAULT_UPSTREAM_REPO = "root-project/root"
DEFAULT_UPSTREAM_SUBDIR = "math/minuit2"
DEFAULT_UPSTREAM_TAG = "v6-36-08"
IN_SCOPE_INVENTORY_ACTIONS = {"implement", "port"}


CPP_KEYWORDS = {
    "if",
    "for",
    "while",
    "switch",
    "catch",
    "return",
    "sizeof",
    "do",
}

ARCHITECTURAL_BASENAMES = {
    "MnApplication",
    "MnPrint",
    "MnPrintImpl",
    "FunctionMinimizer",
    "ModularFunctionMinimizer",
    "MinimumBuilder",
    "MinimumSeedGenerator",
    "MinimumErrorUpdator",
    "GradientCalculator",
    "VariableMetricEDMEstimator",
    "ScanBuilder",
    "ScanMinimizer",
    "CombinedMinimizer",
    "CombinedMinimumBuilder",
    "MnTraceObject",
    "MnConfig",
    "MnTiny",
    "MnEigen",
    "mnteigen",
    "MnVectorTransform",
    "MnParabolaFactory",
}

# Manual file-level mapping to keep the first matrix deterministic.
MANUAL_BASE_TO_RUST = {
    "AnalyticalGradientCalculator": ["src/gradient/analytical.rs"],
    "BFGSErrorUpdator": ["src/migrad/builder.rs"],
    "BasicFunctionGradient": ["src/minimum/gradient.rs"],
    "BasicFunctionMinimum": ["src/minimum/mod.rs", "src/minimum/error.rs", "src/user_parameter_state.rs"],
    "BasicMinimumError": ["src/minimum/error.rs"],
    "BasicMinimumParameters": ["src/minimum/parameters.rs"],
    "BasicMinimumSeed": ["src/minimum/seed.rs"],
    "BasicMinimumState": ["src/minimum/state.rs"],
    "CombinedMinimizer": ["src/minimize/mod.rs"],
    "CombinedMinimumBuilder": ["src/minimize/mod.rs"],
    "ContoursError": ["src/contours/contours_error.rs"],
    "DavidonErrorUpdator": ["src/migrad/builder.rs"],
    "FCNBase": ["src/fcn.rs"],
    "FCNGradientBase": ["src/fcn.rs"],
    "FunctionGradient": ["src/minimum/gradient.rs"],
    "FunctionMinimizer": ["src/migrad/minimizer.rs", "src/simplex/minimizer.rs"],
    "FunctionMinimum": ["src/minimum/mod.rs", "src/minimum/error.rs", "src/user_parameter_state.rs"],
    "GradientCalculator": ["src/gradient/mod.rs"],
    "HessianGradientCalculator": ["src/hesse/gradient.rs"],
    "InitialGradientCalculator": ["src/gradient/initial.rs"],
    "MinimumBuilder": ["src/migrad/builder.rs", "src/simplex/builder.rs", "src/minimize/mod.rs"],
    "MinimumError": ["src/minimum/error.rs"],
    "MinimumErrorUpdator": ["src/migrad/builder.rs"],
    "MinimumParameters": ["src/minimum/parameters.rs"],
    "MinimumSeed": ["src/minimum/seed.rs"],
    "MinimumSeedGenerator": ["src/migrad/seed.rs", "src/simplex/seed.rs"],
    "MinimumState": ["src/minimum/state.rs"],
    "MinosError": ["src/minos/minos_error.rs"],
    "MinuitParameter": ["src/parameter.rs"],
    "MnApplication": ["src/application.rs"],
    "MnConfig": ["src/precision.rs", "src/strategy.rs"],
    "MnContours": ["src/contours/mod.rs"],
    "MnCovarianceSqueeze": ["src/covariance_squeeze.rs"],
    "MnCross": ["src/minos/cross.rs"],
    "MnEigen": ["src/posdef.rs"],
    "MnFcn": ["src/mn_fcn.rs"],
    "MnFunctionCross": ["src/minos/function_cross.rs"],
    "MnGlobalCorrelationCoeff": ["src/global_cc.rs"],
    "MnHesse": ["src/hesse/mod.rs"],
    "MnLineSearch": ["src/linesearch.rs"],
    "MnMachinePrecision": ["src/precision.rs"],
    "MnMigrad": ["src/migrad/mod.rs"],
    "MnMinimize": ["src/minimize/mod.rs"],
    "MnMinos": ["src/minos/mod.rs"],
    "MnParabola": ["src/parabola.rs"],
    "MnParabolaFactory": ["src/parabola.rs"],
    "MnParabolaPoint": ["src/parabola.rs"],
    "MnParameterScan": ["src/scan/mod.rs"],
    "MnPosDef": ["src/posdef.rs"],
    "MnPrint": ["src/print.rs"],
    "MnPrintImpl": ["src/print.rs"],
    "MnScan": ["src/scan/mod.rs"],
    "MnSeedGenerator": ["src/migrad/seed.rs"],
    "MnSimplex": ["src/simplex/mod.rs"],
    "MnStrategy": ["src/strategy.rs"],
    "MnTiny": ["src/precision.rs"],
    "MnTraceObject": [],
    "MnUserCovariance": ["src/user_covariance.rs"],
    "MnUserFcn": ["src/mn_fcn.rs", "src/fcn.rs"],
    "MnUserParameterState": ["src/user_parameter_state.rs", "src/user_parameters.rs", "src/user_transformation.rs"],
    "MnUserParameters": ["src/user_parameters.rs"],
    "MnUserTransformation": ["src/user_transformation.rs", "src/user_parameters.rs"],
    "MnVectorTransform": ["src/user_transformation.rs"],
    "ModularFunctionMinimizer": ["src/migrad/minimizer.rs", "src/simplex/minimizer.rs"],
    "NegativeG2LineSearch": ["src/migrad/builder.rs"],
    "Numerical2PGradientCalculator": ["src/gradient/numerical.rs"],
    "ScanBuilder": ["src/scan/mod.rs"],
    "ScanMinimizer": ["src/scan/mod.rs"],
    "SimplexBuilder": ["src/simplex/builder.rs"],
    "SimplexMinimizer": ["src/simplex/minimizer.rs"],
    "SimplexParameters": ["src/simplex/parameters.rs"],
    "SimplexSeedGenerator": ["src/simplex/seed.rs"],
    "SinParameterTransformation": ["src/transform/sin.rs"],
    "SqrtLowParameterTransformation": ["src/transform/sqrt_low.rs"],
    "SqrtUpParameterTransformation": ["src/transform/sqrt_up.rs"],
    "VariableMetricBuilder": ["src/migrad/builder.rs"],
    "VariableMetricEDMEstimator": ["src/migrad/builder.rs"],
    "VariableMetricMinimizer": ["src/migrad/minimizer.rs"],
    "mnteigen": ["src/posdef.rs"],
}


# Per-basename explicit aliasing for naming drift between C++ and Rust APIs.
SYMBOL_ALIASES: dict[str, dict[str, list[str]]] = {
    "MnStrategy": {
        "GradientNCycles": ["grad_ncycles"],
        "GradientStepTolerance": ["grad_step_tol"],
        "GradientTolerance": ["grad_tol"],
        "HessianNCycles": ["hess_ncycles"],
        "HessianStepTolerance": ["hess_step_tol"],
        "HessianG2Tolerance": ["hess_g2_tol"],
        "HessianGradientNCycles": ["hess_grad_ncycles"],
        "StorageLevel": ["strategy"],
    },
    "BasicFunctionMinimum": {
        "HasAccurateCovar": ["is_accurate"],
        "HasPosDefCovar": ["is_pos_def"],
        "HasCovariance": ["has_covariance"],
        "HasValidCovariance": ["has_covariance"],
        "HesseFailed": ["hesse_failed"],
    },
    "FunctionMinimum": {
        "HasAccurateCovar": ["is_accurate"],
        "HasPosDefCovar": ["is_pos_def"],
        "HasCovariance": ["has_covariance"],
        "HasValidCovariance": ["has_covariance"],
        "HesseFailed": ["hesse_failed"],
    },
}


@dataclass(frozen=True)
class UpstreamFile:
    path: str
    basename: str


@dataclass(frozen=True)
class Symbol:
    name: str
    line: int
    param_count: int


@dataclass
class RustSymbol:
    name: str
    line: int
    file: str
    norm: str


def run(cmd: list[str]) -> str:
    out = subprocess.check_output(cmd, cwd=REPO_ROOT, text=True)
    return out.strip()


def normalize_repo(repo_arg: str) -> str:
    repo_arg = repo_arg.strip()
    m = re.search(r"github\.com[:/]+([^/]+/[^/.]+)", repo_arg)
    if m:
        return m.group(1)
    if "/" in repo_arg and not repo_arg.startswith("http"):
        return repo_arg
    raise ValueError(f"unsupported repo format: {repo_arg}")


def repo_git_url(repo_slug: str) -> str:
    return f"https://github.com/{repo_slug}.git"


def resolve_git_ref(repo_slug: str, ref: str) -> str:
    url = repo_git_url(repo_slug)
    # Resolve annotated tags to the peeled commit first.
    for query in [f"refs/tags/{ref}^{{}}", f"refs/tags/{ref}", f"refs/heads/{ref}", ref]:
        output = run(["git", "ls-remote", url, query])
        if output:
            return output.split()[0]
    raise RuntimeError(f"unable to resolve ref '{ref}' in {repo_slug}")


def parse_inventory_port_files() -> list[UpstreamFile]:
    text = INVENTORY_PATH.read_text()
    files: list[UpstreamFile] = []
    for line in text.splitlines():
        if not line.lstrip().startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 2:
            continue
        if cells[0].lower() in {"file", "------", "---"}:
            continue

        action = cells[1].strip("*` ").lower()
        # Older inventories used a `Port` action; the current inventory uses
        # `Implement` for files that remain in the symbol parity matrix.
        # `Replace` rows document nalgebra substitutions and are tracked in
        # report prose rather than as ROOT symbol-surface rows.
        if action not in IN_SCOPE_INVENTORY_ACTIONS:
            continue

        m = re.search(r"\[([^\]]+)\]", cells[0])
        if not m:
            continue
        rel = m.group(1).strip()
        if not (rel.startswith("inc/") or rel.startswith("src/")):
            continue
        basename = re.sub(r"\.(h|hpp|cxx|cpp|cc)$", "", rel.split("/")[-1])
        files.append(UpstreamFile(path=rel, basename=basename))
    # keep insertion order, drop exact duplicates
    unique: dict[tuple[str, str], UpstreamFile] = {}
    for f in files:
        unique[(f.path, f.basename)] = f
    return list(unique.values())


def download_upstream_file(repo_slug: str, subdir: str, commit: str, rel_path: str) -> Path:
    cache_path = CACHE_DIR / repo_slug.replace("/", "__") / subdir.replace("/", "__") / commit / rel_path
    if cache_path.exists():
        return cache_path
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    prefix = subdir.strip("/")
    full_path = f"{prefix}/{rel_path}" if prefix else rel_path
    url = f"https://raw.githubusercontent.com/{repo_slug}/{commit}/{full_path}"
    try:
        with urllib.request.urlopen(url) as r:  # nosec - fixed github raw host
            data = r.read()
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"failed to download {url}: {e.code}") from e
    cache_path.write_bytes(data)
    return cache_path


def strip_cpp_comments(source: str) -> str:
    source = re.sub(r"/\*.*?\*/", "", source, flags=re.S)
    source = re.sub(r"//.*", "", source)
    return source


def estimate_param_count(params: str) -> int:
    p = params.strip()
    if not p or p == "void":
        return 0
    # simple comma count that works for most Minuit2 signatures
    return p.count(",") + 1


def extract_cpp_symbols(path: Path, basename: str) -> list[Symbol]:
    text = path.read_text(errors="replace")
    stripped = strip_cpp_comments(text)
    lines = stripped.splitlines()
    is_source = path.suffix == ".cxx"

    # Standard function/method declaration/definition with an explicit return type.
    # Anchored to line start to avoid matching regular function calls.
    func_pat = re.compile(
        r"^\s*"
        r"(?:(?:inline|virtual|static|constexpr|friend|explicit|extern)\s+)*"
        r"[A-Za-z_][A-Za-z0-9_:<>\s*&~]*\s+"
        r"(?P<name>(?:~?[A-Za-z_][A-Za-z0-9_:~]*|operator\s*(?:\(\)|\[\]|[^\s(]+(?:\s+[^\s(]+)*)))\s*"
        r"\((?P<params>[^)]*)\)\s*"
        r"(?:(?:const|override|final|noexcept)\s*)*"
        r"(?:=\s*(?:0|default|delete)\s*)?"
        r"(?:;|\{|:|}\s*;?)\s*$"
    )
    # Constructor / destructor declarations (no explicit return type).
    ctor_pat = re.compile(
        r"^\s*(?P<name>[A-Za-z_][A-Za-z0-9_:~]*)\s*"
        r"\((?P<params>[^)]*)\)\s*"
        r"(?:(?:const|override|final|noexcept)\s*)*"
        r"(?:=\s*(?:0|default|delete)\s*)?"
        r"(?:;|\{|:|}\s*;?)\s*$"
    )

    # Build multi-line candidate signatures.
    candidates: list[tuple[int, str]] = []
    i = 0
    n = len(lines)
    depth = 0
    class_depth_stack: list[int] = []
    class_pending = False

    while i < n:
        raw = lines[i]
        line = raw.strip()
        if not line or line.startswith("#"):
            depth += raw.count("{") - raw.count("}")
            while class_depth_stack and depth < class_depth_stack[-1]:
                class_depth_stack.pop()
            i += 1
            continue

        if re.match(r"^\s*(class|struct)\b", raw) and not raw.rstrip().endswith(";"):
            # class definition starts now or on a following line.
            if "{" in raw:
                class_depth_stack.append(depth + raw.count("{"))
            else:
                class_pending = True

        if class_pending and "{" in raw:
            class_depth_stack.append(depth + raw.count("{"))
            class_pending = False

        in_class_scope = bool(class_depth_stack)
        allowed_scope = depth == 0 or in_class_scope

        if allowed_scope and "(" in line:
            # probable signature starts
            if re.match(r"^\s*(?:inline|virtual|static|constexpr|friend|explicit|extern|template|[A-Za-z_~])", raw):
                start = i + 1
                sig = raw.strip()
                # Multi-line declaration with return type above, e.g.
                #   virtual MinimumSeed
                #   operator()(...) const = 0;
                if re.match(r"^\s*operator", raw) and i > 0:
                    prev = lines[i - 1].strip()
                    if prev and "(" not in prev and not re.search(r"[;{}]$", prev):
                        sig = prev + " " + sig
                        start = i
                j = i
                while j + 1 < n and ")" not in sig:
                    j += 1
                    sig += " " + lines[j].strip()
                # consume trailing qualifiers / symbols if needed
                while j + 1 < n and not re.search(r"(;|\{|:|}\s*;?)\s*$", sig):
                    if re.match(r"^\s*(?:const|override|final|noexcept|->|\{|:|;|=)", lines[j + 1]):
                        j += 1
                        sig += " " + lines[j].strip()
                    else:
                        break
                candidates.append((start, sig))
                # advance depth across consumed lines
                for k in range(i, j + 1):
                    depth += lines[k].count("{") - lines[k].count("}")
                    while class_depth_stack and depth < class_depth_stack[-1]:
                        class_depth_stack.pop()
                i = j + 1
                continue

        depth += raw.count("{") - raw.count("}")
        while class_depth_stack and depth < class_depth_stack[-1]:
            class_depth_stack.pop()
        i += 1

    symbols: list[Symbol] = []
    seen_line_name: set[tuple[int, str]] = set()
    for line_no, sig in candidates:
        sig_trim = sig.strip()
        sig_match = re.sub(r"\{.*$", "{", sig_trim)

        if is_source and "{" not in sig_match and ":" not in sig_match:
            # In .cxx files keep definitions (and ctor init-list starts), skip declarations/initializations.
            continue

        m = func_pat.match(sig_match)
        if not m:
            m = ctor_pat.match(sig_match)
        if not m:
            continue
        raw_name = m.group("name").strip()
        short_name = raw_name.split("::")[-1].strip()
        if short_name in CPP_KEYWORDS:
            continue
        if short_name.startswith("~"):
            # Keep destructor only if it belongs to the target class basename.
            if short_name[1:] != basename:
                continue
        if not re.match(r"^[A-Za-z_]|^operator", short_name):
            continue
        # For constructor pattern, keep only if it is really this class ctor/dtor.
        if m.re is ctor_pat and short_name not in {basename, f"~{basename}"}:
            continue

        key = (line_no, short_name)
        if key in seen_line_name:
            continue
        seen_line_name.add(key)
        symbols.append(Symbol(name=short_name, line=line_no, param_count=estimate_param_count(m.group("params"))))

    # Deduplicate per basename by (name, arity), preferring earliest line.
    by_key: dict[tuple[str, int], Symbol] = {}
    for s in symbols:
        key = (s.name, s.param_count)
        old = by_key.get(key)
        if old is None or s.line < old.line:
            by_key[key] = s

    return sorted(by_key.values(), key=lambda s: (s.line, s.name))


def parse_replaces_mirrors_mapping() -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = defaultdict(set)
    for file in sorted((REPO_ROOT / "src").rglob("*.rs")):
        head = "\n".join(file.read_text(errors="replace").splitlines()[:40])
        for m in re.finditer(r"(Replaces|Mirrors)\s+([^\n]+)", head):
            desc = m.group(2)
            for token in re.findall(r"([A-Za-z0-9_./]+\.(?:h|hpp|cxx))", desc):
                base = token.split("/")[-1].rsplit(".", 1)[0]
                mapping[base].add(str(file.relative_to(REPO_ROOT)))
    return mapping


def normalize_symbol(name: str) -> str:
    return "".join(ch for ch in name.lower() if ch.isalnum())


def alternative_upstream_names(name: str, basename: str) -> list[str]:
    n = name.strip()
    alternatives = [n]

    if n.startswith("Get") and len(n) > 3:
        alternatives.append(n[3:])
    if n.startswith("Is") and len(n) > 2:
        alternatives.append(n[2:])
    if n.startswith("Has") and len(n) > 3:
        alternatives.append(n[3:])
    if n == basename:
        alternatives.extend(["new", "default"])
    if n.startswith("~") and n[1:] == basename:
        alternatives.extend(["drop"])
    if n.startswith("operator"):
        alternatives.extend(["call", "value"])

    base_aliases = SYMBOL_ALIASES.get(basename, {})
    if n in base_aliases:
        alternatives.extend(base_aliases[n])

    # preserve order and uniqueness
    out: list[str] = []
    seen: set[str] = set()
    for a in alternatives:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def extract_rust_symbols() -> tuple[list[RustSymbol], dict[str, list[RustSymbol]], dict[str, list[RustSymbol]]]:
    symbols: list[RustSymbol] = []
    by_file: dict[str, list[RustSymbol]] = defaultdict(list)
    by_norm: dict[str, list[RustSymbol]] = defaultdict(list)

    for file in sorted((REPO_ROOT / "src").rglob("*.rs")):
        rel = str(file.relative_to(REPO_ROOT))
        content = file.read_text(errors="replace")
        # Skip in-file test module declarations by cutting at first cfg(test).
        main_part = content.split("#[cfg(test)]", 1)[0]
        pattern = re.compile(r"\bfn\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(")
        for m in pattern.finditer(main_part):
            name = m.group(1)
            line = main_part.count("\n", 0, m.start()) + 1
            norm = normalize_symbol(name)
            rs = RustSymbol(name=name, line=line, file=rel, norm=norm)
            symbols.append(rs)
            by_file[rel].append(rs)
            by_norm[norm].append(rs)
    return symbols, by_file, by_norm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate function-level parity report artifacts")
    parser.add_argument(
        "--repo",
        default=DEFAULT_UPSTREAM_REPO,
        help=f"Upstream GitHub repo slug or URL (default: {DEFAULT_UPSTREAM_REPO})",
    )
    parser.add_argument(
        "--subdir",
        default=DEFAULT_UPSTREAM_SUBDIR,
        help=f"Subdirectory within upstream repo (default: {DEFAULT_UPSTREAM_SUBDIR})",
    )
    parser.add_argument(
        "--tag",
        default=DEFAULT_UPSTREAM_TAG,
        help=f"Upstream tag to resolve when --commit is not provided (default: {DEFAULT_UPSTREAM_TAG})",
    )
    parser.add_argument(
        "--commit",
        help="Upstream commit SHA or ref (overrides --tag)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_slug = normalize_repo(args.repo)
    ref = args.commit or args.tag
    commit = resolve_git_ref(repo_slug, ref)
    subdir = args.subdir.strip("/")
    upstream_files = parse_inventory_port_files()
    if not upstream_files:
        raise RuntimeError(
            f"no in-scope upstream files found in {INVENTORY_PATH.relative_to(REPO_ROOT)}; "
            "expected markdown rows with Implement or Port actions"
        )

    file_groups: dict[str, list[UpstreamFile]] = defaultdict(list)
    for f in upstream_files:
        file_groups[f.basename].append(f)

    replaces_map = parse_replaces_mirrors_mapping()
    _, rust_by_file, rust_by_norm = extract_rust_symbols()

    rows: list[dict[str, str]] = []
    extracted_symbol_count = 0

    for basename in sorted(file_groups):
        # Candidate Rust files from manual map + doc comments.
        candidates = set(MANUAL_BASE_TO_RUST.get(basename, []))
        candidates.update(replaces_map.get(basename, set()))
        candidates = {c for c in candidates if (REPO_ROOT / c).exists()}

        # Extract and dedupe symbols across all upstream files for this basename.
        extracted: dict[tuple[str, int], tuple[Symbol, str]] = {}
        fetch_errors: list[str] = []
        for uf in sorted(file_groups[basename], key=lambda u: u.path):
            try:
                local_path = download_upstream_file(repo_slug, subdir, commit, uf.path)
            except RuntimeError as e:
                fetch_errors.append(f"{uf.path}: {e}")
                continue
            syms = extract_cpp_symbols(local_path, basename)
            for sym in syms:
                key = (sym.name, sym.param_count)
                old = extracted.get(key)
                if old is None:
                    extracted[key] = (sym, uf.path)
                    continue
                old_sym, old_path = old
                # Prefer source definitions over headers if available.
                old_is_src = old_path.startswith("src/")
                new_is_src = uf.path.startswith("src/")
                if new_is_src and not old_is_src:
                    extracted[key] = (sym, uf.path)
                elif new_is_src == old_is_src and sym.line < old_sym.line:
                    extracted[key] = (sym, uf.path)

        extracted_symbol_count += len(extracted)

        if not extracted:
            rows.append(
                {
                    "upstream_repo": repo_slug,
                    "upstream_subdir": subdir,
                    "upstream_ref": ref,
                    "upstream_commit": commit,
                    "upstream_file": ";".join(sorted({f.path for f in file_groups[basename]})),
                    "upstream_symbol": "<no_symbol_extracted>",
                    "upstream_line": "",
                    "rust_file": ";".join(sorted(candidates)),
                    "rust_symbol": "",
                    "rust_line": "",
                    "status": "needs-review",
                    "rationale": (
                        "unable to extract symbols with heuristic parser"
                        if not fetch_errors
                        else f"source unavailable: {fetch_errors[0]}"
                    ),
                }
            )
            continue

        for key in sorted(extracted, key=lambda k: (k[0].lower(), k[1])):
            sym, source_path = extracted[key]
            alts = alternative_upstream_names(sym.name, basename)
            alt_norms = [normalize_symbol(a) for a in alts if a]

            candidate_matches: list[RustSymbol] = []
            # First pass: mapped files only.
            if candidates:
                for rust_file in sorted(candidates):
                    for rs in rust_by_file.get(rust_file, []):
                        if rs.norm in alt_norms:
                            candidate_matches.append(rs)

            # Second pass: global search only when no candidate file mapping exists.
            if not candidate_matches and not candidates:
                seen = set()
                for n in alt_norms:
                    for rs in rust_by_norm.get(n, []):
                        ident = (rs.file, rs.line, rs.name)
                        if ident not in seen:
                            seen.add(ident)
                            candidate_matches.append(rs)

            status = "needs-review"
            rationale = "symbol match ambiguous"
            rust_file = ""
            rust_symbol = ""
            rust_line = ""

            is_ctor = sym.name == basename
            is_dtor = sym.name == f"~{basename}"
            is_operator = sym.name.startswith("operator")

            if candidate_matches and len(candidate_matches) == 1:
                match = candidate_matches[0]
                rust_file = match.file
                rust_symbol = match.name
                rust_line = str(match.line)
                status = "implemented"
                rationale = "symbol name match"
            elif candidate_matches and len(candidate_matches) > 1:
                # Keep one evidence line, but mark for review.
                match = sorted(candidate_matches, key=lambda m: (m.file, m.line, m.name))[0]
                rust_file = match.file
                rust_symbol = match.name
                rust_line = str(match.line)
                status = "needs-review"
                rationale = "multiple Rust symbol candidates"
            else:
                if is_ctor or is_dtor or is_operator:
                    status = "intentionally-skipped"
                    rationale = "constructor/destructor/operator handled idiomatically in Rust"
                    # provide a best-effort constructor evidence if available
                    for rust_file_candidate in sorted(candidates):
                        ctor_candidates = [rs for rs in rust_by_file.get(rust_file_candidate, []) if rs.name in {"new", "default"}]
                        if ctor_candidates:
                            rs = ctor_candidates[0]
                            rust_file = rs.file
                            rust_symbol = rs.name
                            rust_line = str(rs.line)
                            break
                elif basename in ARCHITECTURAL_BASENAMES:
                    status = "needs-review"
                    rationale = "architectural refactor; no 1:1 symbol match"
                elif candidates:
                    status = "missing"
                    rationale = "no symbol match in mapped Rust files"
                else:
                    status = "needs-review"
                    rationale = "no mapped Rust file for upstream basename"

            rows.append(
                {
                    "upstream_repo": repo_slug,
                    "upstream_subdir": subdir,
                    "upstream_ref": ref,
                    "upstream_commit": commit,
                    "upstream_file": source_path,
                    "upstream_symbol": sym.name,
                    "upstream_line": str(sym.line),
                    "rust_file": rust_file,
                    "rust_symbol": rust_symbol,
                    "rust_line": rust_line,
                    "status": status,
                    "rationale": rationale,
                }
            )

    if not rows:
        raise RuntimeError("parity report produced zero rows; refusing to overwrite report artifacts")
    if extracted_symbol_count == 0:
        raise RuntimeError("extracted zero upstream symbols; refusing to overwrite report artifacts")

    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    csv_path = REPORT_DIR / "functions.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            lineterminator="\n",
            fieldnames=[
                "upstream_repo",
                "upstream_subdir",
                "upstream_ref",
                "upstream_commit",
                "upstream_file",
                "upstream_symbol",
                "upstream_line",
                "rust_file",
                "rust_symbol",
                "rust_line",
                "status",
                "rationale",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    missing_csv = REPORT_DIR / "missing.csv"
    with missing_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            lineterminator="\n",
            fieldnames=[
                "upstream_repo",
                "upstream_subdir",
                "upstream_ref",
                "upstream_commit",
                "upstream_file",
                "upstream_symbol",
                "upstream_line",
                "rust_file",
                "rust_symbol",
                "rust_line",
                "status",
                "rationale",
            ],
        )
        writer.writeheader()
        for row in rows:
            if row["status"] == "missing":
                writer.writerow(row)

    review_csv = REPORT_DIR / "needs_review.csv"
    with review_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            lineterminator="\n",
            fieldnames=[
                "upstream_repo",
                "upstream_subdir",
                "upstream_ref",
                "upstream_commit",
                "upstream_file",
                "upstream_symbol",
                "upstream_line",
                "rust_file",
                "rust_symbol",
                "rust_line",
                "status",
                "rationale",
            ],
        )
        writer.writeheader()
        for row in rows:
            if row["status"] == "needs-review":
                writer.writerow(row)

    by_status: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_file_missing: dict[str, int] = defaultdict(int)
    by_file_review: dict[str, int] = defaultdict(int)
    for row in rows:
        by_status[row["status"]].append(row)
        if row["status"] == "missing":
            by_file_missing[row["upstream_file"]] += 1
        if row["status"] == "needs-review":
            by_file_review[row["upstream_file"]] += 1

    total = len(rows)
    implemented = len(by_status["implemented"])
    missing = len(by_status["missing"])
    review = len(by_status["needs-review"])
    skipped = len(by_status["intentionally-skipped"])

    missing_top = sorted(by_file_missing.items(), key=lambda x: (-x[1], x[0]))[:20]
    review_top = sorted(by_file_review.items(), key=lambda x: (-x[1], x[0]))[:20]

    gaps_md = REPORT_DIR / "gaps.md"
    lines = [
        "# Function Parity Gaps",
        "",
        f"Upstream repo: `{repo_slug}`",
        f"Upstream subdir: `{subdir}`",
        f"Upstream ref: `{ref}`",
        f"Upstream commit: `{commit}`",
        "",
        "## Summary",
        "",
        f"- Total upstream symbols in scope: **{total}**",
        f"- `implemented`: **{implemented}**",
        f"- `missing`: **{missing}**",
        f"- `needs-review`: **{review}**",
        f"- `intentionally-skipped`: **{skipped}**",
        "",
        "## Top Files by `missing` Symbols",
        "",
    ]
    if missing_top:
        lines.extend([f"- `{path}`: {count}" for path, count in missing_top])
    else:
        lines.append("- None")

    lines.extend(["", "## Top Files by `needs-review` Symbols", ""])
    if review_top:
        lines.extend([f"- `{path}`: {count}" for path, count in review_top])
    else:
        lines.append("- None")

    lines.extend(["", "## Notes", ""])
    lines.extend(
        textwrap.dedent(
            """
            - Symbol extraction is heuristic (regex-based), not a full C++ parser.
            - `intentionally-skipped` currently captures constructor/destructor/operator-style symbols that map to Rust idioms.
            - `needs-review` includes architectural refactors where strict 1:1 symbol naming is not expected.
            - Use `reports/parity/functions.csv` as the source of truth for triage and manual confirmation.
            """
        )
        .strip()
        .splitlines()
    )
    gaps_md.write_text("\n".join(lines) + "\n")

    print(f"Wrote {csv_path.relative_to(REPO_ROOT)} ({len(rows)} rows)")
    print(f"Wrote {gaps_md.relative_to(REPO_ROOT)}")
    print(f"Wrote {missing_csv.relative_to(REPO_ROOT)}")
    print(f"Wrote {review_csv.relative_to(REPO_ROOT)}")
    print(f"Status counts: implemented={implemented} missing={missing} needs-review={review} intentionally-skipped={skipped}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise
