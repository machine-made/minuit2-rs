#!/usr/bin/env python3
"""
Generate executed-surface mapping and gate summary.

Joins C++ executed functions with the legacy->Rust traceability matrix and
produces:
  - reports/verification/executed_surface_mapping.md
  - reports/verification/executed_surface_gaps.csv
  - reports/verification/executed_surface_manifest.json

Gate rule (strict mode):
  fail if unmapped executed gaps include any P0 or P1 entries.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXECUTED_CSV = REPO_ROOT / "reports" / "verification" / "reference_coverage" / "executed_functions.csv"
DEFAULT_TRACEABILITY_CSV = REPO_ROOT / "reports" / "verification" / "traceability_matrix.csv"
DEFAULT_WAIVER_RULES_CSV = REPO_ROOT / "verification" / "traceability" / "waiver_rules.csv"
DEFAULT_REF_COVERAGE_MANIFEST = REPO_ROOT / "reports" / "verification" / "reference_coverage" / "manifest.json"
DEFAULT_OUT_MD = REPO_ROOT / "reports" / "verification" / "executed_surface_mapping.md"
DEFAULT_OUT_GAPS_CSV = REPO_ROOT / "reports" / "verification" / "executed_surface_gaps.csv"
DEFAULT_OUT_MANIFEST = REPO_ROOT / "reports" / "verification" / "executed_surface_manifest.json"

LOW_PRIORITY_WAIVERS = {
    "intentional",
    "architectural",
    "out-of-scope",
    "upstream-removed",
    "api-shape-drift",
}


@dataclass(frozen=True)
class SymbolInfo:
    symbol: str
    class_name: str
    is_constructor: bool
    is_destructor: bool
    is_operator: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate executed-surface mapping and gap report")
    parser.add_argument("--executed-csv", default=str(DEFAULT_EXECUTED_CSV))
    parser.add_argument("--traceability-csv", default=str(DEFAULT_TRACEABILITY_CSV))
    parser.add_argument("--waiver-rules-csv", default=str(DEFAULT_WAIVER_RULES_CSV))
    parser.add_argument("--reference-manifest", default=str(DEFAULT_REF_COVERAGE_MANIFEST))
    parser.add_argument("--out-md", default=str(DEFAULT_OUT_MD))
    parser.add_argument("--out-gaps-csv", default=str(DEFAULT_OUT_GAPS_CSV))
    parser.add_argument("--out-manifest", default=str(DEFAULT_OUT_MANIFEST))
    parser.add_argument("--strict-gate", action="store_true", help="Fail if any P0/P1 unmapped executed gaps remain")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"missing csv: {path}")
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def normalize_upstream_file(path: str) -> str:
    raw = path.replace("\\", "/").strip()
    marker = "/math/minuit2/"
    if marker in raw:
        return raw.split(marker, 1)[1]
    if raw.startswith("math/minuit2/"):
        return raw[len("math/minuit2/") :]
    return raw


def extract_mangled_name(raw_function: str) -> str:
    idx = raw_function.find("_Z")
    if idx >= 0:
        return raw_function[idx:]
    return raw_function


def demangle_cpp_symbols(symbols: list[str]) -> dict[str, str]:
    unique = list(dict.fromkeys(symbols))
    if not unique:
        return {}
    cxxfilt = shutil.which("c++filt")
    if not cxxfilt:
        return {s: s for s in unique}

    payload = "\n".join(unique) + "\n"
    proc = subprocess.run(
        [cxxfilt, "-n"],
        input=payload,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        return {s: s for s in unique}

    demangled = proc.stdout.splitlines()
    if len(demangled) != len(unique):
        return {s: s for s in unique}
    return dict(zip(unique, demangled))


def strip_templates(text: str) -> str:
    out = text
    while True:
        new = re.sub(r"<[^<>]*>", "", out)
        if new == out:
            return out
        out = new


def canonical_symbol(text: str) -> str:
    out = strip_templates(text).strip().replace(" ", "")
    if out.startswith("operator"):
        if out.startswith("operator()"):
            return "operator()"
        m = re.match(r"operator[^\s(]*", out)
        return m.group(0) if m else "operator"
    return out


def extract_function_prefix(demangled: str) -> str:
    text = demangled.strip()
    if not text:
        return text

    # Find the argument list by matching the final parenthesized segment from the right.
    close_idx = text.rfind(")")
    if close_idx == -1:
        return text

    depth = 0
    open_idx = -1
    for i in range(close_idx, -1, -1):
        ch = text[i]
        if ch == ")":
            depth += 1
        elif ch == "(":
            depth -= 1
            if depth == 0:
                open_idx = i
                break

    if open_idx == -1:
        return text
    return text[:open_idx].strip()


def split_cpp_scope(scope: str) -> list[str]:
    parts: list[str] = []
    chunk: list[str] = []
    depth = 0
    i = 0
    while i < len(scope):
        ch = scope[i]
        if ch == "<":
            depth += 1
        elif ch == ">" and depth > 0:
            depth -= 1

        if ch == ":" and i + 1 < len(scope) and scope[i + 1] == ":" and depth == 0:
            part = "".join(chunk).strip()
            if part:
                parts.append(part)
            chunk = []
            i += 2
            continue

        chunk.append(ch)
        i += 1

    tail = "".join(chunk).strip()
    if tail:
        parts.append(tail)
    return parts


def extract_symbol_info(demangled: str) -> SymbolInfo:
    signature_prefix = extract_function_prefix(demangled)
    parts = split_cpp_scope(signature_prefix)
    if not parts:
        return SymbolInfo(symbol=demangled, class_name="", is_constructor=False, is_destructor=False, is_operator=False)

    symbol = canonical_symbol(parts[-1])
    class_name = canonical_symbol(parts[-2]) if len(parts) >= 2 else ""
    is_operator = symbol.startswith("operator")
    is_constructor = bool(class_name) and symbol == class_name
    is_destructor = bool(class_name) and symbol == f"~{class_name}"
    return SymbolInfo(
        symbol=symbol,
        class_name=class_name,
        is_constructor=is_constructor,
        is_destructor=is_destructor,
        is_operator=is_operator,
    )


def upstream_basename(path: str) -> str:
    if not path:
        return ""
    part = path.replace("\\", "/").strip().split("/")[-1]
    return re.sub(r"\.(h|hpp|cxx|cc|cpp)$", "", part)


def index_traceability(
    rows: list[dict[str, str]]
) -> tuple[
    dict[tuple[str, str], list[dict[str, str]]],
    dict[str, list[dict[str, str]]],
    dict[str, list[dict[str, str]]],
]:
    by_key: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    by_file: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_basename: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        upstream_file = row.get("upstream_file", "").strip()
        upstream_files = [p.strip() for p in upstream_file.split(";") if p.strip()] or [upstream_file]
        upstream_symbol = row.get("upstream_symbol", "").strip().replace(" ", "")
        for upstream_path in upstream_files:
            by_key[(upstream_path, upstream_symbol)].append(row)
            by_file[upstream_path].append(row)
            base = upstream_basename(upstream_path)
            if base:
                by_basename[base].append(row)
    return by_key, by_file, by_basename


def read_waiver_rules(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    out: list[dict[str, str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {
            "raw_status",
            "rationale_contains",
            "upstream_file_regex",
            "upstream_symbol_regex",
            "waiver_type",
            "reason",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"waiver rules file missing columns: {sorted(missing)}")
        for row in reader:
            out.append({k: (row.get(k, "") or "").strip() for k in required})
    return out


def find_unmatched_waiver_rule(
    rules: list[dict[str, str]],
    upstream_file: str,
    upstream_symbol: str,
) -> dict[str, str] | None:
    for rule in rules:
        # Unmatched executed-surface gaps do not have a raw parity status/rationale context.
        # Only dedicated file/symbol rules (with empty raw_status/rationale_contains) apply here.
        if rule.get("raw_status", ""):
            continue
        if rule.get("rationale_contains", ""):
            continue
        file_re = rule.get("upstream_file_regex", "")
        symbol_re = rule.get("upstream_symbol_regex", "")
        if file_re and not re.search(file_re, upstream_file):
            continue
        if symbol_re and not re.search(symbol_re, upstream_symbol):
            continue
        return rule
    return None


def rank_status(rows: list[dict[str, str]]) -> str:
    statuses = [r.get("effective_status", "").strip() for r in rows]
    if "implemented" in statuses:
        return "implemented"
    if "unresolved" in statuses:
        return "unresolved"
    if "waived" in statuses:
        return "waived"
    return "missing"


def file_is_low_priority_waived(file_rows: list[dict[str, str]]) -> bool:
    if not file_rows:
        return False
    for row in file_rows:
        if row.get("effective_status", "").strip() != "waived":
            return False
        waiver_type = row.get("waiver_type", "").strip()
        if waiver_type and waiver_type not in LOW_PRIORITY_WAIVERS:
            return False
    return True


def classify_gap_priority(
    mapping_status: str,
    info: SymbolInfo,
    matched_rows: list[dict[str, str]],
    file_rows: list[dict[str, str]],
) -> str:
    if mapping_status == "implemented":
        return ""
    if mapping_status == "unresolved":
        return "P0"
    if mapping_status == "missing":
        if info.is_constructor or info.is_destructor or info.is_operator:
            return "P2"
        return "P1"

    waiver_types = {r.get("waiver_type", "").strip() for r in matched_rows if r.get("waiver_type", "").strip()}
    if not waiver_types:
        return "P1"
    if waiver_types <= LOW_PRIORITY_WAIVERS:
        return "P2"
    if "tooling" in waiver_types:
        return "P1"
    return "P1"


def load_workloads_from_reference_manifest(path: Path) -> list[str]:
    if not path.exists():
        return []
    data = json.loads(path.read_text())
    workloads = data.get("workloads")
    if isinstance(workloads, list):
        return [str(w) for w in workloads]
    return []


def main() -> int:
    args = parse_args()

    executed_rows = read_csv(Path(args.executed_csv))
    traceability_rows = read_csv(Path(args.traceability_csv))
    waiver_rules = read_waiver_rules(Path(args.waiver_rules_csv))
    trace_by_key, trace_by_file, trace_by_basename = index_traceability(traceability_rows)
    workload_ids = load_workloads_from_reference_manifest(Path(args.reference_manifest))

    mangled_names = [extract_mangled_name(row.get("function", "")) for row in executed_rows]
    demangled_map = demangle_cpp_symbols(mangled_names)

    gaps: list[dict[str, str]] = []
    mapped_count = 0
    priority_counts = Counter()
    file_gap_counts: dict[str, int] = defaultdict(int)

    for row in executed_rows:
        raw_function = row.get("function", "").strip()
        mangled = extract_mangled_name(raw_function)
        demangled = demangled_map.get(mangled, mangled)
        info = extract_symbol_info(demangled)

        upstream_file = normalize_upstream_file(row.get("file", ""))
        file_basename = upstream_basename(upstream_file)
        file_rows = trace_by_file.get(upstream_file, [])
        if not file_rows:
            file_rows = trace_by_basename.get(file_basename, [])
        key = (upstream_file, info.symbol)
        matches = trace_by_key.get(key, [])

        if not matches:
            key_lower = (upstream_file, info.symbol.lower())
            matches = [
                r
                for r in file_rows
                if r.get("upstream_symbol", "").strip().replace(" ", "").lower() == key_lower[1]
            ]
        if not matches:
            base_rows = trace_by_basename.get(file_basename, [])
            matches = [
                r
                for r in base_rows
                if r.get("upstream_symbol", "").strip().replace(" ", "") == info.symbol
            ]
            if not matches:
                matches = [
                    r
                    for r in base_rows
                    if r.get("upstream_symbol", "").strip().replace(" ", "").lower() == info.symbol.lower()
                ]
            if matches:
                file_rows = base_rows

        # For files that are fully waived for low-priority architectural/intentional reasons
        # (e.g. MnPrint), treat unmatched instantiations as waived at file level.
        if not matches and file_is_low_priority_waived(file_rows):
            matches = file_rows
        if not matches:
            rule = find_unmatched_waiver_rule(waiver_rules, upstream_file, info.symbol)
            if rule:
                matches = [
                    {
                        "effective_status": "waived",
                        "waiver_type": rule.get("waiver_type", ""),
                        "rationale": rule.get("reason", ""),
                        "rust_file": "",
                        "rust_symbol": "",
                    }
                ]

        mapping_status = rank_status(matches)
        if mapping_status == "implemented":
            mapped_count += 1
            continue

        priority = classify_gap_priority(mapping_status, info, matches, file_rows)
        if priority:
            priority_counts[priority] += 1
            file_gap_counts[upstream_file] += 1

        waiver_types = sorted(
            {r.get("waiver_type", "").strip() for r in matches if r.get("waiver_type", "").strip()}
        )
        rust_refs = sorted(
            {
                f"{r.get('rust_file', '').strip()}::{r.get('rust_symbol', '').strip()}"
                for r in matches
                if r.get("rust_file", "").strip() and r.get("rust_symbol", "").strip()
            }
        )
        rationale = sorted(
            {r.get("rationale", "").strip() for r in matches if r.get("rationale", "").strip()}
        )

        gaps.append(
            {
                "upstream_file": upstream_file,
                "upstream_symbol": info.symbol,
                "function_mangled": mangled,
                "function_demangled": demangled,
                "call_count": row.get("count", "0"),
                "mapping_status": mapping_status,
                "gap_priority": priority or "P2",
                "waiver_types": ";".join(waiver_types),
                "rust_refs": ";".join(rust_refs),
                "workload_ids": ";".join(workload_ids),
                "notes": " | ".join(rationale),
            }
        )

    gaps.sort(key=lambda r: (r["gap_priority"], r["upstream_file"], r["upstream_symbol"], r["function_mangled"]))

    out_gaps = Path(args.out_gaps_csv)
    out_gaps.parent.mkdir(parents=True, exist_ok=True)
    with out_gaps.open("w", newline="") as f:
        fieldnames = [
            "upstream_file",
            "upstream_symbol",
            "function_mangled",
            "function_demangled",
            "call_count",
            "mapping_status",
            "gap_priority",
            "waiver_types",
            "rust_refs",
            "workload_ids",
            "notes",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(gaps)

    total_executed = len(executed_rows)
    unmapped = len(gaps)
    gate_pass = priority_counts["P0"] == 0 and priority_counts["P1"] == 0
    top_gap_files = sorted(file_gap_counts.items(), key=lambda x: (-x[1], x[0]))[:15]

    manifest = {
        "executed_functions_total": total_executed,
        "mapped_implemented_total": mapped_count,
        "unmapped_total": unmapped,
        "priority_counts": {
            "P0": priority_counts["P0"],
            "P1": priority_counts["P1"],
            "P2": priority_counts["P2"],
        },
        "gate": {
            "rule": "P0 == 0 and P1 == 0",
            "pass": gate_pass,
        },
        "workloads": workload_ids,
        "artifacts": {
            "mapping_md": str(Path(args.out_md)),
            "gaps_csv": str(out_gaps),
        },
    }

    out_manifest = Path(args.out_manifest)
    out_manifest.write_text(json.dumps(manifest, indent=2) + "\n")

    lines: list[str] = []
    lines.append("# Executed Surface Mapping")
    lines.append("")
    lines.append("Join of reference executed C++ functions with traceability matrix mappings.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Executed C++ functions: **{total_executed}**")
    lines.append(f"- Mapped to implemented Rust symbols: **{mapped_count}**")
    lines.append(f"- Unmapped executed functions: **{unmapped}**")
    lines.append(
        f"- Unmapped priority split: P0={priority_counts['P0']}, "
        f"P1={priority_counts['P1']}, P2={priority_counts['P2']}"
    )
    lines.append(f"- Gate (`P0 == 0 and P1 == 0`): **{'PASS' if gate_pass else 'FAIL'}**")
    if workload_ids:
        lines.append(f"- Coverage workloads used: **{len(workload_ids)}**")
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    lines.append("- `reports/verification/executed_surface_mapping.md`")
    lines.append("- `reports/verification/executed_surface_gaps.csv`")
    lines.append("- `reports/verification/executed_surface_manifest.json`")
    lines.append("")
    lines.append("## Top Gap Files")
    lines.append("")
    if top_gap_files:
        for upstream_file, count in top_gap_files:
            lines.append(f"- `{upstream_file}`: {count}")
    else:
        lines.append("- none")
    lines.append("")
    lines.append("## Top P0/P1 Gaps")
    lines.append("")
    lines.append("| Priority | Upstream file | Symbol | Mapping status | Call count |")
    lines.append("|---|---|---|---|---|")
    shown = 0
    for gap in gaps:
        if gap["gap_priority"] not in {"P0", "P1"}:
            continue
        lines.append(
            f"| {gap['gap_priority']} | `{gap['upstream_file']}` | `{gap['upstream_symbol']}` | "
            f"`{gap['mapping_status']}` | {gap['call_count']} |"
        )
        shown += 1
        if shown >= 40:
            break
    if shown == 0:
        lines.append("| - | - | - | - | - |")
    lines.append("")

    out_md = Path(args.out_md)
    out_md.write_text("\n".join(lines) + "\n")

    print(f"Wrote {out_md.relative_to(REPO_ROOT)}")
    print(f"Wrote {out_gaps.relative_to(REPO_ROOT)}")
    print(f"Wrote {out_manifest.relative_to(REPO_ROOT)}")
    print(
        "Gate status: "
        f"P0={priority_counts['P0']} P1={priority_counts['P1']} "
        f"P2={priority_counts['P2']} pass={gate_pass}"
    )

    if args.strict_gate and not gate_pass:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
