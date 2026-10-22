#!/usr/bin/env python3
"""Generate a conservative ROOT-vs-Rust similarity triage report.

This is not a legal conclusion. It is a mechanical screen for provenance review:
exact copied comments/strings, unusual token shingles, and mapped-file lexical
overlap that deserve human inspection.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROOT_REF = REPO_ROOT / "third_party" / "root_ref" / "math" / "minuit2"
DEFAULT_PARITY_CSV = REPO_ROOT / "reports" / "parity" / "functions.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "reports" / "provenance"

SOURCE_SUFFIXES = {".rs", ".h", ".hh", ".hpp", ".hxx", ".c", ".cc", ".cpp", ".cxx"}
RUST_KEYWORDS = {
    "as",
    "break",
    "const",
    "continue",
    "crate",
    "else",
    "enum",
    "extern",
    "false",
    "fn",
    "for",
    "if",
    "impl",
    "in",
    "let",
    "loop",
    "match",
    "mod",
    "move",
    "mut",
    "pub",
    "ref",
    "return",
    "self",
    "Self",
    "static",
    "struct",
    "super",
    "trait",
    "true",
    "type",
    "unsafe",
    "use",
    "where",
    "while",
}
CXX_KEYWORDS = {
    "alignas",
    "alignof",
    "and",
    "auto",
    "bool",
    "break",
    "case",
    "catch",
    "class",
    "const",
    "constexpr",
    "continue",
    "delete",
    "do",
    "double",
    "else",
    "enum",
    "explicit",
    "false",
    "float",
    "for",
    "friend",
    "if",
    "inline",
    "int",
    "long",
    "namespace",
    "new",
    "noexcept",
    "nullptr",
    "operator",
    "private",
    "protected",
    "public",
    "return",
    "short",
    "signed",
    "sizeof",
    "static",
    "struct",
    "switch",
    "template",
    "this",
    "throw",
    "true",
    "try",
    "typedef",
    "typename",
    "using",
    "virtual",
    "void",
    "while",
}
KEYWORDS = RUST_KEYWORDS | CXX_KEYWORDS

TOKEN_RE = re.compile(
    r"[A-Za-z_][A-Za-z0-9_]*"
    r"|\d+\.\d+(?:[eE][+-]?\d+)?"
    r"|\d+(?:[eE][+-]?\d+)?"
    r"|==|!=|<=|>=|&&|\|\||::|->|=>"
    r"|[-+*/%=&|^~!<>?:;.,()[\]{}]"
)
STRING_RE = re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', re.DOTALL)
COMMENT_RE = re.compile(r"//.*?$|/\*.*?\*/", re.DOTALL | re.MULTILINE)
WORD_RE = re.compile(r"[A-Za-z][A-Za-z0-9_]{2,}")


@dataclass(frozen=True)
class FileFeatures:
    path: Path
    rel: str
    text: str
    code_tokens: tuple[str, ...]
    generic_tokens: tuple[str, ...]
    identifiers: frozenset[str]
    strings: frozenset[str]
    comment_words: tuple[str, ...]


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def extract_comments(text: str) -> str:
    return "\n".join(match.group(0) for match in COMMENT_RE.finditer(text))


def strip_comments(text: str) -> str:
    return COMMENT_RE.sub(" ", text)


def strip_strings(text: str) -> str:
    return STRING_RE.sub(" STR ", text)


def normalize_string_literal(token: str) -> str:
    return token.strip("\"'").strip()


def file_features(path: Path, base: Path | None = None) -> FileFeatures:
    text = path.read_text(errors="ignore")
    comments = extract_comments(text)
    code = strip_strings(strip_comments(text))
    raw_tokens = TOKEN_RE.findall(code)
    code_tokens: list[str] = []
    generic_tokens: list[str] = []
    identifiers: set[str] = set()
    for token in raw_tokens:
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", token):
            if token not in KEYWORDS and len(token) > 1:
                identifiers.add(token.lower())
                code_tokens.append(token.lower())
                generic_tokens.append("ID")
            else:
                code_tokens.append(token.lower())
                generic_tokens.append(token.lower())
        elif re.fullmatch(r"\d+\.\d+(?:[eE][+-]?\d+)?|\d+(?:[eE][+-]?\d+)?", token):
            code_tokens.append("NUM")
            generic_tokens.append("NUM")
        else:
            code_tokens.append(token)
            generic_tokens.append(token)

    strings = {
        normalize_string_literal(match.group(0))
        for match in STRING_RE.finditer(text)
        if len(normalize_string_literal(match.group(0))) >= 8
    }
    comment_words = tuple(word.lower() for word in WORD_RE.findall(comments))
    rel_path = str(path.relative_to(base)) if base is not None else rel(path)
    return FileFeatures(
        path=path,
        rel=rel_path,
        text=text,
        code_tokens=tuple(code_tokens),
        generic_tokens=tuple(generic_tokens),
        identifiers=frozenset(identifiers),
        strings=frozenset(strings),
        comment_words=comment_words,
    )


def shingle_counts(tokens: tuple[str, ...], width: int) -> Counter[tuple[str, ...]]:
    if len(tokens) < width:
        return Counter()
    return Counter(tuple(tokens[i : i + width]) for i in range(len(tokens) - width + 1))


def overlap_count(a: Counter[tuple[str, ...]], b: Counter[tuple[str, ...]]) -> int:
    return sum((a & b).values())


def cosine(a: tuple[str, ...], b: tuple[str, ...]) -> float:
    if not a or not b:
        return 0.0
    ca = Counter(a)
    cb = Counter(b)
    dot = sum(count * cb[token] for token, count in ca.items())
    na = math.sqrt(sum(count * count for count in ca.values()))
    nb = math.sqrt(sum(count * count for count in cb.values()))
    if na == 0.0 or nb == 0.0:
        return 0.0
    return dot / (na * nb)


def jaccard(a: frozenset[str], b: frozenset[str]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def common_strings(a: FileFeatures, b: FileFeatures) -> list[str]:
    values = sorted(a.strings & b.strings)
    return values[:5]


def classify(row: dict[str, object]) -> str:
    comment_shingles = int(row["comment_shingles"])
    code_shingles = int(row["code_shingles"])
    generic_shingles = int(row["generic_code_shingles"])
    token_cos = float(row["token_cosine"])
    generic_cos = float(row["generic_token_cosine"])
    string_matches = int(row["shared_string_literals"])
    mapped = row["mapped_by_parity"] == "yes"

    if comment_shingles >= 2 or code_shingles >= 3 or string_matches >= 2:
        return "high"
    if mapped and (generic_shingles >= 25 or token_cos >= 0.22 or generic_cos >= 0.55):
        return "medium"
    if mapped and (generic_shingles >= 8 or token_cos >= 0.12):
        return "low-review"
    return "low"


def load_parity(path: Path) -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rust_files = [item.strip() for item in row.get("rust_file", "").split(";") if item.strip()]
            upstream_file = row.get("upstream_file", "").strip()
            if not upstream_file:
                continue
            for rust_file in rust_files:
                mapping[rust_file].add(upstream_file)
    return mapping


def collect_files(root: Path, suffixes: set[str]) -> list[Path]:
    return sorted(path for path in root.rglob("*") if path.is_file() and path.suffix in suffixes)


def compare_pair(
    rust: FileFeatures,
    upstream: FileFeatures,
    mapped_by_parity: bool,
    code_width: int,
    comment_width: int,
) -> dict[str, object]:
    code_shingles = overlap_count(
        shingle_counts(rust.code_tokens, code_width),
        shingle_counts(upstream.code_tokens, code_width),
    )
    generic_shingles = overlap_count(
        shingle_counts(rust.generic_tokens, code_width),
        shingle_counts(upstream.generic_tokens, code_width),
    )
    comment_shingles = overlap_count(
        shingle_counts(rust.comment_words, comment_width),
        shingle_counts(upstream.comment_words, comment_width),
    )
    strings = common_strings(rust, upstream)
    row: dict[str, object] = {
        "risk": "",
        "rust_file": rust.rel,
        "upstream_file": upstream.rel,
        "mapped_by_parity": "yes" if mapped_by_parity else "no",
        "rust_tokens": len(rust.code_tokens),
        "upstream_tokens": len(upstream.code_tokens),
        "token_cosine": round(cosine(rust.code_tokens, upstream.code_tokens), 4),
        "generic_token_cosine": round(cosine(rust.generic_tokens, upstream.generic_tokens), 4),
        "identifier_jaccard": round(jaccard(rust.identifiers, upstream.identifiers), 4),
        "shared_identifiers": len(rust.identifiers & upstream.identifiers),
        "code_shingles": code_shingles,
        "generic_code_shingles": generic_shingles,
        "comment_shingles": comment_shingles,
        "shared_string_literals": len(strings),
        "example_shared_strings": " | ".join(strings),
    }
    row["risk"] = classify(row)
    return row


def risk_rank(risk: str) -> int:
    return {"high": 0, "medium": 1, "low-review": 2, "low": 3}.get(risk, 9)


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "risk",
        "rust_file",
        "upstream_file",
        "mapped_by_parity",
        "rust_tokens",
        "upstream_tokens",
        "token_cosine",
        "generic_token_cosine",
        "identifier_jaccard",
        "shared_identifiers",
        "code_shingles",
        "generic_code_shingles",
        "comment_shingles",
        "shared_string_literals",
        "example_shared_strings",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(rows: list[dict[str, object]], path: Path, root_ref: Path, parity_csv: Path) -> None:
    counts = Counter(str(row["risk"]) for row in rows)
    mapped_counts = Counter(str(row["risk"]) for row in rows if row["mapped_by_parity"] == "yes")
    top_rows = rows[:25]
    high_rows = [row for row in rows if row["risk"] == "high"]
    medium_rows = [row for row in rows if row["risk"] == "medium"]

    lines: list[str] = []
    lines.append("# Similarity Audit")
    lines.append("")
    lines.append("This is a mechanical provenance triage report, not a legal conclusion.")
    lines.append("")
    lines.append("## Inputs")
    lines.append("")
    lines.append(f"- ROOT reference: `{root_ref.relative_to(REPO_ROOT)}`")
    lines.append(f"- Parity mapping: `{parity_csv.relative_to(REPO_ROOT)}`")
    lines.append("- Rust implementation scope: `src/**/*.rs`")
    lines.append("- ROOT scope: `third_party/root_ref/math/minuit2/**/*.{h,cxx,cpp,cc,c}`")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append("- Strip comments and string literals before code-token comparison.")
    lines.append("- Compare retained lexical tokens and generic tokens, where identifiers become `ID` and numbers become `NUM`.")
    lines.append("- Compare exact code token shingles, exact generic code shingles, exact comment word shingles, and string literals.")
    lines.append("- Mark ROOT/Rust pairs from `reports/parity/functions.csv` as mapped; also scan all unmapped pairs for stronger accidental matches.")
    lines.append("")
    lines.append("Risk labels are triage labels:")
    lines.append("")
    lines.append("- `high`: copied comments/strings or exact retained code-token shingles need immediate human review.")
    lines.append("- `medium`: mapped pair with enough lexical/structural overlap to inspect manually.")
    lines.append("- `low-review`: mapped pair with weak signals worth keeping in the inventory.")
    lines.append("- `low`: no mechanical similarity signal beyond normal algorithm/API vocabulary.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Pairs scanned: **{len(rows)}**")
    lines.append(f"- High risk: **{counts['high']}**")
    lines.append(f"- Medium risk: **{counts['medium']}**")
    lines.append(f"- Low-review: **{counts['low-review']}**")
    lines.append(f"- Low: **{counts['low']}**")
    lines.append(f"- Mapped high risk: **{mapped_counts['high']}**")
    lines.append(f"- Mapped medium risk: **{mapped_counts['medium']}**")
    lines.append("")
    lines.append("## High-Risk Findings")
    lines.append("")
    if high_rows:
        lines.append("| Rust file | ROOT file | Signals |")
        lines.append("|---|---|---|")
        for row in high_rows[:25]:
            signals = (
                f"comments={row['comment_shingles']}, code={row['code_shingles']}, "
                f"strings={row['shared_string_literals']}"
            )
            lines.append(f"| `{row['rust_file']}` | `{row['upstream_file']}` | {signals} |")
    else:
        lines.append("No high-risk exact-copy signals were found by this mechanical pass.")
    lines.append("")
    lines.append("## Medium-Risk Triage Queue")
    lines.append("")
    if medium_rows:
        lines.append("| Rust file | ROOT file | token cosine | generic cosine | generic shingles |")
        lines.append("|---|---|---:|---:|---:|")
        for row in medium_rows[:25]:
            lines.append(
                f"| `{row['rust_file']}` | `{row['upstream_file']}` | "
                f"{row['token_cosine']} | {row['generic_token_cosine']} | {row['generic_code_shingles']} |"
            )
    else:
        lines.append("No medium-risk mapped pairs were found by this mechanical pass.")
    lines.append("")
    lines.append("## Top Mechanical Similarity Rows")
    lines.append("")
    lines.append("| Risk | Rust file | ROOT file | mapped | token cosine | generic shingles | comments | strings |")
    lines.append("|---|---|---|---|---:|---:|---:|---:|")
    for row in top_rows:
        lines.append(
            f"| {row['risk']} | `{row['rust_file']}` | `{row['upstream_file']}` | "
            f"{row['mapped_by_parity']} | {row['token_cosine']} | "
            f"{row['generic_code_shingles']} | {row['comment_shingles']} | {row['shared_string_literals']} |"
        )
    lines.append("")
    lines.append("## Next Manual Review")
    lines.append("")
    lines.append("1. Inspect every `high` row first, if any.")
    lines.append("2. Inspect every `medium` row and decide `keep`, `document`, or `rewrite`.")
    lines.append("3. For algorithm hot spots, compare against papers/manuals rather than ROOT source when rewriting.")
    lines.append("4. Rerun this script after wording cleanup or rewrites.")
    lines.append("")
    lines.append("Current human triage notes, if present: `reports/provenance/manual_review.md`.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/generate_similarity_audit.py")
    lines.append("```")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate ROOT-vs-Rust similarity audit")
    parser.add_argument("--root-ref", default=str(DEFAULT_ROOT_REF))
    parser.add_argument("--parity-csv", default=str(DEFAULT_PARITY_CSV))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--code-width", type=int, default=12)
    parser.add_argument("--comment-width", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root_ref = Path(args.root_ref)
    parity_csv = Path(args.parity_csv)
    out_dir = Path(args.out_dir)
    if not root_ref.exists():
        raise FileNotFoundError(f"ROOT reference checkout not found: {root_ref}")
    if not parity_csv.exists():
        raise FileNotFoundError(f"parity csv not found: {parity_csv}")

    parity = load_parity(parity_csv)
    rust_paths = collect_files(REPO_ROOT / "src", {".rs"})
    upstream_paths = collect_files(root_ref, SOURCE_SUFFIXES - {".rs"})
    rust_features = {rel(path): file_features(path) for path in rust_paths}
    upstream_features = {str(path.relative_to(root_ref)): file_features(path, root_ref) for path in upstream_paths}

    mapped_pairs: set[tuple[str, str]] = set()
    for rust_file, upstream_files in parity.items():
        if rust_file not in rust_features:
            continue
        for upstream_file in upstream_files:
            if upstream_file in upstream_features:
                mapped_pairs.add((rust_file, upstream_file))

    rows: list[dict[str, object]] = []
    for rust_file, rust in rust_features.items():
        mapped_upstream = parity.get(rust_file, set())
        candidates = set(mapped_upstream)
        for upstream_file, upstream in upstream_features.items():
            if len(rust.code_tokens) < 40 or len(upstream.code_tokens) < 40:
                continue
            row = compare_pair(
                rust,
                upstream,
                (rust_file, upstream_file) in mapped_pairs,
                args.code_width,
                args.comment_width,
            )
            if (
                row["mapped_by_parity"] == "yes"
                or row["comment_shingles"]
                or row["code_shingles"]
                or row["shared_string_literals"]
                or row["generic_code_shingles"] >= 8
                or row["token_cosine"] >= 0.12
            ):
                rows.append(row)
            candidates.discard(upstream_file)
        for missing in sorted(candidates):
            rows.append(
                {
                    "risk": "low-review",
                    "rust_file": rust_file,
                    "upstream_file": missing,
                    "mapped_by_parity": "yes",
                    "rust_tokens": len(rust.code_tokens),
                    "upstream_tokens": 0,
                    "token_cosine": 0.0,
                    "generic_token_cosine": 0.0,
                    "identifier_jaccard": 0.0,
                    "shared_identifiers": 0,
                    "code_shingles": 0,
                    "generic_code_shingles": 0,
                    "comment_shingles": 0,
                    "shared_string_literals": 0,
                    "example_shared_strings": "upstream file not present in sparse checkout",
                }
            )

    rows.sort(
        key=lambda row: (
            risk_rank(str(row["risk"])),
            -(int(row["comment_shingles"]) + int(row["code_shingles"]) + int(row["shared_string_literals"])),
            -float(row["token_cosine"]),
            str(row["rust_file"]),
            str(row["upstream_file"]),
        )
    )
    write_csv(rows, out_dir / "similarity_inventory.csv")
    write_markdown(rows, out_dir / "audit.md", root_ref, parity_csv)
    print(f"Wrote {(out_dir / 'similarity_inventory.csv').relative_to(REPO_ROOT)}")
    print(f"Wrote {(out_dir / 'audit.md').relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
