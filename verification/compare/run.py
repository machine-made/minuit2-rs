#!/usr/bin/env python3
"""Write reproducible matched-accuracy measurements from user-side counters."""

import csv
import importlib.util
import json
import math
import os
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SCRATCH = ROOT / "target/compare-scratch"
REPORTS = ROOT / "reports/compare"
GROUPS = ("Hessian", "EDM", "gradient", "bounds")
TOLERANCES = (1.0, 0.1, 0.01, 0.001, 0.0001)
BUNDLED_ROOT = "56be091a18cc07bc12242d65447d63611b4ebda4"


def command(args, **kwargs):
    return subprocess.run(
        list(map(str, args)),
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
        **kwargs,
    )


def uv_command(script, *arguments):
    return [
        "uv",
        "run",
        "--no-project",
        "--with",
        "iminuit==2.33.0",
        "--with",
        "mpmath",
        "--with",
        "numpy",
        "python3",
        "-I",
        script,
        *arguments,
    ]


def load_cases():
    spec = importlib.util.spec_from_file_location("compare_cases", HERE / "cases.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def scenarios(module):
    result = []
    for case in module.cases():
        for strategy in (0, 1, 2, 3):
            result.append(
                (
                    "Hessian",
                    case,
                    dict(case.spec, strategy=strategy, tolerance=0.1, gradient=False),
                )
            )
    by_name = {case.name: case for case in module.cases()}
    for name in ("rosenbrock", "dense2", "dense10", "quadratic2_fixed", "poisson10000"):
        case = by_name[name]
        for tolerance in TOLERANCES:
            result.append(
                (
                    "EDM",
                    case,
                    dict(case.spec, strategy=1, tolerance=tolerance, gradient=False),
                )
            )
    for name in ("quadratic2", "rosenbrock"):
        case = by_name[name]
        result.append(
            (
                "gradient",
                case,
                dict(case.spec, strategy=1, tolerance=0.0001, gradient=True),
            )
        )
    for case in module.bounds_cases():
        result.append(
            ("bounds", case, dict(case.spec, strategy=1, tolerance=0.1, gradient=False))
        )
    return result


def failure_message(error):
    message = (
        error.stderr if isinstance(error, subprocess.CalledProcessError) else str(error)
    )
    message = message.replace(str(ROOT), "<repo>")
    message = re.sub(r"/home/[^\s\"']+", "<path>", message)
    return " ".join(message.split())[-600:]


def rust_run(spec):
    if spec["strategy"] == 3:
        return {"error": "not implemented"}
    try:
        result = command(
            [
                "cargo",
                "run",
                "--release",
                "--locked",
                "--quiet",
                "--example",
                "compare_runner",
            ],
            input=json.dumps(spec),
        )
        return json.loads(result.stdout)
    except (subprocess.CalledProcessError, ValueError) as error:
        return {"error": failure_message(error)}


def root_reference():
    runner = ROOT / "third_party/root_ref_build/ref_runner/ref_runner"
    if not runner.is_file():
        raise RuntimeError("The built ROOT reference runner is missing")
    stored_row = json.loads(
        command([runner, "--workload", "quadratic2_no_g2_migrad"]).stdout
    )
    binary = SCRATCH / "root_counts"
    command(
        [
            "c++",
            "-std=c++17",
            "-O2",
            "-I" + str(ROOT / "third_party/root_ref/math/minuit2/inc"),
            HERE / "root_counts.cpp",
            ROOT / "third_party/root_ref_build/minuit2/lib/libMinuit2.a",
            ROOT / "third_party/root_ref_build/minuit2/lib/libMinuit2Math.a",
            "-lpthread",
            "-o",
            binary,
        ]
    )
    return stored_row, json.loads(command([binary]).stdout)


def covariance_error(case, covariance):
    if covariance is None:
        return None
    n = len(case.free)
    # Rust returns the free submatrix; iminuit includes zero rows for fixed parameters.
    indices = case.free if len(covariance) == len(case.reference) else list(range(n))
    if (
        len(indices) != n
        or len(covariance) not in (n, len(case.reference))
        or any(
            not isinstance(row, list) or len(row) != len(covariance)
            for row in covariance
        )
    ):
        raise ValueError("unexpected covariance dimensions")
    if any(
        value is None or not isinstance(value, (int, float)) or not math.isfinite(value)
        for row in covariance
        for value in row
    ):
        raise ValueError("covariance contains a nonfinite entry")
    reference = case.covariance
    return float(
        max(
            abs(covariance[indices[i]][indices[j]] - reference[i, j])
            / math.sqrt(float(reference[i, i] * reference[j, j]))
            for i in range(n)
            for j in range(n)
        )
    )


def metrics(module, case, group, result):
    mp = module.mp
    phase = "hesse" if group == "Hessian" else "migrad"
    snapshot = result.get(phase)
    if snapshot is None:
        return {
            "status": "not implemented"
            if result.get("error") == "not implemented"
            else "failed",
            "message": result.get("error", "missing result"),
            "value_calls": result.get("value_calls"),
            "gradient_calls": result.get("gradient_calls"),
            "hesse_calls": result.get("hesse_value_calls"),
            "hesse_gradient_calls": result.get("hesse_gradient_calls"),
        }
    valid = snapshot["valid"]
    params = snapshot["params"]
    finite = len(params) == len(case.reference) and all(
        value is not None and math.isfinite(value)
        for value in [*params, snapshot["fval"], snapshot["edm"]]
    )
    covariance_message = ""
    cov_error = None
    try:
        cov_error = covariance_error(case, snapshot["covariance"])
    except ValueError as error:
        covariance_message = str(error)
    distance = (
        float(case.value(list(map(mp.mpf, params)), True) - case.minimum)
        if finite
        else None
    )
    out = {
        "status": "ok" if valid and finite and not covariance_message else "failed",
        "message": covariance_message
        or ("" if valid and finite else "minimum is not valid"),
        "valid": valid,
        "fval": snapshot["fval"],
        "params": params,
        "covariance": snapshot["covariance"],
        "edm": snapshot["edm"],
        "true_distance": distance,
        "value_calls": snapshot["value_calls"],
        "gradient_calls": snapshot["gradient_calls"],
        "hesse_calls": result.get("hesse_value_calls"),
        "hesse_gradient_calls": result.get("hesse_gradient_calls"),
        "total_value_calls": result.get("hesse", {}).get("value_calls"),
        "total_gradient_calls": result.get("hesse", {}).get("gradient_calls"),
        "migrad_value_calls": result["migrad"]["value_calls"],
        "migrad_gradient_calls": result["migrad"]["gradient_calls"],
        "hesse_valid": result.get("hesse", {}).get("valid"),
        "hesse_message": result.get("error", ""),
        "cov_error": cov_error,
    }
    if group == "bounds":
        out["success"] = finite and all(
            abs(params[i] - float(case.reference[i])) < 1e-2 for i in case.free
        )
    return out


def reached(row, library):
    value = row.get(library + "_true_distance")
    return isinstance(value, (float, int)) and math.isfinite(value) and value <= 1e-6


def ratio(a, b):
    if a is None or b is None:
        return None
    if b == 0:
        return 1.0 if a == 0 else math.inf
    return a / b


def medians(rows, group):
    paired = [
        row
        for row in rows
        if row["group"] == group and reached(row, "rust") and reached(row, "iminuit")
    ]
    if group == "Hessian":
        paired = [row for row in paired if row["nonquadratic"]]
    count_key = "hesse_calls" if group == "Hessian" else "value_calls"
    calls = [
        ratio(row.get("rust_" + count_key), row.get("iminuit_" + count_key))
        for row in paired
    ]
    errors = [
        ratio(row.get("rust_cov_error"), row.get("iminuit_cov_error")) for row in paired
    ]
    calls = [value for value in calls if value is not None]
    errors = [value for value in errors if value is not None]
    return (
        statistics.median(calls) if calls else None,
        statistics.median(errors) if errors else None,
        len(paired),
    )


def verdict(rows, group):
    calls, errors, count = medians(rows, group)
    success_differs = group == "bounds" and any(
        row.get("rust_success") != row.get("iminuit_success")
        for row in rows
        if row["group"] == group
    )
    different = (
        (calls is not None and not 0.9 <= calls <= 1.1)
        or (errors is not None and not 0.5 <= errors <= 2)
        or success_differs
    )
    label = "difference" if different else "no difference"
    return (
        f"{label}; median Rust/iminuit calls = {display(calls)}; "
        f"median Rust/iminuit covariance errors = {display(errors)}; paired rows = {count}."
    )


def display(value):
    if value is None:
        return "not reached"
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, (list, dict)):
        return json.dumps(value, separators=(",", ":"), allow_nan=False)
    return str(value)


def paired_display(row, key):
    return " / ".join(
        display(row.get(library + "_" + key)) for library in ("rust", "iminuit")
    )


def matched_counts(rows):
    result = []
    for name in dict.fromkeys(row["case"] for row in rows if row["group"] == "EDM"):
        choices = {}
        for library in ("rust", "iminuit"):
            hit = next(
                (
                    row
                    for row in rows
                    if row["group"] == "EDM"
                    and row["case"] == name
                    and reached(row, library)
                ),
                None,
            )
            choices[library] = (
                (display(hit[library + "_value_calls"]), display(hit["tolerance"]))
                if hit
                else ("not reached", "not reached")
            )
        result.append((name, choices))
    return result


def write_reports(rows, version, root_counts, root_error):
    REPORTS.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (REPORTS / "matched_accuracy.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: display(value) for key, value in row.items()})
    commit = command(["git", "rev-parse", "--short", "HEAD"]).stdout.strip()
    lines = [
        f"iminuit: {version}.",
        f"Bundled Minuit2: ROOT 6.37.01 development snapshot, commit `{BUNDLED_ROOT}`.",
        "ROOT reference tag: `v6-36-08`.",
        f"Rust base commit: `{commit}`. The caller commits these harness files later.",
        "",
        "# Matched accuracy",
        "",
        f"The CSV contains {len(rows)} rows in a fixed order. Values and gradients are counted inside each user function.",
        "Each fit runs Migrad once, then Hesse. Explicit Hesse calls are the difference between the two counter readings.",
        "iminuit uses one Migrad iteration, with wrapper retries and Simplex disabled. Both sides use the same call budgets.",
        "Value calls determine call ratios. Gradient calls are shown separately. Failed fits retain their measurements and a message.",
        "True distance uses the objective evaluated at the returned parameters with 40-digit arithmetic, minus its reference minimum.",
        "Covariance error is the maximum absolute entry error divided by the geometric mean of the corresponding reference diagonal entries.",
        "Reference covariance is twice the error definition times the inverse reference Hessian on free parameters.",
        "The NLL error definition is 0.5. Other cases use 1.0. Fixed parameters are removed before covariance comparisons.",
        "A ratio with two zero errors is one. A positive error divided by zero is infinite. Empty medians say not reached.",
        "Verdicts include rows where both true distances are at most 1e-6. Each median uses those paired rows.",
        "A group says difference if its call median leaves [0.9, 1.1], its covariance median leaves [0.5, 2], or bounds success differs.",
        "Otherwise it says no difference. These measurements describe differences; they do not establish a general speed advantage.",
        "",
        "Gaussian data have seeded noise. Poisson bins contain expected counts, with totals of 10000 and 1000000 events.",
        "The seed is 20261113. All reference minima and Hessians use mpmath at 40 decimal digits.",
        "The one-sided bound distance scale is five units. Double limits also span five units.",
        "Bounds success means every free parameter is within 0.01 of its constrained reference minimum, regardless of the validity flag.",
        "",
        (
            "Version sources: [iminuit source tree](https://github.com/scikit-hep/iminuit/tree/v2.33.0/extern) and "
            f"[bundled ROOT header](https://github.com/root-project/root/blob/{BUNDLED_ROOT}/core/foundation/inc/ROOT/RVersion.hxx)."
        ),
        "",
    ]
    for group in GROUPS:
        lines += [f"## {group}", ""]
        if group == "Hessian":
            lines += [
                "Only non-quadratic cases enter this group's verdict: Rosenbrock, the Gaussian fit, and both Poisson fits.",
                "Strategy 3 rows are iminuit only. Rust strategy 3 is not implemented.",
                "The pinned iminuit setter rejects strategy 3. These rows use its bundled MnStrategy constructor and check the resulting level.",
                "",
            ]
        lines += [
            "Each paired cell gives Rust / iminuit. Full parameters, covariance matrices, and counter readings are in the CSV.",
            "",
            "| case | strategy | tolerance | status | value calls | gradient calls | Hesse calls | EDM | true distance | covariance error | success |",
            "|---|---|---|---|---|---|---|---|---|---|---|",
        ]
        for row in rows:
            if row["group"] != group:
                continue
            if row["case"] == "quadratic2_no_g2_migrad":
                lines += [
                    (
                        f"| quadratic2_no_g2_migrad | ROOT only | | {row['root_status']} | {row.get('root_value_calls', '')} | | | "
                        f"{display(row.get('root_edm'))} | {display(row.get('root_fval'))} | | |"
                    )
                ]
                continue
            cells = [row["case"], display(row["strategy"]), display(row["tolerance"])]
            cells += [
                paired_display(row, key)
                if key != "success" or group == "bounds"
                else ""
                for key in (
                    "status",
                    "value_calls",
                    "gradient_calls",
                    "hesse_calls",
                    "edm",
                    "true_distance",
                    "cov_error",
                    "success",
                )
            ]
            lines += ["| " + " | ".join(cells) + " |"]
        lines += [""]
        if group == "EDM":
            lines += [
                "Matched-accuracy calls use the first tolerance in the listed order that reaches a true distance of at most 1e-6.",
                "",
                "| case | Rust calls | Rust tolerance | iminuit calls | iminuit tolerance |",
                "|---|---|---|---|---|",
            ]
            for name, choices in matched_counts(rows):
                lines += [
                    "| "
                    + " | ".join((name, *choices["rust"], *choices["iminuit"]))
                    + " |"
                ]
            lines += [
                "",
                "The verdict uses all paired tolerance rows, following the same row rule as the other groups.",
                "",
            ]
        if group == "gradient":
            lines += [
                "The gradient rows use tolerance 0.0001 so both libraries can enter the matched-accuracy verdict.",
                "The ROOT-only Hessian row is copied from the existing built runner, using its quadratic2_no_g2_migrad workload.",
                "Its historical count is the runner's reported count. This row is excluded from the paired medians.",
                "",
            ]
        if group == "bounds":
            for name in ("quadratic2_y-1.0", "quadratic2_y-1.2"):
                row = next(row for row in rows if row["case"] == name)
                lines += [
                    (
                        f"{name}: Rust fval {display(row.get('rust_fval'))}, success {display(row.get('rust_success'))}; "
                        f"iminuit fval {display(row.get('iminuit_fval'))}, success {display(row.get('iminuit_success'))}."
                    )
                ]
            bounded = [
                row
                for row in rows
                if row["case"] in ("quadratic2_y-1.0", "quadratic2_y-1.2")
            ]
            if all(
                row.get(library + "_success")
                for row in bounded
                for library in ("rust", "iminuit")
            ):
                lines += [
                    "Both libraries succeed from both starts. The earlier failure is absent at this Rust base commit.",
                    "",
                ]
            else:
                lines += [
                    "The success columns above record whether failure depends on starting exactly on the upper limit.",
                    "",
                ]
        lines += [verdict(rows, group), ""]
    reference = next(
        row
        for row in rows
        if row["group"] == "Hessian"
        and row["case"] == "quadratic2"
        and row["strategy"] == 1
    )
    lines += [
        "## Reference counts",
        "",
        f"Rust measured {reference.get('rust_migrad_value_calls')} Migrad calls and {reference.get('rust_hesse_calls')} additional Hesse calls.",
        "Earlier measurements on cca0b9d were 31 Migrad calls and 13 additional Hesse calls.",
    ]
    if root_counts:
        lines += [
            f"ROOT measured {root_counts['migrad']} Migrad calls and {root_counts['total']} total calls after Hesse.",
            "Earlier ROOT measurements were 28 Migrad calls and 38 total calls after Hesse.",
            "The existing runner lacks this unbounded workload. A small counter probe links its existing ROOT libraries.",
        ]
    else:
        lines += [f"ROOT reference counts failed: {root_error}."]
    failures = [
        (row["group"], row["case"], library, row.get(library + "_message", ""))
        for row in rows
        for library in ("rust", "iminuit")
        if row.get(library + "_status") == "failed"
    ]
    lines += ["", "## Failed fits", ""]
    lines += (
        (
            ["| group | case | library | message |", "|---|---|---|---|"]
            + ["| " + " | ".join(item) + " |" for item in failures]
        )
        if failures
        else ["No paired fit failed."]
    )
    (REPORTS / "matched_accuracy.md").write_text("\n".join(lines) + "\n")


def main():
    started = time.monotonic()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("UV_CACHE_DIR", str(SCRATCH / "uv-cache"))
    if "--inside" not in sys.argv:
        for attempt in range(2):
            process = subprocess.run(
                list(map(str, uv_command(Path(__file__).resolve(), "--inside"))),
                cwd=ROOT,
                check=False,
            )
            if process.returncode == 0:
                return 0
            # Retry dependency acquisition only; completed worker failures are final.
            if process.returncode != 2:
                return process.returncode
        print("RESULT: BLOCKED: dependency acquisition failed twice", file=sys.stderr)
        return 2
    module = load_cases()
    work = scenarios(module)
    specs = [spec for _, _, spec in work]
    command(["cargo", "build", "--release", "--locked", "--example", "compare_runner"])
    python = json.loads(
        command(uv_command(HERE / "iminuit_runner.py"), input=json.dumps(specs)).stdout
    )
    rows = []
    for (group, case, spec), python_result in zip(work, python["results"], strict=True):
        rust = rust_run(spec)
        row = {
            "group": group,
            "case": case.name,
            "strategy": spec["strategy"],
            "tolerance": spec["tolerance"],
            "gradient": spec["gradient"],
            "nonquadratic": case.nonquadratic,
            "scope": "iminuit only" if spec["strategy"] == 3 else "paired",
            "reference_fval": float(case.minimum),
            "reference_params": list(map(float, case.reference)),
        }
        for library, result in (("rust", rust), ("iminuit", python_result)):
            row.update(
                {
                    library + "_" + key: value
                    for key, value in metrics(module, case, group, result).items()
                }
            )
        rows.append(row)
    root_counts, root_error = None, ""
    try:
        root, root_counts = root_reference()
        root_row = {
            "group": "gradient",
            "case": "quadratic2_no_g2_migrad",
            "strategy": 1,
            "tolerance": 0.1,
            "scope": "ROOT only",
            "rust_status": "ROOT only",
            "iminuit_status": "ROOT only",
            "root_status": "ok" if root["valid"] else "failed",
            "root_value_calls": root["nfcn"],
            "root_fval": root["fval"],
            "root_edm": root["edm"],
            "root_params": root["params"],
            "root_covariance": root["covariance"],
            "root_source": "built reference runner",
        }
    except (OSError, RuntimeError, subprocess.CalledProcessError, ValueError) as error:
        root_error = failure_message(error)
        root_row = {
            "group": "gradient",
            "case": "quadratic2_no_g2_migrad",
            "scope": "ROOT only",
            "rust_status": "ROOT only",
            "iminuit_status": "ROOT only",
            "root_status": "failed",
            "root_message": root_error,
        }
    index = next(i for i, row in enumerate(rows) if row["group"] == "bounds")
    rows.insert(index, root_row)
    write_reports(rows, python["version"], root_counts, root_error)
    print(f"rows: {len(rows)}")
    for group in GROUPS:
        print(f"{group}: {verdict(rows, group)}")
    print(f"run seconds: {time.monotonic() - started:.3f}")
    if root_error:
        print(f"RESULT: BLOCKED: ROOT rows incomplete: {root_error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
