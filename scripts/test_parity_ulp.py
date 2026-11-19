"""Exercise binary64 boundaries, field policies, and deterministic reports."""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from typing import cast

from scripts.parity_ulp import (
    ROOT,
    build_report,
    compare_workload,
    main,
    ulp_distance,
)


class UlpTests(unittest.TestCase):
    def test_ulp_equal_is_zero(self) -> None:
        for value in (0.0, -0.0, 1.0, -1.0, sys.float_info.max):
            with self.subTest(value=value):
                self.assertEqual(ulp_distance(value, value), 0)

    def test_ulp_adjacent_is_one(self) -> None:
        for value in (1.0, -1.0, sys.float_info.min, -sys.float_info.min):
            with self.subTest(value=value):
                self.assertEqual(
                    ulp_distance(value, math.nextafter(value, math.inf)), 1
                )

    def test_ulp_signed_zero(self) -> None:
        self.assertEqual(ulp_distance(0.0, -0.0), 1)
        result = compare_workload("sample", {"edm": 0.0}, {"edm": -0.0}, {})[0]
        self.assertEqual(result.sign_mismatch, 1)
        self.assertTrue(result.first_failures)

    def test_ulp_subnormal_pair(self) -> None:
        tiny = math.ulp(0.0)
        self.assertEqual(ulp_distance(tiny, tiny * 2), 1)
        self.assertEqual(ulp_distance(-tiny, tiny), 3)
        result = compare_workload("sample", {"edm": tiny}, {"edm": tiny * 2}, {})[0]
        self.assertEqual(result.subnormal_in_ref, 1)

    def test_ulp_opposite_signs(self) -> None:
        self.assertEqual(ulp_distance(-1.0, 1.0), 2 * 0x3FF0000000000000 + 1)
        self.assertEqual(ulp_distance(1.0, -1.0), ulp_distance(-1.0, 1.0))

    def test_ulp_nan_class(self) -> None:
        self.assertEqual(ulp_distance(math.nan, -math.nan), 0)
        result = compare_workload("sample", {"edm": math.nan}, {"edm": -math.nan}, {})[
            0
        ]
        self.assertFalse(result.first_failures)

    def test_ulp_nan_number(self) -> None:
        self.assertIsNone(ulp_distance(math.nan, 1.0))
        self.assertIsNone(ulp_distance(1.0, math.nan))

    def test_ulp_infinity_largest_finite(self) -> None:
        self.assertIsNone(ulp_distance(math.inf, sys.float_info.max))
        self.assertEqual(ulp_distance(math.inf, math.inf), 0)
        self.assertEqual(ulp_distance(-math.inf, -math.inf), 0)
        self.assertIsNone(ulp_distance(math.inf, -math.inf))


class ComparisonTests(unittest.TestCase):
    def test_covariance_null_on_one_side(self) -> None:
        for reference, candidate in ((None, [[1.0]]), ([[1.0]], None)):
            results = compare_workload(
                "sample", {"covariance": reference}, {"covariance": candidate}, {}
            )
            by_path = {result.path: result for result in results}
            self.assertTrue(by_path["covariance"].structure)
            leaf = by_path["covariance[0][0]"]
            self.assertTrue(leaf.structure)
            self.assertTrue(leaf.first_failures)
            self.assertEqual(leaf.count, 0)

    def test_array_length_and_empty_container(self) -> None:
        results = compare_workload(
            "sample", {"params": [1.0]}, {"params": [1.0, 2.0], "extra": []}, {}
        )
        changed = {result.path for result in results if result.structure}
        self.assertEqual(changed, {"params", "params[1]", "extra"})

    def test_all_numeric_leaves_and_strings_skipped(self) -> None:
        reference = {
            "runner": "reference",
            "minos": {"parameter": 0, "lower": -1.0, "upper": 1.0, "valid": True},
            "params": [1.0],
            "covariance": [[2.0]],
            "errors": [3.0],
            "has_minos": True,
        }
        candidate = dict(reference, runner="candidate")
        results = compare_workload("sample", reference, candidate, {})
        self.assertEqual(sum(result.count for result in results), 8)
        self.assertFalse(any(result.first_failures for result in results))
        self.assertNotIn("runner", {result.path for result in results})

    def test_absolute_tolerance_boundary(self) -> None:
        for value, failed in ((0.25, False), (math.nextafter(0.25, math.inf), True)):
            result = compare_workload(
                "sample", {"edm": 0.0}, {"edm": value}, {"edm_abs": 0.25}
            )[0]
            self.assertEqual(bool(result.first_failures), failed)

    def test_no_tolerance_reports_distance_only(self) -> None:
        result = compare_workload("sample", {"errors": [1.0]}, {"errors": [2.0]}, {})[0]
        self.assertIsNone(result.tolerance)
        self.assertFalse(result.first_failures)
        self.assertEqual(result.max_abs_err, 1.0)
        self.assertEqual(result.max_rel_err, 0.5)

    def test_each_manifest_tolerance(self) -> None:
        for path, key in (
            ("fval", "fval_abs"),
            ("edm", "edm_abs"),
            ("params", "param_abs"),
            ("errors", "error_abs"),
            ("covariance", "cov_abs"),
        ):
            with self.subTest(path=path):
                result = compare_workload(
                    "sample", {path: 1.0}, {path: 2.0}, {key: 0.5}
                )[0]
                self.assertTrue(result.first_failures)
                self.assertIsNotNone(result.tolerance)
        results = compare_workload(
            "sample",
            {"minos": {"lower": -1.0, "upper": 1.0}},
            {"minos": {"lower": -2.0, "upper": 2.0}},
            {"minos_abs": 0.5},
        )
        self.assertTrue(all(result.first_failures for result in results))

    def test_integer_boolean_exactness_and_waiver_preserved(self) -> None:
        results = compare_workload(
            "sample",
            {"nfcn": 10, "valid": True, "parameter": 2**60},
            {"nfcn": 11, "valid": 1, "parameter": 2**60 + 1},
            {"nfcn_rel_warn": 0.5, "nfcn_rel_waiver": "Existing waiver"},
        )
        self.assertTrue(all(result.int_mismatch for result in results))
        counter = next(result for result in results if result.path == "nfcn")
        self.assertIsNotNone(counter.tolerance)
        self.assertEqual(cast("dict[str, object]", counter.tolerance)["rtol"], 0.5)
        self.assertEqual(
            cast("dict[str, object]", counter.tolerance)["waiver"], "Existing waiver"
        )
        integer = next(result for result in results if result.path == "parameter")
        self.assertEqual(integer.max_abs_err, 1)

    def test_nonfinite_metrics_are_json_safe(self) -> None:
        result = compare_workload("sample", {"edm": math.inf}, {"edm": 1.0}, {})[0]
        self.assertEqual(result.nonfinite_mismatch, 1)
        self.assertIsNone(result.max_abs_err)
        self.assertTrue(result.first_failures)
        _ = json.dumps(result.as_json(), allow_nan=False)

    def test_finite_extreme_relative_error(self) -> None:
        result = compare_workload(
            "sample", {"edm": sys.float_info.max}, {"edm": -sys.float_info.max}, {}
        )[0]
        self.assertIsNone(result.max_abs_err)
        self.assertEqual(result.max_rel_err, 2.0)


class ReportTests(unittest.TestCase):
    def test_report_inputs_manifest_scope_and_reproducibility(self) -> None:
        scratch = ROOT / ".tmp"
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as temporary:
            directory = Path(temporary)
            reference, candidate = directory / "ref", directory / "rust"
            reference.mkdir()
            candidate.mkdir()
            payload = '{"edm": 1.0, "nfcn": 10, "valid": true}\n'
            for side in (reference, candidate):
                _ = (side / "sample.json").write_text(payload)
                _ = (side / "ignored.json").write_text("invalid JSON")
            manifest = directory / "workloads.json"
            _ = manifest.write_text(
                json.dumps(
                    {"workloads": [{"id": "sample", "tolerances": {"edm_abs": 0.0}}]}
                )
            )
            report = build_report(reference, candidate, manifest)
            inputs = cast("dict[str, object]", report["inputs"])
            ref_inputs = cast("list[dict[str, str]]", inputs["reference"])
            self.assertEqual(len(ref_inputs), 1)
            self.assertEqual(
                ref_inputs[0]["sha256"], hashlib.sha256(payload.encode()).hexdigest()
            )
            self.assertFalse(Path(ref_inputs[0]["path"]).is_absolute())
            stats = cast("dict[str, object]", report["summary"])
            self.assertEqual(stats["fields"], 3)
            self.assertEqual(stats["values"], 3)
            output = directory / "report.json"
            arguments = [
                "--ref-dir",
                str(reference),
                "--rust-dir",
                str(candidate),
                "--workloads",
                str(manifest),
                "--report",
                str(output),
            ]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(arguments), 0)
                first = output.read_bytes()
                self.assertEqual(main(arguments), 0)
            self.assertEqual(first, output.read_bytes())
            _ = (candidate / "sample.json").write_text('{"edm": 2.0}')
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(arguments), 1)
            _ = (candidate / "sample.json").write_text("invalid JSON")
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(main(arguments), 2)


if __name__ == "__main__":
    _ = unittest.main()
