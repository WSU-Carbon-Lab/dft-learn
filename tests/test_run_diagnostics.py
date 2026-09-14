"""Tests for StoBe run diagnostics and reset helpers."""

from __future__ import annotations

from pathlib import Path

from dftlearn.io.run_diagnostics import (
    expected_stobe_output_path,
    reset_run_directory,
    sniff_stobe_output,
    write_run_diagnostic_report,
)


def test_expected_stobe_output_path_maps_category_tree() -> None:
    run_root = Path("/tmp/als3-001")
    run_file = run_root / "C1" / "C1gnd.run"
    assert expected_stobe_output_path(run_file) == run_root / "GND" / "C1gnd.out"


def test_sniff_stobe_output_detects_loader_error(tmp_path: Path) -> None:
    out = tmp_path / "C1gnd.out"
    out.write_text(
        "StoBe.x: error while loading shared libraries: libmkl_intel_lp64.so.2\n",
        encoding="utf-8",
    )
    sniff = sniff_stobe_output(out)
    assert sniff["exists"] is True
    assert sniff["error_line"] is not None
    assert "shared libraries" in str(sniff["error_line"])


def test_write_run_diagnostic_report_flags_failed_jobs(tmp_path: Path) -> None:
    (tmp_path / "C1").mkdir()
    (tmp_path / "GND").mkdir()
    (tmp_path / "C1" / "C1gnd.run").write_text("#!/bin/bash\n", encoding="utf-8")
    (tmp_path / "GND" / "C1gnd.out").write_text(
        "StoBe.x: error while loading shared libraries\n",
        encoding="utf-8",
    )
    text_path, json_path, report = write_run_diagnostic_report(tmp_path)
    assert text_path.is_file()
    assert json_path.is_file()
    assert report.jobs
    assert report.jobs[0].issues
    assert not report.ready_for_postprocess


def test_reset_run_directory_removes_outputs_keeps_inputs(tmp_path: Path) -> None:
    (tmp_path / "C1").mkdir()
    (tmp_path / "GND").mkdir()
    (tmp_path / "logs").mkdir()
    (tmp_path / "molConfig.py").write_text("x = 1\n", encoding="utf-8")
    (tmp_path / "C1" / "C1gnd.run").write_text("#!/bin/bash\n", encoding="utf-8")
    (tmp_path / "GND" / "C1gnd.out").write_text("fail\n", encoding="utf-8")
    removed = reset_run_directory(tmp_path)
    assert (tmp_path / "molConfig.py").is_file()
    assert (tmp_path / "C1" / "C1gnd.run").is_file()
    assert not (tmp_path / "GND").exists()
    assert "GND/" in removed
