"""
tests/test_two_engine_inputs_rya682.py
======================================
RYA-682 guards for the two-engine driver's inputs.

The driver is the Beta gate. Three ways its inputs could fail quietly, all now
loud:

  A. **Missing generated input.** `data/outputs/` is gitignored, so a clean
     checkout never has the Engine-B synthesis-v2 per-line table. The driver used
     to discover this only AFTER the whole Engine-A leg (GES linelist load, EW
     triage, MOOG baseline). It must now fail before any compute, naming the
     artifact and the regeneration command.

  B. **Present but UNUSABLE generated input.** This is the one that actually bit.
     `ispec/abundances.py:132` assigns a size-1 array into a scalar recarray slot;
     numpy deprecated that in 1.25 and made it an ERROR in 2.3. On numpy >= 2.3
     every element loses its atom code, every line is written `status='failed'`,
     and the run exits 0 having produced a full-length table with nothing in it.
     The frame is not empty, so the RYA-342 empty-set guard does not fire. Reading
     it yields an Engine-A-only result wearing the two-engine floor's name.

  C. **Missing committed input.** The seven dedicated Engine-B measurements are
     tracked in git and were read behind bare `if path.exists()`, so a missing one
     silently shrank the record set.

Plus the `--star` flag, which was parsed and never read: the Engine-B path was
hardwired to solar, so `--star procyon` would have emitted a SOLAR result under
another star's name.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import data_namespace as ns          # noqa: E402
from pipeline import two_engine_inputs as tei      # noqa: E402

OK_ROW = dict(element='Fe', ion='I', wavelength_air_A=5000.0, ew_mA=50.0,
              a_synth=7.5, red_chi2=1.2, status='ok', saturated=False)


def _write(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


# ── the canonical path is single-sourced, not hand-built ─────────────────────

def test_engine_b_path_is_the_rya469_namespaced_location():
    """The driver's input path must come from data_namespace, not a literal.

    Before RYA-682 the driver hand-built `data/outputs/solar/...` while the
    generator wrote via `data_namespace`. They happened to agree; nothing made
    them agree.
    """
    assert tei.engine_b_per_line_path('solar') == \
        ns.output_path('solar', 'per_line_synth_v2.csv', create=False)
    assert tei.engine_b_per_line_path('solar').parent.name == 'solar'
    assert tei.engine_b_per_line_path('solar').name == 'solar_per_line_synth_v2.csv'


def test_engine_b_artifact_is_generated_not_committed():
    """data/outputs/ is gitignored — this input is regenerable by design."""
    gitignore = (ROOT / '.gitignore').read_text().splitlines()
    assert any(l.strip().rstrip('/') == 'data/outputs' for l in gitignore), (
        "data/outputs/ is no longer gitignored — the committed-vs-generated "
        "decision for the Engine-B per-line table has changed and RYA-682's "
        "reasoning needs revisiting")


def test_regen_command_names_the_real_entry_point():
    cmd = tei.regen_command('solar')
    assert 'pipeline.abundances_derive' in cmd and 'synthesis-v2' in cmd
    assert (ROOT / 'pipeline' / 'abundances_derive.py').exists()


# ── A. missing generated input ───────────────────────────────────────────────

def test_missing_artifact_raises_with_the_regeneration_command(tmp_path, monkeypatch):
    monkeypatch.setattr(tei.ns, 'OUTPUTS_ROOT', tmp_path)
    with pytest.raises(SystemExit) as e:
        tei.assert_engine_b_artifact('solar')
    msg = str(e.value)
    assert 'MISSING GENERATED INPUT' in msg
    assert 'solar_per_line_synth_v2.csv' in msg
    assert 'synthesis-v2' in msg          # the regeneration command
    assert 'gitignored' in msg            # says WHY a clean checkout lacks it


# ── B. present-but-unusable generated input (the RYA-682 root cause) ─────────

def test_all_failed_artifact_raises_and_names_the_numpy_cause(tmp_path, monkeypatch):
    """A full-length table with zero status=='ok' rows must not be accepted."""
    monkeypatch.setattr(tei.ns, 'OUTPUTS_ROOT', tmp_path)
    path = tei.engine_b_per_line_path('solar')
    _write(path, [{**OK_ROW, 'status': 'failed', 'a_synth': None} for _ in range(129)])
    with pytest.raises(SystemExit) as e:
        tei.assert_engine_b_artifact('solar')
    msg = str(e.value)
    assert 'EMPTY GENERATED INPUT' in msg
    assert 'numpy' in msg and '2.3' in msg     # names the actual cause
    assert 'Engine-A-only' in msg              # names the downstream consequence


def test_usable_artifact_passes_and_reports_its_row_count(tmp_path, monkeypatch):
    monkeypatch.setattr(tei.ns, 'OUTPUTS_ROOT', tmp_path)
    path = tei.engine_b_per_line_path('solar')
    _write(path, [OK_ROW, {**OK_ROW, 'status': 'failed'}, {**OK_ROW, 'wavelength_air_A': 5001.0}])
    out = tei.assert_engine_b_artifact('solar')
    assert out['usable_rows'] == 2


def test_unparseable_artifact_is_reported_not_swallowed(tmp_path, monkeypatch):
    monkeypatch.setattr(tei.ns, 'OUTPUTS_ROOT', tmp_path)
    path = tei.engine_b_per_line_path('solar')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b'\x00\x01\x02 not,a;csv\x00')
    with pytest.raises(SystemExit):
        tei.assert_engine_b_artifact('solar')


# ── C. missing committed inputs ──────────────────────────────────────────────

def test_missing_committed_input_raises_and_says_restore_not_rerun(tmp_path):
    present = tmp_path / 'there.json'
    present.write_text('{}')
    with pytest.raises(SystemExit) as e:
        tei.assert_committed_inputs({'present': present,
                                     'Sr II synth (RYA-551)': tmp_path / 'gone.json'})
    msg = str(e.value)
    assert 'MISSING COMMITTED INPUT' in msg
    assert 'Sr II synth (RYA-551)' in msg
    assert 'git checkout' in msg          # restore, don't re-run a harness
    assert 'partial record set' in msg


def test_all_committed_inputs_present_passes():
    """The seven dedicated Engine-B inputs the driver declares are really tracked."""
    import subprocess
    paths = {
        'CNO': ROOT / 'data/audit/cno_synthesis/solar_phase_a_cross_arm.json',
        'Mn': ROOT / 'data/audit/mn_hfs_synthesis/solar_mn_hfs_synthesis_rya473.json',
        'CuV': ROOT / 'data/audit/cu_v_hfs_synthesis/solar_cu_v_hfs_synthesis_rya466.json',
        'Sr': ROOT / 'data/results/sr2_synthesis_rya551.json',
        'Zr560': ROOT / 'data/results/zr2_synthesis_rya560.json',
        'Zr585': ROOT / 'data/results/zr2_deblend_rya585.json',
        'Mg592': ROOT / 'data/results/mg_5528_synthesis_rya592.json',
    }
    tei.assert_committed_inputs(paths)                      # present in the tree
    for label, p in paths.items():                          # and actually tracked
        r = subprocess.run(['git', 'ls-files', '--error-unmatch', str(p)],
                           cwd=ROOT, capture_output=True)
        assert r.returncode == 0, f"{label}: {p} is not tracked in git"


# ── the driver wires the preflight in before any compute ─────────────────────

def test_driver_preflights_before_the_engine_a_leg():
    """Preflight must precede _engine_A_perline, or the check is worthless.

    Engine A costs minutes; a missing-input failure after it is a defect even
    though the message is correct.
    """
    src = (ROOT / 'scripts' / 'rya527_two_engine_run.py').read_text()
    i_pre = src.index('assert_engine_b_artifact')
    i_a = src.index('a_pl = _engine_A_perline')
    assert i_pre < i_a, "the Engine-B preflight must run before the Engine-A leg"


def test_driver_rejects_non_solar_star():
    """--star was parsed and never read; a non-solar run produced a solar result."""
    src = (ROOT / 'scripts' / 'rya527_two_engine_run.py').read_text()
    assert "args.star != 'solar'" in src


def test_driver_no_longer_hand_builds_the_output_path():
    src = (ROOT / 'scripts' / 'rya527_two_engine_run.py').read_text()
    assert "'data' / 'outputs'" not in src
    assert 'tei.engine_b_per_line_path' in src


# ── the numpy ceiling that prevents recurrence ───────────────────────────────

def _numpy_specifier(path: Path) -> str:
    for raw in path.read_text().splitlines():
        line = raw.split('#')[0].strip()
        if line.lower().startswith('numpy'):
            return line
    return ''


def test_requirements_pins_numpy_below_the_ispec_break():
    """The declared dependency set must not float onto a numpy that breaks iSpec.

    Ratified 2026-08-08. Without a ceiling the CI venv rebuilds onto the newest
    numpy; from 2.3 on, ispec/abundances.py:132 raises and every synthesis-v2 line
    is written status='failed' while the run exits 0. The guards elsewhere in this
    module make that loud — the ceiling is what stops it happening.
    """
    spec = _numpy_specifier(ROOT / 'requirements.txt')
    assert spec, "numpy is no longer declared in requirements.txt"
    assert '<2.3' in spec.replace(' ', ''), (
        f"requirements.txt numpy specifier is {spec!r} — the <2.3 ceiling is "
        f"load-bearing (RYA-682). Raising it requires fixing iSpec's "
        f"create_free_abundances_structure first.")


def test_numpy_ceiling_is_explained_where_it_is_declared():
    """A bare ceiling invites a future 'why is this pinned?' bump."""
    text = (ROOT / 'requirements.txt').read_text()
    assert 'RYA-682' in text and 'abundances.py' in text


def test_py39_lock_is_below_the_ceiling_too():
    spec = _numpy_specifier(ROOT / 'requirements-lock-py39.txt')
    assert spec.startswith('numpy==1.'), (
        f"the py39 lock pins {spec!r}; it must also sit below the 2.3 break")


# ── the pre-RYA-469 path drift is gone ───────────────────────────────────────

def test_rya305_reads_the_namespaced_per_line_tables():
    """RYA-305 read data/processed/ — the pre-namespacing location nothing writes.

    data/processed/ is gitignored AND no longer a per-line output dir, so a clean
    checkout hit a bare pandas FileNotFoundError while stale worktrees quietly
    served June copies.
    """
    src = (ROOT / 'scripts' / 'rya305_fe2_triage.py').read_text()
    assert "PROC / 'solar_per_line" not in src
    assert "PROC / 'solar_abundances_synth_v2.csv'" not in src
    assert 'engine_b_per_line_path' in src or 'ns.output_path' in src
