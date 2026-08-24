"""RYA-821: the heartbeat must never report progress it has not measured."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("hb", ROOT / "scripts" / "sirius_heartbeat.py")
hb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hb)


def _p(pid=1, cpu=10, rchar=100, wchar=0):
    return dict(pid=pid, comm="interpol_multi_", etimes=60, cpu_s=cpu, pcpu=99.9,
                rss_kb=7000, rchar=rchar, wchar=wchar)


def test_one_sample_is_UNKNOWN_never_a_guess():
    """A process pinned at 99.9% CPU and one spinning on a bad loop are identical in a
    single ps. Reporting either as PROGRESSING would be inventing a measurement."""
    v = hb.verdict(None, dict(procs=[_p()], artifacts={}))
    assert v[0]["state"] == "UNKNOWN"
    assert "two" in v[0]["why"]


def test_burned_cpu_is_progress():
    v = hb.verdict(dict(procs=[_p(cpu=10)]), dict(procs=[_p(cpu=37)], artifacts={}))
    assert v[0]["state"] == "PROGRESSING"
    assert "+27s cpu" in v[0]["why"]


def test_alive_but_moving_nothing_is_STUCK():
    """The failure this exists to catch: alive != working."""
    v = hb.verdict(dict(procs=[_p()]), dict(procs=[_p()], artifacts={}))
    assert v[0]["state"] == "STUCK"


def test_io_alone_counts_as_progress():
    """An I/O-bound job burns no cpu; bytes moved still prove it is alive and working."""
    v = hb.verdict(dict(procs=[_p(rchar=100)]),
                   dict(procs=[_p(rchar=999)], artifacts={}))
    assert v[0]["state"] == "PROGRESSING"


def test_no_process_is_GONE_not_silence():
    v = hb.verdict(dict(procs=[_p()]), dict(procs=[], artifacts={}))
    assert v[0]["state"] == "GONE"


def test_watch_matches_comm_not_the_full_command_line():
    """pgrep -f matched the watcher's OWN ssh command line and hung forever
    (feedback_waiter_self_match). ps -eo comm cannot contain a wrapper."""
    # Check the COMMAND, not prose: the module deliberately documents the pgrep trap,
    # so a naive whole-file grep flags its own explanation.
    assert "comm" in hb._PS, "must select comm"
    assert "pgrep" not in hb._PS, "pgrep -f self-matches; use ps -eo comm"
    code = [l for l in (ROOT / "scripts" / "sirius_heartbeat.py").read_text().splitlines()
            if l.strip() and not l.lstrip().startswith(("#", '"', "*"))]
    assert not [l for l in code if "pgrep" in l and "=" in l and "not in" not in l], \
        "no executable line may build a pgrep command"
