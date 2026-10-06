"""
run_logger.py — Append-only run logs for CLAUD-IA WP1 scripts.

Each script calls start_logging() before doing any work. From then on all
output is mirrored to the terminal and appended to logs/<script>.log, so
every run is kept: the first line of each block is the date and time, and
previous runs are never overwritten.

Log block format:

    2026-10-05 14:32:10 | python scripts/download_tcga.py TCGA-PAAD
    ------------------------------------------------------------------
    ...output of the run...
    ------------------------------------------------------------------
    end 2026-10-05 15:04:51 | elapsed 32m 41s
"""

import atexit
import os
import re
import sys
import time
from datetime import datetime

sys.path.insert(0, os.path.dirname(__file__))
from config import PROJECT_ROOT

LOGS_DIR = os.path.join(PROJECT_ROOT, "logs")
RULE = "-" * 70

# Terminal colour/cursor escapes (ESC [ ... letter). Stripped from the
# log only: a file has no colours, and leaving them in hides messages
# from grep — "ERROR:" written in red does not start with an "E".
ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")


class _Tee:
    """Mirror a stream to the terminal and to the log file.

    Progress bars are collapsed on disk: when a chunk contains a carriage
    return, only the text after the last one is kept, so gdc-client still
    animates in the terminal but the log keeps just each line's final state.
    """

    def __init__(self, stream, logfile):
        self._stream = stream
        self._log = logfile
        self._pending = []

    def write(self, data):
        self._stream.write(data)   # terminal: keep the colours
        self._stream.flush()
        parts = ANSI_RE.sub("", data).split("\n")   # log: plain text
        for part in parts[:-1]:
            self._push(part)
            self._log.write("".join(self._pending) + "\n")
            self._pending = []
        self._push(parts[-1])
        self._log.flush()
        return len(data)

    def _push(self, text):
        if "\r" in text:
            self._pending = [text.rsplit("\r", 1)[1]]
        elif text:
            self._pending.append(text)

    def close_pending(self):
        """Flush a last line that never got its newline."""
        tail = "".join(self._pending).strip()
        self._pending = []
        if tail:
            self._log.write(tail + "\n")
            self._log.flush()

    def flush(self):
        self._stream.flush()
        self._log.flush()

    def isatty(self):
        return self._stream.isatty()


def _format_elapsed(seconds):
    seconds = int(seconds)
    h, m, s = seconds // 3600, (seconds % 3600) // 60, seconds % 60
    if h:
        return f"{h}h {m}m {s}s"
    if m:
        return f"{m}m {s}s"
    return f"{s}s"


def start_logging(name=None):
    """Append this run's output to logs/<script>.log. Returns the log path."""
    script = name or os.path.splitext(os.path.basename(sys.argv[0]))[0] or "run"
    os.makedirs(LOGS_DIR, exist_ok=True)
    path = os.path.join(LOGS_DIR, f"{script}.log")

    logfile = open(path, "a", encoding="utf-8")
    started = time.time()
    cmd = " ".join([f"python scripts/{script}.py", *sys.argv[1:]])
    logfile.write(f"{datetime.now():%Y-%m-%d %H:%M:%S} | {cmd}\n{RULE}\n")
    logfile.flush()

    real_out, real_err = sys.stdout, sys.stderr
    sys.stdout = _Tee(real_out, logfile)
    sys.stderr = _Tee(real_err, logfile)

    def finish():
        sys.stdout.close_pending()
        sys.stderr.close_pending()
        sys.stdout, sys.stderr = real_out, real_err
        logfile.write(f"{RULE}\nend {datetime.now():%Y-%m-%d %H:%M:%S} | "
                      f"elapsed {_format_elapsed(time.time() - started)}\n\n")
        logfile.close()

    atexit.register(finish)
    print(f"  Log: {os.path.relpath(path, PROJECT_ROOT)}")
    return path
