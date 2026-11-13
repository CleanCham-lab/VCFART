"""Logging with MultiQC-style coloured, right-aligned step names (plain ANSI, no deps)."""

import os
import sys

QUIET = False
WARNINGS = []
COLOR = False
TAG_W = 12

PURPLE = "38;5;98"
LILAC = "38;5;141"
GREY = "38;5;245"
GREEN = "38;5;71"
BLUE = "38;5;39"
AMBER = "38;5;214"
RED = "38;5;203"
YELLOW = "38;5;220"
BOLD = "1"
DIM = "2"

# consequence -> colour, shared with the homepage legend
SEVERITY_COLOR = {
    "splice_acceptor_variant": RED, "splice_donor_variant": RED, "stop_gained": RED,
    "frameshift_variant": RED, "start_lost": RED, "stop_lost": RED,
    "missense_variant": AMBER, "inframe_insertion": AMBER, "inframe_deletion": AMBER,
    "protein_altering_variant": AMBER,
    "splice_donor_5th_base_variant": BLUE, "splice_region_variant": BLUE,
    "splice_donor_region_variant": BLUE, "splice_polypyrimidine_tract_variant": BLUE,
    "coding_sequence_variant": BLUE,
    "synonymous_variant": GREEN, "stop_retained_variant": GREEN,
}


def color_wanted(flag_off=False, stream=None):
    stream = stream or sys.stderr
    if flag_off or os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR") or os.environ.get("CLICOLOR_FORCE"):
        return True
    return hasattr(stream, "isatty") and stream.isatty() and os.environ.get("TERM") != "dumb"


def paint(text, *codes):
    if not COLOR or not codes:
        return str(text)
    return f"\033[{';'.join(codes)}m{text}\033[0m"


def _line(tag, msg, tag_color=PURPLE):
    t = tag.rjust(TAG_W)
    return f"{paint(t, tag_color)} {paint('|', GREY)} {msg}"


def log(msg, tag="vcfart"):
    if not QUIET:
        print(_line(tag, msg), file=sys.stderr)


def warn(msg):
    WARNINGS.append(msg)
    print(_line("warning", paint("⚠️  " + msg, YELLOW), YELLOW), file=sys.stderr)


def die(msg, code=1):
    print(_line("error", paint("❌ " + msg, RED), RED), file=sys.stderr)
    sys.exit(code)
