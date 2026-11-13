"""Download UCSC transcript tables (genePred) into a local cache."""

import os
import shutil
import ssl
import urllib.request

from .utils import log

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/{build}/database/{track}.txt.gz"
BUILD_ALIASES = {"grch38": "hg38", "hg38": "hg38", "grch37": "hg19", "hg19": "hg19", "b37": "hg19"}
TRACKS = ("ncbiRefSeq", "refGene", "ncbiRefSeqSelect", "knownGene")


def cache_dir():
    return os.path.expanduser(os.environ.get("VCFART_HOME", "~/.vcfart"))


def ucsc_build(build):
    b = BUILD_ALIASES.get((build or "").lower())
    if not b:
        raise ValueError(f"unknown build '{build}' (use hg19 or hg38)")
    return b


def cached_annotation(build):
    b = ucsc_build(build)
    for track in ("ncbiRefSeq", "refGene", "ncbiRefSeqSelect"):
        p = os.path.join(cache_dir(), f"{b}.{track}.txt.gz")
        if os.path.exists(p):
            return p
    return None


def fetch(build, track="ncbiRefSeq", outdir=None, insecure=False):
    b = ucsc_build(build)
    url = UCSC.format(build=b, track=track)
    outdir = os.path.expanduser(outdir or cache_dir())
    os.makedirs(outdir, exist_ok=True)
    dest = os.path.join(outdir, f"{b}.{track}.txt.gz")
    log(f"🌐 {url}", "fetch")
    ctx = ssl._create_unverified_context() if insecure else None
    tmp = dest + ".part"
    try:
        with urllib.request.urlopen(url, context=ctx, timeout=120) as r, open(tmp, "wb") as out:
            shutil.copyfileobj(r, out)
    except Exception as e:
        if os.path.exists(tmp):
            os.remove(tmp)
        hint = ("\n  behind an SSL-inspecting proxy? retry with --insecure, or download manually:\n"
                f"  wget -O {dest} {url}")
        raise RuntimeError(f"download failed: {e}{hint}")
    os.replace(tmp, dest)
    log(f"💾 {dest} ({os.path.getsize(dest) / 1e6:.1f} MB)", "fetch")
    return dest
