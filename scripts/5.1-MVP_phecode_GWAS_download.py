#!/usr/bin/env python3
"""
Download MVP sumstats from NCBI GIA index.

- Sequential downloads (no concurrency).
- MD5 checksum verification once on download for .tar files when *.md5 sidecar is available.
- Other file types (e.g. .txt) are downloaded without checksum verification.
- Existing files are trusted on subsequent runs; no re-verification is done.
- Missing MD5 for .tar files is treated as "no checksum available", not as a failure.
- Checksum mismatches are logged and recorded but do not abort the run.
- No resume logic; each run cleans partial files per file as needed.
- Exponential backoff retry logic for index fetches and downloads.
- Logging to a timestamped log file in the destination directory.

Usage example:

python scripts/5.1-MVP_phecode_GWAS_download.py --base-url "https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/GIA/" --dest "/mnt/sdg/robert/deepwidegp/ARD-DeepWideGP/data/5-MVP_download"

Dependencies:
  pip install requests tqdm
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import time
from typing import Dict, List, Optional

import requests
from tqdm import tqdm

DEFAULT_PREFIX = "phs002453.MVP_R4.1000G_AGR.GIA.PheCodes_"

# Retry behaviour
MAX_RETRIES = 5
BACKOFF_BASE_SECONDS = 2.0

# Shared HTTP session
session = requests.Session()
session.headers.update({"User-Agent": "mvp-sumstats-downloader/2.3"})


def setup_logging(dest_dir: str, verbose_console: bool = True) -> logging.Logger:
    """Configure logging to file (and optionally to console)."""
    os.makedirs(dest_dir, exist_ok=True)
    logger = logging.getLogger("mvp_downloader")
    logger.setLevel(logging.INFO)

    # Avoid duplicate handlers if setup_logging is called twice.
    if logger.handlers:
        return logger

    ts = time.strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(dest_dir, f"mvp_download_{ts}.log")
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    if verbose_console:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    logger.info("Logging initialised. Log file: %s", log_path)
    return logger


def _sleep_with_backoff(attempt: int, logger: Optional[logging.Logger], context: str) -> None:
    delay = BACKOFF_BASE_SECONDS * (2 ** (attempt - 1))
    if logger:
        logger.info("%s: retrying in %.1f seconds (attempt %d)", context, delay, attempt)
    time.sleep(delay)


def list_index_files(base_url: str, logger: Optional[logging.Logger] = None) -> List[str]:
    """Return all filename URLs from a simple HTTP directory index, with retries."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            if logger:
                logger.info("Fetching index from %s (attempt %d)", base_url, attempt)
            resp = session.get(base_url, timeout=30)
            resp.raise_for_status()
            html = resp.text

            # crude but effective: fetch all hrefs
            hrefs = re.findall(r'href="([^"]+)"', html)
            links = []
            for h in hrefs:
                if h.startswith("#"):
                    continue
                url = requests.compat.urljoin(base_url, h)
                links.append(url)
            if logger:
                logger.info("Found %d links at index.", len(links))
            return links
        except Exception as e:
            if logger:
                logger.error("Error fetching index from %s: %s", base_url, e)
            if attempt == MAX_RETRIES:
                raise
            _sleep_with_backoff(attempt, logger, "Index fetch")


def compute_md5(path: str, chunk_size: int = 8 * 1024 * 1024) -> str:
    """Compute MD5 checksum of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


_md5_cache: Dict[str, Optional[str]] = {}


def fetch_remote_md5(file_url: str, logger: Optional[logging.Logger] = None) -> Optional[str]:
    """
    Fetch remote MD5 for a given data file URL.

    Assumes checksum is at file_url + ".md5" and parses the first token
    of the first line. If the MD5 cannot be obtained or is invalid,
    log a warning and return None (do not abort the run).
    """
    if file_url in _md5_cache:
        return _md5_cache[file_url]

    md5_url = file_url + ".md5"
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            if logger:
                logger.info("Fetching MD5 from %s (attempt %d)", md5_url, attempt)
            resp = session.get(md5_url, timeout=30)
            resp.raise_for_status()
            text = resp.text.strip()

            if not text:
                raise ValueError("Empty MD5 file")

            first_line = text.splitlines()[0]
            md5 = first_line.split()[0]
            if not re.fullmatch(r"[0-9a-fA-F]{32}", md5):
                raise ValueError(f"Invalid MD5 format in {md5_url!r}: {first_line!r}")

            md5 = md5.lower()
            _md5_cache[file_url] = md5
            if logger:
                logger.info("Remote MD5 for %s is %s", file_url, md5)
            return md5
        except Exception as e:
            if logger:
                logger.warning(
                    "Error fetching/parsing MD5 from %s on attempt %d: %s",
                    md5_url,
                    attempt,
                    e,
                )
            if attempt == MAX_RETRIES:
                if logger:
                    logger.warning(
                        "No valid MD5 could be obtained for %s; proceeding without checksum verification.",
                        file_url,
                    )
                _md5_cache[file_url] = None
                return None
            _sleep_with_backoff(attempt, logger, "MD5 fetch")

    # Defensive; loop always returns.
    _md5_cache[file_url] = None
    return None


def download_file_with_retries(
    url: str,
    tmp_path: str,
    logger: Optional[logging.Logger] = None,
) -> float:
    """
    Download a file to tmp_path with retries and a single progress bar.
    Returns elapsed seconds if successful; raises on final failure.
    """
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            if logger:
                logger.info("Starting download of %s to %s (attempt %d)", url, tmp_path, attempt)
            start_time = time.time()
            with session.get(url, stream=True, timeout=60) as r:
                r.raise_for_status()

                total = r.headers.get("Content-Length")
                try:
                    total = int(total) if total is not None else None
                except Exception:
                    total = None

                filename = os.path.basename(tmp_path)
                with open(tmp_path, "wb") as f, tqdm(
                    total=total,
                    unit="B",
                    unit_scale=True,
                    desc=filename,
                ) as pbar:
                    for chunk in r.iter_content(chunk_size=8_388_608):
                        if not chunk:
                            continue
                        f.write(chunk)
                        pbar.update(len(chunk))

            elapsed = time.time() - start_time
            if logger:
                size_bytes = os.path.getsize(tmp_path)
                logger.info(
                    "Completed download of %s (%d bytes) in %.1f seconds",
                    url,
                    size_bytes,
                    elapsed,
                )
            return elapsed
        except Exception as e:
            if logger:
                logger.error("Error downloading %s on attempt %d: %s", url, attempt, e)
            if attempt == MAX_RETRIES:
                raise
            _sleep_with_backoff(attempt, logger, "File download")

    # Defensive; loop always returns or raises.
    raise RuntimeError(f"Failed to download {url}")


def process_one_file(
    entry: Dict,
    dest_dir: str,
    logger: Optional[logging.Logger] = None,
) -> Dict:
    """
    Process a single file:
    - If final file exists, trust it and skip (no checksum re-verification).
    - Otherwise:
        - Remove leftover partial.
        - If filename ends with .tar, fetch remote MD5 (if available).
          For other types, skip checksum logic entirely.
        - Download to .part with retries.
        - Promote to final.
        - For .tar with MD5 available, verify after download; on mismatch, delete and record failure.
    """
    url = entry["url"]
    fname = entry["filename"]
    dest_path = os.path.join(dest_dir, fname)
    tmp_path = dest_path + ".part"

    if logger:
        logger.info("Processing file %s (%s)", fname, url)

    # If final file already exists, trust it and skip any checksum work.
    if os.path.exists(dest_path):
        size_bytes = os.path.getsize(dest_path)
        if logger:
            logger.info(
                "Existing file %s found; trusting existing file and skipping download.",
                dest_path,
            )
        return {
            "filename": fname,
            "url": url,
            "size_bytes": size_bytes,
            "elapsed_seconds": 0.0,
            "status": "skipped_existing_trusted",
        }

    # Clean up any leftover partial file
    if os.path.exists(tmp_path):
        try:
            os.remove(tmp_path)
            if logger:
                logger.info("Removed leftover partial file %s", tmp_path)
        except Exception as e:
            if logger:
                logger.error("Failed to remove leftover partial %s: %s", tmp_path, e)

    # Determine whether we should use checksums
    remote_md5: Optional[str] = None
    if fname.endswith(".tar"):
        remote_md5 = fetch_remote_md5(url, logger=logger)
        if remote_md5 is None and logger:
            logger.info(
                "No remote MD5 available for tar file %s; checksum verification will be skipped after download.",
                fname,
            )
    else:
        if logger:
            logger.info(
                "Non-tar file %s detected; checksum verification will not be attempted.",
                fname,
            )

    # Download fresh copy
    try:
        elapsed = download_file_with_retries(url, tmp_path, logger=logger)
    except Exception as e:
        if logger:
            logger.error("Failed to download %s after retries: %s", fname, e)
        # Wipe temp file on failure
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except Exception as e2:
            if logger:
                logger.error(
                    "Failed to remove temporary file %s after download failure: %s",
                    tmp_path,
                    e2,
                )
        return {
            "filename": fname,
            "url": url,
            "size_bytes": None,
            "elapsed_seconds": None,
            "status": "failed_download",
            "error": str(e),
        }

    # Promote .part to final
    try:
        os.replace(tmp_path, dest_path)
    except Exception as e:
        if logger:
            logger.error(
                "Failed to promote partial file %s to %s: %s", tmp_path, dest_path, e
            )
        # Wipe temp file on failure
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except Exception as e2:
            if logger:
                logger.error(
                    "Failed to remove temporary file %s after promote failure: %s",
                    tmp_path,
                    e2,
                )
        return {
            "filename": fname,
            "url": url,
            "size_bytes": None,
            "elapsed_seconds": elapsed,
            "status": "failed_promote",
            "error": str(e),
        }

    size_bytes = os.path.getsize(dest_path)

    # If this is not a tar file, or no checksum is available, we cannot verify; just record as unverified download.
    if not fname.endswith(".tar") or remote_md5 is None:
        if logger:
            logger.info(
                "Downloaded %s (size=%d) without checksum verification (non-tar or no MD5).",
                fname,
                size_bytes,
            )
        return {
            "filename": fname,
            "url": url,
            "size_bytes": size_bytes,
            "elapsed_seconds": elapsed,
            "status": "downloaded_no_checksum",
        }

    # Verify checksum after download for tar files with MD5.
    try:
        local_md5 = compute_md5(dest_path)
        if local_md5 != remote_md5:
            if logger:
                logger.error(
                    "Checksum mismatch after download for %s (local=%s, remote=%s); deleting file",
                    fname,
                    local_md5,
                    remote_md5,
                )
            try:
                os.remove(dest_path)
            except Exception as e2:
                if logger:
                    logger.error("Failed to delete corrupted file %s: %s", dest_path, e2)
            return {
                "filename": fname,
                "url": url,
                "size_bytes": size_bytes,
                "elapsed_seconds": elapsed,
                "status": "failed_checksum",
            }

        if logger:
            logger.info(
                "Downloaded and verified %s successfully (size=%d, MD5=%s)",
                fname,
                size_bytes,
                local_md5,
            )
        return {
            "filename": fname,
            "url": url,
            "size_bytes": size_bytes,
            "elapsed_seconds": elapsed,
            "status": "downloaded_checksum_ok",
        }
    except Exception as e:
        if logger:
            logger.error("Error verifying checksum for %s: %s", fname, e)
        try:
            os.remove(dest_path)
        except Exception:
            pass
        return {
            "filename": fname,
            "url": url,
            "size_bytes": size_bytes,
            "elapsed_seconds": elapsed,
            "status": "failed_checksum_verification",
            "error": str(e),
        }


def download_mvp_files(
    base_url: str,
    dest_dir: str,
    filename_prefix: str = DEFAULT_PREFIX,
    logger: Optional[logging.Logger] = None,
) -> Dict:
    """
    Discover and download all MVP files under base_url whose filename
    starts with filename_prefix into dest_dir, with one-time checksum
    verification for .tar files when available.

    Downloads are sequential (no concurrency).
    """
    os.makedirs(dest_dir, exist_ok=True)

    links = list_index_files(base_url, logger=logger)
    candidates: List[Dict] = []
    for url in links:
        fname = os.path.basename(url.split("?", 1)[0])
        if fname.startswith(filename_prefix):
            # Skip the .md5 sidecar files themselves; we handle those separately
            if not fname.endswith(".md5"):
                candidates.append({"url": url, "filename": fname})

    if logger:
        logger.info(
            "Found %d candidate files with prefix %r",
            len(candidates),
            filename_prefix,
        )

    if not candidates:
        manifest = {
            "base_url": base_url,
            "dest_dir": os.path.abspath(dest_dir),
            "filename_prefix": filename_prefix,
            "num_candidates": 0,
            "downloads": [],
            "total_elapsed_seconds": 0.0,
            "timestamp": time.time(),
        }
        # Still write manifest for traceability
        manifest_path = os.path.join(dest_dir, "mvp_download_manifest.json")
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        if logger:
            logger.info("No candidates found. Manifest written to %s", manifest_path)
        return manifest

    start_all = time.time()
    results: List[Dict] = []

    for entry in candidates:
        res = process_one_file(entry, dest_dir=dest_dir, logger=logger)
        results.append(res)

    end_all = time.time()

    manifest = {
        "base_url": base_url,
        "dest_dir": os.path.abspath(dest_dir),
        "filename_prefix": filename_prefix,
        "num_candidates": len(candidates),
        "downloads": results,
        "total_elapsed_seconds": end_all - start_all,
        "timestamp": time.time(),
    }

    # Write manifest next to downloads
    manifest_path = os.path.join(dest_dir, "mvp_download_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    if logger:
        logger.info("Manifest written to %s", manifest_path)
    return manifest


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Download MVP sumstats from NCBI GIA index with one-time checksum verification for .tar files."
    )
    parser.add_argument(
        "--base-url",
        required=True,
        help="Base URL of the index, e.g. https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/GIA/",
    )
    parser.add_argument(
        "--dest",
        required=True,
        help="Destination directory for downloads.",
    )
    parser.add_argument(
        "--prefix",
        default=DEFAULT_PREFIX,
        help=f"Filename prefix to match (default: {DEFAULT_PREFIX})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="(Deprecated) Downloads are now sequential; this option is ignored.",
    )
    parser.add_argument(
        "--no-console-log",
        action="store_true",
        help="Disable logging to console; only log to file.",
    )

    args = parser.parse_args(argv)

    logger = setup_logging(args.dest, verbose_console=not args.no_console_log)
    logger.info("Starting MVP download run.")
    logger.info("Base URL: %s", args.base_url)
    logger.info("Destination: %s", args.dest)
    logger.info("Filename prefix: %s", args.prefix)

    manifest = download_mvp_files(
        base_url=args.base_url,
        dest_dir=args.dest,
        filename_prefix=args.prefix,
        logger=logger,
    )

    # Minimal summary to stdout (and log)
    num_downloaded_ok = sum(
        1 for r in manifest["downloads"]
        if r["status"] == "downloaded_checksum_ok"
    )
    num_failed = sum(
        1 for r in manifest["downloads"]
        if r["status"].startswith("failed")
    )
    num_skipped = sum(
        1 for r in manifest["downloads"]
        if r["status"].startswith("skipped")
    )

    summary = {
        "dest_dir": manifest["dest_dir"],
        "num_candidates": manifest["num_candidates"],
        "num_downloaded_ok": num_downloaded_ok,
        "num_skipped": num_skipped,
        "num_failed": num_failed,
        "total_elapsed_seconds": manifest["total_elapsed_seconds"],
    }

    logger.info("Run summary: %s", json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
