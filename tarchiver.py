#!/usr/bin/env python3

"""Simple, fast, high-compression directory archiver.

Creates .tar.zst archives using threaded zstd.

Exit codes:
  0 = full success
  1 = partial success (archive created and readable, but warnings occurred)
  2 = failure (archive missing, empty, or unreadable)
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
from typing import List, Sequence, Tuple


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Archive a directory to .tar.zst with threaded zstd compression "
            "and verification."
        )
    )
    parser.add_argument("directory", help="Directory to archive")
    parser.add_argument(
        "output_directory",
        nargs="?",
        default=".",
        help="Directory where the archive will be written",
    )
    parser.add_argument(
        "--name",
        default=None,
        help="Override archive base name (without .tar.zst)",
    )
    parser.add_argument(
        "--level",
        type=int,
        default=10,
        help="zstd level 1-19 (default: 10)",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Verify an existing archive path supplied as --name or directory",
    )
    parser.add_argument(
        "--max-runtime-minutes",
        type=float,
        default=15.0,
        help="Stop before running if the estimated runtime exceeds this limit",
    )
    parser.add_argument(
        "--no-runtime-guard",
        action="store_true",
        help="Print the estimate but do not stop when it exceeds the limit",
    )
    return parser.parse_args()


def sanitize_name(raw_name: str) -> str:
    """Convert a directory name into a safe archive base name."""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", raw_name)
    safe = re.sub(r"_+", "_", safe).strip("_")
    return safe or "archive"


def tar_help_text() -> str:
    """Return lowercase tar help text, or empty string on failure."""
    try:
        proc = subprocess.run(
            ["tar", "--help"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return proc.stdout.lower()
    except OSError:
        return ""


def supported_tar_flags() -> List[str]:
    """Detect metadata flags supported by local tar."""
    help_text = tar_help_text()
    flags: List[str] = []

    for flag in [
        "--acls",
        "--xattrs",
        "--selinux",
        "--fflags",
        "--numeric-owner",
    ]:
        if flag in help_text:
            flags.append(flag)

    if "--ignore-failed-read" in help_text:
        flags.append("--ignore-failed-read")

    return flags


def estimate_input_stats(source_dir: str) -> Tuple[int, int]:
    """Return (total_bytes, file_count) for a source tree."""
    total_bytes = 0
    file_count = 0

    for root, _, files in os.walk(source_dir):
        for file_name in files:
            file_path = os.path.join(root, file_name)
            try:
                total_bytes += os.path.getsize(file_path)
                file_count += 1
            except OSError:
                continue

    return total_bytes, file_count


def estimate_runtime_minutes(total_bytes: int, file_count: int, level: int) -> float:
    """Estimate archive runtime from tree size and compression level."""
    throughput_by_level = {
        1: 220.0,
        3: 180.0,
        5: 130.0,
        7: 100.0,
        10: 70.0,
        13: 50.0,
        15: 35.0,
        19: 20.0,
    }

    levels = sorted(throughput_by_level)
    if level in throughput_by_level:
        throughput_mb_s = throughput_by_level[level]
    elif level < levels[0]:
        throughput_mb_s = throughput_by_level[levels[0]]
    elif level > levels[-1]:
        throughput_mb_s = throughput_by_level[levels[-1]]
    else:
        lower = max(value for value in levels if value < level)
        upper = min(value for value in levels if value > level)
        fraction = (level - lower) / (upper - lower)
        throughput_mb_s = throughput_by_level[lower] + fraction * (
            throughput_by_level[upper] - throughput_by_level[lower]
        )

    size_mb = total_bytes / (1024.0 * 1024.0)
    tar_overhead_seconds = file_count * 0.002
    runtime_seconds = (size_mb / max(throughput_mb_s, 1.0)) + tar_overhead_seconds
    return runtime_seconds / 60.0


def build_tar_command(
    source_dir: str,
    archive_path: str,
    tar_flags: Sequence[str],
) -> List[str]:
    """Build tar command that writes archive stream to stdout."""
    src_parent = os.path.dirname(source_dir) or "."
    src_base = os.path.basename(source_dir)

    cmd = ["tar", "-cpf", "-"]

    # Avoid self-inclusion when archive path is inside source directory.
    exclude_rel = None
    try:
        exclude_rel = os.path.relpath(archive_path, source_dir)
    except ValueError:
        exclude_rel = None

    if exclude_rel and not exclude_rel.startswith(".."):
        cmd.append(f"--exclude={exclude_rel}")

    cmd.extend(tar_flags)
    cmd.extend(["-C", src_parent, src_base])
    return cmd


def run_archive(
    source_dir: str,
    archive_path: str,
    tar_flags: Sequence[str],
    level: int,
) -> Tuple[int, str]:
    """Run tar | zstd and return (exit_code, combined_stderr)."""
    tar_cmd = build_tar_command(source_dir, archive_path, tar_flags)
    zstd_cmd = ["zstd", "-T0", f"-{level}", "-q", "-o", archive_path]

    with subprocess.Popen(
        tar_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    ) as tar_proc:
        with subprocess.Popen(
            zstd_cmd,
            stdin=tar_proc.stdout,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=False,
        ) as zstd_proc:
            if tar_proc.stdout is not None:
                tar_proc.stdout.close()

            _, tar_err = tar_proc.communicate()
            _, zstd_err = zstd_proc.communicate()

    stderr_bytes = (tar_err or b"") + b"\n" + (zstd_err or b"")
    stderr_text = stderr_bytes.decode("utf-8", errors="replace")

    if tar_proc.returncode != 0:
        return tar_proc.returncode, stderr_text
    if zstd_proc.returncode != 0:
        return zstd_proc.returncode, stderr_text
    return 0, stderr_text


def extract_warning_lines(stderr_text: str) -> List[str]:
    """Extract warning/error lines from stderr."""
    warnings: List[str] = []
    for line in stderr_text.splitlines():
        low = line.lower()
        if not line.strip():
            continue
        if any(
            key in low
            for key in [
                "permission denied",
                "cannot open",
                "cannot stat",
                "warning",
                "file changed as we read it",
                "cannot read",
                "error",
            ]
        ):
            warnings.append(line)
    return warnings


def verify_archive(archive_path: str) -> Tuple[bool, int, str]:
    """Verify archive exists, is non-empty, passes zstd test, and lists."""
    if not os.path.exists(archive_path):
        return False, 0, "Archive file was not created"

    size = os.path.getsize(archive_path)
    if size == 0:
        return False, 0, "Archive file is empty"

    test_proc = subprocess.run(
        ["zstd", "-t", "-q", archive_path],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if test_proc.returncode != 0:
        return False, 0, (test_proc.stderr or "zstd integrity test failed").strip()

    zstd_proc = subprocess.Popen(
        ["zstd", "-dc", archive_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )
    list_proc = subprocess.run(
        ["tar", "-tf", "-"],
        stdin=zstd_proc.stdout,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if zstd_proc.stdout is not None:
        zstd_proc.stdout.close()
    _, zstd_err = zstd_proc.communicate()

    if zstd_proc.returncode != 0:
        msg = (zstd_err or b"").decode("utf-8", errors="replace").strip()
        return False, 0, msg or "zstd decompression failed"

    if list_proc.returncode != 0:
        msg = list_proc.stderr or "Failed to list archive contents"
        return False, 0, msg.strip()

    entries = [line for line in list_proc.stdout.splitlines() if line.strip()]
    return True, len(entries), ""


def main() -> int:
    """CLI entrypoint."""
    args = parse_args()

    if shutil.which("tar") is None:
        print("Error: tar is not available")
        return 2
    if shutil.which("zstd") is None:
        print("Error: zstd is not available")
        return 2

    if args.level < 1 or args.level > 19:
        print("Error: --level must be between 1 and 19")
        return 2

    if args.verify_only:
        archive_candidate = args.name or args.directory
        ok, entry_count, verify_msg = verify_archive(archive_candidate)
        if ok:
            print(f"Archive verification OK: {archive_candidate}")
            print(f"Entries listed: {entry_count}")
            return 0
        print(f"Archive verification failed: {archive_candidate}")
        print(verify_msg)
        return 2

    source_dir = os.path.abspath(args.directory.rstrip("/"))
    if not os.path.isdir(source_dir):
        print(f"Error: '{args.directory}' is not a valid directory")
        return 2

    output_dir = os.path.abspath(args.output_directory)
    if not os.path.isdir(output_dir):
        print(f"Error: '{args.output_directory}' is not a valid directory")
        return 2

    raw_base_name = args.name if args.name else os.path.basename(source_dir)
    base_name = sanitize_name(raw_base_name)
    archive_path = os.path.join(output_dir, f"{base_name}.tar.zst")

    tar_flags = supported_tar_flags()
    if tar_flags:
        print(f"Tar metadata/warning flags enabled: {' '.join(tar_flags)}")

    total_bytes, file_count = estimate_input_stats(source_dir)
    estimated_minutes = estimate_runtime_minutes(total_bytes, file_count, args.level)
    size_gib = total_bytes / (1024.0 ** 3)
    print(f"Input size: {size_gib:.2f} GiB across {file_count} files")
    print(f"Estimated runtime: {estimated_minutes:.1f} minutes")
    if estimated_minutes > args.max_runtime_minutes:
        message = (
            f"Estimated runtime exceeds {args.max_runtime_minutes:.1f} minutes, "
            "which risks HPC job termination"
        )
        if not args.no_runtime_guard:
            print(f"Error: {message}")
            return 2
        print(f"Warning: {message}")

    print("Compression: zstd threaded")
    print(f"zstd level: {args.level}")
    print(f"Creating archive: {archive_path}")

    rc, stderr_text = run_archive(
        source_dir,
        archive_path,
        tar_flags,
        args.level,
    )

    warnings = extract_warning_lines(stderr_text)
    ok, entry_count, verify_msg = verify_archive(archive_path)

    print(f"Archive path: {archive_path}")
    if ok:
        print(f"Verification passed. Entries listed: {entry_count}")
    else:
        print("Verification failed.")
        print(verify_msg)

    if warnings:
        print(f"Warnings/errors captured: {len(warnings)}")
        for line in warnings[:100]:
            print(f"  {line}")
        if len(warnings) > 100:
            print(f"  ... and {len(warnings) - 100} more")
    else:
        print("No warnings captured from tar/zstd.")

    if not ok or entry_count == 0:
        return 2
    if rc != 0 or warnings:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
