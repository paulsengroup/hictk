#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import glob
import hashlib
import pathlib
import subprocess as sp
import sys
import tempfile


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "binaries",
        type=pathlib.Path,
        nargs="+",
        help="Path to one or more binary files to be processed.",
    )
    cli.add_argument(
        "--output-dir",
        type=pathlib.Path,
        required=True,
        help="Path to a folder where to store the output files.",
    )
    cli.add_argument(
        "--prefix",
        type=str,
        required=True,
        help="Path prefix where to look for raw coverage profile data.",
    )
    cli.add_argument(
        "--suffix",
        type=str,
        default=".profraw",
        help="Path suffix used to look for raw coverage profile data.",
    )
    cli.add_argument(
        "--llvm-cov-bin",
        type=pathlib.Path,
        default="llvm-cov",
        help="Path to llvm-cov.",
    )
    cli.add_argument(
        "--llvm-profdata-bin",
        type=pathlib.Path,
        default="llvm-profdata",
        help="Path to llvm-profdata.",
    )
    cli.add_argument(
        "--format",
        type=str,
        choices={"text", "lcov"},
        default="lcov",
        help="Output format.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing file(s).",
    )

    return cli


def hash_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()

    with path.open("rb") as f:
        while True:
            data = f.read(64 << 10)
            if not data:
                return h.hexdigest()
            h.update(data)


def collect_files(dest: pathlib.Path, prefix: str, suffix: str):
    pattern = rf"{prefix}*{suffix}"
    with dest.open("w") as f:
        for prof in sorted(glob.glob(pattern)):
            f.write(f"{prof}\n")

    if dest.stat().st_size == 0:
        raise RuntimeError(f'Unable to collect any coverage profiles using "{pattern}"')


def merge_profiles(llvm_profdata: pathlib.Path, file_list: pathlib.Path, dest: pathlib.Path):
    with dest.open("w") as stdout:
        sp.check_call([llvm_profdata, "merge", "--input-files", file_list], stdout=stdout)


def export_coverage(
    llvm_cov: pathlib.Path,
    binary: pathlib.Path,
    profdata: pathlib.Path,
    dest: pathlib.Path,
    output_format: str,
    force: bool,
):
    if not force and dest.exists():
        raise RuntimeError(f'Refusing to overwrite existing file "{dest}". Pass --force to overwrite')

    with dest.open("w") as stdout:
        sp.check_call(
            [
                llvm_cov,
                "export",
                binary,
                "-instr-profile",
                profdata,
                "--format",
                output_format,
            ],
            stdout=stdout,
        )


def main() -> int:
    args = vars(make_cli().parse_args())

    output_dir = args["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        file_list = tmpdir / "profile_list.txt"
        collect_files(file_list, args["prefix"], args["suffix"])

        profdata_file = tmpdir / "profile.profdata"
        merge_profiles(args["llvm_profdata_bin"], file_list, profdata_file)

        if args["format"] == "text":
            ext = "json"
        elif args["format"] == "lcov":
            ext = "lcov"
        else:
            raise NotImplementedError

        for binary in set(args["binaries"]):
            digest = hash_file(binary)
            dest = output_dir / f"{binary.stem}.{digest}.{ext}"
            export_coverage(
                args["llvm_cov_bin"],
                binary,
                profdata_file,
                dest,
                args["format"],
                args["force"],
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
