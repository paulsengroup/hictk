#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import pathlib
import shutil
import subprocess as sp
import sys
from typing import List


def existing_dir(arg):
    if (path := pathlib.Path(arg)).is_dir():
        return path

    raise FileNotFoundError(f'Path "{arg}" is not an existing directory')


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("--root", type=pathlib.Path, help="Path to the repository root.")
    cli.add_argument("--inplace", action="store_true", default=False, help="Modify CMakeLists.txt files in place.")

    return cli


def infer_root_dir() -> pathlib.Path:
    path = pathlib.Path(sp.check_output(["git", "rev-parse", "--show-toplevel"], encoding="utf-8").strip())

    if not path.is_dir():
        raise RuntimeError("Unable to infer repository root!")

    return path


def parse_target_sources_statement(cmakelists: pathlib.Path, keep_files: bool) -> List[str]:
    lines = []
    target_sources_found = False
    with cmakelists.open() as f:
        for line in f:
            line_stripped = line.strip()

            if line_stripped.startswith("target_sources("):
                target_sources_found = True

            if not target_sources_found:
                continue

            if (not keep_files and line_stripped.startswith("FILES")) or line_stripped == ")":
                lines.append(")")
                break
            lines.append(line.rstrip())

    return lines


def infer_indent(lines: List[str]) -> str:
    for line in lines:
        if not line[0].isspace():
            continue

        indent = ""
        for c in line:
            if c.isspace():
                indent += c
                continue
            return indent

    return ""


def reformat_cmakelists(payload: str, gersemi: pathlib.Path) -> str:
    return sp.check_output([gersemi, "-"], stderr=sp.DEVNULL, input=payload.encode("utf-8")).decode("utf-8")


def process_cmakelists(cmakelists: pathlib.Path, gersemi: pathlib.Path, inplace: bool) -> bool:
    target_sources = parse_target_sources_statement(cmakelists, False)
    if len(target_sources) == 0:
        return False

    assert target_sources[-1].strip() == ")"
    target_sources.pop(-1)

    prefix = cmakelists.parent.resolve().as_posix()
    include_dir = (cmakelists.parent / "include").resolve()
    if not include_dir.is_dir():
        return False

    logging.info("processing %s...", cmakelists)

    header_files = []
    impl_header_files = []

    for f in include_dir.rglob("*.hpp"):
        f = f.as_posix().removeprefix(prefix)
        if "/impl/" in f:
            impl_header_files.append(f)
        else:
            header_files.append(f)

    num_headers = len(header_files) + len(impl_header_files)
    if num_headers == 0:
        logging.warn("unable to find any header files under %s", include_dir)
        return False
    else:
        logging.info("found %d header file(s) under %s", num_headers, include_dir)

    if not inplace:
        print(f"#{cmakelists}")

    buff = "\n".join(target_sources)

    indent = infer_indent(target_sources)
    buff += f"\n{2 * indent}FILES\n"

    for header in sorted(header_files):
        header = header.removeprefix(prefix)
        buff += f'{3 * indent}"${{CMAKE_CURRENT_SOURCE_DIR}}{header}"\n'

    for header in sorted(impl_header_files):
        header = header.removeprefix(prefix)
        buff += f'{3 * indent}"${{CMAKE_CURRENT_SOURCE_DIR}}{header}"\n'

    buff += ")\n"

    if not inplace:
        print(buff)

    old_payload = cmakelists.read_text()
    new_payload = old_payload.replace("\n".join(parse_target_sources_statement(cmakelists, True)), buff)
    new_payload = reformat_cmakelists(new_payload, gersemi)

    if old_payload == new_payload:
        return False

    logging.info("file %s is outdated", cmakelists)
    if inplace:
        logging.info("updating file %s...", cmakelists)
        cmakelists.write_text(new_payload)
        logging.info("DONE!")

    return True


def setup_logger(level: str = "INFO"):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def main() -> int:
    args = vars(make_cli().parse_args())

    root_dir = args["root"]
    if root_dir is None:
        root_dir = infer_root_dir()

    logging.info("found repository root: %s", root_dir)

    gersemi = shutil.which("gersemi")
    if gersemi is None:
        gersemi = shutil.which("gersemi", path=root_dir / "venv/bin")

    if gersemi is None:
        raise RuntimeError("Unable to find gersemi in your PATH!")

    gersemi = pathlib.Path(gersemi)
    logging.info("found gersemi executable: %s", gersemi)

    num_outdated_files = 0
    for cmakelists in root_dir.rglob("src/libhictk/**/CMakeLists.txt"):
        if cmakelists == root_dir / "src/libhictk/CMakeLists.txt":
            continue

        if process_cmakelists(cmakelists, gersemi, args["inplace"]):
            num_outdated_files += 1

    if num_outdated_files == 0:
        logging.info("all CMakeLists.txt are already up-to-date!")
        return 0

    if args["inplace"]:
        logging.info("updated %d CMakeLists.txt file(s)", num_outdated_files)
        return 0

    logging.info("found %d CMakeLists.txt that require updating", num_outdated_files)
    return 1


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
