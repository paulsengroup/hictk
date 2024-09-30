#!/usr/bin/env python3

# Copyright (c) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import os
import pathlib
import re
import subprocess as sp


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "--root-dir",
        type=pathlib.Path,
        required=False,
        help="Path to the root of the doc folder.",
    )
    cli.add_argument(
        "--inplace",
        action="store_true",
        default=False,
        help="Do the replacement in-place.",
    )

    return cli


def infer_root_dir(cwd: pathlib.Path | None = None) -> pathlib.Path:
    if cwd is None:
        cwd = pathlib.Path(__file__).parent.resolve()

    res = sp.check_output(["git", "rev-parse", "--show-toplevel"], cwd=cwd).decode("utf-8").split("\n")[0]

    root_dir = pathlib.Path(res)
    if root_dir.is_dir():
        return root_dir

    if cwd == pathlib.Path(__file__).parent.resolve():
        raise RuntimeError("Unable to infer the root of hictk's repository.")

    return infer_root_dir(pathlib.Path(__file__).parent.resolve())


def patch_index_file(path: pathlib.Path, inplace: bool):
    url = os.getenv("READTHEDOCS_CANONICAL_URL")
    if url is None:
        raise RuntimeError("Unable to read url from the READTHEDOCS_CANONICAL_URL env variable")

    logging.info(f'READTHEDOCS_CANONICAL_URL="{url}"')

    toks = url.removeprefix("https://").rstrip("/").split("/")
    if len(toks) < 2:
        raise RuntimeError("Failed to parse READTHEDOCS_CANONICAL_URL variable")

    tgt_domain = toks[0]
    tgt_branch = toks[-1]

    logging.info(f'new_domain="{tgt_domain}"')
    logging.info(f'new_branch="{tgt_branch}"')

    pdf_pattern = re.compile(r"https://hictk\.readthedocs\.io/_/downloads/en/[\w\-_]+/pdf/")
    html_pattern = re.compile(r"https://hictk\.readthedocs\.io/en/[\w\-_]+/")

    payload = path.read_text()
    payload = pdf_pattern.sub(f"https://{tgt_domain}/_/downloads/en/{tgt_branch}/pdf/", payload)
    payload = html_pattern.sub(f"https://{tgt_domain}/en/{tgt_branch}/", payload)

    if inplace:
        logging.info(f'Updating file "{path}" inplace...')
        path.write_text(payload)
    else:
        print(payload, end="")


def main():
    if "READTHEDOCS" not in os.environ:
        logging.info("Script is not being run by ReadTheDocs. Returning immediately!")
        return

    args = vars(make_cli().parse_args())

    root_dir = args["root_dir"]
    if root_dir is None:
        root_dir = infer_root_dir()

    index_file = root_dir / "docs" / "index.rst"
    if not index_file.exists():
        raise RuntimeError(f'Unable to find file "docs/index.rst" under {root_dir}')

    patch_index_file(index_file, args["inplace"])


if __name__ == "__main__":
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(logging.INFO)
    main()
