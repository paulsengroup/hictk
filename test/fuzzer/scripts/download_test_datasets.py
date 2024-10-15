#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import hashlib
import json
import logging
import multiprocessing as mp
import pathlib
import random
import shutil
import sys
import tempfile
import time
from typing import Any, Dict, List, Tuple
from urllib.request import HTTPError, urlcleanup, urlretrieve


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def existing_file(arg):
        if (path := pathlib.Path(arg)).exists():
            return path

        raise FileNotFoundError(f'File "{arg}" does not exists')

    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    cli.add_argument(
        "json",
        type=existing_file,
        help="Path to a .json file with the list of test files.",
    )
    cli.add_argument(
        "output-dir",
        type=pathlib.Path,
        help="Path where to download file(s).",
    )
    cli.add_argument(
        "--format",
        nargs="+",
        default=["cool"],
        choices={"cool", "hic8", "hic9"},
        help="Format of the file(s) to be downloaded",
    )
    cli.add_argument(
        "--resolution",
        nargs="+",
        type=int,
        default=[50_000],
        help="Resolution(s) of the file(s) to be downloaded",
    )
    cli.add_argument(
        "--dataset",
        nargs="+",
        type=str,
        default=["4DNFIYECESRC"],
        help="Name of the dataset(s) to be downloaded.",
    )
    cli.add_argument(
        "--download-all",
        action="store_true",
        default=False,
        help="Download all files listed in the given .json file",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing file(s).",
    )
    cli.add_argument(
        "--nproc",
        type=positive_int,
        default=1,
        help="Maximum number of parallel downloads.",
    )

    return cli


def generate_download_list_all(config: Dict[str, Any]) -> Dict[str, Tuple[str, str]]:
    downloads = {}

    for format, values in config.items():  # noqa
        for dataset in values["datasets"].values():
            for c in dataset:
                name = c["name"]
                url = c["url"]
                md5sum = c["md5sum"]
                downloads[name] = (url, md5sum)

    logging.info("collected %d URLs", len(downloads))
    return downloads


def generate_download_list(
    config: Dict[str, Any],
    formats: List[str],
    datasets: List[str],
    resolutions: List[int],
) -> Dict[str, Tuple[str, str]]:
    downloads = {}

    for format in formats:  # noqa
        for dataset in datasets:
            for resolution in resolutions:
                found = False
                for c in config[format]["datasets"].get(dataset, []):
                    if c["resolution"] == resolution:
                        name = c["name"]
                        url = c["url"]
                        md5sum = c["md5sum"]
                        downloads[name] = (url, md5sum)
                        found = True
                        break

                if not found:
                    raise RuntimeError(
                        f'Unable to find file corresponding to format="{format}"; dataset="{dataset}"; resolution={resolution}.'
                    )

    logging.info(
        'collected %d URLs for format(s)="%s"; dataset(s)="%s"; resolution(s)=%s',
        len(downloads),
        formats,
        datasets,
        resolutions,
    )
    return downloads


def digest_file(path: pathlib.Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    logging.info('checksumming file "%s"...', path)
    t0 = time.time()
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            data = f.read(chunk_size)
            if not data:
                break
            md5.update(data)

    digest = md5.hexdigest()
    t1 = time.time()

    logging.info('checksumming file "%s" took %fs', path, t1 - t0)
    logging.info('"%s" MD5: %s', path, digest)

    return digest


def download_file(url: str, md5sum: str, dest: pathlib.Path, max_attempts: int = 3) -> pathlib.Path:
    setup_logger(logging.INFO)
    dest.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{dest.parent}/") as tmpdir:
        tmpdest = pathlib.Path(tmpdir) / dest.name

        for attempt in range(1, max_attempts + 1):
            try:
                logging.info('downloading "%s" (attempt %d/%d)...', url, attempt, max_attempts)
                t0 = time.time()
                urlretrieve(url, tmpdest)
                t1 = time.time()
                logging.info('finished downloading "%s": download took %fs.', url, t1 - t0)
                break
            except HTTPError as e:
                logging.warning(f'failed to download file from URL "{url}" (attempt {attempt}/{max_attempts}): {e}')
                if attempt == max_attempts:
                    msg = f'Failed to download file from URL "{url}": {e}'
                    logging.error(msg)
                    raise RuntimeError(msg)
                sleep_time = random.uniform(5, 10)
                logging.debug(
                    'sleeping for %fs before trying to download "%s" one more time...',
                    sleep_time,
                    url,
                )
                time.sleep(sleep_time)
            finally:
                urlcleanup()

        assert tmpdest.exists()

        digest = digest_file(tmpdest)
        if digest != md5sum:
            msg = f'File "{tmpdest}" appears to be corrupted: MD5 checksum failed: expected {md5sum}, found {digest}.'
            logging.error(msg)
            raise RuntimeError(msg)

        shutil.move(tmpdest, dest)
        return dest


def setup_logger(level):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def main() -> int:
    args = vars(make_cli().parse_args())

    with open(args["json"]) as f:
        if args["download_all"]:
            tasks = generate_download_list_all(json.load(f))
        else:
            tasks = generate_download_list(json.load(f), args["format"], args["dataset"], args["resolution"])

    output_dir = args["output-dir"]
    output_dir.mkdir(exist_ok=True)
    force = args["force"]

    dests = [output_dir / name for name in tasks.keys()]
    urls = [url for url, _ in tasks.values()]
    checksums = [checksum for _, checksum in tasks.values()]

    if not force:
        for path in dests:
            if path.exists():
                raise RuntimeError(
                    'Refusing to overwrite existing file "%s". Pass --force to overwrite.',
                    path,
                )

    with mp.Pool(args["nproc"]) as pool:
        pool.starmap(download_file, zip(urls, checksums, dests), chunksize=1)

    return 0


if __name__ == "__main__":
    setup_logger(logging.INFO)
    sys.exit(main())
