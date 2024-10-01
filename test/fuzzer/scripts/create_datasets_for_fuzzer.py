#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import contextlib
import gzip
import logging
import multiprocessing as mp
import pathlib
import re
import shlex
import shutil
import subprocess as sp
import sys
import tempfile
import time
import warnings
from concurrent.futures import ProcessPoolExecutor as Pool
from typing import Dict, List, Tuple, Union

import cooler
import hic2cool
import pandas as pd


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
        "pairs",
        type=existing_file,
        help="Path to a .pairs file to use as input.",
    )
    cli.add_argument(
        "output-prefix",
        type=pathlib.Path,
        help="Path prefix (including parent folder(s)) to use for output.",
    )

    cli.add_argument(
        "--chrom-sizes",
        type=existing_file,
        help="Path to a .chrom.sizes file with the list of chromosomes to use as reference.\n"
        "When not provided, chromosome sizes are extracted from the .pairs file header.",
    )

    cli.add_argument(
        "--resolutions",
        type=positive_int,
        nargs="+",
        help="List of resolutions to be generated.",
    )
    cli.add_argument(
        "--juicer-tools-jar",
        type=existing_file,
        help="Path to the juicer_tools jar file.",
        required=True,
    )
    cli.add_argument(
        "--hic-tools-jar",
        type=existing_file,
        help="Path to the hic_tools jar file.",
        required=True,
    )
    cli.add_argument(
        "--normalization-methods",
        type=str,
        nargs="+",
        help="List of normalization methods to compute.",
        choices={
            "VC",
            "INTER_VC",
            "GW_VC",
            "VC_SQRT",
            "KR",
            "INTER_KR",
            "GW_KR",
            "SCALE",
            "INTER_SCALE",
            "GW_SCALE",
        },
        default=[
            "VC",
            "INTER_VC",
            "GW_VC",
            "VC_SQRT",
            "KR",
            "INTER_KR",
            "GW_KR",
            "SCALE",
            "INTER_SCALE",
            "GW_SCALE",
        ],
    )
    cli.add_argument(
        "--tmpdir",
        type=pathlib.Path,
        default=pathlib.Path(tempfile.gettempdir()),
        help="Path to a folder to use for temporary file(s).",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )
    cli.add_argument(
        "--verbosity",
        type=str,
        choices={"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"},
        default="INFO",
        help="Set log verbosity.",
    )
    cli.add_argument(
        "-Xms",
        type=str,
        default="1024m",
    )
    cli.add_argument(
        "-Xmx",
        type=str,
        default="32g",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    return cli


def setup_logger(level):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def get_compressor(compression_level: int, nproc: int) -> Union[list, None]:
    if cmd := shutil.which("pigz") and nproc > 1:
        return [cmd, f"-{compression_level}", "-p", str(min(16, nproc))]
    if cmd := shutil.which("gzip"):
        return [cmd, f"-{compression_level}"]
    return None


def fetch_pixels(cf: cooler.Cooler, chrom1: str, chrom2: str) -> pd.DataFrame:
    logging.info("[%d] Processing pixels for %s:%s...", cf.binsize, chrom1, chrom2)
    pixels = pd.DataFrame(cf.matrix(balance=False, as_pixels=True, join=True).fetch(chrom1, chrom2))

    pixels["end1"] = 0
    pixels["end2"] = 1
    pixels.rename(columns={"end1": "str1", "end2": "str2"}, inplace=True)

    pixels["frag1"] = pixels["str1"]
    pixels["frag2"] = pixels["str2"]

    return pixels[
        [
            "str1",
            "chrom1",
            "start1",
            "frag1",
            "str2",
            "chrom2",
            "start2",
            "frag2",
            "count",
        ]
    ]


def dump_pixels_helper_(uri: str, dest):
    cf = cooler.Cooler(uri)
    logging.info("[%d] writing pixels to temporary file %s...", cf.binsize, dest)

    num_pixels = 0
    queries = []
    chroms = cf.chromnames
    for i, chrom1 in enumerate(chroms):
        for chrom2 in chroms[i:]:
            queries.append([chrom1, chrom2])

    for chrom1, chrom2 in queries:
        pixels = fetch_pixels(cf, chrom1, chrom2)
        pixels.to_csv(dest, sep="\t", index=False, header=False)
        num_pixels += len(pixels)

    if isinstance(dest, pathlib.Path) or isinstance(dest, str):
        logging.info(
            '[%d] written a total of %s pixels to file "%s".',
            cf.binsize,
            num_pixels,
            dest,
        )
    else:
        logging.info("[%d] written a total of %s pixels.", cf.binsize, num_pixels)


def get_mcool_resolutions(uri: str | pathlib.Path) -> Union[None, List[int]]:
    if not cooler.fileops.is_multires_file(str(uri)):
        return None

    resolutions = []
    for suffix in cooler.fileops.list_coolers(uri):
        resolutions.append(cooler.Cooler(f"{uri}::/{suffix}").binsize)

    return list(sorted(resolutions))


def get_base_resolution_uri(uri: str) -> str:
    if not cooler.fileops.is_multires_file(uri):
        return uri

    base_resolution = get_mcool_resolutions(uri)[0]
    return f"{uri}::/resolutions/{base_resolution}"


def dump_pixels(uri: str, out_prefix: pathlib.Path, compression_level: int, nproc: int) -> pathlib.Path:
    uri = get_base_resolution_uri(uri)

    cmd = get_compressor(compression_level, max(1, nproc - 1))
    if cmd is None:
        dest = pathlib.Path(f"{out_prefix}.pixels.txt")
        with open(dest, "w") as f:
            dump_pixels_helper_(uri, f)
            return dest

    dest = pathlib.Path(f"{out_prefix}.pixels.txt.gz")
    with open(dest, "wb") as f:  # For some reason using shell=True is faster
        with sp.Popen(cmd, stdin=sp.PIPE, stderr=sp.PIPE, stdout=f, shell=True) as compressor:
            dump_pixels_helper_(uri, compressor.stdin)
            compressor.communicate()
            if (code := compressor.returncode) != 0:
                print(compressor.stderr, file=sys.stderr)
                raise RuntimeError(f"{cmd} terminated with code {code}")

    return dest


def dump_chrom_sizes_from_cooler(uri: str, tmpdir: pathlib.Path) -> pathlib.Path:
    uri = get_base_resolution_uri(uri)
    c = cooler.Cooler(uri)
    chrom_sizes = c.chromsizes
    assembly = c.info.get("genome-assembly")
    if assembly is None:
        assembly = c.info.get("assembly")

    if assembly is None:
        dest = tmpdir / "chrom.sizes"
    else:
        dest = tmpdir / f"{assembly}.chrom.sizes"
    chrom_sizes.to_csv(dest, sep="\t", header=False, index=True)
    return dest


def run_juicer_tools_pre(
    path_to_pixels: pathlib.Path,
    path_to_chrom_sizes: pathlib.Path,
    path_to_hic: pathlib.Path,
    path_to_jar: pathlib.Path,
    tmpdir: pathlib.Path,
    nproc: int,
    xms: str,
    xmx: str,
    resolutions: Union[List[int], None],
):
    cmd = [
        shutil.which("java"),
        f"-Xms{xms}",
        f"-Xmx{xmx}",
        "-jar",
        str(path_to_jar),
        "pre",
        "-j",
        str(nproc),
        "--threads",
        str(nproc),
        "-t",
        str(tmpdir),
        "-n",
    ]

    if resolutions is not None:
        cmd.extend(["-r", ",".join((str(res) for res in resolutions))])

    cmd.extend([path_to_pixels, path_to_hic, path_to_chrom_sizes])

    cmd = shlex.split(" ".join(str(tok) for tok in cmd))
    logging.info("running %s...", cmd)

    res = sp.run(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
    for line in res.stdout.decode("utf-8").split("\n"):
        logging.debug("stdout: %s", line)
    for line in res.stderr.decode("utf-8").split("\n"):
        logging.debug("stderr: %s", line)

    if res.returncode != 0:
        raise RuntimeError(f"command {cmd} returned with exit code {res.returncode}")


def file_is_gzipped(path: pathlib.Path) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except:  # noqa
        return False


def try_open_gzipped_file(path: pathlib.Path):
    if file_is_gzipped(path):
        return gzip.open(path, "rt")

    return open(path, "r")


def read_chrom_sizes_from_header(path_to_pairs: pathlib.Path) -> Dict[str, int]:
    chroms = {}
    pattern = re.compile(r"#chromsize:\s+([\w_\-]+)\s+(\d+)$")
    with try_open_gzipped_file(path_to_pairs) as f:
        for line in f:
            if not line.startswith("#"):
                break
            matches = pattern.search(line)
            if matches:
                chroms[matches.group(1)] = int(matches.group(2))

    if len(chroms) == 0:
        raise RuntimeError(f"Unable to read any chromosome size from {path_to_pairs}")
    return chroms


def write_chrom_sizes_to_file(chrom_sizes: Dict[str, int], dest: pathlib.Path):
    with open(dest, "w") as f:
        for chrom, size in chrom_sizes.items():
            f.write(f"{chrom}\t{size}\n")


def run_juicer_tools_add_norm(
    path_to_hic: pathlib.Path,
    path_to_jar: pathlib.Path,
    normalization_methods: List[str],
    resolution: int,
    nproc: int,
    xms: str,
    xmx: str,
):
    cmd = [
        shutil.which("java"),
        f"-Xms{xms}",
        f"-Xmx{xmx}",
        "-jar",
        str(path_to_jar),
        "addNorm",
        "-j",
        str(nproc),
        "-F",
        "-w",
        resolution,
        "-r",
        ",".join(str(resolution) for _ in enumerate(normalization_methods)),
        "-k",
        ",".join(normalization_methods),
        str(path_to_hic),
    ]

    cmd = shlex.split(" ".join(str(tok) for tok in cmd))
    logging.info("running %s...", cmd)

    res = sp.run(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
    for line in res.stdout.decode("utf-8").split("\n"):
        logging.debug("stdout: %s", line)
    for line in res.stderr.decode("utf-8").split("\n"):
        logging.debug("stderr: %s", line)

    if res.returncode != 0:
        raise RuntimeError(f"command {cmd} returned with exit code {res.returncode}")


def check_juicer_tools_is_available(path_to_jar: pathlib.Path):
    java = shutil.which("java")
    if java is None:
        raise RuntimeError("Unable to find java in your PATH!")

    sp.check_output([java, "-jar", str(path_to_jar), "pre", "--help"], stderr=sp.STDOUT)
    sp.check_output([java, "-jar", str(path_to_jar), "addNorm", "--help"], stderr=sp.STDOUT)


def ingest_pairs(
    pairs: pathlib.Path,
    chrom_sizes: pathlib.Path,
    resolution: int,
    tmpdir: pathlib.Path,
) -> pathlib.Path:
    outpath = tmpdir / f"tmpcooler.{resolution}.cool"
    cmd = [
        "cooler",
        "cload",
        "pairs",
        f"{chrom_sizes}:{resolution}",
        pairs,
        outpath,
        "-c1",
        "2",
        "-p1",
        "3",
        "-c2",
        "4",
        "-p2",
        "5",
        "--temp-dir",
        tmpdir,
    ]
    sp.check_call(cmd),
    return outpath


def cooler_zoomify(
    base_clr: pathlib.Path,
    resolutions: List[int],
    tmpdir: pathlib.Path,
) -> pathlib.Path:
    outpath = tmpdir / "tmpcooler.mcool"
    cooler.zoomify_cooler(str(base_clr), str(outpath), resolutions, 50_000_000)
    return outpath


def mcool_to_cool(mclr: pathlib.Path, outprefix: pathlib.Path, force: bool, pool: Pool) -> Dict[int, pathlib.Path]:
    coolers = {}
    results = []
    for res in get_mcool_resolutions(mclr):
        dest = pathlib.Path(f"{outprefix}.{res}.cool")
        if force:
            dest.unlink(missing_ok=True)

        if dest.exists():
            raise RuntimeError(
                f'refusing to overwrite existing file "{dest}". Pass --force to overwrite existing file(s).'
            )

        results.append(pool.submit(cooler.fileops.cp, f"{mclr}::/resolutions/{res}", str(dest), False))
        coolers[res] = dest

    for r in results:
        r.result()

    return coolers


def cool_to_hic(
    path_to_clr: pathlib.Path,
    path_to_chrom_sizes: pathlib.Path,
    out_prefix: pathlib.Path,
    out_suffix: str,
    jar: pathlib.Path,
    tmpdir: pathlib.Path,
    xms: str,
    xmx: str,
    force: bool,
) -> Tuple[int, pathlib.Path]:
    logging.info(f"converting Cooler at {path_to_clr} to .hic format...")
    t0 = time.time()

    with tempfile.TemporaryDirectory(prefix=str(tmpdir)) as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        resolution = cooler.Cooler(str(path_to_clr)).binsize
        pixels = dump_pixels(str(path_to_clr), tmpdir / str(resolution), 1, 1)

        tmpdest = tmpdir / f"{out_prefix.stem}.hic"
        dest = pathlib.Path(f"{out_prefix}.{resolution}.hic{out_suffix}")

        if force:
            dest.unlink(missing_ok=True)

        if dest.exists():
            raise RuntimeError(
                f'refusing to overwrite existing file "{dest}". Pass --force to overwrite existing file(s).'
            )

        run_juicer_tools_pre(pixels, path_to_chrom_sizes, tmpdest, jar, tmpdir, 1, xms, xmx, [resolution])

        tmpdest.rename(dest)

    t1 = time.time()

    logging.info(f"DONE! Converting {path_to_clr} to .hic format took {t1-t0}s!")

    return resolution, dest


def copy_hic_norms(
    path_to_hic: pathlib.Path,
    path_to_cooler: pathlib.Path,
    normalization_methods: List[str],
):
    assert len(normalization_methods) != 0
    logging.info(f"copying weights from {path_to_hic} to {path_to_cooler}...")

    t0 = time.time()
    with contextlib.redirect_stdout(None), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hic2cool.hic2cool_extractnorms(str(path_to_hic), str(path_to_cooler), silent=True)

    columns = cooler.Cooler(str(path_to_cooler)).bins().columns.tolist()
    for norm in normalization_methods:
        if norm not in columns:
            raise RuntimeError(f'failed to copy "{norm}" normalization to Cooler file "{path_to_cooler}"')
    t1 = time.time()

    logging.info(f"DONE! Copying weights from {path_to_hic} to {path_to_cooler}. Copying took {t1 - t0}s.")


def main():
    args = vars(make_cli().parse_args())
    setup_logger(args["verbosity"])

    out_prefix = args["output-prefix"]
    pairs = args["pairs"]
    resolutions = args["resolutions"]
    check_juicer_tools_is_available(args["juicer_tools_jar"])
    check_juicer_tools_is_available(args["hic_tools_jar"])
    jt_jar = args["juicer_tools_jar"]
    ht_jar = args["hic_tools_jar"]
    nproc = args["nproc"]
    force = args["force"]

    with (
        tempfile.TemporaryDirectory(prefix=str(args["tmpdir"] / "hictk-")) as tmpdir,
        Pool(max_workers=nproc, initializer=setup_logger, initargs=(args["verbosity"],)) as pool,
    ):

        if args["chrom_sizes"] is None:
            chrom_sizes = pathlib.Path(tmpdir) / "chrom.sizes"
            write_chrom_sizes_to_file(read_chrom_sizes_from_header(pairs), chrom_sizes)
        else:
            chrom_sizes = args["chrom_sizes"]

        # pairs to cool
        base_clr = ingest_pairs(pairs, chrom_sizes, min(resolutions), pathlib.Path(tmpdir))

        # cool to mcool
        mclr = cooler_zoomify(base_clr, resolutions, pathlib.Path(tmpdir))
        base_clr.unlink()

        out_prefix.parent.mkdir(parents=True, exist_ok=True)

        # mcool to cool(s)
        coolers = mcool_to_cool(mclr, out_prefix, force, pool)
        mclr.unlink()

        # cool(s) to hic (v8 and v9)
        hic_v8_files = {}
        hic_v9_files = {}
        for res, clr in coolers.items():
            hic_v8_files[res] = pool.submit(
                cool_to_hic,
                clr,
                chrom_sizes,
                out_prefix,
                "8",
                jt_jar,
                pathlib.Path(tmpdir),
                args["Xms"],
                args["Xmx"],
                force,
            )

            hic_v9_files[res] = pool.submit(
                cool_to_hic,
                clr,
                chrom_sizes,
                out_prefix,
                "9",
                ht_jar,
                pathlib.Path(tmpdir),
                args["Xms"],
                args["Xmx"],
                force,
            )
        for result in hic_v8_files.values():
            resolution, hf = result.result()
            hic_v8_files[resolution] = hf

        for result in hic_v9_files.values():
            resolution, hf = result.result()
            hic_v9_files[resolution] = hf

        # Balance hic(s)
        norm_methods = [m for m in args["normalization_methods"] if m != "NONE"]
        if len(norm_methods) == 0:
            logging.info("skipping normalization of .hic files")
        else:
            for res, hf in hic_v8_files.items():
                run_juicer_tools_add_norm(
                    hf,
                    jt_jar,
                    norm_methods,
                    res,
                    nproc,
                    args["Xms"],
                    args["Xmx"],
                )

            norm_methods = [m for m in norm_methods if not m.endswith("KR")]
            for res, hf in hic_v9_files.items():
                run_juicer_tools_add_norm(
                    hf,
                    ht_jar,
                    norm_methods,
                    res,
                    nproc,
                    args["Xms"],
                    args["Xmx"],
                )

            # Copy hic weights to cool
            results = []
            for res, hf in hic_v8_files.items():
                results.append(pool.submit(copy_hic_norms, hf, coolers[res], norm_methods))

            for r in results:
                _ = r.result()

        # Balance cool(s)
        for clr in coolers.values():
            logging.info('balancing cooler at "%s"...', clr)
            t0 = time.time()
            cooler.balance_cooler(cooler.Cooler(str(clr)), store=True, map=pool.map)
            t1 = time.time()
            logging.info('DONE! balancing cooler at "%s" took %ss', clr, t1 - t0)


if __name__ == "__main__":
    main()
