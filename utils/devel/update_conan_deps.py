#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import importlib
import inspect
import json
import pathlib
import re
import shutil
import subprocess as sp
import sys
from typing import Any, List, Tuple

from packaging.version import Version


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def existing_file(arg):
        if (path := pathlib.Path(arg)).exists():
            return path

        raise FileNotFoundError(f'File "{arg}" does not exists')

    cli.add_argument(
        "conanfile",
        type=existing_file,
        help="Path to the conanfile.py to use as input.",
    )

    return cli


@functools.cache
def get_conan() -> pathlib.Path:
    conan = shutil.which("conan")
    if not conan:
        raise RuntimeError("Unable to find Conan in your path")

    return pathlib.Path(conan)


def import_from_path(file_path: pathlib.Path, module_name="conanfile"):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def extract_requirements(conanfile: pathlib.Path) -> List[Tuple[str, Any, Any, str]]:
    conanfile = import_from_path(conanfile, conanfile.stem)
    source = inspect.getsource(conanfile.HictkConan.requirements)

    pattern = re.compile(r"^\"([\w\-_]+)/([\w.\-_]+)#(\w+)\"(.*)$")

    requirements = []
    for line in source.split("\n"):
        line = line.strip()
        if not line.startswith("self.requires("):
            continue

        line = line.removeprefix("self.requires(")

        matches = pattern.search(line).groups()
        if len(matches) == 0:
            raise RuntimeError(f'Failed to parse requirements from line "{line}"')

        name = matches[0]
        version = None
        revision = None
        suffix = None
        if len(matches) > 1:
            version = matches[1]

        if len(matches) > 2:
            revision = matches[2]

        if len(matches) > 3:
            suffix = matches[3]

        requirements.append((name, version, revision, suffix))

    return requirements


def get_last_version(package: str, remotes: str) -> Tuple[str, str]:
    query = f"{package}/*"

    res = sp.run([get_conan(), "list", "--remote", remotes, query, "--format", "json"], stdout=sp.PIPE, stderr=sp.PIPE)
    if res.returncode != 0:
        raise RuntimeError(res.stderr.decode("utf-8"))

    version = None
    for _, versions in json.loads(res.stdout.decode("utf-8")).items():
        for pkg in versions.keys():
            found_ver = Version(pkg.partition("/")[-1])
            if version is None:
                version = found_ver
                continue
            version = max(found_ver, version)

    if version is None:
        raise RuntimeError(f'Unable to find any version for "{package}" using "{remotes}" as remote(s)')

    return package, str(version)


def get_last_revision(package: str, version: str, remotes: str) -> Tuple[str, str, str]:
    query = f"{package}/{version}#*"

    res = sp.run([get_conan(), "list", "--remote", remotes, query, "--format", "json"], stdout=sp.PIPE, stderr=sp.PIPE)
    if res.returncode != 0:
        raise RuntimeError(res.stderr.decode("utf-8"))

    revision = None
    timestamp = None
    for _, versions in json.loads(res.stdout.decode("utf-8")).items():
        assert len(versions) == 1
        _, versions = versions.popitem()
        for rev, metadata in versions["revisions"].items():
            if revision is None:
                revision = rev
                timestamp = metadata["timestamp"]
                continue
            if metadata["timestamp"] > timestamp:
                revision = rev
                timestamp = metadata["timestamp"]

    if revision is None:
        raise RuntimeError(f'Unable to find any revision for "{package}/{version}" using "{remotes}" as remote(s)')

    return package, version, revision


def get_last_recipe(package: str, remotes: str) -> Tuple[str, str, str]:
    _, ver = get_last_version(package, remotes)
    _, _, rev = get_last_revision(package, ver, remotes)

    return package, ver, rev


def main():
    args = vars(make_cli().parse_args())

    for package, _, _, suffix in extract_requirements(args["conanfile"]):
        package, version, revision = get_last_recipe(package, "*")
        print(f'self.requires("{package}/{version}#{revision}"{suffix}')


if __name__ == "__main__":
    main()
