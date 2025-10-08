#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import itertools
import os
import pathlib
import platform
import shlex
import shutil
import subprocess as sp
import tarfile
from typing import Any, Dict


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        help="Path where to store the resulting *.cmake files",
    )
    cli.add_argument(
        "--profile",
        nargs="+",
        type=str,
        default=("gcc", "clang"),
        choices={"gcc", "clang", "default"},
        help="Names of the conan profiles to be used.",
    )
    cli.add_argument(
        "--build-type",
        nargs="+",
        type=str,
        default=("Debug", "RelWithDebInfo", "Release"),
        choices={"Debug", "Release", "RelWithDebInfo", "MinSizeRel"},
        help="Conan build types.",
    )
    cli.add_argument(
        "--cppstd",
        type=str,
        default="17",
        choices={"17", "20", "23"},
        help="C++ standard used to compile the dependencies.",
    )
    cli.add_argument(
        "--build-shared-only",
        action="store_true",
        default=False,
        help="Build dependencies as shared libraries only.",
    )
    cli.add_argument(
        "--build-static-only",
        action="store_true",
        default=False,
        help="Build dependencies as static libraries only.",
    )
    cli.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Print the commands that would be executed if --dry-run was not specified, then exit.",
    )

    cli.add_argument(
        "--devel",
        action="store_true",
        default=False,
        help="Install all development deps. When specified, --with-* options are ignored.",
    )

    cli.add_argument(
        "--ci",
        action="store_true",
        default=False,
        help="Install dependencies as required by the CI pipeline.",
    )

    cli.add_argument(
        "--no-update",
        action="store_true",
        default=False,
        help="Do not pass --update to conan install/create.",
    )

    cli.add_argument(
        "--with-cli-tool-deps",
        action="store_true",
        default=False,
        help="Install dependencies required by the CLI tools.",
    )

    cli.add_argument(
        "--with-benchmark-deps",
        action="store_true",
        default=False,
        help="Install dependencies required by the benchmarks.",
    )

    cli.add_argument(
        "--with-arrow",
        action="store_true",
        default=False,
        help="Install arrow.",
    )

    cli.add_argument(
        "--with-eigen",
        action="store_true",
        default=False,
        help="Install eigen.",
    )

    cli.add_argument(
        "--with-telemetry-deps",
        action="store_true",
        default=False,
        help="Install dependencies required by the telemetry.",
    )

    cli.add_argument(
        "--with-unit-testing-deps",
        action="store_true",
        default=False,
        help="Install dependencies required by the unit test suite.",
    )

    cli.add_argument(
        "--with-fuzzy-testing-deps",
        action="store_true",
        default=False,
        help="Install dependencies required by the fuzzy test harness.",
    )

    return cli


def infer_root_dir() -> pathlib.Path:
    path = pathlib.Path(sp.check_output(["git", "rev-parse", "--show-toplevel"], encoding="utf-8").strip())

    if not path.is_dir():
        raise RuntimeError("Unable to infer repository root!")

    return path


def run_or_print(args, env: Dict[str, str], dry_run: bool):
    if dry_run:
        print(shlex.join(str(x) for x in args))
    else:
        sp.check_call(args, env=env)


def run_conan(
    profile: str,
    build_type: str,
    shared: bool,
    cppstd: str,
    output_prefix: pathlib.Path,
    dry_run: bool,
    no_update: bool,
    recipe_options: Dict[str, bool],
):
    output_folder = output_prefix / profile / build_type / ("shared" if shared else "static")

    args = [
        "conan",
        "install",
        infer_root_dir() / "conanfile.py",
        "--build=missing",
        "--profile",
        profile,
        "--settings",
        f"build_type={build_type}",
        "--settings",
        f"compiler.cppstd={cppstd}",
        "--output-folder",
        output_folder,
        "--options",
        f"*/*:shared={shared}",
    ]

    if not no_update:
        args.append("--update")

    for opt, val in recipe_options.items():
        args.extend(("--options", f"hictk/*:{opt}={val}"))

    env = os.environ.copy()
    env["CC"] = profile
    if profile == "gcc":
        env["CXX"] = "g++"
    elif profile == "clang":
        env["CXX"] = "clang++"
    elif profile != "default":
        raise RuntimeError(
            f'Unrecognized compiler "{profile}". Profiles should be either named "gcc", "clang", or "default".'
        )

    run_or_print(args, env=env, dry_run=dry_run)


def collect_recipe_options(args: Dict[str, Any]) -> Dict[str, bool]:
    recipe_options = {k: bool(v) for k, v in args.items() if k.startswith("with_")}
    if args["devel"] or args["ci"]:
        recipe_options = {k: True for k in recipe_options}

    return recipe_options


def run_local(args: Dict[str, Any]):

    profiles = args["profile"]
    build_types = args["build_type"]
    cppstd = args["cppstd"]
    if args["build_shared_only"]:
        shared_build = [True]
    elif args["build_static_only"]:
        shared_build = [False]
    else:
        shared_build = [True, False]

    dry_run = args["dry_run"]

    output_prefix = infer_root_dir() / "conan-envs"
    if not args["dry_run"]:
        if output_prefix.exists():
            shutil.rmtree(output_prefix)

        output_prefix.mkdir(exist_ok=True)

    recipe_options = collect_recipe_options(args)

    for args in itertools.product(profiles, build_types, shared_build):
        run_conan(
            *args,
            cppstd=cppstd,
            output_prefix=output_prefix,
            dry_run=dry_run,
            recipe_options=recipe_options,
        )


def generate_output_folder(build_type: str) -> str:
    if build_type == "Debug":
        return "cmake-prefix-dbg"

    if build_type == "RelWithDebInfo":
        return "cmake-prefix-rwdi"

    if build_type == "Release":
        return "cmake-prefix-rel"

    raise RuntimeError(f'Unknown build type "{build_type}"')


def run_ci_linux(args: Dict[str, Any]) -> str:
    conan_profile = os.environ.get("CONAN_DEFAULT_PROFILE_PATH")
    if conan_profile is None:
        raise RuntimeError("Environment variable CONAN_DEFAULT_PROFILE_PATH is not defined!")

    build_type = args["build_type"]
    assert len(build_type) == 1
    build_type = build_type[0]

    cppstd = args["cppstd"]
    generate_output_folder(build_type),
    conan_args = [
        "conan",
        "install",
        "conanfile.py",
        "--build",
        "missing",
        "--profile:b",
        conan_profile,
        "--profile:h",
        conan_profile,
        "--settings",
        "compiler.libcxx=libstdc++11",
        "--settings",
        f"compiler.cppstd={cppstd}",
        "--settings",
        f"build_type={build_type}",
        "--output-folder",
        generate_output_folder(build_type),
    ]

    for opt, val in collect_recipe_options(args).items():
        conan_args.extend(("--options", f"hictk/*:{opt}={val}"))

    run_or_print(conan_args, os.environ.copy(), False)

    return generate_output_folder(build_type)


def run_ci_macos(args: Dict[str, Any]) -> str:
    build_type = args["build_type"]
    assert len(build_type) == 1
    build_type = build_type[0]

    cppstd = args["cppstd"]
    conan_args = [
        "conan",
        "install",
        "conanfile.py",
        "--build",
        "missing",
        "--profile:b",
        "default",
        "--profile:h",
        "default",
        "--settings",
        "compiler.libcxx=libc++",
        "--settings",
        f"compiler.cppstd={cppstd}",
        "--settings",
        f"build_type={build_type}",
        "--output-folder",
        generate_output_folder(build_type),
    ]

    for opt, val in collect_recipe_options(args).items():
        conan_args.extend(("--options", f"hictk/*:{opt}={val}"))

    run_or_print(conan_args, os.environ.copy(), False)

    return generate_output_folder(build_type)


def run_ci_windows(args: Dict[str, Any]) -> str:
    build_type = args["build_type"]
    assert len(build_type) == 1
    build_type = build_type[0]

    cppstd = args["cppstd"]
    conan_args = [
        "conan",
        "install",
        "conanfile.py",
        "--build",
        "missing",
        "--build",
        "b2/*",
        "--build",
        "catch2/*",
        "--profile:b",
        "default",
        "--profile:h",
        "default",
        "--settings",
        f"compiler.cppstd={cppstd}",
        "--settings",
        f"compiler.runtime_type={build_type}",
        "--settings",
        f"build_type={build_type}",
        "--output-folder",
        generate_output_folder(build_type),
    ]

    for opt, val in collect_recipe_options(args).items():
        conan_args.extend(("--options", f"hictk/*:{opt}={val}"))

    run_or_print(conan_args, os.environ.copy(), False)

    return generate_output_folder(build_type)


def run_ci(args: Dict[str, Any]):
    if platform.system() == "Linux":
        output_folder = run_ci_linux(args)
    elif platform.system() == "Darwin":
        output_folder = run_ci_macos(args)
    elif platform.system() == "Windows":
        output_folder = run_ci_windows(args)
    else:
        raise RuntimeError(f"Platform {platform.system()} is not supported!")

    build_type = args["build_type"]
    assert len(build_type) == 1
    build_type = build_type[0]

    output_file = generate_output_folder(build_type)
    with tarfile.open(f"{output_file}.tar", "w") as tar:
        tar.add(output_folder)

    shutil.rmtree(output_folder)


def main():
    args = vars(make_cli().parse_args())

    if args["ci"]:
        run_ci(args)
    else:
        run_local(args)


if __name__ == "__main__":
    main()
