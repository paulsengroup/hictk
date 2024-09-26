# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib


class URI:
    def __init__(self, path: pathlib.Path | str, resolve: bool = True):
        if isinstance(path, URI):
            self._path = path.path
            self._grp = path.group
        else:
            path, _, grp = str(path).partition("::")
            self._path = pathlib.Path(path)
            if grp:
                self._grp = str(pathlib.Path(grp).as_posix())
            else:
                self._grp = None

        if resolve:
            self._path = self._path.resolve()

    def __str__(self) -> str:
        if self._grp:
            return f"{self._path}::{self._grp}"
        return str(self._path)

    def __repr__(self) -> str:
        if self._grp is not None:
            return f'URI(path="{self._path}", group="{self._grp}")'

        return f'URI(path="{self._path}", group={self._grp})'

    @property
    def path(self) -> pathlib.Path:
        return self._path

    @property
    def group(self) -> str | None:
        return self._grp

    def __eq__(self, other) -> bool:
        return str(self) == str(other)

    def __hash__(self) -> int:
        return hash((self._path, self._grp))

    def __lt__(self, other) -> bool:
        return str(self) < str(other)
