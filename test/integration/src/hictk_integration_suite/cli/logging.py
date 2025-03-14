# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import functools
import platform
import sys
from datetime import datetime, timedelta
from typing import List, Optional

import structlog


class _StructLogPlainStyles(object):
    def __init__(self):
        self.reset = ""
        self.bright = ""
        self.dim = ""

        self.level_critical = ""
        self.level_exception = ""
        self.level_error = ""
        self.level_warn = ""
        self.level_info = ""
        self.level_debug = ""
        self.level_notset = ""

        self.timestamp = ""
        self.logger_name = ""
        self.success = ""
        self.skipped = ""
        self.fail = ""


class _StructLogColorfulStyles(object):
    @staticmethod
    def _try_get_color(key: str):
        """
        Get the requested key (i.e. color) if colorama is available.
        Return "" (i.e. no color) otherwise.
        """
        try:
            import colorama

            return eval(f"colorama.{key}")
        except ImportError:
            return ""

    def __init__(self):
        _try_get_color = _StructLogColorfulStyles._try_get_color
        self.reset = _try_get_color("Style.RESET_ALL")
        self.bright = _try_get_color("Style.BRIGHT")
        self.dim = _try_get_color("Style.DIM")

        self.level_critical = _try_get_color("Fore.RED")
        self.level_exception = _try_get_color("Fore.RED")
        self.level_error = _try_get_color("Fore.RED")
        self.level_warn = _try_get_color("Fore.YELLOW")
        self.level_info = _try_get_color("Fore.GREEN")
        self.level_debug = _try_get_color("Fore.GREEN")
        self.level_notset = _try_get_color("Back.RED")

        self.timestamp = self.dim
        self.logger_name = _try_get_color("Fore.BLUE")
        self.success = self.bright + _try_get_color("Fore.GREEN")
        self.skipped = self.bright + _try_get_color("Fore.YELLOW")
        self.fail = self.bright + _try_get_color("Fore.RED")


@functools.cache
def _map_log_level_to_levelno(level: str) -> int:
    levels = {
        "critical": 50,
        "error": 40,
        "warning": 30,
        "info": 20,
        "debug": 10,
        "notset": 0,
    }

    return levels.get(level.lower(), 0)


def _pretty_format_elapsed_time(
    duration_str: str,
) -> str:
    """
    Format elapsed time between t1 and t0 as a human-readable string.

    Examples:
        123ns
        1.234us
        1.23ms
        1.23s
        1m:2.345s
        1h:2m:3.456s
    """
    microsecond = 1.0e-6
    millisecond = 1.0e-3
    second = 1.0
    minute = 60.0
    hour = 3600.0

    parsed_time = datetime.strptime(duration_str, "%H:%M:%S.%f")

    delta = timedelta(
        hours=parsed_time.hour,
        minutes=parsed_time.minute,
        seconds=parsed_time.second,
        microseconds=parsed_time.microsecond,
    ).total_seconds()

    if delta < microsecond:
        return f"{delta * 1.0e9:.0f}ns"

    if delta < millisecond:
        return f"{delta * 1.0e6:.3f}us"

    if delta < second:
        return f"{delta * 1000:.3f}ms"

    if delta < minute:
        return f"{delta:.3f}s"

    if delta < hour:
        minutes = delta // 60
        seconds = delta - (minutes * 60)
        return f"{minutes:.0f}m:{seconds:.3f}s"

    hours = delta // 3600
    minutes = (delta - (hours * 3600)) // 60
    seconds = delta - (hours * 3600) - (minutes * 60)
    return f"{hours:.0f}h:{minutes:.0f}m:{seconds:.3f}s"


def _format_title(title: str) -> str:
    if not title.startswith("hictk-"):
        return title

    title = title.removeprefix("hictk-")
    if title.endswith("-cli"):
        title = title.removesuffix("-cli")
        title += " (CLI)"

    return title


def _configure_logger_columns(
    colors: bool,
    level_styles: Optional[structlog.dev.Styles] = None,
    event_key: str = "event",
    timestamp_key: str = "timestamp",
    pad_level: bool = True,
) -> List[structlog.dev.Column]:
    """
    The body of this function is an extension of the structlog.dev.ConsoleRenderer.
    In brief, this function configures the columns that will be used to render log messages.

    See the following link for the original implementation:
    https://github.com/hynek/structlog/blob/a60ce7bbb50451ed786ace3c3893fb3a6a01df0a/src/structlog/dev.py#L433
    """

    if level_styles is None:
        level_to_color = structlog.dev.ConsoleRenderer().get_default_level_styles(colors)
    else:
        level_to_color = level_styles

    level_width = 0 if not pad_level else None

    styles: structlog.Styles
    if colors:
        if platform.system() == "Windows":
            # Colorama must be init'd on Windows, but must NOT be
            # init'd on other OSes, because it can break colors.
            import colorama

            colorama.init()

        styles = _StructLogColorfulStyles()
    else:
        styles = _StructLogPlainStyles()

    return [
        structlog.dev.Column(
            timestamp_key,
            structlog.dev.KeyValueColumnFormatter(
                key_style=None,
                value_style=styles.timestamp,
                reset_style=styles.reset,
                value_repr=str,
            ),
        ),
        structlog.dev.Column(
            "level",
            structlog.dev.LogLevelColumnFormatter(
                level_to_color,  # noqa
                reset_style=styles.reset,
                width=level_width,
            ),
        ),
        structlog.dev.Column(
            "title",
            structlog.dev.KeyValueColumnFormatter(
                key_style=None,
                value_repr=_format_title,
                value_style=styles.bright,
                reset_style=styles.reset,
                width=len("rename-chromosomes (CLI)"),
                prefix="[",
                postfix="]",
            ),
        ),
        structlog.dev.Column(
            "id",
            structlog.dev.KeyValueColumnFormatter(
                key_style=None,
                value_style=styles.dim,
                value_repr=lambda x: x[:8],
                reset_style=styles.reset,
                width=8,
                prefix="[id=",
                postfix="]",
            ),
        ),
        structlog.dev.Column(
            "elapsed-time",
            structlog.dev.KeyValueColumnFormatter(
                key_style=None,
                reset_style=styles.reset,
                value_style=styles.dim,
                value_repr=_pretty_format_elapsed_time,
                width=10,
                prefix="[runtime=",
                postfix="]",
            ),
        ),
        structlog.dev.Column(
            "success",
            structlog.dev.KeyValueColumnFormatter(
                key_style=None,
                reset_style=styles.reset,
                value_style=styles.success,
                value_repr=str,
            ),
        ),
        structlog.dev.Column(
            "skipped",
            structlog.dev.KeyValueColumnFormatter(
                key_style=None,
                reset_style=styles.reset,
                value_style=styles.skipped,
                value_repr=str,
            ),
        ),
        structlog.dev.Column(
            "fail",
            structlog.dev.KeyValueColumnFormatter(
                key_style=None,
                reset_style=styles.reset,
                value_style=styles.fail,
                value_repr=str,
            ),
        ),
        structlog.dev.Column(
            "",
            structlog.dev.KeyValueColumnFormatter(
                key_style=styles.bright,
                value_style=styles.bright,
                reset_style=styles.reset,
                value_repr=str,
                postfix=";",
            ),
        ),
        structlog.dev.Column(
            event_key,
            structlog.dev.KeyValueColumnFormatter(
                key_style=None,
                value_style=styles.bright,
                reset_style=styles.reset,
                value_repr=str,
            ),
        ),
    ]


def _warning_handler(message, category, filename, lineno, file=None, line=None):
    from warnings import formatwarning

    structlog.get_logger().warning(
        "\n%s",
        formatwarning(
            message=message,
            category=category,
            filename=filename,
            lineno=lineno,
            line=line,
        ).strip(),
    )


def _install_custom_warning_handler():
    """
    Override the function used to print Python warnings such that warnings are sent to the logger.
    """
    import warnings

    warnings.showwarning = _warning_handler


def setup_logger(log_lvl: str):
    timestamper = structlog.processors.MaybeTimeStamper(fmt="%Y-%m-%d %H:%M:%S.%f")

    def maybe_add_log_level(logger, method_name: str, event_dict):
        if "level" in event_dict:
            return event_dict

        return structlog.processors.add_log_level(logger, method_name, event_dict)

    def preprocess_eventdict(logger, method_name: str, event_dict):
        status = event_dict.pop("status", "")

        if "event" in event_dict and len(event_dict["event"]) == 0:
            event_dict.pop("event")

        keys = ("args", "hictk-runtime", "validation-runtime", "exit-code", "errors", "expect-failure")

        if status == "PASS":
            for k in keys:
                event_dict.pop(k, "")
            event_dict["success"] = status
        elif status == "SKIP":
            for k in keys:
                if k != "args":
                    event_dict.pop(k, "")
            event_dict["skipped"] = status
        elif status == "FAIL":
            event_dict["fail"] = status

        return event_dict

    log_levelno = _map_log_level_to_levelno(log_lvl)

    processors = [
        timestamper,
        maybe_add_log_level,
        preprocess_eventdict,
        structlog.processors.StackInfoRenderer(),
        structlog.dev.ConsoleRenderer(
            columns=_configure_logger_columns(colors=sys.stderr.isatty()),
        ),
    ]

    structlog.configure(
        cache_logger_on_first_use=True,
        wrapper_class=structlog.make_filtering_bound_logger(log_levelno),
        processors=processors,
        logger_factory=structlog.PrintLoggerFactory(file=sys.stderr),
    )

    _install_custom_warning_handler()
