from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Mapping, Sequence


def setup_logger(log_file: Path | None = None, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger("hawcquant")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    stream = logging.StreamHandler()
    stream.setFormatter(fmt)
    stream.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(stream)
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file)
        fh.setFormatter(fmt)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)
    return logger


def require_tools(tools: Iterable[str]) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise RuntimeError("Missing required command(s): " + ", ".join(missing))


def run_cmd(
    cmd: Sequence[str] | str,
    logger: logging.Logger | None = None,
    shell: bool = False,
    cwd: Path | None = None,
    log_path: Path | None = None,
    timeout: int | None = None,
    stream: bool = True,
) -> None:
    """Run an external command.

    Parameters
    ----------
    stream:
        If true, stdout/stderr are streamed line-by-line to the logger and optional log file.
        If false, stdout/stderr are written directly to the optional log file, which is safer
        for tools that print progress using carriage returns.
    """
    cmd_display = cmd if isinstance(cmd, str) else " ".join(map(str, cmd))
    if logger:
        logger.info("$ %s", cmd_display)
        if log_path:
            logger.info("Command log: %s", log_path)

    if not stream:
        stdout = stderr = None
        fh = None
        try:
            if log_path:
                log_path.parent.mkdir(parents=True, exist_ok=True)
                fh = open(log_path, "a")
                fh.write(f"\n$ {cmd_display}\n")
                fh.flush()
                stdout = fh
                stderr = subprocess.STDOUT
            proc = subprocess.run(
                cmd,
                shell=shell,
                cwd=str(cwd) if cwd else None,
                stdout=stdout,
                stderr=stderr,
                timeout=timeout,
            )
            if proc.returncode != 0:
                raise RuntimeError(f"Command failed with exit code {proc.returncode}: {cmd_display}")
        finally:
            if fh:
                fh.close()
        return

    log_handle = None
    try:
        if log_path:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_handle = open(log_path, "a")
            log_handle.write(f"\n$ {cmd_display}\n")
            log_handle.flush()

        proc = subprocess.Popen(
            cmd,
            shell=shell,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        try:
            for line in proc.stdout or []:
                line = line.rstrip("\n")
                if log_handle:
                    log_handle.write(line + "\n")
                    log_handle.flush()
                if logger and line.strip():
                    logger.info(line)
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            raise RuntimeError(f"Command timed out after {timeout} seconds: {cmd_display}")

        if rc != 0:
            raise RuntimeError(f"Command failed with exit code {rc}: {cmd_display}")
    finally:
        if log_handle:
            log_handle.close()


def pathify(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


def mkdir(path: str | Path) -> Path:
    p = pathify(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_list(path: str | Path) -> list[str]:
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip() and not line.startswith("#")]


def write_tsv(rows: list[Mapping[str, object]], out: Path, columns: list[str] | None = None) -> None:
    import pandas as pd
    df = pd.DataFrame(rows)
    if columns:
        df = df[columns]
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
