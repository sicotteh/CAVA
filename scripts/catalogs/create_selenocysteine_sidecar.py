#!/usr/bin/env python3
"""Run CAVA's create_selenodata.py and install the result as .cesis."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--refseqs", type=Path)
    group.add_argument("--ensembl-refseq", type=Path)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    script = args.repo.resolve() / "cava" / "ensembldb" / "create_selenodata.py"
    if not script.is_file():
        raise SystemExit("Missing {}".format(script))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="cava-seleno-") as temporary:
        work = Path(temporary)
        command = [sys.executable, str(script), "-D", str(work), "-t", args.tag]
        if args.refseqs:
            command.extend(["--refseqs", str(args.refseqs.resolve())])
        else:
            command.extend(["--ensembl_refseq", str(args.ensembl_refseq.resolve())])
        subprocess.run(command, cwd=args.repo.resolve(), check=True)
        expected = work / "SECIS_in_refseq_pos.homo_sapiens.{}.txt".format(args.tag)
        if not expected.is_file() or expected.stat().st_size == 0:
            found = sorted(path.name for path in work.iterdir())
            raise SystemExit(
                "create_selenodata.py did not create {}. Found: {}".format(expected, found)
            )
        shutil.copy2(expected, args.output)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
