"""List and install versioned CAVA transcript catalogs stored with Git LFS.

The small ``cava_catalogs.tsv`` manifest is ordinary Git content. Catalog
payloads are immutable Git LFS objects. This keeps catalog publication
independent from CAVA code releases while allowing one catalog to be reused by
many CAVA versions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath
from urllib.error import URLError
from urllib.parse import urlparse
from urllib.request import Request, url2pathname, urlopen

DEFAULT_REPOSITORY = os.environ.get(
    "CAVA_CATALOG_REPOSITORY", "https://github.com/sicotteh/CAVA.git"
)
DEFAULT_MANIFEST_URL = (
    "https://raw.githubusercontent.com/sicotteh/CAVA/master/cava_catalogs.tsv"
)
DEFAULT_REF = os.environ.get("CAVA_CATALOG_REF", "master")
DEFAULT_INSTALL_DIR = os.environ.get("CAVA_DATA_DIR", "./cava_data")
USER_AGENT = "cava_data/2.0.15"

REQUIRED_COLUMNS = {
    "catalog_name",
    "file_name",
    "Build",
    "source",
    "source_version",
    "transcript_ID_type",
    "repository",
    "ref",
    "lfs_path",
    "index_file_name",
    "index_path",
    "transcript_map_file_name",
    "transcript_map_path",
    "selenocysteine_file_name",
    "selenocysteine_path",
    "sha256",
    "index_sha256",
    "transcript_map_sha256",
    "selenocysteine_sha256",
    "lfs_oid_sha256",
    "lfs_size",
    "index_lfs_oid_sha256",
    "index_lfs_size",
    "transcript_map_lfs_oid_sha256",
    "transcript_map_lfs_size",
    "selenocysteine_lfs_oid_sha256",
    "selenocysteine_lfs_size",
}

FILE_GROUPS = (
    ("", "file_name", "lfs_path", "sha256", "lfs_oid_sha256", "lfs_size"),
    (
        "index_",
        "index_file_name",
        "index_path",
        "index_sha256",
        "index_lfs_oid_sha256",
        "index_lfs_size",
    ),
    (
        "transcript_map_",
        "transcript_map_file_name",
        "transcript_map_path",
        "transcript_map_sha256",
        "transcript_map_lfs_oid_sha256",
        "transcript_map_lfs_size",
    ),
    (
        "selenocysteine_",
        "selenocysteine_file_name",
        "selenocysteine_path",
        "selenocysteine_sha256",
        "selenocysteine_lfs_oid_sha256",
        "selenocysteine_lfs_size",
    ),
)


class CavaDataError(RuntimeError):
    """A user-facing catalog-management error."""


class Catalog:
    def __init__(self, row: dict[str, str]):
        self.row = row

    @property
    def name(self) -> str:
        return self.row["catalog_name"].strip()

    @property
    def file_name(self) -> str:
        return _safe_output_name(self.row["file_name"])

    @property
    def repository(self) -> str:
        return self.row.get("repository", "").strip() or DEFAULT_REPOSITORY

    @property
    def ref(self) -> str:
        return self.row.get("ref", "").strip() or DEFAULT_REF

    @property
    def install_files(self) -> list[tuple[str, str, str, str, str]]:
        """Return repository path, installed name, SHA, LFS OID and LFS size."""
        result: list[tuple[str, str, str, str, str]] = []
        for _, name_column, path_column, sha_column, oid_column, size_column in FILE_GROUPS:
            output_name = _safe_output_name(self.row.get(name_column, "").strip())
            repository_path = _safe_repo_path(self.row.get(path_column, "").strip())
            result.append(
                (
                    repository_path,
                    output_name,
                    self.row.get(sha_column, "").strip(),
                    self.row.get(oid_column, "").strip(),
                    self.row.get(size_column, "").strip(),
                )
            )
        return result


# ---------------------------------------------------------------------------
# Manifest handling


def _safe_output_name(value: str) -> str:
    path = PurePosixPath(value)
    if not value or path.is_absolute() or len(path.parts) != 1 or value in (".", ".."):
        raise CavaDataError(
            "Unsafe or empty installed file name in manifest: {!r}".format(value)
        )
    return value


def _safe_repo_path(value: str) -> str:
    path = PurePosixPath(value)
    if not value or path.is_absolute() or ".." in path.parts:
        raise CavaDataError(
            "Unsafe or empty repository path in manifest: {!r}".format(value)
        )
    return path.as_posix()


def _run(
    command: list[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    capture: bool = False,
) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            command,
            cwd=str(cwd) if cwd else None,
            env=env,
            check=True,
            text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE if capture else None,
        )
    except OSError as exc:
        raise CavaDataError(
            "Required executable was not found: {}".format(command[0])
        ) from exc
    except subprocess.CalledProcessError as exc:
        detail = ""
        if capture:
            detail = (exc.stderr or exc.stdout or "").strip()
        suffix = "\n" + detail if detail else ""
        raise CavaDataError(
            "Command failed with exit code {}: {}{}".format(
                exc.returncode, " ".join(command), suffix
            )
        ) from exc


def _require_git_lfs() -> None:
    _run(["git", "--version"], capture=True)
    _run(["git", "lfs", "version"], capture=True)


def _read_text(location: str) -> str:
    parsed = urlparse(location)
    if parsed.scheme in ("http", "https"):
        request = Request(location, headers={"User-Agent": USER_AGENT})
        try:
            with urlopen(request, timeout=60) as response:
                return response.read().decode("utf-8")
        except (URLError, TimeoutError, UnicodeDecodeError) as exc:
            raise CavaDataError(
                "Could not read catalog manifest {}: {}".format(location, exc)
            ) from exc
    if parsed.scheme == "file":
        try:
            return Path(url2pathname(parsed.path)).read_text(encoding="utf-8")
        except OSError as exc:
            raise CavaDataError(
                "Could not read catalog manifest {}: {}".format(location, exc)
            ) from exc
    try:
        return Path(location).expanduser().read_text(encoding="utf-8")
    except OSError as exc:
        raise CavaDataError(
            "Could not read catalog manifest {}: {}".format(location, exc)
        ) from exc


def _bundled_manifest_candidates() -> list[Path]:
    return [
        Path.cwd() / "cava_catalogs.tsv",
        Path(__file__).resolve().parent / "cava_catalogs.tsv",
        Path(sys.prefix) / "share" / "cava" / "cava_catalogs.tsv",
    ]


def _manifest_text(explicit_location: str | None) -> tuple[str, str]:
    environment_location = os.environ.get("CAVA_CATALOG_MANIFEST")
    if explicit_location:
        return _read_text(explicit_location), explicit_location
    if environment_location:
        return _read_text(environment_location), environment_location

    remote_error: CavaDataError | None = None
    try:
        return _read_text(DEFAULT_MANIFEST_URL), DEFAULT_MANIFEST_URL
    except CavaDataError as exc:
        remote_error = exc

    checked: list[str] = []
    for candidate in _bundled_manifest_candidates():
        checked.append(str(candidate))
        if candidate.is_file():
            return candidate.read_text(encoding="utf-8"), str(candidate)
    raise CavaDataError(
        "The current catalog manifest could not be downloaded and no local snapshot "
        "was found. Checked: {}. Original error: {}".format(
            ", ".join(checked), remote_error
        )
    )


def read_manifest(location: str | None = None) -> list[Catalog]:
    text, resolved_location = _manifest_text(location)
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    columns = set(reader.fieldnames or [])
    missing = sorted(REQUIRED_COLUMNS - columns)
    if missing:
        raise CavaDataError(
            "Catalog manifest is missing required columns: " + ", ".join(missing)
        )

    catalogs: list[Catalog] = []
    for row in reader:
        normalized = {key: value or "" for key, value in row.items()}
        catalog = Catalog(normalized)
        if not catalog.name:
            raise CavaDataError("Catalog manifest contains an empty catalog_name")
        catalog.install_files  # Validate every path before Git is invoked.
        catalogs.append(catalog)

    if not catalogs:
        raise CavaDataError(
            "Catalog manifest contains no catalog rows: {}".format(resolved_location)
        )
    names = [catalog.name for catalog in catalogs]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise CavaDataError(
            "Catalog names must be unique; duplicates: " + ", ".join(duplicates)
        )
    return catalogs


def _select(catalogs: list[Catalog], name: str) -> Catalog:
    for catalog in catalogs:
        if catalog.name == name:
            return catalog
    raise CavaDataError(
        "Unknown catalog {!r}. Run 'cava_data list' to see available catalogs.".format(
            name
        )
    )


# ---------------------------------------------------------------------------
# Git LFS materialization and validation


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalize_hash(value: str) -> str:
    value = value.strip().lower()
    return value[7:] if value.startswith("sha256:") else value


def _validate_materialized_file(
    path: Path, expected_sha: str, expected_oid: str, expected_size: str
) -> None:
    if not path.is_file():
        raise CavaDataError("Git LFS did not materialize expected file: {}".format(path))
    actual_size = path.stat().st_size
    if expected_size.lower() not in ("", "pending", "unknown"):
        try:
            wanted_size = int(expected_size)
        except ValueError as exc:
            raise CavaDataError(
                "Invalid Git LFS size in manifest: {!r}".format(expected_size)
            ) from exc
        if actual_size != wanted_size:
            raise CavaDataError(
                "Size mismatch for {}: expected {}, found {}".format(
                    path.name, wanted_size, actual_size
                )
            )

    actual_sha = _sha256(path)
    for label, expected in (("SHA-256", expected_sha), ("Git LFS OID", expected_oid)):
        normalized = _normalize_hash(expected)
        if normalized not in ("", "pending", "unknown") and actual_sha != normalized:
            raise CavaDataError(
                "{} mismatch for {}: expected {}, found {}".format(
                    label, path.name, normalized, actual_sha
                )
            )


def _clone_pointer_checkout(catalog: Catalog, checkout: Path) -> None:
    environment = dict(os.environ)
    environment["GIT_LFS_SKIP_SMUDGE"] = "1"
    _run(["git", "init", "--quiet", str(checkout)], env=environment)
    _run(
        ["git", "-C", str(checkout), "remote", "add", "origin", catalog.repository],
        env=environment,
    )
    _run(
        [
            "git",
            "-C",
            str(checkout),
            "fetch",
            "--quiet",
            "--depth",
            "1",
            "origin",
            catalog.ref,
        ],
        env=environment,
    )
    _run(
        [
            "git",
            "-C",
            str(checkout),
            "checkout",
            "--quiet",
            "--detach",
            "FETCH_HEAD",
        ],
        env=environment,
    )


def _validate_installed_catalog(catalog: Catalog, destination: Path) -> list[Path]:
    targets: list[Path] = []
    for _, output_name, expected_sha, expected_oid, expected_size in catalog.install_files:
        target = destination / output_name
        _validate_materialized_file(target, expected_sha, expected_oid, expected_size)
        targets.append(target)
    return targets


def install_catalog(
    catalog: Catalog, destination: Path, *, force: bool = False, quiet: bool = False
) -> list[Path]:
    destination.mkdir(parents=True, exist_ok=True)
    targets = [destination / item[1] for item in catalog.install_files]
    existing = [target for target in targets if target.exists()]
    if existing and not force:
        if len(existing) == len(targets):
            try:
                valid = _validate_installed_catalog(catalog, destination)
            except CavaDataError as exc:
                raise CavaDataError(
                    "Existing catalog installation is invalid; use --force: {}".format(
                        exc
                    )
                ) from exc
            if not quiet:
                for target in valid:
                    print(target)
            return valid
        raise CavaDataError(
            "Partial catalog installation exists; use --force to repair: {}".format(
                ", ".join(str(path) for path in existing)
            )
        )

    _require_git_lfs()
    with tempfile.TemporaryDirectory(prefix="cava-data-") as temporary:
        checkout = Path(temporary) / "repo"
        _clone_pointer_checkout(catalog, checkout)
        include_paths = ",".join(item[0] for item in catalog.install_files)
        _run(
            [
                "git",
                "-C",
                str(checkout),
                "lfs",
                "pull",
                "--include={}".format(include_paths),
                "--exclude=",
                "origin",
            ]
        )

        for repository_path, output_name, expected_sha, expected_oid, expected_size in catalog.install_files:
            source = checkout / repository_path
            _validate_materialized_file(
                source, expected_sha, expected_oid, expected_size
            )
            target = destination / output_name
            partial = target.with_name(target.name + ".part")
            shutil.copy2(source, partial)
            os.replace(partial, target)
            if not quiet:
                print(target)
    return targets


# ---------------------------------------------------------------------------
# Commands


def _normalize_build(value: str) -> str:
    return "grch37" if value.lower() == "hg19" else value.lower()


def _format_rows(catalogs: list[Catalog]) -> str:
    headers = [
        "catalog_name",
        "Build",
        "source",
        "source_version",
        "transcript_ID_type",
        "status",
        "file_name",
    ]
    rows = [headers] + [
        [catalog.row.get(column, "") for column in headers] for catalog in catalogs
    ]
    widths = [max(len(row[index]) for row in rows) for index in range(len(headers))]
    return "\n".join(
        "  ".join(value.ljust(widths[index]) for index, value in enumerate(row)).rstrip()
        for row in rows
    )


def command_list(args: argparse.Namespace) -> int:
    catalogs = read_manifest(args.manifest)
    if args.build:
        wanted = _normalize_build(args.build)
        catalogs = [
            catalog
            for catalog in catalogs
            if _normalize_build(catalog.row.get("Build", "")) == wanted
        ]
    if args.source:
        catalogs = [
            catalog
            for catalog in catalogs
            if catalog.row.get("source", "").lower() == args.source.lower()
        ]
    if not args.all_statuses:
        catalogs = [
            catalog
            for catalog in catalogs
            if catalog.row.get("status", "published").lower() == "published"
        ]

    if args.tsv:
        fieldnames = list(catalogs[0].row) if catalogs else sorted(REQUIRED_COLUMNS)
        writer = csv.DictWriter(
            sys.stdout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for catalog in catalogs:
            writer.writerow(catalog.row)
    elif catalogs:
        print(_format_rows(catalogs))
    else:
        print("No matching published catalogs were found.", file=sys.stderr)
    return 0


def _require_publishable(catalog: Catalog, allow_unpublished: bool) -> None:
    status = catalog.row.get("status", "published").lower()
    if status != "published" and not allow_unpublished:
        raise CavaDataError(
            "Catalog {!r} has status {!r}; use --allow-unpublished to use it".format(
                catalog.name, status
            )
        )


def command_install(args: argparse.Namespace) -> int:
    catalog = _select(read_manifest(args.manifest), args.catalog_name)
    _require_publishable(catalog, args.allow_unpublished)
    install_catalog(
        catalog, Path(args.path).expanduser().resolve(), force=args.force
    )
    return 0


def _candidate_templates(explicit: str | None):
    if explicit:
        yield Path(explicit).expanduser()
    environment_template = os.environ.get("CAVA_CONFIG_TEMPLATE")
    if environment_template:
        yield Path(environment_template).expanduser()
    yield Path.cwd() / "config_template.txt"
    yield Path(__file__).resolve().parents[1] / "config_template.txt"
    yield Path(sys.prefix) / "share" / "cava" / "config_template.txt"


def _find_template(explicit: str | None) -> Path:
    checked: list[str] = []
    for candidate in _candidate_templates(explicit):
        checked.append(str(candidate))
        if candidate.is_file():
            return candidate.resolve()
    raise CavaDataError(
        "Could not locate config_template.txt. Checked: " + ", ".join(checked)
    )


def _render_config(template: str, database: Path) -> str:
    pattern = re.compile(
        r"^(?P<prefix>\s*@ensembl\s*=\s*)(?P<value>.*?)(?P<eol>\r?\n)?$",
        re.IGNORECASE,
    )
    output: list[str] = []
    replacements = 0
    for line in template.splitlines(True):
        match = pattern.match(line)
        if match:
            output.append(
                "{}{}{}".format(
                    match.group("prefix"), database, match.group("eol") or ""
                )
            )
            replacements += 1
        else:
            output.append(line)
    if replacements != 1:
        raise CavaDataError(
            "config_template.txt must contain exactly one active '@ensembl = ...' "
            "setting; found {}".format(replacements)
        )
    return "".join(output)


def command_config(args: argparse.Namespace) -> int:
    catalog = _select(read_manifest(args.manifest), args.catalog_name)
    _require_publishable(catalog, args.allow_unpublished)
    destination = Path(args.path).expanduser().resolve()

    try:
        _validate_installed_catalog(catalog, destination)
        installed = True
    except CavaDataError:
        installed = False
    if args.force_install or not installed:
        install_catalog(catalog, destination, force=True, quiet=True)
    _validate_installed_catalog(catalog, destination)

    template = _find_template(args.template)
    output = (
        Path(args.output).expanduser().resolve()
        if args.output
        else destination / "config.{}.txt".format(catalog.name)
    )
    if output.exists() and not args.force:
        raise CavaDataError("Config file already exists: {}; use --force".format(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    rendered = _render_config(
        template.read_text(encoding="utf-8"),
        (destination / catalog.file_name).resolve(),
    )
    partial = output.with_name(output.name + ".part")
    partial.write_text(rendered, encoding="utf-8")
    os.replace(partial, output)
    print(output)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cava_data",
        description="List and install CAVA transcript catalogs stored with Git LFS.",
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help=(
            "Local path, file:// URL or HTTP(S) URL for cava_catalogs.tsv. "
            "Default: CAVA_CATALOG_MANIFEST, then the current CAVA master "
            "manifest, then the installed snapshot."
        ),
    )
    subparsers = parser.add_subparsers(dest="command")

    list_parser = subparsers.add_parser("list", help="List available catalogs")
    list_parser.add_argument("--build", choices=("GRCh38", "GRCh37", "hg19"))
    list_parser.add_argument("--source")
    list_parser.add_argument("--all-statuses", action="store_true")
    list_parser.add_argument("--tsv", action="store_true")
    list_parser.set_defaults(func=command_list)

    install_parser = subparsers.add_parser("install", help="Install one catalog")
    install_parser.add_argument("catalog_name")
    install_parser.add_argument("path", nargs="?", default=DEFAULT_INSTALL_DIR)
    install_parser.add_argument("--force", action="store_true")
    install_parser.add_argument("--allow-unpublished", action="store_true")
    install_parser.set_defaults(func=command_install)

    config_parser = subparsers.add_parser(
        "config", help="Install one catalog if needed and create a CAVA config"
    )
    config_parser.add_argument("catalog_name")
    config_parser.add_argument("path", nargs="?", default=DEFAULT_INSTALL_DIR)
    config_parser.add_argument("--output")
    config_parser.add_argument("--template")
    config_parser.add_argument("--force", action="store_true")
    config_parser.add_argument("--force-install", action="store_true")
    config_parser.add_argument("--allow-unpublished", action="store_true")
    config_parser.set_defaults(func=command_config)
    return parser


def main(argv: list[str] | None = None) -> int:
    try:
        parser = build_parser()
        args = parser.parse_args(argv)
        if not hasattr(args, "func"):
            parser.print_help(sys.stderr)
            return 2
        return int(args.func(args))
    except CavaDataError as exc:
        print("cava_data: error: {}".format(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
