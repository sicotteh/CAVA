from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().parent
_CAVA_PACKAGE_DIR = _PACKAGE_DIR.parent / "cava" / "ensembldb"

__path__ = [str(_PACKAGE_DIR), str(_CAVA_PACKAGE_DIR)]
__all__ = []
