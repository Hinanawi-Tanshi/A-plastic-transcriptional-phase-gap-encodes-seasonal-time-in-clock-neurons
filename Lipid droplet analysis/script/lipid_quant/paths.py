"""Shared path discovery for single- and multi-directory workflows."""
from pathlib import Path
from typing import Sequence

PathInput = str | Path | Sequence[str | Path]


def normalize_paths(paths: PathInput) -> list[Path]:
    values = [paths] if isinstance(paths, (str, Path)) else list(paths)
    if not values:
        raise ValueError("Provide at least one input path")
    unique: dict[str, Path] = {}
    for value in values:
        path = Path(value).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Input path does not exist: {path}")
        unique.setdefault(str(path).casefold(), path)
    return list(unique.values())


def discover_ims_files(inputs: PathInput, pattern: str = "*.ims", *, recursive: bool = False) -> list[Path]:
    """Expand each input, deduplicating overlapping directories/files by path."""
    found: dict[str, Path] = {}
    for path in normalize_paths(inputs):
        if path.is_file():
            if path.suffix.lower() != ".ims":
                raise ValueError(f"Expected an IMS file: {path}")
            matches = [path]
        else:
            matches = sorted(p for p in (path.rglob(pattern) if recursive else path.glob(pattern))
                             if p.is_file() and p.suffix.lower() == ".ims")
            if not matches:
                raise FileNotFoundError(f"No IMS files matching {pattern!r} in {path}")
        for item in matches:
            found.setdefault(str(item.resolve()).casefold(), item.resolve())
    return sorted(found.values(), key=lambda p: str(p).casefold())
