"""Small atomic writers for deterministic research artifacts."""

from __future__ import annotations

import json
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, TextIO


@contextmanager
def atomic_text_writer(path: Path) -> Iterator[TextIO]:
    """Yield a same-directory temporary stream, then atomically replace path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            yield stream
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            temporary.unlink(missing_ok=True)
        finally:
            raise


def atomic_write_json(path: Path, value: Any) -> None:
    with atomic_text_writer(path) as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
