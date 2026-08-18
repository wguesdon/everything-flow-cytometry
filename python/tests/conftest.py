"""Put the scripts folder on the import path.

The analysis modules live in `scripts/` next to the scripts that run them, and
the Python project is not installed as a package, so the tests add that folder
rather than importing from an installed distribution.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
