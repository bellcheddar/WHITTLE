"""Record a QC result that cannot be measured from Python. Spec Section 10.

Two of the ten QC items in spec Section 10 can only be checked in a real
browser: item 7, that the Mol* WebGL context count stays bounded after 30
consecutive drawer opens, and item 8, that there is no horizontal page scroll
at phone width. Playwright's Chromium is explicitly not to be installed
(spec 0.2), so those are run against a real browser and the outcome is
recorded here, with the measurement that decided it, rather than asserted.

    python tools/record_qc.py molstar_contexts pass \
        "30 opens, context count peaked at 1, viewer alive"
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.common import connect, log  # noqa: E402

KNOWN = {
    "molstar_contexts": "Spec 10 item 7: the Mol* WebGL context count stays "
                        "bounded after 30 consecutive drawer opens.",
    "phone_width": "Spec 10 item 8: no horizontal page scroll at phone width; "
                   "the explorer collapses and the drawer stacks.",
    "drawer_deep_link": "A deep link opens the drawer directly.",
}


def main() -> int:
    if len(sys.argv) < 4:
        print(__doc__)
        print("known keys: " + ", ".join(sorted(KNOWN)))
        return 2
    key, verdict, detail = sys.argv[1], sys.argv[2].lower(), " ".join(sys.argv[3:])
    if verdict not in ("pass", "fail"):
        print("the verdict must be pass or fail")
        return 2
    note = KNOWN.get(key, "")
    conn = connect()
    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 (key, f"{detail}" + (f" | {note}" if note else ""),
                  1 if verdict == "pass" else 0))
    conn.commit()
    conn.close()
    log(f"qc: recorded {key} = {verdict.upper()} ({detail})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
