"""Keep the test suite out of the real build log.

`pipeline.common.log` appends every line to `BUILD_LOG.md`, which is a
deliverable: spec Section 12 lists it as part of the build record. It already
honours `WHITTLE_LOG` so that a test run can point it elsewhere, and nothing
was setting it, so any test that exercised a code path with a `log` call wrote
into the real record.

That stayed invisible until a test was added that deliberately triggers
`score_external`'s comparability refusal. It logged
`EGFR specificity is not comparable and stays empty` against candidate ids of
the form EGFR-0001, which is indistinguishable from a real pipeline event, and
wrote one such line on every run of the suite. Eight reached BUILD_LOG.md
before this fixture existed and were removed by hand; the note in the log says
so.

Session-scoped and autouse, so it covers collection-time logging too.
"""
from __future__ import annotations

import os
import pytest


@pytest.fixture(scope="session", autouse=True)
def _log_to_a_scratch_file(tmp_path_factory):
    previous = os.environ.get("WHITTLE_LOG")
    os.environ["WHITTLE_LOG"] = str(
        tmp_path_factory.mktemp("logs") / "BUILD_LOG.md")
    yield
    if previous is None:
        os.environ.pop("WHITTLE_LOG", None)
    else:
        os.environ["WHITTLE_LOG"] = previous
