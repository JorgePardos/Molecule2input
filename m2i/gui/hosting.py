"""Whether the interface is served to other people rather than run on your own machine.

Set ``M2I_HOSTED=1`` (the Docker image does). A few things only make sense
locally and are left out when hosted: a folder of the server to save into,
where the files were written on its disk, and ``m2i setup`` instructions
nobody visiting the page can follow.
"""

from __future__ import annotations

import os

HOSTED = os.environ.get("M2I_HOSTED", "").strip().lower() in ("1", "true", "yes")
