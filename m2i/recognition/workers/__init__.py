"""Standalone scripts executed inside the per-backend virtual environments.

Nothing here may import m2i: these run under a different interpreter with a
different, deliberately incompatible set of dependencies. They talk to m2i only
through the JSON request/response files described in ``_subprocess.py``.
"""
