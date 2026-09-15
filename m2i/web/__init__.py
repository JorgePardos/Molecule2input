"""The web interface: a FastAPI app serving ``/api`` and a static front end.

Launch it with ``python -m m2i.cli gui``, or ``uvicorn m2i.web.main:app``.
Nothing in here is chemistry: every route is a thin wrapper around the engine
(``m2i.pipeline``, ``m2i.crystal``, ``m2i.organometallic``), and anything that
looks like a decision about a molecule belongs there instead.
"""
