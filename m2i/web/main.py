"""The web application: ``/api`` routes and the static front end, on one port.

    uvicorn m2i.web.main:app --host 0.0.0.0 --port 7860
"""

from __future__ import annotations

from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

from .. import __version__, config
from ..gui.hosting import HOSTED
from ..recognition.manual import STRUCTURE_FILE_SUFFIXES
from ..recognition.registry import describe_backends
from . import complexes, crystal, organic, recipes, warmup

STATIC = Path(__file__).parent / "static"
VIEWER = Path(__file__).parent.parent / "gui" / "viewer"


@asynccontextmanager
async def lifespan(app: FastAPI):
    warmup.start()  # the photo model loads while the first visitor chooses a file
    yield


app = FastAPI(title="m2i", version=__version__, lifespan=lifespan, docs_url=None, redoc_url=None)


@app.exception_handler(HTTPException)
async def engine_message(request: Request, exc: HTTPException):
    """Errors as ``{"error": message}``: the engine's own words, shown verbatim."""
    return JSONResponse({"error": exc.detail}, status_code=exc.status_code)


@app.middleware("http")
async def upload_limit(request: Request, call_next):
    """Refuse an oversized upload before it is read, not after it fills the disk."""
    length = request.headers.get("content-length")
    limit = organic.max_upload_bytes() + 2**20  # room for the multipart envelope
    if request.method == "POST" and length and length.isdigit() and int(length) > limit:
        return JSONResponse(
            {"error": f"That file is larger than {organic.max_upload_bytes() // 2**20} MB."}, status_code=413
        )
    return await call_next(request)


@app.middleware("http")
async def security_headers(request: Request, call_next):
    response = await call_next(request)
    if not request.url.path.startswith("/api/"):
        # Revalidate the front end on every load (cheap, with ETags), so an
        # update to the server reaches browsers without a hard refresh.
        response.headers.setdefault("Cache-Control", "no-cache")
    response.headers.setdefault("X-Content-Type-Options", "nosniff")
    response.headers.setdefault("Referrer-Policy", "same-origin")
    response.headers.setdefault("X-Frame-Options", "SAMEORIGIN")
    response.headers.setdefault(
        "Content-Security-Policy",
        # Everything is served from here, fonts included; 3Dmol needs eval-free
        # WebGL only, and inline styles are used for computed widths.
        "default-src 'self'; img-src 'self' data: blob:; style-src 'self' 'unsafe-inline'; "
        "script-src 'self'; connect-src 'self'; frame-ancestors 'self'",
    )
    return response


app.include_router(organic.router)
app.include_router(crystal.router)
app.include_router(complexes.router)


@app.get("/api/health")
def health():
    """Always 200 while the server runs: Docker and the tunnel wait on this."""
    return warmup.status()


@app.get("/api/about")
def about():
    return {
        "version": __version__,
        "hosted": HOSTED,
        "max_upload_mb": organic.max_upload_bytes() // 2**20,
        "types": {
            "picture": [s.lstrip(".") for s in organic.IMAGE_SUFFIXES],
            "file": [s.lstrip(".") for s in STRUCTURE_FILE_SUFFIXES],
            "cif": ["cif"],
        },
        "backends": [row for row in describe_backends() if row["name"] != "manual"],
    }


@app.get("/api/profiles")
def profiles(program: str = "gaussian"):
    return {
        "program": program,
        "qm": program in recipes.QM_PROGRAMS,
        "recipes": recipes.recipes(program) if program in recipes.QM_PROGRAMS else [],
        "dispersion": recipes.dispersion_choices(program),
    }


@app.post("/api/profile/validate")
def validate_profile(settings: recipes.Settings):
    try:
        profile = recipes.build(settings)
    except config.ProfileError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {"ok": True, "jobs": " ".join(profile.jobs) or "single point"}



app.mount("/viewer", StaticFiles(directory=VIEWER), name="viewer")
app.mount("/", StaticFiles(directory=STATIC, html=True), name="web")
