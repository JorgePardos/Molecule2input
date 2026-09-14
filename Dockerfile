# The browser interface as a container. deploy/lab/ runs it on a machine of
# your own behind a Cloudflare tunnel; it runs the same on any Docker host
# (docker build -t m2i . && docker run -p 7860:7860 m2i).
FROM python:3.11-slim

# RDKit draws the 2D check image through libXrender; OpenCV, which DECIMER
# reads pictures with, needs libGL and GLib. The slim image has none of them.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libxrender1 libxext6 libgl1 libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# An unprivileged user (uid 1000, as Spaces expect); owning the home keeps the
# per-session temporary folders writable, and puts DECIMER's weights (which
# it keeps in ~/.data) inside the image.
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    M2I_HOSTED=1 \
    STREAMLIT_SERVER_PORT=7860 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_MAX_UPLOAD_SIZE=20 \
    STREAMLIT_SERVER_ENABLE_CORS=true \
    STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION=true
# Served on its own address (a tunnel, a reverse proxy), Streamlit's XSRF
# protection works and stays on. Only where the page is shown inside an iframe
# on another domain -- a Hugging Face Space -- does its cookie never arrive;
# there both settings have to be false, in the Space's variables.

WORKDIR $HOME/app
# Layers from slowest to change to fastest, so that a code change rebuilds in
# seconds instead of reinstalling TensorFlow.
COPY --chown=user requirements-web.txt .
RUN pip install --no-cache-dir --user -r requirements-web.txt

# DECIMER, for hand-drawn pictures, in its own environment exactly as
# `m2i setup decimer` makes it locally: the packages, the weights, and a real
# reading of a probe drawing before it is marked ready. Only the modules the
# installer needs are copied here.
COPY --chown=user requirements-decimer.txt .
COPY --chown=user m2i/__init__.py m2i/backends.py m2i/types.py ./m2i/
COPY --chown=user m2i/recognition ./m2i/recognition
RUN PIP_CONSTRAINT=$HOME/app/requirements-decimer.txt PIP_NO_CACHE_DIR=1 \
    python -c "from m2i.backends import install; install('decimer')"

COPY --chown=user m2i ./m2i

EXPOSE 7860
HEALTHCHECK --interval=60s --timeout=5s --start-period=120s \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:7860/_stcore/health')"
CMD ["streamlit", "run", "m2i/gui/app.py"]
