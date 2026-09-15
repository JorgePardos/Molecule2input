# The web interface as a container. deploy/lab/ runs it on a machine of your
# own behind a Cloudflare tunnel; it runs the same on any Docker host
# (docker build -t m2i . && docker run -p 7860:7860 m2i).
FROM python:3.11-slim

# RDKit draws the 2D depictions through libXrender; OpenCV, which DECIMER
# reads pictures with, needs libGL and GLib. The slim image has none of them.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libxrender1 libxext6 libgl1 libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# An unprivileged user; owning the home directory keeps the per-session
# temporary folders writable, and puts DECIMER's weights (which it keeps in
# ~/.data) inside the image.
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    M2I_HOSTED=1 \
    M2I_MAX_UPLOAD_MB=20

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
# Healthy as soon as the server answers: the photo model keeps loading in the
# background, and the page shows that on its own.
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:7860/api/health')"
# One worker: the sessions and the loaded model live in that process's memory.
# --proxy-headers: behind the tunnel the visitor's connection is HTTPS, and the
# session cookie is then marked secure.
CMD ["uvicorn", "m2i.web.main:app", "--host", "0.0.0.0", "--port", "7860", "--workers", "1", \
     "--proxy-headers", "--forwarded-allow-ips", "*"]
