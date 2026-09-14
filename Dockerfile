# The browser interface as a container: what Hugging Face Spaces builds, and
# what runs the same anywhere else (docker build -t m2i . && docker run -p 7860:7860 m2i).
FROM python:3.11-slim

# RDKit draws the 2D check image through libXrender, which the slim image lacks.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Spaces run the container as uid 1000; owning the home directory keeps the
# per-session temporary folders writable.
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
    STREAMLIT_SERVER_ENABLE_CORS=false \
    STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION=false
# XSRF protection is off because a Space shows the app inside an iframe on
# another domain, where Streamlit's XSRF cookie never arrives and every upload
# would be refused. Nothing here is behind a login for such a request to abuse.

WORKDIR $HOME/app
# Dependencies first, pinned: a code change then rebuilds in seconds, and the
# image is the combination the test suite passed with.
COPY --chown=user requirements-web.txt .
RUN pip install --no-cache-dir --user -r requirements-web.txt

COPY --chown=user m2i ./m2i

EXPOSE 7860
HEALTHCHECK --interval=60s --timeout=5s --start-period=60s \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:7860/_stcore/health')"
CMD ["streamlit", "run", "m2i/gui/app.py"]
