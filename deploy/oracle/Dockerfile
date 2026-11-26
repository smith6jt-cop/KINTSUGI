# KINTSUGI MCP server, streamable HTTP transport.
#
# Multi-architecture: builds natively on Oracle Cloud Ampere A1 (linux/arm64)
# and on x86_64. Build from the repository root:
#
#   docker build -f deploy/oracle/Dockerfile -t kintsugi-mcp .
#
# The container serves one KINTSUGI project mounted at /project and needs
# KINTSUGI_MCP_TOKEN (>= 32 characters) and KINTSUGI_MCP_PUBLIC_URL (the HTTPS
# origin of the reverse proxy). See docs/hubmap_deployment.md.
FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# libvips for pyvips; libglib for OpenCV; curl for the health check.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libvips42 libglib2.0-0 curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/kintsugi
COPY pyproject.toml constraints.txt README.md License.txt ./
COPY src ./src
RUN pip install ".[mcp]" -c constraints.txt \
    && kintsugi mcp tools > /dev/null

RUN useradd --create-home --uid 1000 kintsugi \
    && mkdir -p /project /home/kintsugi/.kintsugi \
    && chown kintsugi:kintsugi /project /home/kintsugi/.kintsugi
USER kintsugi

# The root filesystem is read-only in compose.yaml: caches go to the /tmp tmpfs;
# the learning database to the ~/.kintsugi volume.
ENV KINTSUGI_PROJECT_DIR=/project \
    MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp/.cache
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
    CMD curl -fsS http://127.0.0.1:8000/healthz || exit 1

ENTRYPOINT ["kintsugi", "mcp", "start", "--transport", "streamable-http", \
            "--host", "0.0.0.0", "--port", "8000", "--project-dir", "/project"]
