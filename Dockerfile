FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=8501

WORKDIR /app

COPY requirements.deploy.txt ./
RUN python -m pip install --upgrade pip && \
    python -m pip install -r requirements.deploy.txt

COPY . .

RUN mkdir -p /app/outputs /app/outputs_mouse /app/outputs/cache /app/outputs/history
RUN chmod +x /app/docker-entrypoint.sh

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=5s --start-period=90s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=3)" || exit 1

ENTRYPOINT ["/app/docker-entrypoint.sh"]
