FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libffi-dev \
    libssl-dev \
    libsqlite3-dev \
    libxml2-dev \
    libxslt1-dev \
    && pip install --upgrade pip \
    && pip install boldigger3==2.0.0 \
    && pip install "lxml[html_clean]" \
    && apt-get remove -y gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set a writable directory (Singularity will bind-mount over this)
RUN mkdir -p /data && cp -r /usr/local/lib/python3.11/site-packages/boldigger3 /data/boldigger3

WORKDIR /data

CMD ["python", "-m", "boldigger3"]
