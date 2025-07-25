FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive
#WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libffi-dev \
    libssl-dev \
    libsqlite3-dev \
    libxml2-dev \
    libxslt1-dev \
    build-essential

 RUN pip install --upgrade pip \
    && pip install boldigger3==2.0.1 \
    && pip install "lxml[html_clean]" \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data

# Copy patched scripts from local repo into container
# Make sure these are in the same GitHub repo where the Dockerfile lives
COPY modifications/metadata_download.py /usr/local/lib/python3.11/site-packages/boldigger3/metadata_download.py
COPY modifications/__main__.py /usr/local/lib/python3.11/site-packages/boldigger3/__main__.py
COPY modifications/add_metadata.py /usr/local/lib/python3.11/site-packages/boldigger3/add_metadata.py

# Default command
ENTRYPOINT ["python", "-m", "boldigger3"]
