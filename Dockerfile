FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libffi-dev \
    libssl-dev \
    libsqlite3-dev \
    libxml2-dev \
    libxslt1-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip \
    && pip install boldigger3==2.0.0 \
    && pip install "lxml[html_clean]"

RUN mkdir -p /data
WORKDIR /data

COPY modifications/metadata_download.py /tmp/
COPY modifications/__main__.py /tmp/

# Overwrite the installed package files with your patched versions
RUN set -eux; \
    pkg_path=$(python -c "import boldigger3, os; print(os.path.dirname(boldigger3.__file__))"); \
    cp /tmp/metadata_download.py "$pkg_path/metadata_download.py"; \
    cp /tmp/__main__.py "$pkg_path/__main__.py"

# Default command
CMD ["python", "-m", "boldigger3"]
