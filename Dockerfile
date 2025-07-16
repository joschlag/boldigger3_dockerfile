FROM mcr.microsoft.com/playwright/python:v1.53.0

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PLAYWRIGHT_BROWSERS_PATH=/ms-playwright

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    git \
    zlib1g-dev \
    libhdf5-serial-dev \
    sudo && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip setuptools packaging tables httpx

WORKDIR /home

RUN git clone https://github.com/DominikBuchner/BOLDigger3.git

WORKDIR /home/BOLDigger3

RUN pip install .

RUN playwright install --with-deps
