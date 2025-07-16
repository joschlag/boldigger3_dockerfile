FROM mcr.microsoft.com/playwright/python:v1.53.0

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    git \
    zlib1g-dev \
    libhdf5-serial-dev \
    sudo && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip && \
    pip install packaging tables httpx

WORKDIR /home

RUN git clone https://github.com/DominikBuchner/BOLDigger3.git

WORKDIR /home/BOLDigger3

RUN python setup.py install
