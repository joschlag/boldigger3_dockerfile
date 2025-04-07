FROM mcr.microsoft.com/playwright/python:v1.51.0-noble

RUN apt update && \
	apt upgrade -y && \
	apt install -y build-essential wget zlib1g zlib1g-dev python3 pip && \
	apt clean
RUN apt-get install -y git
RUN apt-get install -y libhdf5-serial-dev
RUN pip install packaging tables httpx
RUN cd /home && \
	git clone https://github.com/DominikBuchner/BOLDigger3.git
WORKDIR /home/BOLDigger3
RUN python3 setup.py install
RUN apt install -y sudo
