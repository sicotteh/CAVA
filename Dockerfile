FROM ubuntu:latest
RUN apt-get update
RUN apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common \
    python3 \
    python3-distutils

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3 get-pip.py