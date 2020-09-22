FROM ubuntu:latest

RUN export DEBIAN_FRONTEND=noninteractive; \
    export DEBCONF_NONINTERACTIVE_SEEN=true; \
    echo 'tzdata tzdata/Areas select Etc' | debconf-set-selections; \
    echo 'tzdata tzdata/Zones/Etc select UTC' | debconf-set-selections; \
    apt-get update -qqy \
 && apt-get install -qqy --no-install-recommends \
        tzdata \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update

RUN apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    build-essential \
    autoconf \
    libtool \
    pkg-config \
    gnupg-agent \
    software-properties-common \
    python3 \
    python3-dev \
    python3-distutils \
    git \
    libcurl4-nss-dev \
    zlib1g-dev


RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3 get-pip.py
RUN pip install appdirs Cython flake8 packaging py pylint pyparsing pysam pytest pytest-cov radon six sphinx gevent CrossMap

RUN git clone https://github.com/Steven-N-Hart/CAVA.git
RUN cd CAVA && python3 setup.py install && cd /
