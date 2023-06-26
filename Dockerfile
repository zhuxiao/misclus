# Use Ubuntu as the base image
FROM ubuntu

# Install necessary packages
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    g++ \
    gcc \
    git \
    libbz2-dev\
    libcurl4-gnutls-dev \
    liblzma-dev \
    libssl-dev \
    libtool \
    make \
    mlpack-bin \
    perl \
    python3 \
    zlib1g-dev

# Download and install htslib from source
RUN git clone https://github.com/samtools/htslib.git && \
    cd htslib && \
    git submodule update --init --recursive && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install

# Download and compile misclas
RUN git clone https://github.com/zhuxiao/misclas.git && \
    cd misclas && \
    bash autogen.sh


RUN ln -s /misclas/bin/misclas_bin /usr/local/bin/
RUN ln -s /misclas/misclas.py /usr/local/bin/

# Set the entry point or run any desired command
ENTRYPOINT ["/usr/local/bin/misclas.py"]
