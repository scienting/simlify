FROM ubuntu:latest

LABEL org.opencontainers.image.title="simlify"
LABEL org.opencontainers.image.description="Environment for testing simlify Python package."
LABEL org.opencontainers.image.authors="OASCI <us@oasci.org>"
LABEL org.opencontainers.image.source="https://gitlab.com/oasci/software/simlify"

RUN apt-get update
RUN apt-get install build-essential -y
RUN apt-get install wget -y

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR
ENV PATH=$CONDA_DIR/bin:$PATH

WORKDIR /simlify
COPY Makefile .
COPY conda-lock.yml .
COPY poetry.lock .
COPY pyproject.toml .

RUN make conda-create
RUN make conda-setup
RUN make from-conda-lock
RUN conda run -n simlify-dev poetry install --no-root
