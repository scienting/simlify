FROM continuumio/miniconda3:latest

LABEL org.opencontainers.image.title="simlify-env"
LABEL org.opencontainers.image.description="Environment for testing simlify Python package."
LABEL org.opencontainers.image.authors="OASCI <us@oasci.org>"
LABEL org.opencontainers.image.source="https://gitlab.com/oasci/software/simlify"

RUN apt-get update
RUN apt-get install build-essential -y
RUN apt-get install curl -y

WORKDIR /simlify
COPY Makefile .
COPY conda-lock.yml .
COPY poetry.lock .
COPY pyproject.toml .

RUN make conda-create
RUN make conda-setup
RUN make from-conda-lock
RUN conda run -n simlify-dev poetry install --no-root
