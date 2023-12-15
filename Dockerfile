FROM ubuntu:latest

LABEL org.opencontainers.image.title="simlify"
LABEL org.opencontainers.image.description="Environment for testing simlify Python package."
LABEL org.opencontainers.image.authors="OASCI <us@oasci.org>"
LABEL org.opencontainers.image.source="https://gitlab.com/oasci/software/simlify"

RUN apt-get update \
    && apt-get install --no-install-recommends g++ gcc make ca-certificates wget git -y \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh \
    && /bin/sh ~/miniconda.sh -b -p $CONDA_DIR \
    && rm -f ~/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

WORKDIR /simlify

RUN conda create -n simlify-dev
ENV CONDA_RUN conda run -n simlify-dev
ENV CONDA_INSTALL $CONDA_RUN conda install -y
RUN $CONDA_INSTALL -y python=3.11 --freeze-installed \
    && $CONDA_INSTALL -c conda-forge --freeze-installed poetry conda-poetry-liaison ambertools \
    && conda clean -afy \
    && $CONDA_RUN cpl-clean --env_name simlify-dev

ENV POETRY_DYNAMIC_VERSIONING_BYPASS 0.0.0
COPY pyproject.toml .
RUN $CONDA_RUN poetry install --no-interaction --no-root \
    && $CONDA_RUN poetry cache clear pypi --all
