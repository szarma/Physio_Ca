FROM ideonate/containds-all-scipy:latest

USER root
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential \
    git \
    default-jdk
RUN mkdir -p /opt/islets/Physio_Ca_framework
RUN mkdir /data
RUN /usr/local/bin/fix-permissions /opt/islets/
RUN /usr/local/bin/fix-permissions /data
RUN chown -R ${NB_UID}:users ${HOME}
RUN chown -R ${NB_UID}:users /data

USER ${NB_UID}
RUN curl -sSL https://install.python-poetry.org | python3 -

COPY pyproject.toml /opt/islets/Physio_Ca_framework/pyproject.toml
COPY poetry.lock /opt/islets/Physio_Ca_framework/poetry.lock
COPY README.rst /opt/islets/Physio_Ca_framework/README.rst
COPY islets /opt/islets/Physio_Ca_framework/islets
COPY scripts /opt/islets/Physio_Ca_framework/scripts

WORKDIR /opt/islets/Physio_Ca_framework/
ENV PATH="${HOME}/.local/bin:${PATH}"
RUN poetry config virtualenvs.create false
RUN poetry install

WORKDIR ${HOME}
