FROM ideonate/containds-all-scipy:latest

USER root
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential \
    git \
    default-jdk
RUN mkdir /opt/islets/
RUN /usr/local/bin/fix-permissions /opt/islets/
RUN chown -R ${NB_UID}:users ${HOME}

USER ${NB_UID}
RUN curl -sSL https://install.python-poetry.org | python3 -

WORKDIR /opt/islets/
RUN git clone https://github.com/Hannnsen/Physio_Ca_framework.git

WORKDIR /opt/islets/Physio_Ca_framework/
ENV PATH="${HOME}/.local/bin:${PATH}"
# RUN poetry config virtualenvs.create false
RUN poetry build --format wheel
RUN pip install --no-cache dist/*.whl
