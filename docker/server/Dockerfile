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
RUN pip install --no-cache \
        ffmpeg \
        tifffile

WORKDIR /opt/islets/
RUN git clone https://github.com/szarma/Physio_Ca.git

WORKDIR /opt/islets/Physio_Ca/
RUN pip install .

WORKDIR ${HOME}
RUN mkdir -p ${HOME}/work/app/
RUN cp /opt/islets/Physio_Ca/app/explorer.py ${HOME}/work/app/explorer.py
