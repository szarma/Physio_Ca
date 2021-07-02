FROM jupyter/scipy-notebook:latest

# update package manager and create physio environment:
RUN python3 -m pip install --upgrade pip
RUN conda create \
          -n Physio \
          -c conda-forge \
          python \
          cython \
          caiman \
          jupyter_contrib_nbextensions \ 
          dash \
          cdsdashboards-singleuser \
          ffmpeg

# copy dash-specific config files:
RUN rm -rf /home/jovyan/work
ARG SERVER_DOCKER_PATH
ARG JUPYTER_NOTEBOOK_CONFIG_EXTRA="$SERVER_DOCKER_PATH/jupyter_notebook_config_extra.py"
USER $NB_UID
COPY $SERVER_DOCKER_PATH/jupyter_notebook_config_extra.py /etc/jupyter/
COPY $SERVER_DOCKER_PATH/voila.json /etc/jupyter/
COPY $SERVER_DOCKER_PATH/github-tokens.sh /usr/local/bin/before-notebook.d/github-tokens.sh
RUN cat /etc/jupyter/jupyter_notebook_config_extra.py >> /etc/jupyter/jupyter_notebook_config.py
RUN rm /etc/jupyter/jupyter_notebook_config_extra.py
RUN conda init bash

# prepare installation of islets:
USER root
COPY ./islets/ /opt/islets/
COPY ./app/ /opt/app/
COPY ${SERVER_DOCKER_PATH}/setup.py /opt/setup.py
COPY requirements.txt /opt/requirements.txt
COPY README.md /opt/README.md
COPY VERSION /opt/VERSION

# install jvm for bioformats and fix permissions:
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential \
    default-jdk
RUN ln -s /opt/app/ /home/jovyan/app
RUN /usr/local/bin/fix-permissions "/opt/"

# install islets and enable kernel:
USER $NB_UID
RUN /opt/conda/envs/Physio/bin/pip install -e /opt/
RUN /opt/conda/envs/Physio/bin/python3 -m ipykernel install --user --name Physio --display-name "Physio"
