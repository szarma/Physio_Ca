ARG HUB_VERSION
FROM jupyter/scipy-notebook:latest

ARG SERVER_DOCKER_PATH
ARG JUPYTER_NOTEBOOK_CONFIG_EXTRA="${SERVER_DOCKER_PATH}/jupyter_notebook_config_extra.py"

# install jvm for bioformats
USER root
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential \
    default-jdk

USER ${NB_UID}
WORKDIR /tmp
RUN mkdir ${HOME}/islets/
RUN pip3 install --upgrade --no-cache setuptools && \
    python3 -m pip install --upgrade --no-cache pip
RUN pip3 install --no-cache \
         dash \
         cdsdashboards[user] \
         ffmpeg \
         jupyter_contrib_nbextensions

COPY ${SERVER_DOCKER_PATH}/jupyter_notebook_config_extra.py /etc/jupyter/
COPY ${SERVER_DOCKER_PATH}/voila.json /etc/jupyter/
COPY ${SERVER_DOCKER_PATH}/github-tokens.sh /usr/local/bin/before-notebook.d/github-tokens.sh
RUN cat /etc/jupyter/jupyter_notebook_config_extra.py >> /etc/jupyter/jupyter_notebook_config.py
RUN rm /etc/jupyter/jupyter_notebook_config_extra.py

# prepare installation of islets:
COPY ./islets/ ${HOME}/islets/
COPY ./app/ ${HOME}/app/
COPY setup.py ${HOME}/islets/setup.py
COPY requirements.txt ${HOME}/islets/requirements.txt
COPY README.md ${HOME}/islets/README.md
COPY VERSION ${HOME}/islets/VERSION

# install islets and enable kernel:
WORKDIR ${HOME}/islets/
RUN /usr/local/bin/fix-permissions "${HOME}/islets/"
RUN python3 setup.py bdist_wheel
RUN pip3 install "dist/islets-0.5.0-py3-none-any.whl"

WORKDIR "${HOME}"
