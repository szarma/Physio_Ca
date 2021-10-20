ARG HUB_VERSION
FROM jupyter/scipy-notebook:latest

ARG SERVER_DOCKER_PATH
ARG ISLETS_VERSION
ARG JUPYTER_NOTEBOOK_CONFIG_EXTRA="${SERVER_DOCKER_PATH}/jupyter_notebook_config_extra.py"

# install jvm for bioformats
USER root
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential \
    default-jdk

USER ${NB_UID}
RUN mkdir ${HOME}/islets/
RUN pip install --no-cache \
        dash \
        cdsdashboards[user] \
        ffmpeg \
        tifffile \
        jupyter_contrib_nbextensions

# copy additional files for cdsdashboards:
COPY ${SERVER_DOCKER_PATH}/jupyter_notebook_config_extra.py /etc/jupyter/
COPY ${SERVER_DOCKER_PATH}/voila.json /etc/jupyter/
COPY ${SERVER_DOCKER_PATH}/github-tokens.sh /usr/local/bin/before-notebook.d/github-tokens.sh
RUN cat /etc/jupyter/jupyter_notebook_config_extra.py >> /etc/jupyter/jupyter_notebook_config.py
RUN rm /etc/jupyter/jupyter_notebook_config_extra.py

# prepare installation of islets:
COPY ./islets/ ${CONDA_DIR}/lib/python3.9/site-packages/islets/
COPY ./app/ ${HOME}/app/
COPY setup.py ${CONDA_DIR}/lib/python3.9/site-packages/islets/setup.py
COPY requirements.txt ${CONDA_DIR}/lib/python3.9/site-packages/islets/requirements.txt
COPY README.md ${CONDA_DIR}/lib/python3.9/site-packages/islets/README.md
COPY VERSION ${CONDA_DIR}/lib/python3.9/site-packages/islets/VERSION

# install islets and enable kernel:
WORKDIR ${CONDA_DIR}/lib/python3.9/site-packages/islets/
USER root
RUN python3 setup.py bdist_wheel
RUN pip3 install "dist/islets-${ISLETS_VERSION}-py3-none-any.whl"
# RUN python -m ipykernel install --display-name "Physio"
RUN /usr/local/bin/fix-permissions "${CONDA_DIR}/lib/python3.9/site-packages/islets/"

USER ${NB_UID}
WORKDIR ${HOME}
