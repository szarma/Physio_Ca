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

USER $NB_UID
RUN conda install --yes \
          -c conda-forge \
          python \
          cython \
          # caiman \
          jupyter_contrib_nbextensions \ 
          dash \
          cdsdashboards-singleuser && \
    conda clean --all -f -y

COPY ${SERVER_DOCKER_PATH}/jupyter_notebook_config_extra.py /etc/jupyter/
COPY ${SERVER_DOCKER_PATH}/voila.json /etc/jupyter/
COPY ${SERVER_DOCKER_PATH}/github-tokens.sh /usr/local/bin/before-notebook.d/github-tokens.sh
RUN cat /etc/jupyter/jupyter_notebook_config_extra.py >> /etc/jupyter/jupyter_notebook_config.py
RUN rm /etc/jupyter/jupyter_notebook_config_extra.py

# prepare installation of islets:
COPY ./islets/ /home/jovyan/.islets/
COPY ./app/ /opt/app/
COPY ${SERVER_DOCKER_PATH}/setup.py /home/jovyan/.islets/setup.py
COPY requirements.txt /home/jovyan/.islets/requirements.txt
COPY README.md /home/jovyan/.islets/README.md
COPY VERSION /home/jovyan/.islets/VERSION

# install islets and enable kernel:
RUN ln -s /opt/app/ /home/jovyan/app
USER root
RUN /opt/conda/bin/pip install -e /home/jovyan/.islets/
USER $NB_UID
# RUN /opt/conda/envs/Physio/bin/python3 -m ipykernel install --user --name Physio --display-name "Physio"
