FROM jupyter/scipy-notebook:hub-3.1.1

# Switch to root user to install tools to build islets and to run bioformats JVM
USER root

# Install binaries to install and execute the islets module
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential \
    git \
    default-jdk

# Create a directory to copy the islets source into it.
# The directory will have fixed permissions.
RUN mkdir /opt/islets/
RUN /usr/local/bin/fix-permissions /opt/islets/
RUN chown -R ${NB_UID}:users ${HOME}

# Switch back to the default user as we do not need root rights anymore and don't want to give them to the user.
USER ${NB_UID}

# Upgrade pip and install ffmpeg and tiff dependencies
RUN pip install --upgrade --no-cache pip && \
    pip install --upgrade --no-cache ffmpeg tifffile

# Install JupyterDash, the corresponding proxy to stream webapps, and needed labextension
RUN pip install --upgrade --no-cache dash jupyter-dash jupyter-server-proxy
RUN jupyter labextension install jupyterlab-plotly
RUN jupyter lab build

# Install the islets module
WORKDIR /opt/islets/
RUN git clone --branch=docker \
              --depth=1 \
              --single-branch \
              https://github.com/szarma/Physio_Ca.git
WORKDIR /opt/islets/Physio_Ca/
RUN pip install .

# Link to data and set the home directory as default work directory
RUN ln -s /data ${HOME}/work/local_data
WORKDIR ${HOME}/work
