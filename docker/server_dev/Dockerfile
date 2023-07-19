# Use a separate build stage to clone github content
FROM alpine:latest AS git-puller
WORKDIR /opt/islets/
RUN apk --no-cache add git
RUN git clone --depth=1 https://github.com/szarma/Physio_Ca.git

FROM jupyter/datascience-notebook:hub-3.1.1

# Switch to root user to install tools to build islets and to run bioformats JVM
USER root

RUN apt update --yes && \
    apt upgrade --yes

# Install bottom as resource monitor in the terminal
RUN curl -LO https://github.com/ClementTsang/bottom/releases/download/0.8.0/bottom_0.8.0_amd64.deb
RUN dpkg -i bottom_0.8.0_amd64.deb
RUN rm bottom_0.8.0_amd64.deb

# Install bat
RUN apt update --yes && \
    apt install --yes --no-install-recommends \
    bat

# Install binaries to install and execute the islets module
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential \
    git \
    default-jdk && \
    apt-get clean

# Create a directory to copy the islets source into it.
# The directory will have fixed permissions.
RUN mkdir /opt/islets/
RUN mkdir /data
RUN /usr/local/bin/fix-permissions /opt/islets/
RUN /usr/local/bin/fix-permissions /data
RUN chown -R ${NB_UID}:users ${HOME}

# Switch back to the default user as we do not need root rights anymore and don't want to give them to the user.
USER ${NB_UID}

# Upgrade pip and install ffmpeg and tiff dependencies
RUN pip install --upgrade --no-cache pip && \
    pip install --upgrade --no-cache ffmpeg tifffile black isort

# Install graph-tool
RUN conda install -c conda-forge graph-tool 

# Make bat accessible to the user
RUN mkdir -p ~/.local/bin
RUN ln -s /usr/bin/batcat ~/.local/bin/bat

# Install JupyterDash, the corresponding proxy to stream webapps, and needed labextension
RUN pip install --upgrade --no-cache \
        dash jupyter-dash jupyter-server-proxy \
        jupyterlab-drawio jupyter-archive jupyter-resource-usage jupyterlab-spreadsheet-editor \
        jupyterlab-code-formatter jupyterlab-system-monitor jupyterlab_templates jupyterlab_hdf hdf5plugin
RUN jupyter labextension install jupyterlab-plotly @jupyterlab/hdf5
RUN jupyter lab build
RUN echo "c.JupyterLabTemplates.template_dirs = ['/data/useful/templates']" >> ${HOME}/.jupyter/jupyter_notebook_config.py

# Install the islets module
WORKDIR /opt/islets/
COPY --from=git-puller /opt/islets/ ./
WORKDIR /opt/islets/Physio_Ca/
RUN pip install --editable .
RUN chmod +x /opt/islets/Physio_Ca/scripts/full_process.py

# Link to data and set the home directory as default work directory
WORKDIR ${HOME}/work
RUN ln -s /data ${HOME}/work/local_data
RUN ln -s /opt/islets/Physio_Ca ${HOME}/work/Physio_Ca
