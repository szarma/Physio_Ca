# Do not forget to pin down the version
FROM jupyterhub/jupyterhub:latest

ARG HUB_DOCKER_PATH

# Copy the JupyterHub configuration in the container
COPY $HUB_DOCKER_PATH/jupyterhub_config.py .

# Install dependencies (for culling and spawning)
RUN pip install \
    dockerspawner \
    oauthenticator \
    jupyterhub-idle-culler \
    cdsdashboards

RUN useradd --create-home --shell /bin/bash admin --password "$(openssl passwd -1 admin)"
