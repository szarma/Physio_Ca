FROM jupyterhub/jupyterhub:3.1

# Copy the JupyterHub configuration in the container
COPY jupyterhub/jupyterhub_config.py .

RUN useradd --create-home \
            --shell /bin/bash admin \
            --password "$(openssl passwd -1 admin)"


# Install dependencies (for culling and spawning)
USER $NB_UID
RUN pip install --no-cache \
    dockerspawner \
    oauthenticator \
    jupyterhub-idle-culler \
    cdsdashboards

RUN mkdir /data
