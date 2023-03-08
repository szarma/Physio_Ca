FROM jupyterhub/jupyterhub:3.1.1

# Install dependencies (for culling and spawning)
RUN pip install --upgrade pip && \
    pip install --no-cache \
    dockerspawner \
    oauthenticator \
    jupyter-server-proxy \
    jupyterhub-idle-culler

# Copy the JupyterHub configuration in the container
COPY jupyterhub/jupyterhub_config.py .

RUN mkdir /data
