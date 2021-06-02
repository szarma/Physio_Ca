FROM jupyter/scipy-notebook

# install basic: python environment
RUN python3 -m pip install --upgrade pip
RUN conda install python=3.7
RUN conda install -c conda-forge caiman

# Install CDS components according to https://github.com/ideonate/cdsdashboards/blob/master/docker-images/singleuser-example/containds-all-example/Dockerfile
RUN rm -rf /home/jovyan/work
# ARG JHSINGLENATIVEPROXY_LINE="jhsingle-native-proxy>=0.6.1"
# RUN pip install $JHSINGLENATIVEPROXY_LINE plotlydash-tornado-cmd>=0.0.4 bokeh-root-cmd>=0.0.5 rshiny-server-cmd>=0.0.2 voila-materialstream>=0.2.6
ARG JHSINGLENATIVEPROXY_LINE="jhsingle-native-proxy"
RUN pip install $JHSINGLENATIVEPROXY_LINE plotlydash-tornado-cmd bokeh-root-cmd rshiny-server-cmd voila-materialstream
ARG FRAMEWORKS_LINE="voila streamlit dash bokeh panel holoviews"
RUN pip install $FRAMEWORKS_LINE
# for github repo setting the argument is needed to let docker know where the conent is located
ARG SERVER_DOCKER_PATH
ARG JUPYTER_NOTEBOOK_CONFIG_EXTRA="$SERVER_DOCKER_PATH/jupyter_notebook_config_extra.py"
USER $NB_UID
COPY $SERVER_DOCKER_PATH/jupyter_notebook_config_extra.py /etc/jupyter/
RUN cat /etc/jupyter/jupyter_notebook_config_extra.py >> /etc/jupyter/jupyter_notebook_config.py
RUN rm /etc/jupyter/jupyter_notebook_config_extra.py
COPY ${SERVER_DOCKER_PATH}/voila.json /etc/jupyter/
RUN conda init bash
COPY ${SERVER_DOCKER_PATH}/github-tokens.sh /usr/local/bin/before-notebook.d/github-tokens.sh
RUN conda install -c conda-forge jupyter_contrib_nbextensions

# install islets (not finished yet):
COPY ./islets/ /opt/islets/
