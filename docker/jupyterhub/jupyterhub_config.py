import os
import sys

from cdsdashboards.app import CDS_TEMPLATE_PATHS
from cdsdashboards.hubextension import cds_extra_handlers

c = get_config()

c.JupyterHub.log_level = 'DEBUG'
c.JupyterHub.admin_access = True
c.JupyterHub.hub_ip = '0.0.0.0'
c.JupyterHub.allow_named_servers = True

c.JupyterHub.template_paths = CDS_TEMPLATE_PATHS
c.JupyterHub.extra_handlers = cds_extra_handlers
c.CDSDashboardsConfig.builder_class = 'cdsdashboards.builder.dockerbuilder.DockerBuilder'

c.Authenticator.admin_users = {'admin'}
c.LocalAuthenticator.create_system_users = True

#c.JupyterHub.spawner_class = 'dockerspawner.DockerSpawner'
c.JupyterHub.spawner_class = 'cdsdashboards.hubextension.spawners.variabledocker.VariableDockerSpawner'
c.DockerSpawner.name_template = "{prefix}-{username}-{servername}"
c.DockerSpawner.debug = True
c.DockerSpawner.image = os.environ['DOCKER_JUPYTER_CONTAINER']
c.DockerSpawner.remove = True
c.DockerSpawner.network_name = os.environ['DOCKER_NETWORK_NAME']

c.JupyterHub.services = [
    {
        'name': 'idle-culler',
        'admin': True,
        'command': [
            sys.executable,
            '-m', 'jupyterhub_idle_culler',
            '--timeout=3600'
        ],
    }
]
