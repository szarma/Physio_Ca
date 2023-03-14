import os
import sys

import dockerspawner
import tomli

c = get_config()

def get_allowed_images(spawner):
    users = []
    images = {}

    with open('/etc/jupyterhub/user_config.toml', 'rb') as f:
        toml = tomli.load(f)
        users = list(toml['users'].keys())
        if spawner.user in users:
            if toml['users'][spawner.user]['admin']:
                images = {
                    'stable': os.environ['DOCKER_JUPYTER_CONTAINER'] + ':' + os.environ['ISLETS_VERSION'],
                    'development': os.environ['DOCKER_JUPYTER_CONTAINER'] + '_dev:' + os.environ['ISLETS_VERSION'],
                }

    return images

class CustomDockerSpawner(dockerspawner.DockerSpawner):
    def create_object(self):
        self.extra_create_kwargs['ports'] = {
            "%i/tcp" % self.port: None,
            "%i/tcp" % 8050: None
        }
        return super().create_object()


c.JupyterHub.log_level = 'INFO'
c.JupyterHub.logging_config = {
    'handlers': {
        'file': {
            'class': 'logging.FileHandler',
            'level': 'INFO',
            'filename': '/var/log/jupyterhub.log',
        }
    },
    'loggers': {
        'jupyterhub_log': {
            'level': 'INFO',
            'handlers': ['console', 'file'],
        }
    }
}

c.JupyterHub.admin_access = True
c.JupyterHub.hub_ip = '0.0.0.0'
c.JupyterHub.hub_connect_ip = os.environ.get('HUB_IP') or 'jupyterhub'
c.JupyterHub.allow_named_servers = False
c.JupyterHub.spawner_class = CustomDockerSpawner

c.Authenticator.admin_users = {'johannes'}
c.LocalAuthenticator.create_system_users = True

notebook_dir = os.environ.get('DOCKER_NOTEBOOK_DIR', default='/home/jovyan/work')
c.DockerSpawner.notebook_dir = notebook_dir
c.DockerSpawner.environment = {'URL': os.getenv('HOST', default='ctn2.physiologie.meduniwien.ac.at')}
c.DockerSpawner.volumes = {
    'jupyterhub-user-{username}': notebook_dir,
    '/data': {
        'bind': '/data',
        'mode': 'rw',
    }
}
c.DockerSpawner.name_template = "{prefix}-{username}"
c.DockerSpawner.debug = True
c.DockerSpawner.allowed_images = get_allowed_images
# c.DockerSpawner.image = os.environ['DOCKER_JUPYTER_CONTAINER'] + ':' + os.environ['ISLETS_VERSION']
c.DockerSpawner.remove = True
c.DockerSpawner.network_name = os.environ['DOCKER_NETWORK_NAME']
c.DockerSpawner.use_internal_ip = True

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
