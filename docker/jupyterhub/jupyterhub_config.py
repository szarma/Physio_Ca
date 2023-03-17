import os
import sys

import dockerspawner
import tomli

c = get_config()

def get_allowed_images(spawner):
    """
    Gets the images, which the user is allowed to spawn. In general, normal user will get the compiled stable image.
    For developers, there will be an additional unstable image, which has an editable version of islets installed.
    User information is read from the user config at `/etc/jupyterhub/user_config.toml` (read-only).

    Parameters
    ----------
    spawner
        The spawner object. It is used to retrieve the logged in user.

    """

    images = {}
    with open('/etc/jupyterhub/user_config.toml', 'rb') as f:
        toml = tomli.load(f)
        users = [user['name'] for user in toml['users']]
        if spawner.user.name in users:
            user = toml['users'][users.index(spawner.user.name)]
            if user['developer'] == True:
                images = {
                    'Stable (normal productive use)': os.environ['DOCKER_JUPYTER_CONTAINER'] + ':' + os.environ['ISLETS_VERSION'],
                    'Development (Development of new features, UNSTABLE!)': os.environ['DOCKER_JUPYTER_CONTAINER'] + '_dev:' + os.environ['ISLETS_VERSION'],
                }
            else:
                images = {
                    'Stable (normal productive use)': os.environ['DOCKER_JUPYTER_CONTAINER'] + ':' + os.environ['ISLETS_VERSION'],
                }

    return images

def get_admin_users():
    """
    Gets the verified admin users.
    User information is read from the user config at `/etc/jupyterhub/user_config.toml` (read-only).

    """

    admin_users = {}
    with open('/etc/jupyterhub/user_config.toml', 'rb') as f:
        toml = tomli.load(f)
        admin_users = { user['name'] for user in toml['users'] if user['admin'] == True }

    return admin_users

def get_allowed_users():
    """
    Gets the users, which are allowed to sign in.
    User information is read from the user config at `/etc/jupyterhub/user_config.toml` (read-only).

    """

    allowed_users = {}
    with open('/etc/jupyterhub/user_config.toml', 'rb') as f:
        toml = tomli.load(f)
        allowed_users = { user['name'] for user in toml['users'] }

    return allowed_users


class CustomDockerSpawner(dockerspawner.DockerSpawner):
    """
    Customized DockerSpawner to expose port 8050 of the spawned image in addition to the default port 8888.

    Methods
    -------
    create_object()
        Hook, which gets called when the spawner gets created.

    """
    def create_object(self):
        self.extra_create_kwargs['ports'] = {
            "%i/tcp" % self.port: None,
            "%i/tcp" % 8050: None
        }
        return super().create_object()


# Add log support for the jupyterhub.
# It will log information and more severe stuff to "/var/log/jupyterhub.log"
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

c.Authenticator.allowed_users = get_allowed_users()
c.Authenticator.admin_users = get_admin_users()
# c.LocalAuthenticator.create_system_users = True

notebook_dir = os.environ.get('DOCKER_NOTEBOOK_DIR', default='/home/jovyan/work')
c.DockerSpawner.notebook_dir = notebook_dir
c.DockerSpawner.environment = {'URL': os.getenv('HOST', default='ctn2.physiologie.meduniwien.ac.at')}
c.DockerSpawner.volumes = {
    '/data': {
        'bind': '/data',
        'mode': 'z',
    },
    'jupyterhub-user-{username}': notebook_dir,
}
c.DockerSpawner.name_template = "{prefix}-{username}"
#c.DockerSpawner.debug = True
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
