## Docker Deployment (still in development and to be changed in near future)
JupyterHub can be deployed as a docker image giving users individual sandboxes to work with and easily upgradable containers to developers.

### Prerequisites
To install the CTN-JupyterHub make sure to have:
1. administrative rights
2. internet connection
3. Installed dependencies:
   - docker engine version 19.03.0 or higher
   - docker-compose

### Preparation
To install CTN make sure to clone the repository with its content to one of your folders on the server.
The essential parts for the docker installation beside the notebooks and the islets module are:
- Compose file letting docker know what services to have and how to connect: ./docker-compose.yml
- Environment file defining specific necessary settings: ./.env
- Folder containing the dockerfiles and helpers for image creation: ./docker/

Before installation you need to edit the environment file!
Inside the file there are different parameters to be set:
- COMPOSE_PROJECT_NAME - Pproject name (the created default network will use this value)
- HOST - URL to be used for referring with traefik. **For productive systems change to correct URL!**
- SERVER_IMAGE_NAME - This value defines how the image will be called. This is especially important for the DockerSpawner.
- HUB_CONTAINER_NAME - This value defines how the image of the jupyterhub will be called. This is important for communication between the single-user servers and the hub.
- COMPOSE_PROJECT_NAME

### Important considerations
As you will notice, we use the host PAM authentication, because it is easy to use. For this it is needed to give the JupyterHub container access to the following files:
```
/etc/passwd
/etc/groups
/etc/shadow
```

### SSL configuration
As we wanted to be open to add services to the setup, we use traefik as reverse-proxy. SSL-encryption will be added on the level of the proxy.
Thus, it is needed to provide traefik access to `/etc/ssl` so it can use your custom TLS certificate.
To let traefik know which certificate shall be used, we need to provide access to `/etc/traefik/dynamic` as well.
In this directory you can store a dynamic configuration `certificates.yml` containing:
```
tls:
  stores:
    default:
      defaultCertificate:
        certFile: /path/to/cert.for.yourdomain.org.crt
        keyFile: /path/to/cert.for.yourdomain.org.key
  certificates:
    - certFile: /path/to/cert.for.yourdomain.org.crt
      keyFile: /path/to/cert.for.yourdomain.org.key
      stores:
        - default
```

### Installation
Installation of the images is rather easy.
Once navigated into the cloned directory we need to build the images using the following:
```
docker compose build
```
This will create 3 different images:
- traefik (Webserver to dispatch requests from the global network and between containers)
- jupyterhub (The hub itself for user management)
- jupyterhub_server (The image for a single person calculating server)

To start up the system we just need to type:
```
docker compose up
# optional: add -d as parameter to run it as detached process
```
After this the server should be ready to work.
